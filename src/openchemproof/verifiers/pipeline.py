"""Deterministic status precedence with an auditable, explicitly bounded interpretation."""

import os
from pathlib import Path

from openchemproof import __version__
from openchemproof.schema.export import make_bundle, persist_bundle
from openchemproof.schema.models import (
    AuditPolicy,
    AuditVerdict,
    Modality,
    ModalityResult,
    Status,
    VerificationRecord,
)
from openchemproof.verifiers.common import flag
from openchemproof.verifiers.ms_audit import audit_lcms
from openchemproof.verifiers.nmr_audit import audit_nmr
from openchemproof.verifiers.optical_audit import audit_optical
from openchemproof.verifiers.structure_audit import audit_structure


def run_verification(
    record: VerificationRecord,
    policy: AuditPolicy | None = None,
    audit_dir: Path | None = None,
    modalities: tuple[Modality, ...] | None = None,
) -> AuditVerdict:
    """Verify supplied evidence and persist the full input/policy/result snapshot.

    VERIFIED means no implemented consistency check failed on supplied evidence.
    It does not establish synthesis success, biological activity or molecular identity.
    """
    # Revalidate to catch mutation of nested lists/dicts in otherwise frozen models.
    record = VerificationRecord.model_validate_json(record.model_dump_json())
    policy = policy or AuditPolicy()
    requested = modalities or ("structure", "lcms", "optical", "nmr")
    results: list[ModalityResult] = []
    if "structure" in requested:
        results.append(audit_structure(record.compound))
    ev = record.evidence
    if "lcms" in requested and ev.lcms is not None:
        results.append(audit_lcms(record.compound, ev.lcms, policy))
    if "optical" in requested and ev.optical is not None:
        results.append(audit_optical(ev.optical, policy))
    if "nmr" in requested and ev.nmr is not None:
        results.append(audit_nmr(record.compound, ev.nmr, policy))
    flags = [item for result in results for item in result.flags]
    metrics = {
        f"{result.modality}.{key}": value for result in results for key, value in result.corrected_metrics.items()
    }
    checks = [f"{result.modality}.{check}" for result in results for check in result.checks_performed]
    critical = any(item.severity == "CRITICAL_ARTIFACT" for item in flags)
    incomplete = any(result.inconclusive for result in results)
    has_spectrum = any(result.modality != "structure" for result in results)
    status: Status
    negative = record.negative_result
    if critical:
        status = "ARTIFACT_DETECTED"
    elif negative is not None:
        if negative.controls_passed and not incomplete:
            status = "NEGATIVE_RESULT"
            flags.append(
                flag(
                    "DOCUMENTED_NEGATIVE",
                    "INFO",
                    "provenance",
                    "Failure is documented with protocol, evidence references and detection limits; "
                    "the external evidence has not been independently authenticated",
                )
            )
        else:
            status = "INCONCLUSIVE"
            flags.append(
                flag(
                    "NEGATIVE_CONTROLS_UNRESOLVED",
                    "WARNING",
                    "provenance",
                    "Controls or analytical consistency are unresolved for the negative result",
                )
            )
    elif incomplete or not has_spectrum:
        status = "INCONCLUSIVE"
        if not has_spectrum:
            flags.append(
                flag(
                    "NO_SPECTRAL_EVIDENCE",
                    "WARNING",
                    "provenance",
                    "Structure metadata alone cannot establish an analytical verification",
                )
            )
    else:
        status = "VERIFIED"
    # Declared rule certainty, not empirical/calibrated confidence in chemical truth.
    score = {
        "ARTIFACT_DETECTED": 1.0,
        "NEGATIVE_RESULT": 0.8,
        "INCONCLUSIVE": 0.0,
        "VERIFIED": 0.8,
    }[status]
    verdict = AuditVerdict(
        record_id=record.record_id,
        status=status,
        confidence_score=score,
        flags=flags,
        corrected_metrics=metrics,
        checks_performed=checks,
        deterministic_engine_version=__version__,
    )
    bundle = make_bundle(record, policy, verdict)
    directory = audit_dir or Path(os.environ.get("OPENCHEMPROOF_AUDIT_DIR", ".openchemproof/audits"))
    persist_bundle(bundle, directory)
    return bundle.verdict
