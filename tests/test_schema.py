import json
from concurrent.futures import ThreadPoolExecutor

import pytest
from pydantic import ValidationError

from openchemproof.schema.export import persist_bundle, read_parquet, verify_bundle, write_parquet
from openchemproof.schema.models import AuditBundle, LCMSEvidence, NMREvidence, OpticalEvidence, VerificationRecord
from openchemproof.verifiers.pipeline import run_verification


@pytest.mark.parametrize("bad", ["123.4", float("nan"), float("inf"), -2.0, 0.0, True])
def test_strict_mass_values(bad):
    with pytest.raises(ValidationError):
        LCMSEvidence(observed_mz=bad, ionization_mode="positive")


def test_extra_frozen_and_lists(ethanol):
    with pytest.raises(ValidationError):
        ethanol.record_id = "edited"
    data = ethanol.model_dump()
    data["unknown"] = 1
    with pytest.raises(ValidationError):
        VerificationRecord.model_validate(data)
    with pytest.raises(ValidationError):
        NMREvidence(solvent="CDCl3", proton_shifts_ppm=[1.0], proton_integrals=[])
    with pytest.raises(ValidationError):
        LCMSEvidence(observed_mz=47.0, ionization_mode="positive", relative_isotope_intensities=[1.0, 0.1])


def test_bad_qy_is_representable():
    ev = OpticalEvidence(
        excitation_nm=400.0,
        emission_nm=500.0,
        absorbance_at_ex=0.01,
        absorbance_at_em=0.01,
        observed_quantum_yield=-0.1,
    )
    assert ev.observed_quantum_yield == -0.1


def test_parquet_lossless(ethanol, tmp_path):
    path = tmp_path / "records.parquet"
    write_parquet([ethanol], path)
    assert read_parquet(path, VerificationRecord) == [ethanol]


def test_audit_idempotent_tamper_evident(ethanol, tmp_path):
    verdict = run_verification(ethanol, audit_dir=tmp_path)
    second = run_verification(ethanol, audit_dir=tmp_path)
    assert verdict == second
    path = tmp_path / f"{verdict.audit_sha256}.json"
    bundle = AuditBundle.model_validate_json(path.read_text())
    assert verify_bundle(bundle)
    with ThreadPoolExecutor(max_workers=4) as executor:
        assert len(list(executor.map(lambda _: persist_bundle(bundle, tmp_path), range(16)))) == 16
    changed = json.loads(path.read_text())
    changed["input"]["record_id"] = "tampered"
    assert not verify_bundle(AuditBundle.model_validate(changed))
    path.chmod(0o644)
    path.write_text("modified")
    with pytest.raises(ValueError, match="modified"):
        persist_bundle(bundle, tmp_path)


def test_nested_mutation_revalidated(ethanol):
    ethanol.source_references.append(123)
    with pytest.warns(UserWarning, match="serializer warnings"), pytest.raises(ValidationError):
        run_verification(ethanol)
