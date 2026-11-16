"""Geometry-qualified fluorescence consistency checks and comparative IFE correction."""

import math

from scipy.stats import spearmanr

from openchemproof.schema.models import AuditPolicy, ModalityResult, OpticalEvidence
from openchemproof.verifiers.common import flag


def ife_factor(absorbance_ex: float, absorbance_em: float) -> float:
    if not all(math.isfinite(a) and a >= 0 for a in (absorbance_ex, absorbance_em)):
        raise ValueError("Absorbances must be finite and nonnegative")
    exponent = (absorbance_ex + absorbance_em) / 2
    if exponent > 300:
        raise ValueError("Correction exceeds numerical range")
    return float(10.0**exponent)


def audit_optical(evidence: OpticalEvidence, policy: AuditPolicy | None = None) -> ModalityResult:
    policy = policy or AuditPolicy()
    flags = []
    metrics: dict[str, float | str] = {}
    checks = ["stokes_shift", "quantum_yield", "scatter_overlap"]
    inconclusive = False
    shift_nm = evidence.emission_nm - evidence.excitation_nm
    shift_cm = 1e7 / evidence.excitation_nm - 1e7 / evidence.emission_nm
    metrics.update(stokes_shift_nm=shift_nm, stokes_shift_cm=shift_cm)
    if shift_nm <= 0 and not evidence.upconversion_declared:
        flags.append(
            flag(
                "ANTI_STOKES_UNEXPLAINED",
                "CRITICAL_ARTIFACT",
                "optical",
                "Anti-Stokes signal is inconsistent with the ordinary one-photon assay model; "
                "requires independent mechanism evidence",
                stokes_shift_nm=shift_nm,
            )
        )
    if abs(shift_nm) < 15:
        flags.append(
            flag(
                "RAYLEIGH_OVERLAP",
                "CRITICAL_ARTIFACT",
                "optical",
                "Emission within the configured 15 nm excitation exclusion window; "
                "spectral blanks and filter bandwidths are needed",
                separation_nm=abs(shift_nm),
            )
        )
    if evidence.solvent_raman_shift_cm is not None and abs(shift_cm - evidence.solvent_raman_shift_cm) < 150:
        flags.append(
            flag(
                "RAMAN_OVERLAP",
                "WARNING",
                "optical",
                "Emission overlaps a supplied solvent Raman shift; confirm with solvent blank",
            )
        )
        inconclusive = True
    if evidence.observed_quantum_yield is not None and not 0 <= evidence.observed_quantum_yield <= 1:
        flags.append(
            flag(
                "UNPHYSICAL_QUANTUM_YIELD",
                "CRITICAL_ARTIFACT",
                "optical",
                "Quantum yield is outside the 0–1 ordinary fluorescence model",
                observed=evidence.observed_quantum_yield,
            )
        )
    aex, aem = evidence.absorbance_at_ex, evidence.absorbance_at_em
    bex, bem = evidence.baseline_absorbance_at_ex, evidence.baseline_absorbance_at_em
    densities = [aex, aem] + ([bex, bem] if bex is not None and bem is not None else [])
    suitable = evidence.geometry == "right_angle_center" and max(densities) <= policy.ife_max_absorbance
    if not suitable:
        inconclusive = True
        flags.append(
            flag(
                "IFE_OUTSIDE_VALIDITY",
                "WARNING",
                "optical",
                "Center-view right-angle correction requires matched path geometry and moderate OD; "
                "dilution or instrument-specific correction is needed",
            )
        )
    else:
        factor = ife_factor(aex, aem)
        metrics["ife_correction_factor"] = factor
        checks.append("ife_correction")
        if evidence.apparent_fold_change is not None:
            if bex is None or bem is None:
                inconclusive = True
                flags.append(
                    flag(
                        "MISSING_REFERENCE_ABSORBANCE",
                        "WARNING",
                        "optical",
                        "Fold correction requires absorbances of both sample and reference",
                    )
                )
            else:
                corrected = evidence.apparent_fold_change * factor / ife_factor(bex, bem)
                metrics["corrected_fold_change"] = corrected
                apparent = evidence.apparent_fold_change
                high_od = max(densities) > policy.ife_absorbance_trigger
                apparent_response = abs(apparent - 1.0) > policy.ife_baseline_tolerance
                restored = abs(corrected - 1.0) <= policy.ife_baseline_tolerance
                if high_od and apparent_response and restored:
                    flags.append(
                        flag(
                            "IFE_RESPONSE_ARTIFACT",
                            "CRITICAL_ARTIFACT",
                            "optical",
                            "Apparent response is explained by differential inner filtering within "
                            "the declared geometry model",
                            apparent_fold_change=apparent,
                            corrected_fold_change=corrected,
                        )
                    )
    if max(densities) > policy.ife_absorbance_trigger:
        flags.append(
            flag(
                "HIGH_ABSORBANCE",
                "WARNING",
                "optical",
                "Absorbance exceeds the conservative IFE screening threshold",
                max_absorbance=max(densities),
            )
        )
    if evidence.expected_solvatochromism is not None:
        series = evidence.solvent_series
        if len(series) < 3 or len({item.polarity_index for item in series}) < 3:
            inconclusive = True
            flags.append(
                flag(
                    "INSUFFICIENT_SOLVENT_SERIES",
                    "WARNING",
                    "optical",
                    "At least three distinct solvent polarities are required",
                )
            )
        else:
            values = [item.emission_nm for item in series]
            rho = (
                0.0
                if len(set(values)) == 1
                else float(spearmanr([item.polarity_index for item in series], values).statistic)
            )
            metrics["solvatochromic_spearman_rho"] = rho
            checks.append("solvent_trend")
            if (evidence.expected_solvatochromism == "positive" and rho < 0.5) or (
                evidence.expected_solvatochromism == "negative" and rho > -0.5
            ):
                inconclusive = True
                flags.append(
                    flag(
                        "SOLVATOCHROMIC_TREND",
                        "WARNING",
                        "optical",
                        "Observed solvent series does not support the declared trend; "
                        "specific solvation and state changes remain alternatives",
                        rho=rho,
                    )
                )
    return ModalityResult(
        modality="optical",
        flags=flags,
        corrected_metrics=metrics,
        checks_performed=checks,
        inconclusive=inconclusive,
    )
