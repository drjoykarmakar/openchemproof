"""Calibrated integral budgets and tentative, solvent-specific contamination alerts."""

from rdkit import Chem

from openchemproof.schema.models import AuditPolicy, CompoundSpec, ModalityResult, NMREvidence
from openchemproof.verifiers.common import flag, prepare

# Approximate ppm at ambient temperature. Assignments require orthogonal peak evidence.
RESIDUAL_PROTON = {"CDCl3": 7.26, "DMSO-d6": 2.50, "CD3OD": 3.31, "D2O": 4.79, "Acetone-d6": 2.05}
RESIDUAL_CARBON = {
    "CDCl3": [77.16],
    "DMSO-d6": [39.52],
    "CD3OD": [49.0],
    "D2O": [],
    "Acetone-d6": [29.84, 206.26],
}
CONTAMINANTS = {
    "CDCl3": {
        "ethyl_acetate": [4.12, 2.05, 1.26],
        "dichloromethane": [5.30],
        "water": [1.56],
        "grease": [1.26, 0.86],
    },
    "DMSO-d6": {"ethyl_acetate": [4.03, 1.99, 1.17], "dichloromethane": [5.76], "water": [3.33]},
    "CD3OD": {"ethyl_acetate": [4.09, 2.01, 1.24], "dichloromethane": [5.49], "water": [4.87]},
    "D2O": {"water": [4.79]},
    "Acetone-d6": {"ethyl_acetate": [4.05, 1.97, 1.20], "dichloromethane": [5.63], "water": [2.84]},
}


def proton_budget(smiles: str) -> tuple[int, int]:
    mol = Chem.AddHs(prepare(smiles).parent)
    total = sum(a.GetAtomicNum() == 1 and a.GetIsotope() not in (2, 3) for a in mol.GetAtoms())
    exchangeable = sum(
        1
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 1
        and a.GetIsotope() not in (2, 3)
        and any(n.GetAtomicNum() in (7, 8, 16) for n in a.GetNeighbors())
    )
    return total - exchangeable, exchangeable


def audit_nmr(compound: CompoundSpec, evidence: NMREvidence, policy: AuditPolicy | None = None) -> ModalityResult:
    policy = policy or AuditPolicy()
    flags = []
    metrics: dict[str, float | str] = {}
    checks = ["exchangeable_proton_budget", "residual_solvent_candidates"]
    try:
        prepared = prepare(compound.proposed_smiles)
        nonexchangeable, exchangeable = proton_budget(compound.proposed_smiles)
    except ValueError:
        return ModalityResult(
            modality="nmr",
            inconclusive=True,
            flags=[flag("NMR_INVALID_STRUCTURE", "WARNING", "nmr", "Cannot calculate a proton budget")],
        )
    if prepared.unresolved_fragments:
        return ModalityResult(
            modality="nmr",
            inconclusive=True,
            flags=[
                flag(
                    "NMR_MIXTURE",
                    "WARNING",
                    "nmr",
                    "Integrals require a defined molecular component",
                )
            ],
        )
    lower = nonexchangeable
    upper = nonexchangeable + exchangeable
    if evidence.solvent in ("D2O", "CD3OD") or evidence.exchangeable_protons_observed is False:
        upper = nonexchangeable
    elif evidence.exchangeable_protons_observed is True:
        lower = upper
    # DMSO commonly preserves OH/NH, but moisture, temperature and exchange kinetics matter.
    metrics.update(
        nonexchangeable_h=float(nonexchangeable),
        exchangeable_h=float(exchangeable),
        expected_integral_lower_h=float(lower),
        expected_integral_upper_h=float(upper),
    )
    candidates: list[str] = []
    for name, positions in CONTAMINANTS[evidence.solvent].items():
        if all(any(abs(peak - position) <= 0.035 for peak in evidence.proton_shifts_ppm) for position in positions):
            candidates.append(name)
    residual = RESIDUAL_PROTON[evidence.solvent]
    if any(abs(peak - residual) <= 0.035 for peak in evidence.proton_shifts_ppm):
        candidates.append("residual_deuterated_solvent")
    for candidate in candidates:
        flags.append(
            flag(
                "SOLVENT_PEAK_CANDIDATE",
                "WARNING",
                "nmr",
                "Peak positions overlap a common solvent/impurity; coincidence alone is not identification",
                candidate=candidate,
                solvent=evidence.solvent,
            )
        )
    inconclusive = evidence.integral_scale is None or not evidence.integration_complete
    if evidence.integral_scale is None:
        flags.append(
            flag(
                "UNCALIBRATED_INTEGRALS",
                "WARNING",
                "nmr",
                "Supply independent integral units per proton; total normalization would hide artifacts",
            )
        )
    if not evidence.integration_complete:
        flags.append(
            flag(
                "INCOMPLETE_INTEGRATION",
                "WARNING",
                "nmr",
                "A partial peak list cannot support a total proton balance",
            )
        )
    if evidence.integral_scale is not None and evidence.integration_complete:
        observed = sum(evidence.proton_integrals) / evidence.integral_scale
        tolerance = max(policy.nmr_absolute_h_tolerance, upper * policy.nmr_relative_h_tolerance)
        metrics["observed_integral_h"] = observed
        checks.append("calibrated_integral_balance")
        if not lower - tolerance <= observed <= upper + tolerance:
            flags.append(
                flag(
                    "PROTON_INTEGRAL_MISMATCH",
                    "CRITICAL_ARTIFACT",
                    "nmr",
                    "Calibrated total integration conflicts with the solvent-specific proton budget",
                    observed_h=observed,
                    lower_h=float(lower),
                    upper_h=float(upper),
                    tolerance_h=tolerance,
                )
            )
            if candidates and observed > upper + tolerance:
                flags.append(
                    flag(
                        "SOLVENT_INTEGRAL_INFLATION",
                        "CRITICAL_ARTIFACT",
                        "nmr",
                        "Excess integration coincides with candidate contaminant peaks; "
                        "contaminant identity remains tentative",
                        candidates=", ".join(candidates),
                    )
                )
    if evidence.carbon_shifts_ppm:
        solvent_carbons = RESIDUAL_CARBON[evidence.solvent]
        analyte_peaks = [
            x
            for x in evidence.carbon_shifts_ppm
            if not any(abs(x - solvent_peak) < 0.6 for solvent_peak in solvent_carbons)
        ]
        # Deduplicate numerical repeats; symmetry/overlap can reduce observed resonances.
        distinct: list[float] = []
        for peak in sorted(analyte_peaks):
            if not distinct or abs(peak - distinct[-1]) > 0.05:
                distinct.append(peak)
        carbons = sum(a.GetAtomicNum() == 6 for a in prepared.parent.GetAtoms())
        metrics.update(
            expected_carbon_atoms=float(carbons),
            observed_distinct_carbon_peaks=float(len(distinct)),
        )
        checks.append("carbon_peak_upper_bound")
        if len(distinct) > carbons:
            flags.append(
                flag(
                    "EXCESS_CARBON_RESONANCES",
                    "WARNING",
                    "nmr",
                    "More carbon resonances than parent carbon atoms; mixtures, rotamers or contaminants possible",
                )
            )
            inconclusive = True
    return ModalityResult(
        modality="nmr",
        flags=flags,
        corrected_metrics=metrics,
        checks_performed=checks,
        inconclusive=inconclusive,
    )
