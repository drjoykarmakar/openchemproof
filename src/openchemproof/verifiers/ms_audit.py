"""Singly charged ESI mass/adduct checks and nominal M/M+1/M+2 isotope envelopes."""

import json
from collections import Counter
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from openchemproof.schema.models import AuditPolicy, CompoundSpec, LCMSEvidence, ModalityResult
from openchemproof.verifiers.common import flag, prepare


@dataclass(frozen=True)
class Adduct:
    mode: str
    multiplier: int
    mass_delta: float
    composition_delta: tuple[tuple[str, int], ...]


ADDUCTS = {
    "[M+H]+": Adduct("positive", 1, 1.007276, (("H", 1),)),
    "[M+NH4]+": Adduct("positive", 1, 18.033826, (("N", 1), ("H", 4))),
    "[M+Na]+": Adduct("positive", 1, 22.989218, (("Na", 1),)),
    "[M+K]+": Adduct("positive", 1, 38.963158, (("K", 1),)),
    "[2M+H]+": Adduct("positive", 2, 1.007276, (("H", 1),)),
    "[M-H]-": Adduct("negative", 1, -1.007276, (("H", -1),)),
    "[M+Cl]-": Adduct("negative", 1, 34.969402, (("Cl", 1),)),
    "[M+FA-H]-": Adduct("negative", 1, 44.998203, (("C", 1), ("H", 1), ("O", 2))),
}

# Representative terrestrial abundances, not sample-specific isotope enrichment.
# Offsets are nominal neutron increments, relative to each element's light isotope.
ISOTOPES: dict[str, tuple[float, float, float]] = {
    "C": (0.9893, 0.0107, 0.0),
    "H": (0.999885, 0.000115, 0.0),
    "N": (0.99636, 0.00364, 0.0),
    "O": (0.99757, 0.00038, 0.00205),
    "S": (0.9499, 0.0075, 0.0425),
    "Cl": (0.7576, 0.0, 0.2424),
    "Br": (0.5069, 0.0, 0.4931),
    "F": (1.0, 0.0, 0.0),
    "P": (1.0, 0.0, 0.0),
    "I": (1.0, 0.0, 0.0),
    "Na": (1.0, 0.0, 0.0),
    "K": (0.932581, 0.000117, 0.067302),
}


def ppm_error(observed: float, theoretical: float) -> float:
    if theoretical <= 0:
        raise ValueError("Theoretical m/z must be positive")
    return (observed - theoretical) / theoretical * 1e6


def theoretical_mz(smiles: str, adduct: str) -> float:
    item = ADDUCTS[adduct]
    mol = prepare(smiles).parent
    return float(rdMolDescriptors.CalcExactMolWt(mol)) * item.multiplier + item.mass_delta


def isotope_envelope(smiles: str, adduct: str = "[M+H]+") -> list[float]:
    """Exact truncated multinomial convolution, normalized to M=1 (not to sum=1).

    Includes all cross terms, dimer stoichiometry, and adduct atoms. Explicitly
    isotope-labeled atoms remain fixed; their natural distribution is not reapplied.
    This is a nominal envelope, not a resolved isotope fine-structure spectrum.
    """
    item = ADDUCTS[adduct]
    mol = Chem.AddHs(prepare(smiles).parent)
    counts = Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetIsotope() == 0)
    counts = Counter({k: v * item.multiplier for k, v in counts.items()})
    for element, number in item.composition_delta:
        counts[element] += number
    if any(number < 0 for number in counts.values()):
        raise ValueError("Adduct removes an atom absent from the unlabeled parent")
    envelope = [1.0, 0.0, 0.0]
    for element, number in sorted(counts.items()):
        if element not in ISOTOPES:
            raise ValueError(f"Natural isotope model unavailable for {element}")
        p0, p1, p2 = ISOTOPES[element]
        r1, r2 = p1 / p0, p2 / p0
        for _ in range(number):
            envelope = [1.0, envelope[1] + r1, envelope[2] + envelope[1] * r1 + r2]
    return envelope


def audit_lcms(compound: CompoundSpec, evidence: LCMSEvidence, policy: AuditPolicy | None = None) -> ModalityResult:
    policy = policy or AuditPolicy()
    flags = []
    metrics: dict[str, float | str] = {}
    checks = ["adduct_mass_matching"]
    inconclusive = False
    try:
        prepared = prepare(compound.proposed_smiles)
    except ValueError:
        return ModalityResult(
            modality="lcms",
            inconclusive=True,
            flags=[
                flag(
                    "MS_INVALID_STRUCTURE",
                    "WARNING",
                    "lcms",
                    "Cannot calculate masses for invalid structure",
                )
            ],
        )
    if prepared.unresolved_fragments or Chem.GetFormalCharge(prepared.parent) != 0:
        return ModalityResult(
            modality="lcms",
            inconclusive=True,
            flags=[
                flag(
                    "MS_UNSUPPORTED_PARENT",
                    "WARNING",
                    "lcms",
                    "Need a unique neutral parent for v0.1 adducts",
                )
            ],
        )
    possible = {
        name: theoretical_mz(compound.proposed_smiles, name)
        for name, item in ADDUCTS.items()
        if item.mode == evidence.ionization_mode
    }
    possible = {name: mass for name, mass in possible.items() if mass > 0}
    if not possible:
        return ModalityResult(modality="lcms", inconclusive=True)
    best = min(possible, key=lambda name: abs(ppm_error(evidence.observed_mz, possible[name])))
    claimed = evidence.observed_adduct
    if claimed is not None and claimed not in ADDUCTS:
        flags.append(flag("UNKNOWN_ADDUCT", "WARNING", "lcms", "Unsupported reported adduct", adduct=claimed))
        inconclusive = True
        claimed = None
    elif claimed is not None and ADDUCTS[claimed].mode != evidence.ionization_mode:
        flags.append(
            flag(
                "ADDUCT_POLARITY_MISMATCH",
                "CRITICAL_ARTIFACT",
                "lcms",
                "Adduct charge contradicts ionization polarity",
            )
        )
        claimed = None
    selected = claimed or best
    theoretical = possible[selected]
    ppm = ppm_error(evidence.observed_mz, theoretical)
    metrics.update(theoretical_mz=theoretical, ppm_error=ppm, matched_adduct=best, selected_adduct=selected)
    if abs(ppm) > policy.hrms_ppm_tolerance + 1e-8:
        flags.append(
            flag(
                "HRMS_MASS_MISMATCH",
                "CRITICAL_ARTIFACT",
                "lcms",
                "Observed m/z exceeds the declared HRMS tolerance",
                ppm_error=ppm,
                tolerance_ppm=policy.hrms_ppm_tolerance,
            )
        )
    if (
        claimed is not None
        and best != claimed
        and abs(ppm_error(evidence.observed_mz, possible[best])) <= policy.alternative_ppm_tolerance
    ):
        flags.append(
            flag(
                "ADDUCT_MISASSIGNMENT",
                "CRITICAL_ARTIFACT",
                "lcms",
                "Observed peak matches another parent adduct",
                claimed=claimed,
                alternative=best,
            )
        )
    alternatives: list[dict[str, float | str]] = []
    for material in evidence.starting_materials:
        try:
            candidate = prepare(material.proposed_smiles)
        except ValueError:
            flags.append(
                flag(
                    "INVALID_STARTING_MATERIAL",
                    "WARNING",
                    "lcms",
                    "Starting-material structure cannot be parsed",
                )
            )
            inconclusive = True
            continue
        if candidate.unresolved_fragments or Chem.GetFormalCharge(candidate.parent) != 0:
            inconclusive = True
            continue
        for name, item in ADDUCTS.items():
            if item.mode != evidence.ionization_mode:
                continue
            value = theoretical_mz(material.proposed_smiles, name)
            if value <= 0:
                continue
            error = ppm_error(evidence.observed_mz, value)
            if abs(error) <= policy.alternative_ppm_tolerance:
                alternatives.append(dict(compound_id=material.compound_id, adduct=name, ppm=error))
    if alternatives:
        code = "STARTING_MATERIAL_MISASSIGNMENT" if abs(ppm) > policy.hrms_ppm_tolerance else "ISOBARIC_AMBIGUITY"
        severity = "CRITICAL_ARTIFACT" if code.startswith("STARTING") else "WARNING"
        flags.append(
            flag(
                code,
                severity,
                "lcms",
                "Peak also matches a supplied starting-material adduct; identity needs orthogonal evidence",
                candidates=json.dumps(alternatives, sort_keys=True),
            )
        )
        if severity == "WARNING":
            inconclusive = True
    if evidence.relative_isotope_intensities:
        try:
            expected = isotope_envelope(compound.proposed_smiles, selected)
        except ValueError as exc:
            flags.append(flag("ISOTOPE_MODEL_UNAVAILABLE", "WARNING", "lcms", str(exc)))
            inconclusive = True
        else:
            observed = [x / evidence.relative_isotope_intensities[0] for x in evidence.relative_isotope_intensities]
            checks.append("isotope_envelope")
            for i in (1, 2):
                metrics[f"isotope_M+{i}_expected"] = expected[i]
                metrics[f"isotope_M+{i}_observed"] = observed[i]
                tolerance = max(
                    policy.isotope_absolute_tolerance,
                    expected[i] * policy.isotope_relative_tolerance,
                )
                if abs(observed[i] - expected[i]) > tolerance:
                    flags.append(
                        flag(
                            "ISOTOPE_PATTERN_MISMATCH",
                            "CRITICAL_ARTIFACT",
                            "lcms",
                            "Nominal isotope envelope conflicts with the proposed ion composition",
                            isotope=f"M+{i}",
                            expected=expected[i],
                            observed=observed[i],
                            tolerance=tolerance,
                        )
                    )
    else:
        inconclusive = True
        flags.append(
            flag(
                "MISSING_ISOTOPE_EVIDENCE",
                "WARNING",
                "lcms",
                "No MS1 M/M+1/M+2 envelope supplied; MS2 fragment peaks cannot substitute",
            )
        )
    return ModalityResult(
        modality="lcms",
        flags=flags,
        corrected_metrics=metrics,
        checks_performed=checks,
        inconclusive=inconclusive,
    )
