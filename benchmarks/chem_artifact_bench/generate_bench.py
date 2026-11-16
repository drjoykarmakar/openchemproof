"""Seeded integration fixtures with construction-time labels independent of predictions.

This is NOT an experimental validation cohort. Public-source acquisitions live
separately in external_data/ and retain their measured-versus-derived provenance.
"""

import argparse
import hashlib
import json
import random
from collections import Counter
from pathlib import Path
from typing import Literal

from rdkit.Chem import rdMolDescriptors

from benchmarks.chem_artifact_bench.models import BenchmarkCase, Category
from openchemproof.schema.export import canonical_json, write_parquet
from openchemproof.schema.models import (
    AssayContext,
    CompoundSpec,
    LCMSEvidence,
    NegativeResultRecord,
    NMREvidence,
    OpticalEvidence,
    SpectralEvidence,
    Status,
    VerificationRecord,
)
from openchemproof.verifiers.common import compound_from_smiles, prepare
from openchemproof.verifiers.ms_audit import isotope_envelope, theoretical_mz
from openchemproof.verifiers.nmr_audit import proton_budget

DEFAULT_DATA = Path("benchmarks/chem_artifact_bench/data")


def clean_ms(compound: CompoundSpec, rng: random.Random, adduct: str = "[M+H]+") -> LCMSEvidence:
    return LCMSEvidence(
        observed_mz=theoretical_mz(compound.proposed_smiles, adduct) * (1 + rng.uniform(-1.5, 1.5) * 1e-6),
        ionization_mode="positive",
        observed_adduct=adduct,
        relative_isotope_intensities=isotope_envelope(compound.proposed_smiles, adduct),
    )


def generate_cases(seed: int = 42) -> list[BenchmarkCase]:
    rng = random.Random(seed)
    cases: list[BenchmarkCase] = []

    def add(
        category: Category,
        subtype: str,
        compound: CompoundSpec,
        evidence: SpectralEvidence,
        status: Status = "ARTIFACT_DETECTED",
        codes: list[str] | None = None,
        negative: NegativeResultRecord | None = None,
    ) -> None:
        token = hashlib.sha256(rng.randbytes(32)).hexdigest()[:20]
        assay_type: Literal["fluorescence_probe", "synthesis_qc"] = (
            "fluorescence_probe" if evidence.optical is not None else "synthesis_qc"
        )
        record = VerificationRecord(
            record_id=f"rec-{token}",
            compound=compound,
            assay=AssayContext(assay_type=assay_type, solvent="DMSO", temperature_c=25.0),
            evidence=evidence,
            negative_result=negative,
        )
        cases.append(
            BenchmarkCase(
                case_id=f"case-{token}",
                category=category,
                subtype=subtype,
                expected_status=status,
                expected_flag_codes=codes or [],
                record=record,
            )
        )

    # 40 inversions, 20 missing stereo, 20 E/Z swaps, 20 explicit salt-mass assertions.
    stereo_pairs = [
        ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"),
        ("C[C@H](O)c1ccccc1", "C[C@@H](O)c1ccccc1"),
        ("C[C@H](O)C(=O)O", "C[C@@H](O)C(=O)O"),
        ("CC[C@H](O)C", "CC[C@@H](O)C"),
    ]
    for i in range(100):
        if i < 40:
            actual, reference = stereo_pairs[i % 4]
            compound = compound_from_smiles(f"m{i}", actual).model_copy(update={"reference_smiles": reference})
            subtype, code = "stereo_inversion", "STEREO_MISMATCH"
        elif i < 60:
            compound = compound_from_smiles(f"m{i}", "CC(O)C(=O)O").model_copy(update={"expected_stereocenters": 1})
            subtype, code = "missing_stereo", "MISSING_STEREOCHEMISTRY"
        elif i < 80:
            compound = compound_from_smiles(f"m{i}", "C/C=C/C").model_copy(update={"reference_smiles": "C/C=C\\C"})
            subtype, code = "ez_inversion", "STEREO_MISMATCH"
        else:
            compound = compound_from_smiles(f"m{i}", "C[NH2+]C.[Cl-]")
            salt_mass = float(rdMolDescriptors.CalcExactMolWt(prepare(compound.proposed_smiles).original))
            compound = compound.model_copy(update={"monoisotopic_mass": salt_mass, "salt_form": "hydrochloride"})
            subtype, code = "salt_parent_mass", "SALT_MASS_MISMATCH"
        add("STEREO_TAUTOMER_TRAP", subtype, compound, SpectralEvidence(), codes=[code])

    ms_smiles = ["CCO", "CCOC(=O)C", "CC(=O)Nc1ccccc1", "COc1ccccc1", "CCc1ccccc1"]
    for i in range(100):
        compound = compound_from_smiles(f"m{i}", ms_smiles[i % len(ms_smiles)])
        if i < 60:
            actual_adduct = "[M+Na]+" if i < 40 else "[2M+H]+"
            ms = clean_ms(compound, rng, actual_adduct).model_copy(update={"observed_adduct": "[M+H]+"})
            subtype, code = "adduct" if i < 40 else "dimer", "ADDUCT_MISASSIGNMENT"
        elif i < 80:
            compound = compound_from_smiles(f"m{i}", "Clc1ccccc1" if i % 2 else "Brc1ccccc1")
            ms = clean_ms(compound, rng)
            ms = ms.model_copy(update={"relative_isotope_intensities": [1.0, ms.relative_isotope_intensities[1], 0.01]})
            subtype, code = "missing_halogen_signature", "ISOTOPE_PATTERN_MISMATCH"
        else:
            starting = compound_from_smiles("starting", "CCO")
            compound = compound_from_smiles(f"m{i}", "CC(=O)Nc1ccccc1")
            ms = clean_ms(starting, rng, "[2M+H]+" if i % 2 else "[M+Na]+").model_copy(
                update={"observed_adduct": "[M+H]+", "starting_materials": [starting]}
            )
            subtype, code = "starting_material", "STARTING_MATERIAL_MISASSIGNMENT"
        add("LCMS_ADDUCT_DIMER_TRAP", subtype, compound, SpectralEvidence(lcms=ms), codes=[code])

    for i in range(100):
        compound = compound_from_smiles(f"m{i}", "O=c1ccc2ccccc2o1")
        ex = rng.uniform(350, 500)
        aex, aem = rng.uniform(0.3, 0.8), rng.uniform(0.1, 0.3)
        bex, bem = 0.01, 0.01
        if 35 <= i < 70:
            aex, aem, bex, bem = bex, bem, aex, aem
        # Independent algebraic synthetic attenuation; no engine calls generate labels.
        fold = 10 ** (-(aex + aem - bex - bem) / 2)
        optical = OpticalEvidence(
            excitation_nm=ex,
            emission_nm=ex + 70,
            absorbance_at_ex=aex,
            absorbance_at_em=aem,
            baseline_absorbance_at_ex=bex,
            baseline_absorbance_at_em=bem,
            apparent_fold_change=fold,
            observed_quantum_yield=0.5,
        )
        subtype, code = "ife_turn_off" if i < 35 else "ife_turn_on", "IFE_RESPONSE_ARTIFACT"
        if 70 <= i < 90:
            optical = optical.model_copy(update={"emission_nm": ex + rng.uniform(1, 10), "apparent_fold_change": None})
            subtype, code = "scatter_overlap", "RAYLEIGH_OVERLAP"
        elif i >= 90:
            optical = optical.model_copy(
                update={
                    "observed_quantum_yield": rng.uniform(1.1, 2.0),
                    "apparent_fold_change": None,
                }
            )
            subtype, code = "invalid_qy", "UNPHYSICAL_QUANTUM_YIELD"
        add(
            "OPTICAL_IFE_ARTIFACT",
            subtype,
            compound,
            SpectralEvidence(optical=optical),
            codes=[code],
        )

    for i in range(100):
        compound = compound_from_smiles(f"m{i}", "CCO" if i % 2 else "OCCO")
        nonex, exchangeable = proton_budget(compound.proposed_smiles)
        if i < 40:
            nmr = NMREvidence(
                solvent="D2O" if i % 2 else "CD3OD",
                proton_shifts_ppm=[3.6],
                proton_integrals=[float(nonex + exchangeable)],
                integral_scale=1.0,
            )
            subtype, code = "exchange_rules", "PROTON_INTEGRAL_MISMATCH"
        else:
            shifts, extra = (
                ([4.12, 2.05, 1.26], [2.0, 3.0, 3.0]) if i < 70 else (([5.30], [2.0]) if i % 2 else ([1.56], [2.0]))
            )
            nmr = NMREvidence(
                solvent="CDCl3",
                proton_shifts_ppm=[3.65, *shifts],
                proton_integrals=[float(nonex + exchangeable), *extra],
                integral_scale=1.0,
            )
            subtype, code = "residual_contaminant", "SOLVENT_INTEGRAL_INFLATION"
        add("NMR_SOLVENT_INTEGRAL_TRAP", subtype, compound, SpectralEvidence(nmr=nmr), codes=[code])

    controls = [
        "CCO",
        "CCOC(=O)C",
        "CC(=O)Nc1ccccc1",
        "COc1ccccc1",
        "Clc1ccccc1",
        "Brc1ccccc1",
        "CCN",
        "CCOC",
        "OCCO",
        "CCCCO",
        "CC(=O)O",
        "C1CCCCC1",
    ]
    for i in range(100):
        compound = compound_from_smiles(f"m{i}", controls[i % len(controls)])
        if i < 60:
            ev = SpectralEvidence(lcms=clean_ms(compound, rng))
            if i < 10:
                compound = compound_from_smiles(f"m{i}", "Oc1ccccn1").model_copy(
                    update={"reference_smiles": "O=c1cccc[nH]1"}
                )
                ev = SpectralEvidence(lcms=clean_ms(compound, rng))
            add("VERIFIED_TRUE_AND_NEGATIVE", "clean_consistency_control", compound, ev, "VERIFIED")
        else:
            negative = NegativeResultRecord(
                experiment_id=f"exp-{i}",
                outcome="failed_synthesis",
                observation="No product above the stated analytical detection limit.",
                protocol_reference="synthetic-fixture://protocol/v1",
                evidence_references=["synthetic-fixture://control-and-chromatogram"],
                detection_limit_description="Synthetic LCMS target detection limit 0.1 micromolar",
                controls_passed=True,
                measured_yield_percent=0.0,
            )
            add(
                "VERIFIED_TRUE_AND_NEGATIVE",
                "documented_synthetic_failure",
                compound,
                SpectralEvidence(),
                "NEGATIVE_RESULT",
                negative=negative,
            )
    rng.shuffle(cases)
    return cases


def write_dataset(destination: Path = DEFAULT_DATA, seed: int = 42) -> list[BenchmarkCase]:
    cases = generate_cases(seed)
    destination.mkdir(parents=True, exist_ok=True)
    path = destination / "chem_artifact_bench_v1.jsonl"
    path.write_text("".join(canonical_json(case.model_dump(mode="json")) + "\n" for case in cases))
    write_parquet(cases, destination / "chem_artifact_bench_v1.parquet")
    manifest = {
        "seed": seed,
        "count": len(cases),
        "origin": "synthetic_unit_fixture",
        "categories": dict(Counter(case.category for case in cases)),
        "sha256_jsonl": hashlib.sha256(path.read_bytes()).hexdigest(),
        "warning": "Internal consistency fixtures share chemistry primitives with the engine; "
        "100% agreement does not establish external validity or LLM performance.",
    }
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return cases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    cases = write_dataset(args.output_dir, args.seed)
    print(f"Generated {len(cases)} synthetic integration fixtures in {args.output_dir}")


if __name__ == "__main__":
    main()
