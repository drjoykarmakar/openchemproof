"""Curate and freeze a source-linked pilot BEFORE calling the verification engine.

This script intentionally imports no OpenChemProof verifier or schema modules.
Expected labels are curator-authored, not copied from engine output.
"""

# Keep the frozen curator descriptions and source strings intact on one line.
# ruff: noqa: E501

import hashlib
import json
import re
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASES = []
OPTICAL = []


def read(name):
    return [json.loads(line) for line in (ROOT / "external_data" / name).read_text().splitlines()]


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source(row, family):
    path = HERE / "sources" / f"{row['record_id']}.json"
    write(path, row)
    if family == "optical":
        url = row["source"]["reference_url"]
        license_ = row["source"]["license"]
    elif family == "lcms":
        url = row["source"]["url"]
        license_ = row["source"]["license_literal"]
    elif family == "nmr":
        url = row["source"]["record_url"]
        license_ = "NMRShiftDB2 database license; see included license"
    else:
        url = row["source_url"]
        license_ = row.get("source_license", row.get("license"))
    return {
        "record_id": row["record_id"],
        "family": family,
        "url": url,
        "license": license_,
        "snapshot": path.relative_to(HERE).as_posix(),
        "snapshot_sha256": sha(path),
    }


def record(cid, smiles, solvent="not reported", assay="synthesis_qc"):
    return {
        "schema_version": "1.0",
        "record_id": cid,
        "compound": {"compound_id": cid, "proposed_smiles": smiles},
        "assay": {"assay_type": assay, "solvent": solvent},
        "evidence": {},
        "source_references": [],
    }


def add(
    rec,
    category,
    title,
    sources,
    kind,
    original,
    modification,
    expected,
    rationale,
    flags=(),
    limitations="",
    assumptions=None,
):
    rec["source_references"] = [s["url"] for s in sources]
    path = HERE / "cases" / f"{rec['record_id']}.json"
    write(path, rec)
    CASES.append(
        {
            "case_id": rec["record_id"],
            "category": category,
            "title": title,
            "evidence_kind": kind,
            "sources": sources,
            "original_evidence": original,
            "deliberate_modification": modification,
            "expected_status": expected,
            "expected_reason": rationale,
            "required_flags": list(flags),
            "limitations": limitations,
            "assumptions": assumptions or {},
            "input_file": path.relative_to(HERE).as_posix(),
            "input_sha256": sha(path),
            "human_review_status": "Pending",
            "human_review_notes": "",
        }
    )


def ord_smiles(row):
    return next(v["value"] for v in row["product_identifiers"] if v["type"] == "SMILES")


def nmr_text(row):
    return row["identity_analyses"]["NMR of product_1"]["data"]["peaks"]["string_value"]


def parse_nmr(row):
    # Preserve source counts. No rescaling to the molecular formula.
    txt = nmr_text(row)
    proton = re.split(r"13\s*C", txt, maxsplit=1)[0]
    assert "CDCl3" in proton
    peaks = []
    pattern = r"(-?\d+\.\d+)(?:\s*[-–]\s*(-?\d+\.\d+))?\s*\(([^()]*)\)"
    for match in re.finditer(pattern, proton):
        count = re.search(r"\b(\d+(?:\.\d+)?)\s*H\b", match[3])
        if not count:
            continue
        ppm = (float(match[1]) + float(match[2])) / 2 if match[2] else float(match[1])
        mult = match[3].split(",")[0].strip()
        peaks.append((ppm, float(count[1]), mult if mult in ["s", "d", "t", "q", "m", "br"] else "unknown"))
    assert peaks and len(peaks) == len(re.findall(r"\b\d+(?:\.\d+)?\s*H\b", proton)) - 1
    return {
        "solvent": "CDCl3",
        "proton_shifts_ppm": [p[0] for p in peaks],
        "proton_integrals": [p[1] for p in peaks],
        "multiplicities": [p[2] for p in peaks],
        "integral_scale": 1.0,
        "integration_complete": True,
    }


def main():
    if (HERE / "expectations.lock.json").exists():
        raise SystemExit("Frozen pilot exists. Create a new version instead of overwriting its expectations.")
    mass = read("lcms/massbank_100.jsonl")
    dyes = read("optical/records.jsonl")
    nmrs = read("nmr/nmrshiftdb2_100.jsonl")
    stereo = read("structure/records.jsonl")
    ords = read("ord/records.jsonl")

    cat = "LCMS_ADDUCT_DIMER_TRAP"
    src = mass[0]
    s = src["compound"]["smiles_canonical"]
    p = min(src["spectrum"]["peaks"], key=lambda p: abs(p["mz"] - src["acquisition"]["precursor_mz_source_reported"]))
    assert abs(p["mz"] - src["acquisition"]["precursor_mz_source_reported"]) < 0.005
    r = record("MS01", s, "water / methanol gradient")
    r["evidence"]["lcms"] = {"observed_mz": p["mz"], "ionization_mode": "positive", "observed_adduct": "[M+H]+"}
    add(
        r,
        cat,
        "Measured surviving precursor; no MS1 envelope",
        [source(src, "lcms")],
        "source_only",
        f"Experimental MS2 peak {p['mz']} m/z; source adduct [M+H]+. Precursor metadata {src['acquisition']['precursor_mz_source_reported']} is not independent measured mass.",
        "None. Use the measured surviving precursor peak; omit missing MS1 isotope evidence.",
        "INCONCLUSIVE",
        "One accurate surviving precursor peak cannot replace a measured MS1 isotope envelope.",
        ["MISSING_ISOTOPE_EVIDENCE"],
        "The peak is from processed MS2 data. Fragment matching and chromatography are outside this pilot.",
    )
    for cid, idx, ion, offset, factor in [("MS02", 1, "[M+Na]+", 22.989218, 1), ("MS03", 2, "[2M+H]+", 1.007276, 2)]:
        src = mass[idx]
        s = src["compound"]["smiles_canonical"]
        mz = factor * Descriptors.ExactMolWt(Chem.MolFromSmiles(s)) + offset
        r = record(cid, s)
        r["evidence"]["lcms"] = {"observed_mz": mz, "ionization_mode": "positive", "observed_adduct": "[M+H]+"}
        add(
            r,
            cat,
            f"Simulated {ion} mislabeled [M+H]+",
            [source(src, "lcms")],
            "source_based_simulation",
            f"Source structure {s}; experimental source spectrum retained in snapshot.",
            f"Replace m/z with formula-derived {ion} value {mz:.8f}; claim [M+H]+.",
            "ARTIFACT_DETECTED",
            "The asserted protonated monomer cannot have this m/z. The deliberately generated alternative adduct explains the discrepancy.",
            ["ADDUCT_MISASSIGNMENT"],
            "Simulated ion, not an observed sodium-adduct or dimer spectrum.",
            {"ion_mass_offset_Da": offset, "molecule_multiplier": factor, "simulated_mz": mz},
        )
    src = next(x for x in mass if x["record_id"] == "MSBNK-Eawag-EA006901")
    s = src["compound"]["smiles_canonical"]
    mol = Chem.MolFromSmiles(s)
    r = record("MS04", s)
    r["evidence"]["lcms"] = {
        "observed_mz": Descriptors.ExactMolWt(mol) + 1.007276,
        "ionization_mode": "positive",
        "observed_adduct": "[M+H]+",
        "relative_isotope_intensities": [1.0, 0.011 * sum(a.GetSymbol() == "C" for a in mol.GetAtoms()), 0.0],
    }
    add(
        r,
        cat,
        "Chlorinated molecule with simulated missing M+2",
        [source(src, "lcms")],
        "source_based_simulation",
        "Source structure contains one Cl and one S; no measured MS1 isotope envelope is available.",
        "Construct a protonated exact mass and artificial [M,M+1,M+2] intensities; force M+2 to zero.",
        "ARTIFACT_DETECTED",
        "A natural-abundance chlorinated ion cannot have a truly absent M+2 peak at the stipulated measurement precision.",
        ["ISOTOPE_PATTERN_MISMATCH"],
        "Noise-free simulated envelope. Does not establish experimental isotope sensitivity.",
    )
    a = next(x for x in mass if x["record_id"] == "MSBNK-Eawag-EA009401")
    b = next(x for x in mass if x["record_id"] == "MSBNK-Eawag-EA014401")
    s = a["compound"]["smiles_canonical"]
    r = record("MS05", s)
    r["evidence"]["lcms"] = {
        "observed_mz": Descriptors.ExactMolWt(Chem.MolFromSmiles(s)) + 1.007276,
        "ionization_mode": "positive",
        "observed_adduct": "[M+H]+",
        "starting_materials": [{"compound_id": b["record_id"], "proposed_smiles": b["compound"]["smiles_canonical"]}],
    }
    add(
        r,
        cat,
        "Real constitutional isomers share exact mass",
        [source(a, "lcms"), source(b, "lcms")],
        "source_based_simulation",
        "Two MassBank structures have formula C16H16N2O4 and distinct connectivity.",
        "Use calculated [M+H]+ mass and put the second isomer in the candidate field (named starting_materials by the schema).",
        "INCONCLUSIVE",
        "Same-formula isomers share exact mass and nominal isotope composition. Mass evidence alone cannot identify which isomer is present.",
        ["ISOBARIC_AMBIGUITY"],
        "No real synthetic starting-material relationship is asserted. No MS2 discrimination is evaluated.",
    )

    cat = "OPTICAL_IFE_ARTIFACT"
    dye_rows = [x for x in dyes if x["emission_max_nm"] - x["absorption_max_nm"] > 20][:5]
    designs = [
        ("OP01", None, None, "INCONCLUSIVE"),
        ("OP02", 0.6, 0.02, "ARTIFACT_DETECTED"),
        ("OP03", 0.6, 0.6, "VERIFIED"),
        ("OP04", 1.5, 0.02, "INCONCLUSIVE"),
        ("OP05", 0.02, 0.02, "VERIFIED"),
    ]
    for src, (cid, aex, refex, expected) in zip(dye_rows, designs, strict=True):
        r = record(cid, src["chromophore_smiles"], src["solvent_name"], "fluorescence_probe")
        original = f"Absorption maximum {src['absorption_max_nm']} nm; emission {src['emission_max_nm']} nm; QY {src['quantum_yield']}; log10 epsilon {src['log10_epsilon_at_absorption_max']} (M^-1 cm^-1)."
        assumptions = {}
        if aex is None:
            title = "Source lacks excitation setting and absorbances"
            change = "None. Optical input cannot represent these missing required values, so the source data remain in the snapshot."
            why = "Reported maxima and QY alone do not quantify cuvette inner-filter effects."
            flags = ["NO_SPECTRAL_EVIDENCE"]
        else:
            aem, refem = 0.1 * aex, 0.1 * refex
            factor = 10 ** ((aex + aem - refex - refem) / 2)
            apparent = 1 / factor
            concentration = aex / src["epsilon_at_absorption_max_M_inv_cm_inv"]
            assumptions = {
                "path_length_cm": 1.0,
                "excitation_equals_absorption_max": True,
                "epsilon_em_over_epsilon_ex": 0.1,
                "sample_concentration_M": concentration,
                "reference_concentration_M": refex / src["epsilon_at_absorption_max_M_inv_cm_inv"],
                "intrinsic_normalized_response_ratio": 1.0,
                "apparent_response_ratio": apparent,
                "sample_Aex": aex,
                "sample_Aem": aem,
                "reference_Aex": refex,
                "reference_Aem": refem,
                "response_normalization": "Per concentration/absorbed-photon response; not raw intensity across unequal concentrations.",
            }
            r["assay"]["concentration_um"] = concentration * 1e6
            r["evidence"]["optical"] = {
                "excitation_nm": src["absorption_max_nm"],
                "emission_nm": src["emission_max_nm"],
                "absorbance_at_ex": aex,
                "absorbance_at_em": aem,
                "baseline_absorbance_at_ex": refex,
                "baseline_absorbance_at_em": refem,
                "observed_quantum_yield": src["quantum_yield"],
                "apparent_fold_change": apparent,
                "geometry": "right_angle_center",
            }
            title = {
                "OP02": "Moderate OD simulates apparent quenching",
                "OP03": "Matched absorbance control",
                "OP04": "High OD exceeds correction policy",
                "OP05": "Dilute simulated control",
            }[cid]
            change = f"Assume 1 cm right-angle cuvette, excitation at absorption maximum and epsilon_em=0.1*epsilon_ex. Generate sample Aex={aex}, reference Aex={refex}, normalized apparent ratio={apparent:.6f}."
            why = {
                "OP02": "The stipulated unchanged intrinsic response is restored by the differential IFE correction; apparent quenching is an optical artifact.",
                "OP03": "Sample and reference have equal attenuation, so their normalized ratio remains 1. High absorbance alone is not a response artifact.",
                "OP04": "Aex=1.5 exceeds the engine's declared A<=1 correction domain. Abstention is required under that policy, not a universal physical cutoff.",
                "OP05": "Dilute matched sample/reference with a physically possible QY pass the implemented consistency checks.",
            }[cid]
            flags = {"OP02": ["IFE_RESPONSE_ARTIFACT"], "OP03": [], "OP04": ["IFE_OUTSIDE_VALIDITY"], "OP05": []}[cid]
        add(
            r,
            cat,
            title,
            [source(src, "optical")],
            "source_only" if aex is None else "source_based_simulation",
            original,
            change,
            expected,
            why,
            flags,
            "Absorption maximum is not a measured excitation setting. OD, path length, epsilon_em, concentrations and response ratio are missing in the source; any supplied values are simulated. QY is held fixed as a reference property, not a high-concentration measurement.",
            assumptions,
        )
        OPTICAL.append({"case_id": cid, "source": src, "assumptions": assumptions})

    cat = "NMR_SOLVENT_INTEGRAL_TRAP"
    for cid, idx, mode in [
        ("NM01", 0, "control"),
        ("NM02", 3, "ethyl_acetate"),
        ("NM03", 5, "no_scale"),
        ("NM04", 7, "dcm"),
    ]:
        src = ords[idx]
        ev = parse_nmr(src)
        r = record(cid, ord_smiles(src), "CDCl3")
        orig_total = sum(ev["proton_integrals"])
        if mode == "ethyl_acetate":
            ev["proton_shifts_ppm"] += [2.05, 4.12, 1.26]
            ev["proton_integrals"] += [3.0, 2.0, 3.0]
            ev["multiplicities"] += ["s", "q", "t"]
        elif mode == "dcm":
            ev["proton_shifts_ppm"] += [5.30]
            ev["proton_integrals"] += [2.0]
            ev["multiplicities"] += ["s"]
        elif mode == "no_scale":
            ev["integral_scale"] = None
        r["evidence"]["nmr"] = ev
        status = "VERIFIED" if mode == "control" else "INCONCLUSIVE" if mode == "no_scale" else "ARTIFACT_DETECTED"
        change = {
            "control": "No peak/count modification. Source nH assignments are represented in proton equivalents (scale=1).",
            "ethyl_acetate": "Append synthetic EtOAc signals: 2.05 s 3H, 4.12 q 2H, 1.26 t 3H; wrongly count them as product.",
            "dcm": "Append synthetic DCM signal: 5.30 s 2H; wrongly count it as product.",
            "no_scale": "Remove the reported-proton-equivalent scale; retain source peak positions and numerical integrals.",
        }[mode]
        flags = (
            []
            if mode == "control"
            else ["UNCALIBRATED_INTEGRALS"]
            if mode == "no_scale"
            else ["SOLVENT_INTEGRAL_INFLATION"]
        )
        add(
            r,
            cat,
            {
                "control": "Source NMR proton-count control",
                "ethyl_acetate": "EtOAc falsely included in product integral",
                "dcm": "DCM falsely included in product integral",
                "no_scale": "Unknown integral calibration",
            }[mode],
            [source(src, "ord")],
            "source_only" if mode == "control" else "source_plus_perturbation",
            nmr_text(src),
            change,
            status,
            "Source peak-count consistency is testable without assigning molecular identity."
            if mode == "control"
            else "Without a supplied integral scale the proton total is not calibrated."
            if mode == "no_scale"
            else "The deliberately added contaminant increases the complete product integral above its molecular proton budget.",
            flags,
            "Reported nH values are author assignments, not re-integrated raw spectra or an independently calibrated qNMR measurement. Range shifts use midpoints; complex multiplicities map to unknown. 13C evidence is not tested.",
            {"original_reported_H_total": orig_total, "impurity_reference": "https://doi.org/10.1021/jo971176v"},
        )
    src = next(x for x in nmrs if x["solvent"] == "CD3OD")
    r = record("NM05", src["structure"]["canonical_smiles"], "CD3OD")
    r["evidence"]["nmr"] = {
        "solvent": "CD3OD",
        "proton_shifts_ppm": [],
        "proton_integrals": [],
        "integration_complete": False,
    }
    add(
        r,
        cat,
        "Experimental assignments without measured integrals",
        [source(src, "nmr")],
        "source_only",
        "Experimental 1H/13C atom-assigned shifts in CD3OD; measured_integral is null for each assignment.",
        "None. Preserve assigned shifts in the source snapshot; submit an explicitly incomplete integration list because paired arrays cannot encode null integrals.",
        "INCONCLUSIVE",
        "Atom assignments are not measured integrated intensities. A complete calibrated proton sum is unavailable.",
        ["INCOMPLETE_INTEGRATION"],
        "This tests missing-data handling, not CD3OD exchange or chemical-shift assignment accuracy.",
    )

    cat = "STEREO_TAUTOMER_TRAP"
    for number in range(1, 6):
        src = stereo[(number - 1) * 2]
        mate = stereo[(number - 1) * 2 + 1]
        assert src["pair_id"] == mate["pair_id"]
        if number == 5:
            enumerator = rdMolStandardize.TautomerEnumerator()
            src = next(
                x
                for x in stereo[10:]
                if len(enumerator.Enumerate(Chem.MolFromSmiles(x["measurement"]["canonical_smiles"]))) > 1
            )
        original = src["measurement"]["canonical_smiles"]
        r = record(f"ST{number:02}", original, assay="binding_assay")
        r["compound"]["reference_smiles"] = original
        refs = [source(src, "structure")]
        status = "INCONCLUSIVE"
        flags = []
        change = "None. Retain the source stereochemical structure."
        title = "Matching structure without analytical evidence"
        why = "Matching structural metadata alone cannot verify an experimental sample or biological activity."
        if number == 2:
            r["compound"]["reference_smiles"] = mate["measurement"]["canonical_smiles"]
            refs.append(source(mate, "structure"))
            status = "ARTIFACT_DETECTED"
            flags = ["STEREO_MISMATCH"]
            title = "Swap reference to the paired stereoisomer"
            change = "Replace reference_smiles with the real paired stereoisomer from the same assay."
            why = "The proposed and asserted reference structures have different stereochemistry."
        elif number == 3:
            r["compound"]["proposed_smiles"] = src["structure_identity"]["connectivity_smiles"]
            r["compound"]["expected_stereocenters"] = src["structure_identity"]["specified_tetrahedral_centers"]
            status = "ARTIFACT_DETECTED"
            flags = ["MISSING_STEREOCHEMISTRY"]
            title = "Remove a required stereochemical assignment"
            change = "Strip the stereochemical annotations; retain the source expectation for assigned centers."
            why = "The modified structure omits stereochemistry required by the source identity assertion."
        elif number == 4:
            r["compound"]["inchikey"] = mate["structure_identity"]["rdkit_inchikey"]
            refs.append(source(mate, "structure"))
            status = "ARTIFACT_DETECTED"
            flags = ["INCHI_LAYER_MISMATCH"]
            title = "Swap the full InChIKey stereo layer"
            change = "Attach the paired stereoisomer's derived full InChIKey to the source SMILES."
            why = "Matching connectivity blocks do not make different full stereochemical identifiers interchangeable."
        elif number == 5:
            enumerator = rdMolStandardize.TautomerEnumerator()
            canonical = Chem.MolToSmiles(Chem.MolFromSmiles(original))
            alternatives = [
                Chem.MolToSmiles(x)
                for x in enumerator.Enumerate(Chem.MolFromSmiles(original))
                if Chem.MolToSmiles(x) != canonical
            ]
            assert alternatives
            r["compound"]["reference_smiles"] = alternatives[0]
            title = "Generated tautomer retains parent connectivity"
            change = "Replace reference with a deterministic RDKit-enumerated tautomer, not an experimentally observed tautomer."
            why = "A tautomer-equivalent representation is not sufficient evidence of an identity error; no spectrum supports verification."
            flags = ["TAUTOMER_EQUIVALENT"]
        add(
            r,
            cat,
            title,
            refs,
            "source_only" if number == 1 else "source_plus_perturbation",
            f"ChEMBL {src['measurement']['molecule_chembl_id']}; {src['measurement']['standard_type']}={src['measurement']['standard_value']} {src['measurement']['standard_units']}; assay {src['measurement']['assay_chembl_id']}; SMILES {original}.",
            change,
            status,
            why,
            flags,
            "This tests structure identity assertions. Activity differences do not establish a statistically significant activity flip. Aggregation and experimental tautomer populations are not evaluated.",
        )

    cat = "VERIFIED_TRUE_AND_NEGATIVE"
    for cid, idx in [("TR01", 2), ("TR02", 6), ("TR03", 8), ("TR04", 60), ("TR05", 61)]:
        src = ords[idx]
        r = record(cid, ord_smiles(src), "CDCl3" if idx < 60 else "see source procedure")
        r["ord_reaction_id"] = src["record_id"]
        status = "VERIFIED"
        flags = []
        modification = "None. Extract source-reported NMR peak positions and nH assignments."
        title = "Source-characterized product: NMR consistency"
        reason = "Complete source-reported proton equivalents agree with the molecular proton budget. VERIFIED is limited to the supplied consistency checks."
        if idx < 60:
            r["evidence"]["nmr"] = parse_nmr(src)
            original = nmr_text(src) + f" Recorded yield: {src['recorded_yield_percentages']}%."
            if cid == "TR03":
                mz = src["identity_analyses"]["MS of product_1"]["data"]["found"]["float_value"]
                r["evidence"]["lcms"] = {"observed_mz": mz, "ionization_mode": "positive"}
                status = "INCONCLUSIVE"
                title = "Real HRMS value with unspecified adduct"
                modification = "None to the measured m/z. Leave absent adduct and isotope envelope missing. Positive polarity is an explicit adapter assumption."
                reason = "The source excerpt does not specify the ion/adduct and has no isotope envelope. The observed value cannot be rejected by silently treating it as [M+H]+."
                original += f" HRMS found {mz} m/z; adduct absent from retained source excerpt."
        else:
            title = "Source reports zero yield; controls unresolved"
            r["negative_result"] = {
                "experiment_id": src["record_id"],
                "outcome": "failed_synthesis",
                "observation": src["procedure_details"],
                "protocol_reference": src["source_url"],
                "evidence_references": [src["source_url"]],
                "detection_limit_description": "Not reported in the source record.",
                "controls_passed": False,
                "measured_yield_percent": 0.0,
            }
            status = "INCONCLUSIVE"
            flags = ["NEGATIVE_CONTROLS_UNRESOLVED"]
            original = src["procedure_details"] + " Recorded yield: 0%."
            modification = "None. controls_passed=False conservatively represents not established, not a documented failed control. Detection-limit text explicitly states missingness."
            reason = "A reported failure and zero yield are real observations, but absent detection limits and control validation do not support a qualified NEGATIVE_RESULT."
        add(
            r,
            cat,
            title,
            [source(src, "ord")],
            "source_only",
            original,
            modification,
            status,
            reason,
            flags,
            "Yield is provenance, not independently verified by this engine. NMR uses author-reported nH assignments. No qualified negative control is available in this source subset.",
        )

    assert len(CASES) == 25 and set(Counter(c["category"] for c in CASES).values()) == {5}
    assert len({c["case_id"] for c in CASES}) == 25
    write(HERE / "expected_cases.json", CASES)
    write(HERE / "optical_calculations.json", OPTICAL)
    files = [
        HERE / "expected_cases.json",
        HERE / "optical_calculations.json",
        *sorted((HERE / "cases").glob("*.json")),
        *sorted((HERE / "sources").glob("*.json")),
    ]
    write(
        HERE / "expectations.lock.json",
        {
            "created_utc": datetime.now(UTC).isoformat(),
            "protocol": "Curator expectations written before the first pilot engine run; no independent/blinded review claimed.",
            "expected_counts": dict(Counter(c["expected_status"] for c in CASES)),
            "files": {p.relative_to(HERE).as_posix(): sha(p) for p in files},
        },
    )
    print(
        json.dumps(
            {
                "cases": len(CASES),
                "categories": dict(Counter(c["category"] for c in CASES)),
                "expected_counts": dict(Counter(c["expected_status"] for c in CASES)),
                "lock_sha256": sha(HERE / "expectations.lock.json"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
