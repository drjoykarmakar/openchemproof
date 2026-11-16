#!/usr/bin/env python3
"""Reacquire pinned, real MassBank LC-MS/MS records; never synthesize observations.

Run from the repository with its RDKit environment:
  python external_data/lcms/acquire_massbank.py --count 100
The raw records are licensed separately from OpenChemProof code (see manifest).
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import urllib.request

from rdkit import Chem, rdBase

COMMIT = "72980b3af504cfd1e89ac54d7698d3da45d85a58"
RELEASE = "2026.11"
TREE_SHA = "1cd69cbf8796c73124237e928e1fc432a953f8ed"
REPO = "https://github.com/MassBank/MassBank-data"
ROOT = Path(__file__).resolve().parent


def get(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "OpenChemProof/0.1 reproducible research acquisition"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=40) as response:
                return response.read()
        except Exception:
            if attempt == 2:
                raise
    raise AssertionError("unreachable")


def raw_fields(text: str) -> dict[str, list[str]]:
    fields: dict[str, list[str]] = defaultdict(list)
    for line in text.splitlines():
        if line and not line.startswith(" ") and ": " in line:
            key, value = line.split(": ", 1)
            fields[key].append(value)
    return dict(fields)


def normalize(data: bytes, filename: str) -> dict:
    text = data.decode("utf-8")
    fields = raw_fields(text)
    def one(key):
        return fields.get(key, [None])[0]
    def sub(key, name):
        return next((v[len(name)+1:] for v in fields.get(key, []) if v.startswith(name + " ")), None)
    assert one("AC$INSTRUMENT_TYPE").startswith("LC-"), "not LC"
    assert sub("AC$MASS_SPECTROMETRY", "MS_TYPE") == "MS2", "not MS2"
    assert one("LICENSE") == "CC BY", "unexpected license"
    assert "CONFIDENCE standard compound" in fields.get("COMMENT", []), "not standard compound"
    original_smiles = one("CH$SMILES")
    mol = Chem.MolFromSmiles(original_smiles)
    assert mol is not None, "invalid source SMILES"
    smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    precursor = sub("MS$FOCUSED_ION", "PRECURSOR_M/Z")
    adduct = sub("MS$FOCUSED_ION", "PRECURSOR_TYPE")
    mode = sub("AC$MASS_SPECTROMETRY", "ION_MODE")
    assert precursor and adduct and mode in {"POSITIVE", "NEGATIVE"}, "missing precursor fields"
    assert math.isfinite(float(precursor)) and float(precursor) > 0
    peaks = []
    reading = False
    for line in text.splitlines():
        if line.startswith("PK$PEAK:"):
            reading = True
        elif reading and line.startswith("  "):
            values = line.split()
            assert len(values) == 3
            mz, intensity, relative = map(float, values)
            assert all(math.isfinite(v) for v in (mz, intensity, relative))
            assert mz > 0 and intensity >= 0 and relative >= 0
            peaks.append({"mz": mz, "intensity": intensity, "relative_intensity": relative})
        elif reading:
            reading = False
    assert len(peaks) == int(one("PK$NUM_PEAK")) and peaks, "peak count mismatch"
    accession = one("ACCESSION")
    path = "Eawag/" + filename
    return {
        "schema_version": "openchemproof.real_lcms.v1",
        "record_id": accession,
        "data_origin": "public_experimental_library_record",
        "source": {"database": "MassBank", "release": RELEASE, "commit": COMMIT,
                   "path": path, "url": f"{REPO}/blob/{COMMIT}/{path}",
                   "raw_url": f"https://raw.githubusercontent.com/MassBank/MassBank-data/{COMMIT}/{path}",
                   "sha256": hashlib.sha256(data).hexdigest(),
                   "raw_file": f"raw/{filename}", "authors": one("AUTHORS"),
                   "copyright": one("COPYRIGHT"), "license_literal": one("LICENSE"),
                   "license_version": None, "publication": fields.get("PUBLICATION", []),
                   "record_date": one("DATE")},
        "compound": {"name": one("CH$NAME"), "names": fields.get("CH$NAME", []),
                     "smiles_original": original_smiles, "smiles_canonical": smiles,
                     "canonicalization": f"RDKit {rdBase.rdkitVersion}; canonical isomeric SMILES; no tautomer/salt/charge standardization",
                     "formula_source_reported": one("CH$FORMULA"),
                     "exact_mass_source_reported": float(one("CH$EXACT_MASS")),
                     "structure_status": "source_asserted_standard_compound; not independently reverified"},
        "acquisition": {"instrument": one("AC$INSTRUMENT"), "instrument_type": one("AC$INSTRUMENT_TYPE"),
                        "ms_level": 2, "ion_mode": mode.lower(), "adduct": adduct,
                        "precursor_mz_source_reported": float(precursor),
                        "precursor_mz_original_text": precursor,
                        "precursor_mz_status": "source_reported_metadata; MassBank_prefers_calculated_exact_mass; not_independent_measured_mass",
                        "collision_energy_original": sub("AC$MASS_SPECTROMETRY", "COLLISION_ENERGY"),
                        "retention_time_original": sub("AC$CHROMATOGRAPHY", "RETENTION_TIME")},
        "spectrum": {"kind": "experimental_processed_ms2_peak_list", "peaks": peaks,
                     "splash": one("PK$SPLASH"),
                     "processing_original": fields.get("MS$DATA_PROCESSING", []),
                     "raw_instrument_data_available_in_this_bundle": False},
        "isotope_envelope": {"status": "unavailable_in_ms2_record", "M": None, "M+1": None, "M+2": None,
                             "reason": "MS2 fragment peaks cannot establish an MS1 precursor isotope envelope."},
        "benchmark_labels": {"structure_match": "source_asserted", "independent_positive_verification": False,
                             "synthetic_perturbation": False, "expected_verifier_verdict": None},
        "original_fields": fields,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=100)
    parser.add_argument("--offline", action="store_true", help="Rebuild from cached source records without network.")
    args = parser.parse_args()
    assert args.count > 0
    (ROOT / "cache").mkdir(exist_ok=True)
    (ROOT / "raw").mkdir(exist_ok=True)
    tree_path = ROOT / "cache" / "eawag_tree.json"
    tree_url = f"https://api.github.com/repos/MassBank/MassBank-data/git/trees/{TREE_SHA}"
    if not tree_path.exists():
        assert not args.offline, "tree cache missing"
        tree_path.write_bytes(get(tree_url))
    tree = json.loads(tree_path.read_bytes())
    assert not tree.get("truncated", False)
    # One candidate from each compound-number/polarity block, then interleave
    # positive/negative acquisition blocks. Canonical SMILES de-duplicates later.
    groups = {}
    for entry in sorted(tree["tree"], key=lambda x: x["path"]):
        match = re.fullmatch(r"MSBNK-Eawag-EA(\d{4})(\d{2})\.txt", entry["path"])
        if match:
            groups.setdefault((match[1], int(match[2]) >= 50), entry["path"])
    candidates = sorted(groups.values(), key=lambda x: (int(re.search(r"EA(\d{4})", x)[1]), x))
    failures = []
    selected = []
    seen = set()
    def acquire(filename):
        cache = ROOT / "cache" / filename
        try:
            if not cache.exists():
                assert not args.offline, "source cache missing"
                cache.write_bytes(get(f"https://raw.githubusercontent.com/MassBank/MassBank-data/{COMMIT}/Eawag/{filename}"))
            data = cache.read_bytes()
            record = normalize(data, filename)
            return filename, data, record, None
        except Exception as exc:
            return filename, None, None, f"{type(exc).__name__}: {exc}"
    for offset in range(0, len(candidates), 24):
        batch = candidates[offset:offset+24]
        with ThreadPoolExecutor(max_workers=8) as pool:
            results = list(pool.map(acquire, batch))
        # Prefer alternating polarity when both modes exist for the next molecule.
        by_compound = defaultdict(list)
        for result in results:
            by_compound[re.search(r"EA(\d{4})", result[0])[1]].append(result)
        for compound_results in by_compound.values():
            desired = "negative" if len(selected) % 2 else "positive"
            compound_results.sort(key=lambda r: bool(r[2] and r[2]["acquisition"]["ion_mode"] != desired))
            for filename, data, record, error in compound_results:
                if error:
                    failures.append({"path": filename, "reason": error})
                    continue
                smiles = record["compound"]["smiles_canonical"]
                if smiles in seen:
                    continue
                seen.add(smiles)
                selected.append(record)
                (ROOT / "raw" / filename).write_bytes(data)
                if len(selected) == args.count:
                    break
            if len(selected) == args.count:
                break
        print(f"selected {len(selected)} unique structures; scanned {offset+len(batch)} candidate records", flush=True)
        if len(selected) == args.count:
            break
    assert len(selected) == args.count, f"only acquired {len(selected)} qualified unique records"
    output = ROOT / "massbank_100.jsonl"
    output.write_text("".join(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n" for row in selected), encoding="utf-8")
    manifest = {
        "dataset": "OpenChemProof v0.1 public experimental LC-MS/MS acquisition subset",
        "created_utc": datetime.now(timezone.utc).isoformat(), "requested_count": args.count,
        "record_count": len(selected), "unique_canonical_smiles": len(seen),
        "source_repository": REPO, "source_release": RELEASE, "source_commit": COMMIT,
        "source_release_doi": "https://doi.org/10.5281/zenodo.24082491",
        "source_tree_url": tree_url, "source_tree_sha256": hashlib.sha256(tree_path.read_bytes()).hexdigest(),
        "normalized_file": output.name, "normalized_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "rdkit_version": rdBase.rdkitVersion,
        "licenses": dict(Counter(r["source"]["license_literal"] for r in selected)),
        "license_note": "Each original record says CC BY without a version. Preserve attribution, record URLs and exact source license; do not relicense these records under the code license or infer a version.",
        "ion_modes": dict(Counter(r["acquisition"]["ion_mode"] for r in selected)),
        "adducts": dict(Counter(r["acquisition"]["adduct"] for r in selected)),
        "instrument_types": dict(Counter(r["acquisition"]["instrument_type"] for r in selected)),
        "experimental_ms2_peak_lists": len(selected),
        "measured_ms1_isotope_envelopes": 0, "synthetic_records": 0,
        "selection": "Pinned Eawag EA series; ascending compound IDs; first available record per polarity block; prefer alternating polarity; only source-standard compounds with LC instrument, valid SMILES, precursor fields and MS2 peaks; unique canonical SMILES.",
        "limitations": [
            "Convenience subset of one contributor/series; not a representative or blind evaluation.",
            "Records are processed experimental MS2 reference spectra, not original raw instrument files.",
            "All M/M+1/M+2 isotope envelope fields are unavailable; MS2 fragments are not substituted.",
            "Reported precursor m/z is metadata; MassBank format explicitly prefers calculated exact precursor mass. Do not use this field as independent measured mass accuracy or PASS truth.",
            "Reference identity is contributor asserted for standard compounds; no engine-derived truth labels are assigned.",
            "Canonical SMILES parsing does not validate chemistry, stereochemistry or source assignment.",
            "No synthetic perturbations and no expected verifier verdict labels are included."
        ],
        "rejected_candidates": failures,
    }
    (ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("record_count", "unique_canonical_smiles", "ion_modes", "adducts", "experimental_ms2_peak_lists", "measured_ms1_isotope_envelopes")}, indent=2))


if __name__ == "__main__":
    main()
