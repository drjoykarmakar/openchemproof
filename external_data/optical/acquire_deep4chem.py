#!/usr/bin/env python3
"""Reproduce the 100-record experimental optical subset from official Figshare v2.

Run with the project's environment (RDKit required). Cached raw inputs are verified;
pass --refresh to fetch them again. No spectra or optical artifacts are synthesized.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import urllib.request
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit import __version__ as rdkit_version

BASE = Path(__file__).resolve().parent
API_URL = "https://api.figshare.com/v2/articles/12045567/versions/2"
DOWNLOAD_URL = "https://ndownloader.figshare.com/files/23637518"
DATASET_DOI = "10.6084/m9.figshare.12045567.v2"
EXPECTED_MD5 = "6619ca296ebbf14e8c80dda084b43ac5"
SOLVENT_NAMES = {
    "O": "water",
    "CO": "methanol",
    "CCO": "ethanol",
    "CC#N": "acetonitrile",
    "ClCCl": "dichloromethane",
    "ClC(Cl)Cl": "chloroform",
    "Cc1ccccc1": "toluene",
    "C1CCOC1": "tetrahydrofuran",
    "CS(C)=O": "dimethyl sulfoxide",
    "CN(C)C=O": "N,N-dimethylformamide",
    "C1CCCCC1": "cyclohexane",
    "CCOC(C)=O": "ethyl acetate",
    "C1COCCO1": "1,4-dioxane",
    "CCCCCC": "n-hexane",
    "CC(C)=O": "acetone",
    "c1ccccc1": "benzene",
    "CCOCC": "diethyl ether",
    "CC(C)O": "2-propanol",
    "OCC(F)(F)F": "2,2,2-trifluoroethanol",
    "CCCCCCCCO": "1-octanol",
    "ClC(Cl)(Cl)Cl": "carbon tetrachloride",
    "CN1CCCC1=O": "N-methyl-2-pyrrolidone",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def number(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def fetch_cached(path: Path, url: str, refresh: bool) -> bytes:
    if refresh or not path.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "OpenChemProof/0.1"})
        with urllib.request.urlopen(request, timeout=90) as response:
            data = response.read()
        path.write_bytes(data)
    return path.read_bytes()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    raw_dir = BASE / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    metadata_bytes = fetch_cached(raw_dir / "figshare_v2_metadata.json", API_URL, args.refresh)
    metadata = json.loads(metadata_bytes)
    if metadata["license"]["name"] != "CC BY 4.0":
        raise ValueError("Unexpected source license; inspect metadata before redistribution")
    blob = fetch_cached(raw_dir / "deep4chem_v2.csv", DOWNLOAD_URL, args.refresh)
    if hashlib.md5(blob).hexdigest() != EXPECTED_MD5:  # noqa: S324 (source integrity, not security)
        raise ValueError("Figshare v2 CSV checksum differs from pinned publication")
    rows = list(csv.DictReader(io.StringIO(blob.decode("utf-8-sig"))))
    counters: Counter[str] = Counter()
    eligible = []
    RDLogger.DisableLog("rdApp.error")
    for csv_line, row in enumerate(rows, 2):
        absorption = number(row["Absorption max (nm)"])
        emission = number(row["Emission max (nm)"])
        qy = number(row["Quantum yield"])
        log_epsilon = number(row["log(e/mol-1 dm3 cm-1)"])
        if any(x is None for x in (absorption, emission, qy, log_epsilon)):
            counters["missing_required_optical_property"] += 1
            continue
        if not (absorption > 0 and emission > 0 and 0 <= qy <= 1 and -10 < log_epsilon < 10):
            counters["invalid_numeric_range"] += 1
            continue
        doi_match = re.search(r"10\.\d{4,9}/[^\s]+", row["Reference"])
        if not doi_match:
            counters["missing_parseable_reference_doi"] += 1
            continue
        molecule = Chem.MolFromSmiles(row["Chromophore"])
        solvent = Chem.MolFromSmiles(row["Solvent"])
        if molecule is None or solvent is None:
            counters["unparseable_structure_or_solvent"] += 1
            continue
        if "." in row["Solvent"]:
            counters["mixture_solvent_excluded"] += 1
            continue
        eligible.append(
            (
                row,
                csv_line,
                Chem.MolToSmiles(molecule),
                Chem.MolToSmiles(solvent),
                doi_match.group(0),
                absorption,
                emission,
                qy,
                log_epsilon,
            )
        )

    # Stable hash order avoids selecting a contiguous block from one paper. One
    # canonical chromophore and one original literature reference per subset row.
    eligible.sort(key=lambda item: sha256(f"openchemproof-optical-v1:{item[0]['Tag']}".encode()))
    selected = []
    seen_molecules: set[str] = set()
    seen_references: set[str] = set()
    for item in eligible:
        row, csv_line, smiles, solvent_smiles, doi, absorption, emission, qy, log_epsilon = item
        if smiles in seen_molecules or doi.lower() in seen_references:
            continue
        record = {
            "record_id": f"deep4chem-v2-tag-{row['Tag']}",
            "record_kind": "literature_compiled_experimental",
            "is_simulated": False,
            "chromophore_smiles": smiles,
            "chromophore_smiles_original": row["Chromophore"],
            "solvent_smiles": solvent_smiles,
            "solvent_smiles_original": row["Solvent"],
            "solvent_name": SOLVENT_NAMES.get(row["Solvent"]),
            "absorption_max_nm": absorption,
            "emission_max_nm": emission,
            "quantum_yield": qy,
            "log10_epsilon_at_absorption_max": log_epsilon,
            "epsilon_at_absorption_max_M_inv_cm_inv": 10**log_epsilon,
            "epsilon_derivation": "10 ** source_log10_epsilon; unit conversion only",
            "fluorescence_lifetime_ns": number(row["Lifetime (ns)"]),
            "absorption_fwhm_nm": number(row["abs FWHM (nm)"]),
            "emission_fwhm_nm": number(row["emi FWHM (nm)"]),
            "absorption_fwhm_cm_inv": number(row["abs FWHM (cm-1)"]),
            "emission_fwhm_cm_inv": number(row["emi FWHM (cm-1)"]),
            "molecular_weight_g_mol": number(row["Molecular weight (g mol-1)"]),
            "excitation_wavelength_nm": None,
            "epsilon_at_excitation_M_inv_cm_inv": None,
            "epsilon_at_emission_M_inv_cm_inv": None,
            "absorbance_at_excitation": None,
            "absorbance_at_emission": None,
            "concentration_M": None,
            "path_length_cm": None,
            "temperature_K": None,
            "pH": None,
            "absorption_spectrum": None,
            "emission_spectrum": None,
            "source": {
                "dataset": "Deep4Chem / DB for chromophore",
                "dataset_doi": DATASET_DOI,
                "dataset_version": 2,
                "source_tag": row["Tag"],
                "csv_line": csv_line,
                "download_url": DOWNLOAD_URL,
                "source_file_sha256": sha256(blob),
                "reference_original": row["Reference"],
                "reference_doi": doi,
                "reference_url": f"https://doi.org/{doi}",
                "license": "CC-BY-4.0",
            },
            "provenance_limitations": [
                "Compiled experimental values; original article not independently re-extracted.",
                "Absorption maximum is not a recorded excitation wavelength.",
                "Epsilon is reported at absorption maximum, not excitation or emission wavelength.",
                "No measured OD, concentration, path length or spectral trace in this source CSV.",
            ],
        }
        selected.append(record)
        seen_molecules.add(smiles)
        seen_references.add(doi.lower())
        if len(selected) == 100:
            break
    if len(selected) != 100:
        raise ValueError(f"Expected 100 records, found {len(selected)}")
    selected.sort(key=lambda item: int(item["source"]["source_tag"]))
    output = "".join(
        json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n" for record in selected
    )
    (BASE / "records.jsonl").write_text(output, encoding="utf-8")
    selected_tags = {record["source"]["source_tag"] for record in selected}
    with (BASE / "source_rows.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(row for row in rows if row["Tag"] in selected_tags)
    manifest = {
        "schema_version": "openchemproof.external.optical.v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "source_dataset_doi": DATASET_DOI,
        "source_publication_doi": "10.1038/s41597-020-00634-8",
        "source_metadata_url": API_URL,
        "source_download_url": DOWNLOAD_URL,
        "source_license": metadata["license"],
        "source_attribution": metadata["citation"],
        "source_csv_bytes": len(blob),
        "source_csv_md5": EXPECTED_MD5,
        "source_csv_sha256": sha256(blob),
        "metadata_sha256": sha256(metadata_bytes),
        "records_sha256": sha256(output.encode()),
        "source_record_count": len(rows),
        "eligible_record_count": len(eligible),
        "excluded_by_first_reason": dict(counters),
        "selected_record_count": len(selected),
        "unique_canonical_chromophores": len(seen_molecules),
        "unique_original_literature_dois": len(seen_references),
        "unique_canonical_solvents": len({r["solvent_smiles"] for r in selected}),
        "selection": (
            "SHA256(openchemproof-optical-v1:<Tag>) order; first 100 with distinct RDKit "
            "canonical chromophore and distinct reference DOI. Require positive abs/em max, "
            "QY in [0,1], finite log10 epsilon in (-10,10), parseable chromophore and "
            "single-component solvent SMILES, DOI. Sort output by numeric Tag."
        ),
        "selection_bias": (
            "Complete-property, unique-dye, unique-reference convenience subset; "
            "not random or representative of the full dataset."
        ),
        "selected_tags": [record["source"]["source_tag"] for record in selected],
        "normalization": (
            "Preserve source SMILES; add RDKit canonical SMILES; numeric strings to JSON "
            "numbers; NaN to null; epsilon=10**log10_epsilon; parse reference DOI. No imputation."
        ),
        "rdkit_version": rdkit_version,
        "experimental_records": 100,
        "simulated_records": 0,
        "independently_reextracted_original_articles": 0,
        "missingness": {
            field: sum(record[field] is None for record in selected)
            for field in [
                "absorption_max_nm",
                "emission_max_nm",
                "quantum_yield",
                "epsilon_at_absorption_max_M_inv_cm_inv",
                "fluorescence_lifetime_ns",
                "excitation_wavelength_nm",
                "absorbance_at_excitation",
                "absorbance_at_emission",
                "absorption_spectrum",
                "concentration_M",
                "path_length_cm",
            ]
        },
    }
    (BASE / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: manifest[key]
                for key in [
                    "selected_record_count",
                    "unique_canonical_chromophores",
                    "unique_original_literature_dois",
                    "unique_canonical_solvents",
                    "records_sha256",
                ]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
