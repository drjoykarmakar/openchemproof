"""Build a uniform 500-record provenance index without inventing analytical evidence."""

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Literal

from benchmarks.chem_artifact_bench.models import Category
from openchemproof.schema.export import canonical_json, write_parquet
from openchemproof.schema.models import StrictModel


class SourceRecord(StrictModel):
    record_id: str
    category: Category
    source_name: str
    source_url: str
    source_license: str
    normalized_file: str
    line_number: int
    canonical_row_sha256: str
    proposed_smiles: str
    experimental_evidence: list[str]
    missing_for_requested_verification: list[str]
    eligibility: Literal[
        "partial_analytical_evidence", "identity_comparison_only", "source_characterized", "source_reported_negative"
    ]
    independently_verified: Literal[False] = False
    expected_verdict: None = None
    synthetic_perturbation: Literal[False] = False


FILES = {
    "lcms": "lcms/massbank_100.jsonl",
    "optical": "optical/records.jsonl",
    "nmr": "nmr/nmrshiftdb2_100.jsonl",
    "structure": "structure/records.jsonl",
    "ord": "ord/records.jsonl",
}


def normalize_index_row(source: str, row: dict[str, Any], line: int) -> SourceRecord:
    fields: dict[str, Any] = dict(
        record_id=row["record_id"],
        normalized_file=FILES[source],
        line_number=line,
        canonical_row_sha256=hashlib.sha256(canonical_json(row).encode()).hexdigest(),
    )
    if source == "lcms":
        fields.update(
            category="LCMS_ADDUCT_DIMER_TRAP",
            source_name="MassBank",
            source_url=row["source"]["url"],
            source_license=row["source"]["license_literal"],
            proposed_smiles=row["compound"]["smiles_canonical"],
            experimental_evidence=["processed MS2 fragment spectrum", "source asserted standard compound"],
            missing_for_requested_verification=[
                "independently measured MS1 precursor m/z",
                "MS1 isotope envelope",
                "sodium/dimer examples in selected subset",
            ],
            eligibility="partial_analytical_evidence",
        )
    elif source == "optical":
        fields.update(
            category="OPTICAL_IFE_ARTIFACT",
            source_name="Deep4Chem",
            source_url=row["source"]["reference_url"],
            source_license=row["source"]["license"],
            proposed_smiles=row["chromophore_smiles"],
            experimental_evidence=[
                "compiled absorption maximum",
                "compiled emission maximum",
                "quantum yield",
                "log epsilon at absorption maximum",
            ],
            missing_for_requested_verification=[
                "measured excitation wavelength",
                "sample and reference absorbance at excitation and emission",
                "cuvette geometry and concentration",
                "independent measured fluorescence response",
            ],
            eligibility="partial_analytical_evidence",
        )
    elif source == "nmr":
        fields.update(
            category="NMR_SOLVENT_INTEGRAL_TRAP",
            source_name="NMRShiftDB2",
            source_url=row["source"]["record_url"],
            source_license=row["source"]["license"],
            proposed_smiles=row["structure"]["canonical_smiles"],
            experimental_evidence=["assigned 1H shifts", "assigned 13C shifts", "solvent metadata"],
            missing_for_requested_verification=[
                "measured proton integrals",
                "independent integral calibration",
                "sample-level pairing confirmation",
            ],
            eligibility="partial_analytical_evidence",
        )
    elif source == "structure":
        fields.update(
            category="STEREO_TAUTOMER_TRAP",
            source_name="ChEMBL",
            source_url=row["source_url"],
            source_license=row["source_license"],
            proposed_smiles=row["structure_identity"]["canonical_isomeric_smiles"],
            experimental_evidence=["source stereo structure", "numeric activity in an assay-matched pair"],
            missing_for_requested_verification=[
                "uncertainty for statistical activity flip",
                "tautomer subset",
                "experimentally confirmed aggregation controls",
            ],
            eligibility="identity_comparison_only",
        )
    else:
        negative = row["cohort"] != "source_characterized_success"
        smiles = next(item["value"] for item in row["product_identifiers"] if item["type"] == "SMILES")
        fields.update(
            category="VERIFIED_TRUE_AND_NEGATIVE",
            source_name="ORD",
            source_url=row["source_url"],
            source_license=row["license"],
            proposed_smiles=smiles,
            experimental_evidence=["source zero yield", "failure narrative"]
            if negative
            else ["source yield >=50%", "product-linked NMR text", "observed HRMS m/z"],
            missing_for_requested_verification=["validated detection limit", "documented successful controls"]
            if negative
            else ["independently matched raw spectra", "chromatographic purity and LCMS evidence"],
            eligibility="source_reported_negative" if negative else "source_characterized",
        )
    return SourceRecord.model_validate(fields)


def build_index(root: Path) -> list[SourceRecord]:
    records: list[SourceRecord] = []
    for source, filename in FILES.items():
        rows = [json.loads(line) for line in (root / filename).read_text().splitlines() if line.strip()]
        if len(rows) != 100 or len({row["record_id"] for row in rows}) != 100:
            raise ValueError(f"Expected 100 unique {source} records")
        records.extend(normalize_index_row(source, row, i) for i, row in enumerate(rows, start=1))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("external_data"))
    args = parser.parse_args()
    records = build_index(args.root)
    (args.root / "source_index.jsonl").write_text("".join(canonical_json(r.model_dump()) + "\n" for r in records))
    write_parquet(records, args.root / "source_index.parquet")
    print(f"Indexed {len(records)} public-source records; no experimental ground-truth verdicts invented")


if __name__ == "__main__":
    main()
