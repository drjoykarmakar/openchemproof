import json
from collections import Counter
from pathlib import Path

from rdkit import Chem

from benchmarks.chem_artifact_bench.source_index import FILES, build_index

ROOT = Path(__file__).resolve().parents[1] / "external_data"


def test_500_source_rows_preserve_evidence_boundaries():
    records = build_index(ROOT)
    assert len(records) == 500
    assert set(Counter(row.category for row in records).values()) == {100}
    assert all(row.expected_verdict is None and not row.synthetic_perturbation for row in records)
    assert all(row.missing_for_requested_verification for row in records)
    assert all(Chem.MolFromSmiles(row.proposed_smiles) is not None for row in records)


def test_missing_measurements_not_imputed():
    data = {name: [json.loads(line) for line in (ROOT / path).read_text().splitlines()] for name, path in FILES.items()}
    assert all(row["isotope_envelope"]["M+2"] is None for row in data["lcms"])
    assert all(row["excitation_wavelength_nm"] is None for row in data["optical"])
    assert all(
        peak["measured_integral"] is None
        for row in data["nmr"]
        for spec in row["spectra"]
        for peak in spec["assigned_shifts"]
    )
    assert len({row["pair_id"] for row in data["structure"]}) == 50
    assert sum(row["cohort"] == "source_characterized_success" for row in data["ord"]) == 60
    assert not any(row["strict_negative_eligible"] for row in data["ord"])
