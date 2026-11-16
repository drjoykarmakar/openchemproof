from collections import Counter

from benchmarks.chem_artifact_bench.eval_runner import run_oracle, sample_cases
from benchmarks.chem_artifact_bench.generate_bench import generate_cases
from benchmarks.chem_artifact_bench.report import calculate_metrics


def test_determinism_balance_and_hidden_labels():
    cases = generate_cases()
    assert cases == generate_cases()
    assert cases != generate_cases(43)
    assert len({case.record.record_id for case in cases}) == 500
    assert set(Counter(case.category for case in cases).values()) == {100}
    assert Counter(case.expected_status for case in cases) == {
        "ARTIFACT_DETECTED": 400,
        "VERIFIED": 60,
        "NEGATIVE_RESULT": 40,
    }
    for case in cases:
        record = case.record.model_dump_json()
        assert case.category not in record
        assert "expected_status" not in record
        assert case.evidence_origin == "synthetic_unit_fixture"
    assert list(Counter(case.category for case in sample_cases(cases, 25, 42)).values()) == [5] * 5


def test_entire_oracle_without_relabeling(tmp_path):
    metrics = run_oracle(generate_cases(), tmp_path)
    assert metrics["accuracy"] == 1.0
    assert metrics["expected_flag_failures"] == []
    assert metrics["artifact_confusion"] == {"tp": 400, "tn": 100, "fp": 0, "fn": 0}


def test_all_verified_catches_metric_direction():
    cases = generate_cases()
    predictions = [{"case_id": case.case_id, "status": "VERIFIED"} for case in cases]
    metrics = calculate_metrics(cases, predictions)
    assert metrics["accuracy"] == 0.12
    assert metrics["artifact_recall"] == 0
    assert metrics["false_positive_acceptance_rate"] == 1
    assert metrics["negative_accepted_as_hit_rate"] == 1
