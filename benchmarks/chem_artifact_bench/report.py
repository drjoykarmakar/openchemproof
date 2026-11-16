"""Explicit multiclass and artifact-detection metrics with full denominators."""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from benchmarks.chem_artifact_bench.models import BenchmarkCase


def calculate_metrics(cases: Sequence[BenchmarkCase], predictions: Sequence[dict[str, Any]]) -> dict[str, Any]:
    by_id = {row["case_id"]: row for row in predictions}
    if len(by_id) != len(predictions):
        raise ValueError("Duplicate predictions")
    if set(by_id) != {case.case_id for case in cases}:
        raise ValueError("Predictions must cover exactly the evaluated cases")
    tp = fp = fn = tn = accepted_artifacts = accepted_negatives = correct = errors = 0
    confusion: dict[str, dict[str, int]] = {}
    category_counts: dict[str, dict[str, int]] = {}
    for case in cases:
        prediction = by_id[case.case_id]
        predicted = prediction["status"]
        actual = case.expected_status
        error = bool(prediction.get("error"))
        errors += error
        correct += predicted == actual and not error
        confusion.setdefault(actual, {})[predicted] = confusion.setdefault(actual, {}).get(predicted, 0) + 1
        group = category_counts.setdefault(case.category, {"correct": 0, "total": 0})
        group["total"] += 1
        group["correct"] += predicted == actual and not error
        positive = actual == "ARTIFACT_DETECTED"
        called = predicted == "ARTIFACT_DETECTED" and not error
        tp += positive and called
        fp += not positive and called
        fn += positive and not called
        tn += not positive and not called
        accepted_artifacts += positive and predicted == "VERIFIED" and not error
        accepted_negatives += actual == "NEGATIVE_RESULT" and predicted == "VERIFIED" and not error
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    negatives = sum(case.expected_status == "NEGATIVE_RESULT" for case in cases)
    return {
        "count": len(cases),
        "correct": correct,
        "accuracy": correct / len(cases) if cases else 0.0,
        "artifact_precision": precision,
        "artifact_recall": recall,
        "artifact_f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
        "false_positive_acceptance_rate": accepted_artifacts / (tp + fn) if tp + fn else None,
        "negative_accepted_as_hit_rate": accepted_negatives / negatives if negatives else None,
        "artifact_confusion": {"tp": tp, "fp": fp, "fn": fn, "tn": tn},
        "errors": errors,
        "confusion_matrix": confusion,
        "per_category": category_counts,
        "definitions": {
            "positive_class": "ARTIFACT_DETECTED",
            "false_positive_acceptance_rate": "artifact cases predicted VERIFIED / all true artifact cases",
            "undefined_precision_recall": "0.0 when denominator is zero",
            "errors": "Recorded as INCONCLUSIVE; counted incorrect, never removed from denominator",
        },
    }


def write_reports(results: dict[str, dict[str, Any]], destination: Path, metadata: dict[str, Any]) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    payload = {"metadata": metadata, "results": results}
    (destination / "leaderboard.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    lines = [
        "# chem-artifact-bench",
        "",
        str(metadata["interpretation"]),
        "",
        "| Mode | Cases | Accuracy | Artifact precision | Recall | F1 | Artifact acceptance | Errors |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for mode, row in results.items():
        acceptance = row["false_positive_acceptance_rate"]
        value = "n/a" if acceptance is None else f"{acceptance:.4f}"
        lines.append(
            f"| {mode} | {row['count']} | {row['accuracy']:.4f} | "
            f"{row['artifact_precision']:.4f} | {row['artifact_recall']:.4f} | "
            f"{row['artifact_f1']:.4f} | {value} | {row['errors']} |"
        )
    lines.extend(
        [
            "",
            "The positive class for precision/recall/F1 is ARTIFACT_DETECTED.",
            "Artifact acceptance is the fraction of true artifact cases accepted as VERIFIED.",
            "Negative-result acceptance as a hit is reported separately in JSON.",
            "",
        ]
    )
    (destination / "leaderboard.md").write_text("\n".join(lines))
