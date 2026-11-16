"""Run the frozen pilot without network access or an API key.

From the repository root: uv run python validation/pilot25/run_pilot.py
The exit code reports execution/integrity errors, not scientific disagreement.
Use --fail-on-disagreement for CI regression checks.
"""

import argparse
import csv
import hashlib
import json
import platform
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from openchemproof.schema.models import VerificationRecord
from openchemproof.verifiers.pipeline import run_verification

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "results")
    parser.add_argument("--fail-on-disagreement", action="store_true")
    args = parser.parse_args()
    out = args.output_dir
    lock = json.loads((HERE / "expectations.lock.json").read_text())
    for filename, expected_hash in lock["files"].items():
        if sha(HERE / filename) != expected_hash:
            raise SystemExit(f"Frozen file changed: {filename}. Restore it or create a new pilot version.")
    cases = json.loads((HERE / "expected_cases.json").read_text())
    results = []
    for case in cases:
        row = {
            "case_id": case["case_id"],
            "category": case["category"],
            "expected_status": case["expected_status"],
            "evidence_kind": case["evidence_kind"],
        }
        try:
            record = VerificationRecord.model_validate_json((HERE / case["input_file"]).read_text())
            verdict = run_verification(record, audit_dir=out / "audits")
            data = verdict.model_dump(mode="json")
            write(out / "verdicts" / f"{case['case_id']}.json", data)
            flags = [f["code"] for f in data["flags"]]
            row.update(
                actual_status=verdict.status,
                flags=flags,
                status_agrees=verdict.status == case["expected_status"],
                missing_required_flags=sorted(set(case["required_flags"]) - set(flags)),
                audit_sha256=verdict.audit_sha256,
                corrected_metrics=data["corrected_metrics"],
                error=None,
            )
        except Exception as exc:
            row.update(
                actual_status="EXECUTION_ERROR",
                flags=[],
                status_agrees=False,
                missing_required_flags=case["required_flags"],
                error=f"{type(exc).__name__}: {exc}",
            )
        row["diagnostic_agrees"] = row["status_agrees"] and not row["missing_required_flags"]
        row["review_status"] = "Pending human review"
        results.append(row)
    root = HERE.parents[1]
    engine_files = sorted((root / "src" / "openchemproof").rglob("*.py")) + [root / "pyproject.toml", root / "uv.lock"]
    file_hashes = {p.relative_to(root).as_posix(): sha(p) for p in engine_files}
    categories = list(dict.fromkeys(c["category"] for c in cases))
    statuses = ["VERIFIED", "ARTIFACT_DETECTED", "INCONCLUSIVE", "NEGATIVE_RESULT", "EXECUTION_ERROR"]
    report = {
        "run_utc": datetime.now(UTC).isoformat(),
        "python": platform.python_version(),
        "expectations_lock_sha256": sha(HERE / "expectations.lock.json"),
        "lock_created_utc": lock["created_utc"],
        "engine_file_sha256": file_hashes,
        "engine_fingerprint_sha256": hashlib.sha256(json.dumps(file_hashes, sort_keys=True).encode()).hexdigest(),
        "total": len(results),
        "status_agreement_count": sum(r["status_agrees"] for r in results),
        "status_and_diagnostic_agreement_count": sum(r["diagnostic_agrees"] for r in results),
        "execution_error_count": sum(r["error"] is not None for r in results),
        "expected_counts": dict(Counter(r["expected_status"] for r in results)),
        "actual_counts": dict(Counter(r["actual_status"] for r in results)),
        "evidence_kind_counts": dict(Counter(r["evidence_kind"] for r in results)),
        "by_category": {
            c: {
                "total": sum(r["category"] == c for r in results),
                "status_agreement": sum(r["category"] == c and r["status_agrees"] for r in results),
                "diagnostic_agreement": sum(r["category"] == c and r["diagnostic_agrees"] for r in results),
            }
            for c in categories
        },
        "confusion_matrix_expected_rows_actual_columns": {
            e: {a: sum(r["expected_status"] == e and r["actual_status"] == a for r in results) for a in statuses}
            for e in statuses[:-1]
        },
        "independent_blinded_validation": False,
        "human_review_complete": False,
        "results": results,
    }
    write(out / "report.json", report)
    with (out / "case_results.csv").open("w", newline="") as f:
        fields = [
            "case_id",
            "category",
            "title",
            "evidence_kind",
            "source_urls",
            "original_evidence",
            "deliberate_modification",
            "expected_status",
            "expected_reason",
            "actual_status",
            "status_agrees",
            "required_flags",
            "actual_flags",
            "missing_required_flags",
            "review_status",
            "human_review_notes",
            "limitations",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for c, r in zip(cases, results, strict=True):
            w.writerow(
                {
                    "case_id": c["case_id"],
                    "category": c["category"],
                    "title": c["title"],
                    "evidence_kind": c["evidence_kind"],
                    "source_urls": "\n".join(s["url"] for s in c["sources"]),
                    "original_evidence": c["original_evidence"],
                    "deliberate_modification": c["deliberate_modification"],
                    "expected_status": c["expected_status"],
                    "expected_reason": c["expected_reason"],
                    "actual_status": r["actual_status"],
                    "status_agrees": r["status_agrees"],
                    "required_flags": ", ".join(c["required_flags"]),
                    "actual_flags": ", ".join(r["flags"]),
                    "missing_required_flags": ", ".join(r["missing_required_flags"]),
                    "review_status": r["review_status"],
                    "human_review_notes": "",
                    "limitations": c["limitations"],
                }
            )
    print(
        json.dumps(
            {
                k: report[k]
                for k in [
                    "total",
                    "status_agreement_count",
                    "status_and_diagnostic_agreement_count",
                    "execution_error_count",
                    "expected_counts",
                    "actual_counts",
                    "by_category",
                ]
            },
            indent=2,
        )
    )
    for r in results:
        if not r["diagnostic_agrees"]:
            print(json.dumps(r, indent=2))
    raise SystemExit(
        1
        if report["execution_error_count"]
        or args.fail_on_disagreement
        and report["status_and_diagnostic_agreement_count"] != len(results)
        else 0
    )


if __name__ == "__main__":
    main()
