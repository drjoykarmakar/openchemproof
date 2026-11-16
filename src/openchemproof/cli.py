"""Local JSON auditing, schema export and audit-integrity verification."""

import argparse
import json
from pathlib import Path

from openchemproof.schema.export import json_schema, verify_bundle
from openchemproof.schema.models import AuditBundle, VerificationRecord
from openchemproof.verifiers.pipeline import run_verification


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    audit = commands.add_parser("verify")
    audit.add_argument("input", type=Path)
    audit.add_argument("--audit-dir", type=Path)
    commands.add_parser("schema")
    integrity = commands.add_parser("check-audit")
    integrity.add_argument("input", type=Path)
    args = parser.parse_args()
    if args.command == "schema":
        print(json.dumps(json_schema(), indent=2))
    elif args.command == "check-audit":
        bundle = AuditBundle.model_validate_json(args.input.read_text())
        if not verify_bundle(bundle):
            raise SystemExit("Audit integrity check failed")
        print(bundle.sha256)
    else:
        record = VerificationRecord.model_validate_json(args.input.read_text())
        verdict = run_verification(record, audit_dir=args.audit_dir)
        print(verdict.model_dump_json(indent=2))


if __name__ == "__main__":
    main()
