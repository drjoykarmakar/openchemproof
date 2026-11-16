"""Canonical JSON, lossless Parquet payloads, and append-only content-addressed audits."""

import hashlib
import importlib.metadata
import json
import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
from pydantic import BaseModel

from openchemproof.schema.models import AuditBundle, AuditPolicy, AuditVerdict, VerificationRecord


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def json_schema() -> dict[str, Any]:
    return AuditBundle.model_json_schema()


def write_parquet(records: Sequence[BaseModel], path: Path) -> None:
    """JSON column preserves strict nested types; schema is embedded in Arrow metadata."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payloads = [canonical_json(item.model_dump(mode="json")) for item in records]
    metadata = {
        b"openchemproof_schema": canonical_json(records[0].model_json_schema() if records else {}).encode(),
        b"format": b"canonical-json-envelope-v1",
    }
    schema = pa.schema([("payload_json", pa.string()), ("sha256", pa.string())], metadata=metadata)
    table = pa.Table.from_pydict(
        {
            "payload_json": payloads,
            "sha256": [hashlib.sha256(item.encode()).hexdigest() for item in payloads],
        },
        schema=schema,
    )
    pq.write_table(table, path, compression="zstd")


def read_parquet[T: BaseModel](path: Path, model: type[T]) -> list[T]:
    table = pq.read_table(path).to_pydict()
    result = []
    for payload, digest in zip(table["payload_json"], table["sha256"], strict=True):
        if hashlib.sha256(payload.encode()).hexdigest() != digest:
            raise ValueError("Parquet row integrity failure")
        result.append(model.model_validate_json(payload))
    return result


def _bundle_payload(bundle: AuditBundle) -> dict[str, Any]:
    payload = bundle.model_dump(mode="json", exclude={"sha256"})
    payload["verdict"]["audit_sha256"] = ""
    return payload


def verify_bundle(bundle: AuditBundle) -> bool:
    calculated = hashlib.sha256(canonical_json(_bundle_payload(bundle)).encode()).hexdigest()
    return calculated == bundle.sha256 == bundle.verdict.audit_sha256


def make_bundle(record: VerificationRecord, policy: AuditPolicy, verdict: AuditVerdict) -> AuditBundle:
    versions = {
        name: importlib.metadata.version(name)
        for name in ("rdkit", "pydantic", "fastmcp", "openchemproof", "numpy", "scipy")
    }
    provisional = AuditBundle(input=record, policy=policy, verdict=verdict, runtime_versions=versions, sha256="0" * 64)
    digest = hashlib.sha256(canonical_json(_bundle_payload(provisional)).encode()).hexdigest()
    sealed_verdict = AuditVerdict.model_validate({**verdict.model_dump(), "audit_sha256": digest})
    return AuditBundle(
        input=record,
        policy=policy,
        verdict=sealed_verdict,
        runtime_versions=versions,
        sha256=digest,
    )


def persist_bundle(bundle: AuditBundle, directory: Path) -> Path:
    """Exclusive-create snapshots; repeated identical inputs are idempotent.

    Hashes detect modification. An OS owner can still delete/replace files; for
    regulated immutability use external WORM storage and a signed manifest.
    """
    if not verify_bundle(bundle):
        raise ValueError("Refusing an invalid audit digest")
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / f"{bundle.sha256}.json"
    contents = canonical_json(bundle.model_dump(mode="json")) + "\n"
    # Write complete temporary bytes, then atomically link without replacing any destination.
    import tempfile

    descriptor, temporary = tempfile.mkstemp(prefix=".audit-", dir=directory)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(contents)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, 0o444)
        try:
            os.link(temporary, destination)
        except FileExistsError:
            if destination.read_text(encoding="utf-8") != contents:
                raise ValueError("Existing audit was modified; refusing overwrite") from None
    finally:
        Path(temporary).unlink(missing_ok=True)
    return destination
