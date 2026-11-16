from pathlib import Path

import pytest

from openchemproof.schema.models import AssayContext, VerificationRecord
from openchemproof.verifiers.common import compound_from_smiles


@pytest.fixture(autouse=True)
def isolated_audits(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("OPENCHEMPROOF_AUDIT_DIR", str(tmp_path / "audits"))


@pytest.fixture
def ethanol() -> VerificationRecord:
    return VerificationRecord(
        record_id="test-ethanol",
        compound=compound_from_smiles("ethanol", "CCO"),
        assay=AssayContext(assay_type="synthesis_qc", solvent="water"),
    )
