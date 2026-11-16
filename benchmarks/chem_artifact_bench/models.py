"""Labels live outside the model-visible observations."""

from typing import Literal

from pydantic import Field

from openchemproof.schema.models import Status, StrictModel, VerificationRecord

Category = Literal[
    "STEREO_TAUTOMER_TRAP",
    "LCMS_ADDUCT_DIMER_TRAP",
    "OPTICAL_IFE_ARTIFACT",
    "NMR_SOLVENT_INTEGRAL_TRAP",
    "VERIFIED_TRUE_AND_NEGATIVE",
]
CATEGORIES: tuple[Category, ...] = (
    "STEREO_TAUTOMER_TRAP",
    "LCMS_ADDUCT_DIMER_TRAP",
    "OPTICAL_IFE_ARTIFACT",
    "NMR_SOLVENT_INTEGRAL_TRAP",
    "VERIFIED_TRUE_AND_NEGATIVE",
)


class BenchmarkCase(StrictModel):
    case_id: str
    category: Category
    subtype: str
    expected_status: Status
    expected_flag_codes: list[str] = Field(default_factory=list)
    evidence_origin: Literal["synthetic_unit_fixture"] = "synthetic_unit_fixture"
    record: VerificationRecord
