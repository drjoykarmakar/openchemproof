"""Versioned, strict analytical observations. Reported assertions are never trusted calculations."""

from typing import Annotated, Literal, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

Finite = Annotated[float, Field(allow_inf_nan=False)]
NonNegative = Annotated[Finite, Field(ge=0)]
Positive = Annotated[Finite, Field(gt=0)]
Status = Literal["VERIFIED", "ARTIFACT_DETECTED", "INCONCLUSIVE", "NEGATIVE_RESULT"]
Modality = Literal["structure", "lcms", "optical", "nmr", "provenance"]


class StrictModel(BaseModel):
    model_config = ConfigDict(strict=True, extra="forbid", frozen=True, validate_default=True)


class CompoundSpec(StrictModel):
    compound_id: str = Field(min_length=1)
    proposed_smiles: str = Field(min_length=1, max_length=10000)
    # Optional assertions are checked independently; omitted values are calculated.
    canonical_smiles: str | None = None
    inchikey: str | None = Field(default=None, pattern=r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
    inchi_connectivity_layer: str | None = Field(default=None, pattern=r"^[A-Z]{14}$")
    molecular_formula: str | None = None
    monoisotopic_mass: Positive | None = None
    salt_form: str | None = None
    mass_basis: Literal["neutral_parent", "as_drawn"] = "neutral_parent"
    expected_stereocenters: int = Field(default=0, ge=0)
    reference_smiles: str | None = None
    ord_compound_identifier: str | None = None


class AssayContext(StrictModel):
    assay_type: Literal["synthesis_qc", "fluorescence_probe", "binding_assay", "reaction_screen"]
    solvent: str = Field(min_length=1)
    ph: Finite | None = Field(default=None, ge=-2, le=16)
    temperature_c: Finite = Field(default=25.0, gt=-273.15)
    concentration_um: Positive | None = None


class LCMSEvidence(StrictModel):
    observed_mz: Positive
    ionization_mode: Literal["positive", "negative"]
    observed_adduct: str | None = None
    relative_isotope_intensities: list[NonNegative] = Field(default_factory=list, max_length=3)
    retention_time_min: NonNegative | None = None
    purity_uv_percent: Finite | None = Field(default=None, ge=0, le=100)
    starting_materials: list[CompoundSpec] = Field(default_factory=list, max_length=20)

    @model_validator(mode="after")
    def isotope_triplet(self) -> Self:
        if len(self.relative_isotope_intensities) not in (0, 3):
            raise ValueError("Supply M, M+1, M+2 together, or omit isotope evidence")
        if self.relative_isotope_intensities and self.relative_isotope_intensities[0] <= 0:
            raise ValueError("M must be positive to normalize isotope ratios")
        return self


class SolventObservation(StrictModel):
    solvent: str
    polarity_index: Finite
    emission_nm: Positive


class OpticalEvidence(StrictModel):
    excitation_nm: Positive
    emission_nm: Positive
    absorbance_at_ex: NonNegative
    absorbance_at_em: NonNegative
    # Keep physically invalid QYs representable so the engine can audit them.
    observed_quantum_yield: Finite | None = None
    apparent_fold_change: NonNegative | None = None
    baseline_absorbance_at_ex: NonNegative | None = None
    baseline_absorbance_at_em: NonNegative | None = None
    geometry: Literal["right_angle_center", "other"] = "right_angle_center"
    upconversion_declared: bool = False
    solvent_raman_shift_cm: Positive | None = None
    solvent_series: list[SolventObservation] = Field(default_factory=list)
    expected_solvatochromism: Literal["positive", "negative"] | None = None

    @model_validator(mode="after")
    def paired_reference(self) -> Self:
        if (self.baseline_absorbance_at_ex is None) != (self.baseline_absorbance_at_em is None):
            raise ValueError("Both reference absorbances are required together")
        return self


class NMREvidence(StrictModel):
    solvent: Literal["CDCl3", "DMSO-d6", "CD3OD", "D2O", "Acetone-d6"]
    proton_shifts_ppm: list[Finite] = Field(max_length=2000)
    proton_integrals: list[NonNegative] = Field(max_length=2000)
    carbon_shifts_ppm: list[Finite] = Field(default_factory=list, max_length=2000)
    integral_scale: Positive | None = None  # integral units per proton; never infer from total
    multiplicities: list[Literal["s", "d", "t", "q", "m", "br", "unknown"]] = Field(default_factory=list)
    integration_complete: bool = True
    exchangeable_protons_observed: bool | None = None

    @model_validator(mode="after")
    def aligned_peaks(self) -> Self:
        if len(self.proton_shifts_ppm) != len(self.proton_integrals):
            raise ValueError("Proton shifts and integrals must have equal length")
        if self.multiplicities and len(self.multiplicities) != len(self.proton_shifts_ppm):
            raise ValueError("Multiplicities must align with proton peaks")
        return self


class SpectralEvidence(StrictModel):
    lcms: LCMSEvidence | None = None
    optical: OpticalEvidence | None = None
    nmr: NMREvidence | None = None


class NegativeResultRecord(StrictModel):
    experiment_id: str = Field(min_length=1)
    outcome: Literal["no_product", "decomposition", "no_response", "failed_synthesis"]
    observation: str = Field(min_length=10)
    protocol_reference: str = Field(min_length=1)
    evidence_references: list[str] = Field(min_length=1)
    detection_limit_description: str = Field(min_length=5)
    controls_passed: bool
    # Zero yield is not evidence of absence unless the detection limit and controls are given.
    measured_yield_percent: Finite | None = Field(default=None, ge=0, le=100)


class VerificationRecord(StrictModel):
    schema_version: Literal["1.0"] = "1.0"
    record_id: str = Field(min_length=1)
    compound: CompoundSpec
    assay: AssayContext
    evidence: SpectralEvidence = Field(default_factory=SpectralEvidence)
    negative_result: NegativeResultRecord | None = None
    ord_reaction_id: str | None = Field(default=None, pattern=r"^ord-[0-9a-f]{32}$")
    source_references: list[str] = Field(default_factory=list)


class AuditFlag(StrictModel):
    code: str
    severity: Literal["INFO", "WARNING", "CRITICAL_ARTIFACT"]
    modality: Modality
    message: str
    quantitative_evidence: dict[str, Finite | str] = Field(default_factory=dict)


class ModalityResult(StrictModel):
    modality: Modality
    flags: list[AuditFlag] = Field(default_factory=list)
    corrected_metrics: dict[str, Finite | str] = Field(default_factory=dict)
    checks_performed: list[str] = Field(default_factory=list)
    inconclusive: bool = False


class AuditPolicy(StrictModel):
    policy_version: Literal["1.0"] = "1.0"
    hrms_ppm_tolerance: Positive = 5.0
    alternative_ppm_tolerance: Positive = 3.0
    isotope_absolute_tolerance: Positive = 0.05
    isotope_relative_tolerance: Positive = 0.25
    ife_absorbance_trigger: Positive = 0.15
    ife_max_absorbance: Positive = 1.0
    ife_baseline_tolerance: Positive = 0.15
    nmr_absolute_h_tolerance: Positive = 0.35
    nmr_relative_h_tolerance: Positive = 0.05


class AuditVerdict(StrictModel):
    record_id: str
    status: Status
    confidence_score: Finite = Field(ge=0, le=1)
    flags: list[AuditFlag]
    corrected_metrics: dict[str, Finite | str]
    deterministic_engine_version: str
    checks_performed: list[str] = Field(default_factory=list)
    confidence_basis: str = "Deterministic rule coverage; not a calibrated probability of identity"
    audit_sha256: str = Field(default="", pattern=r"^([a-f0-9]{64})?$")


class AuditBundle(StrictModel):
    schema_version: Literal["1.0"] = "1.0"
    input: VerificationRecord
    policy: AuditPolicy
    verdict: AuditVerdict
    runtime_versions: dict[str, str]
    sha256: str = Field(pattern=r"^[a-f0-9]{64}$")
