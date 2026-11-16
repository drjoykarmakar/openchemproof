"""Typed MCP tools; all calls use the same deterministic, audited Python pipeline."""

import argparse
import json
from dataclasses import asdict

from fastmcp import FastMCP

from openchemproof.schema.export import json_schema
from openchemproof.schema.models import AuditPolicy, AuditVerdict, VerificationRecord
from openchemproof.verifiers.ms_audit import ADDUCTS
from openchemproof.verifiers.nmr_audit import CONTAMINANTS, RESIDUAL_CARBON, RESIDUAL_PROTON
from openchemproof.verifiers.pipeline import run_verification

mcp = FastMCP(
    "OpenChemProof",
    instructions=(
        "Deterministic checks of supplied analytical evidence. VERIFIED denotes consistency only. "
        "Warnings are risks, not confirmed false positives. Negative results require explicit evidence. "
        "Every tool call persists a content-addressed audit locally."
    ),
)


@mcp.tool()
def verify_structure_identity(record: VerificationRecord) -> AuditVerdict:
    """Audit parent/salt mass, asserted identity, stereo/tautomer relations and structural risk alerts."""
    return run_verification(record, modalities=("structure",))


@mcp.tool()
def audit_lcms_spectrum(record: VerificationRecord) -> AuditVerdict:
    """Audit measured MS1 precursor/adduct/isotope evidence; MS2 fragments are not isotope peaks."""
    if record.evidence.lcms is None:
        raise ValueError("LCMS evidence is required")
    return run_verification(record, modalities=("structure", "lcms"))


@mcp.tool()
def audit_optical_spectrum(record: VerificationRecord) -> AuditVerdict:
    """Audit fluorescence physics and differential IFE using matched reference absorbances."""
    if record.evidence.optical is None:
        raise ValueError("Optical evidence is required")
    return run_verification(record, modalities=("structure", "optical"))


@mcp.tool()
def audit_nmr_consistency(record: VerificationRecord) -> AuditVerdict:
    """Reconcile independently calibrated integrals and solvent exchange rules; flag candidate impurities."""
    if record.evidence.nmr is None:
        raise ValueError("NMR evidence is required")
    return run_verification(record, modalities=("structure", "nmr"))


@mcp.tool()
def run_full_verification_pipeline(record: VerificationRecord) -> AuditVerdict:
    """Run all supplied modalities, preserve documented negative outcomes, and seal the audit snapshot."""
    return run_verification(record)


@mcp.resource("openchemproof://schema/v1")
def verification_schema() -> str:
    """JSON Schema for complete version-1 input, policy, verdict and audit bundle."""
    return json.dumps(json_schema(), sort_keys=True)


@mcp.resource("openchemproof://reference/contaminants")
def reference_tables() -> str:
    """Approximate solvent shifts, singly charged adduct definitions and default audit thresholds."""
    return json.dumps(
        {
            "adducts": {key: asdict(value) for key, value in ADDUCTS.items()},
            "residual_proton_ppm": RESIDUAL_PROTON,
            "residual_carbon_ppm": RESIDUAL_CARBON,
            "contaminants_ppm": CONTAMINANTS,
            "policy": AuditPolicy().model_dump(),
        },
        sort_keys=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transport", choices=["stdio", "http"], default="stdio")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    if args.transport == "http":
        mcp.run(transport="http", host=args.host, port=args.port, show_banner=False)
    else:
        mcp.run(transport="stdio", show_banner=False)


if __name__ == "__main__":
    main()
