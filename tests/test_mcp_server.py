import json

import pytest
from fastmcp import Client
from fastmcp.exceptions import ToolError

from openchemproof.mcp_server.server import mcp
from openchemproof.schema.models import LCMSEvidence, NMREvidence, OpticalEvidence, SpectralEvidence
from openchemproof.verifiers.ms_audit import isotope_envelope, theoretical_mz


async def test_every_tool_and_resource(ethanol):
    evidence = SpectralEvidence(
        lcms=LCMSEvidence(
            observed_mz=theoretical_mz("CCO", "[M+H]+"),
            observed_adduct="[M+H]+",
            ionization_mode="positive",
            relative_isotope_intensities=isotope_envelope("CCO"),
        ),
        optical=OpticalEvidence(excitation_nm=400.0, emission_nm=500.0, absorbance_at_ex=0.01, absorbance_at_em=0.01),
        nmr=NMREvidence(solvent="D2O", proton_shifts_ppm=[1.2, 3.6], proton_integrals=[3.0, 2.0], integral_scale=1.0),
    )
    record = ethanol.model_copy(update={"evidence": evidence})
    async with Client(mcp) as client:
        tools = await client.list_tools()
        assert len(tools) == 5
        for tool in tools:
            result = await client.call_tool(tool.name, {"record": record.model_dump(mode="json")})
            assert result.structured_content["record_id"] == record.record_id
            assert len(result.structured_content["audit_sha256"]) == 64
        resources = await client.list_resources()
        assert len(resources) == 2
        for resource in resources:
            result = await client.read_resource(str(resource.uri))
            assert json.loads(result[0].text)


async def test_strict_mcp_rejects_coercion(ethanol):
    record = ethanol.model_dump(mode="json")
    record["assay"]["temperature_c"] = "25.0"
    async with Client(mcp) as client:
        with pytest.raises(ToolError):
            await client.call_tool("run_full_verification_pipeline", {"record": record})


async def test_missing_modality_is_explicit_error(ethanol):
    async with Client(mcp) as client:
        with pytest.raises(ToolError):
            await client.call_tool("audit_lcms_spectrum", {"record": ethanol.model_dump(mode="json")})
