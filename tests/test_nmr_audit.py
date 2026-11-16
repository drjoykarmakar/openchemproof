import pytest

from openchemproof.schema.models import NMREvidence
from openchemproof.verifiers.common import compound_from_smiles
from openchemproof.verifiers.nmr_audit import audit_nmr, proton_budget


@pytest.mark.parametrize(
    "smiles,budget", [("CCO", (5, 1)), ("CCN", (5, 2)), ("CCS", (5, 1)), ("OCCO", (4, 2)), ("C[2H]", (3, 0))]
)
def test_hydrogen_counts(smiles, budget):
    assert proton_budget(smiles) == budget


@pytest.mark.parametrize(
    "solvent,total,critical",
    [
        ("D2O", 5.0, False),
        ("D2O", 6.0, True),
        ("CD3OD", 5.0, False),
        ("DMSO-d6", 6.0, False),
        ("DMSO-d6", 5.0, False),
        ("CDCl3", 5.5, False),
    ],
)
def test_exchange_rules(solvent, total, critical):
    result = audit_nmr(
        compound_from_smiles("x", "CCO"),
        NMREvidence(solvent=solvent, proton_shifts_ppm=[3.6], proton_integrals=[total], integral_scale=1.0),
    )
    assert any(f.severity == "CRITICAL_ARTIFACT" for f in result.flags) is critical


def test_uncalibrated_abstains():
    result = audit_nmr(
        compound_from_smiles("x", "CCO"), NMREvidence(solvent="CDCl3", proton_shifts_ppm=[3.6], proton_integrals=[1.0])
    )
    assert result.inconclusive


def test_candidate_solvent_is_not_automatic_contamination():
    compound = compound_from_smiles("ethyl_acetate", "CCOC(C)=O")
    ev = NMREvidence(
        solvent="CDCl3", proton_shifts_ppm=[4.12, 2.05, 1.26], proton_integrals=[2.0, 3.0, 3.0], integral_scale=1.0
    )
    result = audit_nmr(compound, ev)
    assert "SOLVENT_PEAK_CANDIDATE" in {f.code for f in result.flags}
    assert not any(f.severity == "CRITICAL_ARTIFACT" for f in result.flags)


def test_carbon_symmetry_and_rotamer_limit():
    compound = compound_from_smiles("benzene", "c1ccccc1")
    ev = NMREvidence(
        solvent="CDCl3",
        proton_shifts_ppm=[7.3],
        proton_integrals=[6.0],
        integral_scale=1.0,
        carbon_shifts_ppm=[128.0, 77.16],
    )
    assert not audit_nmr(compound, ev).inconclusive
    result = audit_nmr(compound, ev.model_copy(update={"carbon_shifts_ppm": [float(x) for x in range(100, 108)]}))
    assert result.inconclusive
    assert "EXCESS_CARBON_RESONANCES" in {f.code for f in result.flags}
