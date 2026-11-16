import pytest
from hypothesis import given
from hypothesis import strategies as st
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from openchemproof.schema.models import CompoundSpec
from openchemproof.verifiers.common import compound_from_smiles, exact_mass
from openchemproof.verifiers.structure_audit import audit_structure


def codes(result):
    return {f.code for f in result.flags}


@pytest.mark.parametrize("smiles", ["", "not-smiles", "C(C)(C)(C)(C)C"])
def test_invalid_structure(smiles):
    if not smiles:
        smiles = " "
    assert "INVALID_SMILES" in codes(audit_structure(CompoundSpec(compound_id="x", proposed_smiles=smiles)))


@pytest.mark.parametrize(
    "pair", [("N[C@H](C)C(=O)O", "N[C@@H](C)C(=O)O"), ("C/C=C/C", "C/C=C\\C"), ("CC(O)C(=O)O", "C[C@H](O)C(=O)O")]
)
def test_stereo(pair):
    spec = compound_from_smiles("x", pair[0]).model_copy(update={"reference_smiles": pair[1]})
    assert "STEREO_MISMATCH" in codes(audit_structure(spec))


def test_tautomer_equivalence():
    spec = compound_from_smiles("x", "Oc1ccccn1").model_copy(update={"reference_smiles": "O=c1cccc[nH]1"})
    result = audit_structure(spec)
    assert "TAUTOMER_EQUIVALENT" in codes(result)
    assert not any(f.severity == "CRITICAL_ARTIFACT" for f in result.flags)


def test_salt_mass_and_fragment_tracking():
    spec = compound_from_smiles("x", "C[NH2+]C.[Cl-]")
    result = audit_structure(spec)
    assert "SALT_STRIPPED" in codes(result)
    assert result.corrected_metrics["molecular_formula"] == "C2H7N"
    spec = spec.model_copy(update={"monoisotopic_mass": result.corrected_metrics["as_drawn_monoisotopic_mass"]})
    assert "SALT_MASS_MISMATCH" in codes(audit_structure(spec))


def test_inchi_connectivity_and_other_layers():
    a = compound_from_smiles("a", "N[C@H](C)C(=O)O")
    b = compound_from_smiles("b", "N[C@@H](C)C(=O)O")
    assert a.inchi_connectivity_layer == b.inchi_connectivity_layer
    assert "INCHI_LAYER_MISMATCH" in codes(audit_structure(a.model_copy(update={"inchikey": b.inchikey})))


def test_alerts_are_not_confirmed_artifacts():
    spec = compound_from_smiles("risk", "Oc1ccccc1O")
    result = audit_structure(spec)
    assert "STRUCTURAL_ALERT" in codes(result)
    assert all(f.severity != "CRITICAL_ARTIFACT" for f in result.flags)
    aromatic = compound_from_smiles("risk", "c1ccc(-c2ccc(-c3ccccc3)cc2)cc1")
    assert "AGGREGATION_RISK" in codes(audit_structure(aromatic))


@given(st.integers(min_value=1, max_value=25), st.sampled_from(["O", "N", "Cl", "Br", "S", "F"]))
def test_exact_mass_matches_rdkit(chain, end):
    smiles = "C" * chain + end
    assert exact_mass(smiles) == pytest.approx(rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(smiles)))


def test_unknown_mixture_abstains():
    spec = compound_from_smiles("mix", "c1ccccc1.CCCCCC")
    assert audit_structure(spec).inconclusive
