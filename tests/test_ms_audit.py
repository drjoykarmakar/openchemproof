import pytest
from hypothesis import given
from hypothesis import strategies as st

from openchemproof.schema.models import LCMSEvidence
from openchemproof.verifiers.common import compound_from_smiles
from openchemproof.verifiers.ms_audit import ADDUCTS, audit_lcms, isotope_envelope, ppm_error, theoretical_mz


@pytest.mark.parametrize("adduct", list(ADDUCTS))
def test_supported_adducts(adduct):
    compound = compound_from_smiles("x", "CCO")
    evidence = LCMSEvidence(
        observed_mz=theoretical_mz("CCO", adduct),
        ionization_mode=ADDUCTS[adduct].mode,
        observed_adduct=adduct,
        relative_isotope_intensities=isotope_envelope("CCO", adduct),
    )
    result = audit_lcms(compound, evidence)
    assert not result.inconclusive
    assert not any(f.severity == "CRITICAL_ARTIFACT" for f in result.flags)


@given(st.floats(min_value=20, max_value=2000), st.floats(min_value=-50, max_value=50))
def test_ppm_formula(mass, ppm):
    assert ppm_error(mass * (1 + ppm / 1e6), mass) == pytest.approx(ppm, abs=1e-8)


@pytest.mark.parametrize("ppm,failed", [(4.999, False), (5.0, False), (5.001, True), (-5.001, True)])
def test_gate_boundary(ppm, failed):
    mz = theoretical_mz("CCO", "[M+H]+")
    result = audit_lcms(
        compound_from_smiles("x", "CCO"),
        LCMSEvidence(
            observed_mz=mz * (1 + ppm * 1e-6),
            ionization_mode="positive",
            observed_adduct="[M+H]+",
            relative_isotope_intensities=isotope_envelope("CCO"),
        ),
    )
    assert ("HRMS_MASS_MISMATCH" in {f.code for f in result.flags}) is failed


def test_isotope_known_ratios_and_cross_terms():
    # Cl2 M+2/M = 2*p37/p35, with no hydrogen/carbon cross terms for [M+Na]+.
    assert isotope_envelope("ClCl", "[M+Na]+")[2] == pytest.approx(2 * 0.2424 / 0.7576)
    assert isotope_envelope("BrBr", "[M+Na]+")[2] == pytest.approx(2 * 0.4931 / 0.5069)
    mono = isotope_envelope("c1ccccc1", "[M+Na]+")
    assert mono[2] > 15 * (0.0107 / 0.9893) ** 2  # Includes H and mixed terms.
    assert isotope_envelope("[13CH4]", "[M+Na]+")[1] < 0.001


def test_starting_material_alternative_and_isobaric_abstention():
    starting = compound_from_smiles("a", "CCO")
    ev = LCMSEvidence(
        observed_mz=theoretical_mz("CCO", "[2M+H]+"),
        ionization_mode="positive",
        observed_adduct="[M+H]+",
        starting_materials=[starting],
        relative_isotope_intensities=[1.0, 0.04, 0.01],
    )
    result = audit_lcms(compound_from_smiles("target", "CCCCCCO"), ev)
    assert "STARTING_MATERIAL_MISASSIGNMENT" in {f.code for f in result.flags}
    ev = ev.model_copy(
        update={"observed_mz": theoretical_mz("CCO", "[M+H]+"), "relative_isotope_intensities": isotope_envelope("CCO")}
    )
    result = audit_lcms(compound_from_smiles("target", "COC"), ev)
    assert result.inconclusive
    assert "ISOBARIC_AMBIGUITY" in {f.code for f in result.flags}


def test_missing_isotope_evidence_abstains():
    ev = LCMSEvidence(observed_mz=47.049141, ionization_mode="positive")
    assert audit_lcms(compound_from_smiles("x", "CCO"), ev).inconclusive
