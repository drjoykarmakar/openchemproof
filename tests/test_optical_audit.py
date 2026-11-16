import pytest
from hypothesis import given
from hypothesis import strategies as st

from openchemproof.schema.models import OpticalEvidence, SolventObservation
from openchemproof.verifiers.optical_audit import audit_optical, ife_factor


def optical(**updates):
    return OpticalEvidence(
        excitation_nm=450.0,
        emission_nm=520.0,
        absorbance_at_ex=0.5,
        absorbance_at_em=0.3,
        baseline_absorbance_at_ex=0.0,
        baseline_absorbance_at_em=0.0,
        apparent_fold_change=10**-0.4,
    ).model_copy(update=updates)


@given(
    st.floats(min_value=0, max_value=1), st.floats(min_value=0, max_value=1), st.floats(min_value=1e-3, max_value=1e6)
)
def test_ife_restores_intensity(aex, aem, true_intensity):
    observed = true_intensity * 10 ** (-(aex + aem) / 2)
    assert observed * ife_factor(aex, aem) == pytest.approx(true_intensity)


def test_known_correction_and_reference():
    assert ife_factor(0.2, 0.4) == pytest.approx(1.9952623149688795)
    result = audit_optical(optical())
    assert result.corrected_metrics["corrected_fold_change"] == pytest.approx(1.0)
    assert "IFE_RESPONSE_ARTIFACT" in {f.code for f in result.flags}
    result = audit_optical(
        optical(baseline_absorbance_at_ex=0.5, baseline_absorbance_at_em=0.3, apparent_fold_change=0.5)
    )
    assert result.corrected_metrics["corrected_fold_change"] == pytest.approx(0.5)
    assert "IFE_RESPONSE_ARTIFACT" not in {f.code for f in result.flags}


@pytest.mark.parametrize(
    "updates,code",
    [
        ({"observed_quantum_yield": -0.1}, "UNPHYSICAL_QUANTUM_YIELD"),
        ({"observed_quantum_yield": 1.1}, "UNPHYSICAL_QUANTUM_YIELD"),
        ({"emission_nm": 440.0}, "ANTI_STOKES_UNEXPLAINED"),
        ({"emission_nm": 455.0}, "RAYLEIGH_OVERLAP"),
        ({"absorbance_at_ex": 1000.0}, "IFE_OUTSIDE_VALIDITY"),
        ({"baseline_absorbance_at_ex": None, "baseline_absorbance_at_em": None}, "MISSING_REFERENCE_ABSORBANCE"),
    ],
)
def test_optical_flags(updates, code):
    assert code in {f.code for f in audit_optical(optical(**updates)).flags}


def test_solvatochromic_trend():
    series = [
        SolventObservation(solvent=str(i), polarity_index=float(i), emission_nm=float(600 - i * 30)) for i in range(3)
    ]
    result = audit_optical(optical(solvent_series=series, expected_solvatochromism="positive"))
    assert "SOLVATOCHROMIC_TREND" in {f.code for f in result.flags}


@pytest.mark.parametrize("value", [-0.1, float("inf"), float("nan"), 1000.0])
def test_numeric_guard(value):
    with pytest.raises(ValueError):
        ife_factor(value, value)
