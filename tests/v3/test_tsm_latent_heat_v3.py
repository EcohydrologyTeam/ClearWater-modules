"""Standalone unit tests for v3 ``Temperature.latent_heat_vaporization``.

Purpose: pin the latent-heat-of-vaporization formula's behavior at known
reference temperatures so that any regression of the unit convention
inside the function (e.g. someone re-introduces a Kelvin offset) fails
with a specific, readable error rather than rippling out into the
integration tests as opaque downstream temperature drifts.

Background: v3's ``Temperature.latent_heat_vaporization`` takes water
temperature directly in **degrees Celsius**. The polynomial coefficients
(``2,499,999 J/kg`` intercept, ``-2385.74 J/kg/K`` slope) are calibrated
for Celsius. v2's pre-fix form applied the polynomial to Kelvin (or
equivalently, ``2.5e6 - 2385.74 * (T_C + 273.16)``), producing
Lv ~ 1.80 MJ/kg at 20 C instead of the correct ~ 2.45 MJ/kg
(roughly 27 percent underestimate across typical surface-water
temperatures). Symptom: simulated water temperatures biased warm
because evaporative cooling was systematically too small.

This is the v3 port of ``tests/test_tsm_latent_heat.py`` (v1). The v1
test passes Kelvin and the kernel converts internally; the v3 test
passes Celsius directly.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.temperature import Temperature


# Textbook reference values (J/kg) for the latent heat of vaporization
# of water. The simplified linear formula used in the kernel matches
# these to better than 0.5% across the 0-100 C range. Sources: CRC
# Handbook 92nd ed., NIST Webbook for water saturation states.
REFERENCE_LV = {
    0.0:   2_500_900.0,
    10.0:  2_477_300.0,
    20.0:  2_453_500.0,
    25.0:  2_441_700.0,
    30.0:  2_429_800.0,
    100.0: 2_257_000.0,
}

# Tolerance: 0.5% covers the simplified linear formula's deviation from
# the textbook reference at any point in 0-100 C; the largest gap is at
# 100 C (~0.2%).
RTOL = 0.005


@pytest.fixture(scope="module")
def temperature_module() -> Temperature:
    """Stub Temperature instance.

    The latent-heat-of-vaporization method does not depend on wind or
    sediment parameters, so any valid stub values for the constructor
    are fine.
    """
    return Temperature(0.3, 1.5, 3.0)


def _lv(temp: Temperature, T_celsius: float) -> float:
    """Helper: invoke the formula on a scalar Celsius input and return
    the float result."""
    return float(temp.latent_heat_vaporization(T_celsius))


@pytest.mark.parametrize("T_C, expected", list(REFERENCE_LV.items()))
def test_lv_matches_reference_at_known_temps(
    temperature_module: Temperature, T_C: float, expected: float
):
    """Lv at standard temperatures must match published reference values
    within 0.5%. This is the primary regression guard against the
    Kelvin-vs-Celsius unit bug: if someone re-introduces a Kelvin
    offset, the 0 C case alone (``Lv(273.15) ~ 1.85 MJ/kg``) fails by
    ~26%, far outside the 0.5% tolerance.
    """
    assert _lv(temperature_module, T_C) == pytest.approx(expected, rel=RTOL)


def test_lv_celsius_input_no_kelvin_offset(temperature_module: Temperature):
    """Direct sanity check that the input is treated as Celsius (no
    internal Kelvin offset).

    The formula ``2,499,999 - 2385.74 * T_C`` evaluated at T_C = 0 must
    give ~2.5 MJ/kg. If a Kelvin offset is re-introduced (i.e. the
    formula is mistakenly applied to ``T_C + 273.15`` or
    ``T_C + 273.16``), the result at 0 C is ~1.85 MJ/kg or
    ~1.80 MJ/kg respectively -- those are the regressions we want to
    catch.
    """
    lv_at_zero_C = _lv(temperature_module, 0.0)
    assert lv_at_zero_C == pytest.approx(2_499_999.0, abs=1.0), (
        f"Lv at 0 C = {lv_at_zero_C:.0f} J/kg; expected ~2.50e6 J/kg. "
        "If this is ~1.85e6 J/kg, a Kelvin offset has been "
        "re-introduced (T+273.15); if ~1.80e6 J/kg, the v2 pre-fix "
        "form (T+273.16) has regressed."
    )


def test_lv_not_prefix_v2_value(temperature_module: Temperature):
    """Anti-regression: the result at 20 C must NOT match the pre-fix
    v2 form ``2.5e6 - 2385.74 * (T_C + 273.16)``, which gives
    ~1.80 MJ/kg.
    """
    prefix_v2_at_20 = 2.5e6 - 2385.74 * (20.0 + 273.16)
    lv_at_20 = _lv(temperature_module, 20.0)
    assert abs(lv_at_20 - prefix_v2_at_20) > 1.0e5, (
        f"Lv at 20 C = {lv_at_20:.0f} J/kg matches the pre-fix v2 "
        f"value {prefix_v2_at_20:.0f} J/kg; the Kelvin-offset bug "
        "has regressed."
    )
    # And it should match the documented expected value ~2.45 MJ/kg.
    assert lv_at_20 == pytest.approx(2_452_284.0, abs=1.0)


def test_lv_documented_expected_values(temperature_module: Temperature):
    """Pin the exact expected values quoted in the v3 source docstring.

    Per ``Temperature.latent_heat_vaporization`` docstring:
      Lv(0 C)  = 2,499,999 J/kg
      Lv(20 C) = 2,452,284 J/kg  (== 2499999 - 2385.74 * 20)
      Lv(25 C) = 2,440,357 J/kg  (== 2499999 - 2385.74 * 25, rounded)
    """
    assert _lv(temperature_module, 0.0) == pytest.approx(2_499_999.0, abs=1.0)
    assert _lv(temperature_module, 20.0) == pytest.approx(2_452_284.0, abs=1.0)
    # 2499999 - 2385.74 * 25 = 2440355.5 J/kg; the docstring rounds to
    # 2,440,357. Use a small absolute tolerance to cover the rounding.
    assert _lv(temperature_module, 25.0) == pytest.approx(2_440_357.0, abs=2.0)


def test_lv_monotonic_decreasing_with_temperature(
    temperature_module: Temperature,
):
    """Lv decreases monotonically with temperature across the typical
    range. This catches sign-flip regressions in the slope coefficient.
    """
    temps_C = np.linspace(0.0, 50.0, 11)
    lvs = [_lv(temperature_module, T) for T in temps_C]
    for i in range(1, len(lvs)):
        assert lvs[i] < lvs[i - 1], (
            f"Lv not monotonically decreasing: "
            f"Lv({temps_C[i-1]:.1f}C)={lvs[i-1]:.0f} < "
            f"Lv({temps_C[i]:.1f}C)={lvs[i]:.0f}"
        )


def test_lv_arraylike_input_works(temperature_module: Temperature):
    """Function must accept array-shaped Celsius temperatures and return
    array-shaped Lv values, since the kernel calls it with cell-shaped
    DataArrays. Covers scalar, ndarray, and xr.DataArray inputs.

    Bare Python lists are not part of the ``ArrayLike`` contract: scalar
    broadcasting (``2.499999e6 - 2385.74 * x``) requires a numeric type
    that supports scalar multiplication, which lists don't. v1's
    ``mf_latent_heat_vaporization`` has the same limitation. The kernel
    always calls this method with xarray DataArrays or numpy arrays, so
    list input isn't a needed contract.
    """
    # Scalar
    out_scalar = temperature_module.latent_heat_vaporization(20.0)
    assert float(out_scalar) == pytest.approx(REFERENCE_LV[20.0], rel=RTOL)

    # numpy ndarray
    T_C_arr = np.array([0.0, 20.0, 40.0])
    out_np = np.asarray(temperature_module.latent_heat_vaporization(T_C_arr))
    assert out_np.shape == (3,)
    assert out_np[0] == pytest.approx(REFERENCE_LV[0.0], rel=RTOL)
    assert out_np[1] == pytest.approx(REFERENCE_LV[20.0], rel=RTOL)
    assert out_np[2] < out_np[1]

    # xarray DataArray
    T_C_da = xr.DataArray(np.array([0.0, 20.0, 40.0]))
    out_da = temperature_module.latent_heat_vaporization(T_C_da)
    out_da_np = np.asarray(out_da)
    assert out_da_np.shape == (3,)
    assert out_da_np[0] == pytest.approx(REFERENCE_LV[0.0], rel=RTOL)
    assert out_da_np[1] == pytest.approx(REFERENCE_LV[20.0], rel=RTOL)
    assert out_da_np[2] < out_da_np[1]
