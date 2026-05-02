"""Standalone unit tests for `mf_latent_heat_vaporization`.

Purpose: pin the latent-heat-of-vaporization formula's behavior at known
reference temperatures so that any regression of the unit conversion
inside the function (e.g. someone reverts the K→C step) fails with a
specific, readable error rather than rippling out into the integration
tests as opaque downstream temperature drifts.

Background (2026-05-01 fix): the polynomial coefficients in the formula
(`2,499,999 J/kg` intercept, `-2385.74 J/kg/K` slope) are calibrated for
water temperature in degrees CELSIUS. The function's input, however, is
registered as Kelvin (per `tsm/dynamic_variables.py` `water_temp_k`).
Before the fix the formula was applied directly to the Kelvin value,
producing Lv ≈ 1.80 MJ/kg at 20°C instead of the correct ≈ 2.45 MJ/kg
(roughly 27% underestimate across typical surface-water temperatures).
Symptom: simulated water temperatures biased warm because evaporative
cooling was systematically too small. See
`design/tsm_stability_thin_water.md` for the broader TSM stability work
landed alongside this fix.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.tsm.processes import mf_latent_heat_vaporization


# Textbook reference values (J/kg) for the latent heat of vaporization
# of water. The simplified linear formula used in the kernel matches
# these to better than 0.2% across the 0–100 °C range. Sources: CRC
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
# the textbook reference at any point in 0–100 °C; the largest gap is at
# 100 °C (~0.2%).
RTOL = 0.005


def _lv(T_celsius: float) -> float:
    """Helper: invoke the formula on a scalar Kelvin DataArray and return
    the float result."""
    T_kelvin = T_celsius + 273.15
    return float(mf_latent_heat_vaporization(xr.DataArray(T_kelvin)))


@pytest.mark.parametrize("T_C, expected", list(REFERENCE_LV.items()))
def test_lv_matches_reference_at_known_temps(T_C: float, expected: float):
    """Lv at standard temperatures must match published reference values
    within 0.5%. This is the primary regression guard against the
    Kelvin-vs-Celsius unit bug: if someone reverts the K→C conversion,
    the 0 °C case alone (`Lv(273.15K) ≈ 1.85 MJ/kg`) fails by ~26%, far
    outside the 0.5% tolerance.
    """
    assert _lv(T_C) == pytest.approx(expected, rel=RTOL)


def test_lv_celsius_conversion_present():
    """Direct sanity check that the K→C conversion is in place.

    The formula `2,499,999 - 2385.74 * T_C` evaluated at T_C = 0 (i.e.
    T_K = 273.15) must give ~2.5 MJ/kg. If the conversion is missing, the
    formula is mistakenly applied to T_K = 273.15 directly, giving
    ~1.85 MJ/kg — that is the regression we want to catch.
    """
    lv_at_zero_C = _lv(0.0)
    assert lv_at_zero_C == pytest.approx(2_499_999.0, abs=1.0), (
        f"Lv at 0 °C = {lv_at_zero_C:.0f} J/kg; expected ~2.50e6 J/kg. "
        "If this is ~1.85e6 J/kg, the K→C conversion has regressed."
    )


def test_lv_monotonic_decreasing_with_temperature():
    """Lv decreases monotonically with temperature across the typical
    range. This catches sign-flip regressions in the slope coefficient."""
    temps_C = np.linspace(0.0, 50.0, 11)
    lvs = [_lv(T) for T in temps_C]
    for i in range(1, len(lvs)):
        assert lvs[i] < lvs[i - 1], (
            f"Lv not monotonically decreasing: "
            f"Lv({temps_C[i-1]:.1f}°C)={lvs[i-1]:.0f} < "
            f"Lv({temps_C[i]:.1f}°C)={lvs[i]:.0f}"
        )


def test_lv_arraylike_input_works():
    """Function must accept an array of Kelvin temperatures and return
    an array of Lv values, since the kernel calls it with cell-shaped
    DataArrays."""
    T_K_arr = xr.DataArray(np.array([273.15, 293.15, 313.15]))
    out = mf_latent_heat_vaporization(T_K_arr)
    expected = np.array([
        REFERENCE_LV[0.0],
        REFERENCE_LV[20.0],
        2_429_800.0 - 23.857 * 1000,  # ~Lv(40°C); rough check, loose tol
    ])
    out_np = np.asarray(out)
    assert out_np.shape == (3,)
    # Check first two against the table (tight); third is a rough monotonic check.
    assert out_np[0] == pytest.approx(REFERENCE_LV[0.0], rel=RTOL)
    assert out_np[1] == pytest.approx(REFERENCE_LV[20.0], rel=RTOL)
    assert out_np[2] < out_np[1]
