"""Standalone unit tests for v3 ``Temperature`` minor-finding cleanups.

These tests pin behavior introduced or preserved by the MINOR
review-finding cleanups applied to
``src/clearwater_modules_v3/processes/temperature.py`` on 2026-05-04.
The corresponding findings are documented in
``design/clearwater_modules_v3_review_findings.md`` (m1, m2, m4, m6,
m9, m19) and overlap-finding m7 (which collapses with m1).

Each test pins a single behavior so that any regression fails with a
specific, readable message rather than as a noisy aggregate.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.temperature import Temperature


@pytest.fixture(scope="module")
def temperature_module() -> Temperature:
    """Stub Temperature instance.

    Wind parameters are arbitrary stub values; none of the methods
    exercised in this file depend on the specific wind coefficients.
    """
    return Temperature(0.3, 1.5, 3.0)


# ---------- m1 / m7: simplified ``mixing_ratio_air`` form ----------


def test_mixing_ratio_air_simplified_form_preserves_c4_cases(
    temperature_module: Temperature,
) -> None:
    """Pin that the simplified single-predicate form of
    ``mixing_ratio_air`` still satisfies the four C4 cases (normal,
    equality, negative-denominator, vectorized-mixed).

    The cleanup collapsed the prior nested ``xr.where`` into one
    divide-safe denominator and one outer guard. Behavior must be
    identical to the prior nested form.
    """
    # 1. Normal case: e_air = 20 mb, P_air = 1013 mb returns the
    # standard formula value to floating-point precision.
    e_normal = 20.0
    P_normal = 1013.0
    expected_normal = 0.622 * e_normal / (P_normal - e_normal)
    result_normal = float(
        temperature_module.mixing_ratio_air(e_normal, P_normal)
    )
    assert result_normal == pytest.approx(
        expected_normal, rel=1e-12, abs=1e-15
    ), (
        f"normal-case regression: mixing_ratio_air({e_normal}, {P_normal}) "
        f"= {result_normal}; expected {expected_normal}"
    )

    # 2. Equality case: e_air == P_air -> denom == 0, guard returns 0.0.
    result_eq = float(temperature_module.mixing_ratio_air(1013.0, 1013.0))
    assert result_eq == 0.0, (
        f"equality-case guard regression: got {result_eq}; expected 0.0"
    )

    # 3. Negative-denominator case (C4 fix): e_air > P_air -> guard
    # returns 0.0 rather than a negative mixing ratio.
    result_neg = float(temperature_module.mixing_ratio_air(1100.0, 1013.0))
    assert result_neg == 0.0, (
        f"C4 negative-denominator guard regression: got {result_neg}; "
        "expected 0.0. A negative result (~ -7.86) indicates the C4 "
        "fix has regressed."
    )
    assert result_neg >= 0.0, (
        "C4 guard regression: mixing_ratio_air returned a negative value "
        f"({result_neg}); the simplified form must still return 0.0 here."
    )

    # 4. Vectorized mixed case: a DataArray with cells in normal /
    # equality / e>P regimes returns [normal, 0.0, 0.0].
    e_vec = xr.DataArray(np.array([20.0, 1013.0, 1100.0]))
    P_vec = xr.DataArray(np.array([1013.0, 1013.0, 1013.0]))
    out = temperature_module.mixing_ratio_air(e_vec, P_vec)
    out_np = np.asarray(out)
    expected_vec = np.array([expected_normal, 0.0, 0.0])
    np.testing.assert_allclose(
        out_np,
        expected_vec,
        rtol=1e-12,
        atol=1e-15,
        err_msg=(
            "vectorized mixing_ratio_air result does not match "
            "[normal, 0.0, 0.0] under the simplified form; check that "
            "the divide-safe denominator and outer ``> 0.0`` guard are "
            "consistent."
        ),
    )
    assert out_np[2] >= 0.0, (
        "C4 guard regression in vectorized form: e>P cell produced a "
        "negative value."
    )


# ---------- m4: NaN-safe ``water_specific_heat`` ----------


def test_water_specific_heat_propagates_nan(
    temperature_module: Temperature,
) -> None:
    """Pin that ``water_specific_heat(NaN)`` returns NaN (m4 fix).

    Without the fix, ``np.select`` evaluates every comparison against
    NaN as False and silently returns the ``default=4178.0`` value.
    The NaN-safe wrapper must override that and let NaN propagate so
    missing-data defects surface visibly.
    """
    result_scalar = temperature_module.water_specific_heat(np.nan)
    assert np.isnan(result_scalar), (
        "m4 regression: water_specific_heat(NaN) returned "
        f"{result_scalar!r}; expected NaN. A finite return value "
        "(typically 4178.0) indicates the NaN guard has been removed."
    )

    # Vectorized case: NaN cells must remain NaN, finite cells must
    # return the calibrated lookup value.
    temperatures = xr.DataArray(np.array([np.nan, 12.0, 22.0]))
    result_vec = temperature_module.water_specific_heat(temperatures)
    result_np = np.asarray(result_vec)
    assert np.isnan(result_np[0]), (
        "m4 vectorized regression: NaN cell returned a finite value."
    )
    assert result_np[1] == pytest.approx(4186.0)
    assert result_np[2] == pytest.approx(4180.0)


# ---------- m6: ``flux_sediment`` returns DataArray when disabled ----------


def test_flux_sediment_disabled_returns_dataarray_zeros(
    temperature_module: Temperature,
) -> None:
    """Pin that ``flux_sediment`` with ``use_sediment_temperature=False``
    returns an ``xr.DataArray`` of zeros with the same shape as
    ``water_temperature`` (m6 fix).

    Previously the disabled branch returned the Python scalar ``0.0``,
    which broadcasts but produces inconsistent dtypes in dask graphs
    and frustrates static type stubs.
    """
    disabled = Temperature(0.3, 1.5, 3.0, use_sediment_temperature=False)
    water_temperature = xr.DataArray(np.array([20.0, 21.0, 22.0]))
    sediment_temperature = xr.DataArray(np.array([15.0, 15.0, 15.0]))
    sediment_thickness = xr.DataArray(np.array([0.1, 0.1, 0.1]))

    result = disabled.flux_sediment(
        water_temperature=water_temperature,
        sediment_temperature=sediment_temperature,
        sediment_thickness=sediment_thickness,
    )

    assert isinstance(result, xr.DataArray), (
        f"m6 regression: flux_sediment(disabled) returned "
        f"{type(result).__name__}, expected xr.DataArray. The disabled "
        "branch must return ``xr.zeros_like(water_temperature)`` rather "
        "than the Python float 0.0."
    )
    assert result.shape == water_temperature.shape, (
        f"m6 shape regression: got {result.shape}, expected "
        f"{water_temperature.shape}."
    )
    np.testing.assert_array_equal(
        np.asarray(result),
        np.zeros_like(np.asarray(water_temperature)),
        err_msg="m6 value regression: disabled branch must return zeros.",
    )


# ---------- m19: per-process volume>0 guard removed ----------


def test_temperature_change_no_longer_self_guards_zero_volume() -> None:
    """Pin that ``Temperature.run`` no longer applies its own
    ``volume > 0`` guard.

    The pre-condition is that ``temperature_change`` itself produces
    a non-finite (NaN or inf) delta when ``volume == 0`` (division by
    zero in the ``flux_net * surface_area * dt / (V * rho * cp)``
    update). The orchestration-layer wet-mask in
    ``Model.__apply_wet_mask`` handles the dry-cell masking at the
    registry level instead. With the per-process guard restored, the
    middle delta would be 0.0 rather than NaN/inf, so this test is the
    direct regression guard for the m19 / C5 cleanup.
    """
    process = Temperature(0.3, 1.5, 3.0, use_sediment_temperature=False)

    water_temperature = xr.DataArray(np.array([20.0, 20.0, 20.0]))
    surface_area = xr.DataArray(np.array([10.0, 10.0, 10.0]))
    volume = xr.DataArray(np.array([100.0, 0.0, 100.0]))
    cloudiness = xr.DataArray(np.array([0.0, 0.0, 0.0]))
    air_temperature = xr.DataArray(np.array([25.0, 25.0, 25.0]))
    solar_flux = xr.DataArray(np.array([200.0, 200.0, 200.0]))
    wind_speed = xr.DataArray(np.array([3.0, 3.0, 3.0]))
    sediment_temperature = xr.DataArray(np.array([15.0, 15.0, 15.0]))
    sediment_thickness = xr.DataArray(np.array([0.1, 0.1, 0.1]))
    atmospheric_pressure = xr.DataArray(np.array([1013.0, 1013.0, 1013.0]))
    atmospheric_vapor_pressure = xr.DataArray(np.array([20.0, 20.0, 20.0]))

    # Suppress the divide-by-zero warning emitted by the unguarded
    # update; it is the expected and documented Phase 3 behavior.
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = process.temperature_change(
            water_temperature=water_temperature,
            surface_area=surface_area,
            volume=volume,
            cloudiness=cloudiness,
            air_temperature=air_temperature,
            solar_flux=solar_flux,
            wind_speed=wind_speed,
            sediment_temperature=sediment_temperature,
            sediment_thickness=sediment_thickness,
            atmospheric_pressure=atmospheric_pressure,
            atmospheric_vapor_pressure=atmospheric_vapor_pressure,
        )

    delta_np = np.asarray(delta)

    # Outer cells (volume = 100) must produce finite, nonzero deltas.
    assert np.isfinite(delta_np[0]), (
        f"unexpected non-finite delta on wet cell 0: {delta_np[0]!r}."
    )
    assert np.isfinite(delta_np[2]), (
        f"unexpected non-finite delta on wet cell 2: {delta_np[2]!r}."
    )

    # Middle cell (volume = 0) must produce a non-finite delta. With
    # the per-process guard restored, this would be 0.0, which is the
    # specific regression this test guards.
    assert not np.isfinite(delta_np[1]), (
        "m19 / C5 regression: temperature_change returned a finite "
        f"delta ({delta_np[1]!r}) for volume == 0. The per-process "
        "``xr.where(volume > 0, ...)`` guard appears to have been "
        "reintroduced; the dry-cell pre-condition required by "
        "``Model.__apply_wet_mask`` is now broken."
    )
