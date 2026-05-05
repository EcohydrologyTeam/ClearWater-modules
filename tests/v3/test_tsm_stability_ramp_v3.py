"""Tests for v3's thin-water stability fix (depth ramp + rate cap).

This is the v3 port of ``tests/test_tsm_stability_ramp.py``. The v1 test
guards the regularizations baked into ``dTdt_water_c``; this v3 port
guards the same regularizations in
``Temperature.temperature_change`` (see
``src/clearwater_modules_v3/processes/temperature.py``, ~lines 295-376).

Key v1 -> v3 conventions to keep straight:

- v1's ``dTdt_water_c`` takes a precomputed ``q_net`` that has *already*
  been multiplied by ``86400 * dt_days`` (so it is energy-per-substep
  per m^2, J/m^2). v3's ``temperature_change`` instead computes the net
  flux internally as ``flux_net`` (W/m^2) and multiplies by
  ``self.time_step_seconds`` inside the method. The math at the call
  site is therefore different even though the regularization logic is
  identical.

- v3's per-substep delta T base form is
  ``flux_net * SA * time_step_seconds / (V * rho * cp)``.
  For a 5-minute substep ``time_step_seconds = 300`` and the rate cap
  ``dTdt_max_per_hour * dt_hours = 5.0 * (5/60) = 0.4167 K``.

To isolate the regularization logic from the full energy budget, T1-T5
monkeypatch ``Temperature.flux_net`` to return a controlled W/m^2 value.
T6 calls ``temperature_change`` directly with multi-cell inputs, also
with ``flux_net`` patched, to verify per-cell vector behavior of ramp +
cap on a heterogeneous depth distribution.

Design memo: ``design/tsm_stability_thin_water.md``.
"""
from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.temperature import Temperature


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _arr(value, n: int = 1) -> xr.DataArray:
    """Wrap a scalar (or scalar-broadcast) as a 1-D DataArray of length n."""
    return xr.DataArray(np.full(n, value, dtype=float), dims=["nface"])


def _stub_temperature_change_inputs(
    *,
    water_temperature: float = 20.0,
    surface_area: float = 1.0,
    volume: float = 1.0,
    n: int = 1,
) -> dict:
    """Build the kwargs ``temperature_change`` expects.

    Only ``water_temperature``, ``surface_area``, and ``volume`` matter
    once ``flux_net`` is patched to return a constant; the remaining
    fields just need to be present so the dispatch doesn't error.
    """
    return {
        "water_temperature": _arr(water_temperature, n),
        "surface_area": _arr(surface_area, n),
        "volume": _arr(volume, n),
        "cloudiness": _arr(0.1, n),
        "air_temperature": _arr(15.0, n),
        "solar_flux": _arr(0.0, n),
        "wind_speed": _arr(3.0, n),
        "sediment_temperature": _arr(15.0, n),
        "sediment_thickness": _arr(0.1, n),
        "atmospheric_pressure": _arr(1013.0, n),
        "atmospheric_vapor_pressure": _arr(10.0, n),
    }


def _patch_flux_net(temp: Temperature, monkeypatch, flux_w_per_m2):
    """Replace ``flux_net`` with a stub that returns a known W/m^2 value
    broadcast over the input cells. ``flux_w_per_m2`` may be a scalar
    or a 1-D array matching the per-cell shape.
    """

    def _stub(self, **kwargs):
        # Use water_temperature's shape as the broadcast template so the
        # output has matching dims for the downstream xr.where calls.
        wt = kwargs["water_temperature"]
        return xr.zeros_like(wt) + np.asarray(flux_w_per_m2)

    monkeypatch.setattr(Temperature, "flux_net", _stub, raising=True)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def temp_module() -> Temperature:
    """Default Temperature module: 5-min substep, ramp_ref=0.3, cap=5 K/hr."""
    return Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
    )


# ---------------------------------------------------------------------------
# T1 - parity with prior kernel when both regularisations are disabled
# ---------------------------------------------------------------------------


def test_t1_disabled_matches_legacy_arithmetic(monkeypatch):
    """ramp=disabled, cap=disabled -> bit-exact prior formula.

    The base form is ``flux_net * SA * dt_seconds / (V * rho * cp)``. With
    both regularizations off, ``temperature_change`` must reduce to that
    expression exactly.
    """
    temp = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.0,
        dTdt_max_per_hour=float("inf"),
    )

    flux_w = 100.0  # W/m^2
    _patch_flux_net(temp, monkeypatch, flux_w)

    water_temp_c = 20.0
    surface_area = 1.0
    depth = 1.0
    volume = depth * surface_area
    inputs = _stub_temperature_change_inputs(
        water_temperature=water_temp_c,
        surface_area=surface_area,
        volume=volume,
    )

    out = temp.temperature_change(**inputs).values[0]

    rho = float(temp.water_density(water_temp_c))
    cp = float(temp.water_specific_heat(water_temp_c))
    dt_seconds = temp.time_step_seconds
    expected = flux_w * surface_area * dt_seconds / (volume * rho * cp)

    np.testing.assert_allclose(out, expected, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# T2 - ramp activates on shallow cell
# ---------------------------------------------------------------------------


def test_t2_ramp_activates_on_shallow_cell(monkeypatch):
    """Shallow cell (depth < D_ref): ramp = depth/D_ref damps the flux."""
    temp = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.3,
        dTdt_max_per_hour=float("inf"),  # isolate ramp from cap
    )

    flux_w = 600.0
    _patch_flux_net(temp, monkeypatch, flux_w)

    water_temp_c = 20.0
    surface_area = 1.0
    depth = 0.05
    volume = depth * surface_area
    inputs = _stub_temperature_change_inputs(
        water_temperature=water_temp_c,
        surface_area=surface_area,
        volume=volume,
    )

    out = temp.temperature_change(**inputs).values[0]

    rho = float(temp.water_density(water_temp_c))
    cp = float(temp.water_specific_heat(water_temp_c))
    dt_seconds = temp.time_step_seconds
    unramped = flux_w * surface_area * dt_seconds / (volume * rho * cp)
    expected_ramp = depth / 0.3
    expected = unramped * expected_ramp

    np.testing.assert_allclose(out, expected, rtol=1e-9, atol=0.0)
    # Sanity: ramped result is strictly less than half of the unramped.
    assert abs(out) < abs(unramped) * 0.5


# ---------------------------------------------------------------------------
# T3 - ramp inactive on deep cell
# ---------------------------------------------------------------------------


def test_t3_ramp_inactive_on_deep_cell(monkeypatch):
    """Deep cell (depth > D_ref): ramp == 1.0; result equals legacy form."""
    temp = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.3,
        dTdt_max_per_hour=float("inf"),
    )

    flux_w = 300.0
    _patch_flux_net(temp, monkeypatch, flux_w)

    water_temp_c = 20.0
    surface_area = 1.0
    depth = 1.0
    volume = depth * surface_area
    inputs = _stub_temperature_change_inputs(
        water_temperature=water_temp_c,
        surface_area=surface_area,
        volume=volume,
    )

    out = temp.temperature_change(**inputs).values[0]

    rho = float(temp.water_density(water_temp_c))
    cp = float(temp.water_specific_heat(water_temp_c))
    dt_seconds = temp.time_step_seconds
    expected = flux_w * surface_area * dt_seconds / (volume * rho * cp)

    np.testing.assert_allclose(out, expected, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# T4 - cap activates on extreme flux
# ---------------------------------------------------------------------------


def test_t4_cap_activates_on_extreme_flux(monkeypatch):
    """Cap clips |delta_T| to dTdt_max_per_hour * dt_hours.

    For a 5-min substep and cap=5 K/hr, the cap is 5/12 ~= 0.4167 K.
    """
    temp = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.0,  # isolate cap from ramp
        dTdt_max_per_hour=5.0,
    )

    flux_w = 50_000.0  # unphysical, force the cap
    _patch_flux_net(temp, monkeypatch, flux_w)

    water_temp_c = 20.0
    surface_area = 1.0
    depth = 1.0
    volume = depth * surface_area
    inputs = _stub_temperature_change_inputs(
        water_temperature=water_temp_c,
        surface_area=surface_area,
        volume=volume,
    )

    out = temp.temperature_change(**inputs).values[0]

    dt_hours = temp.time_step_seconds / 3600.0
    cap_value = 5.0 * dt_hours
    np.testing.assert_allclose(cap_value, 5.0 / 12.0, rtol=1e-12)
    assert abs(out) <= cap_value + 1e-12
    # And it should saturate, not be smaller:
    np.testing.assert_allclose(out, cap_value, rtol=1e-12, atol=0.0)


def test_t4b_cap_activates_on_extreme_negative_flux(monkeypatch):
    """Cap is symmetric: clips |delta_T| for large negative flux too."""
    temp = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.0,
        dTdt_max_per_hour=5.0,
    )

    flux_w = -50_000.0
    _patch_flux_net(temp, monkeypatch, flux_w)

    inputs = _stub_temperature_change_inputs(
        water_temperature=20.0,
        surface_area=1.0,
        volume=1.0,
    )
    out = temp.temperature_change(**inputs).values[0]

    dt_hours = temp.time_step_seconds / 3600.0
    cap_value = 5.0 * dt_hours
    np.testing.assert_allclose(out, -cap_value, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# T5 - cap disabled (np.inf) passes extreme flux through unchanged
# ---------------------------------------------------------------------------


def test_t5_cap_disabled_passes_extreme_flux(monkeypatch):
    """dTdt_max_per_hour = +inf -> no clipping, raw value preserved."""
    temp = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.0,
        dTdt_max_per_hour=float("inf"),
    )

    flux_w = 50_000.0
    _patch_flux_net(temp, monkeypatch, flux_w)

    water_temp_c = 20.0
    surface_area = 1.0
    depth = 1.0
    volume = depth * surface_area
    inputs = _stub_temperature_change_inputs(
        water_temperature=water_temp_c,
        surface_area=surface_area,
        volume=volume,
    )

    out = temp.temperature_change(**inputs).values[0]

    rho = float(temp.water_density(water_temp_c))
    cp = float(temp.water_specific_heat(water_temp_c))
    dt_seconds = temp.time_step_seconds
    expected_raw = flux_w * surface_area * dt_seconds / (volume * rho * cp)
    # Should be a non-trivial number (>> the 0.4167 K cap)
    assert expected_raw > 1.0
    np.testing.assert_allclose(out, expected_raw, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# T6 - multi-cell round trip with mixed depths
# ---------------------------------------------------------------------------


def test_t6_mixed_depths_round_trip(monkeypatch):
    """Multi-cell DataArray with mixed depths: ramped/unramped cells coexist.

    Build a 5-cell input. Cell 0 is shallow (depth=0.05 m); cells 1-4 are
    deep (depth=1.0 m). With ``q_net_depth_ramp_ref=0.3``, cell 0's ramp
    factor is 0.05/0.3 and cells 1-4's ramp factor is 1.0. With
    ``dTdt_max_per_hour=float('inf')`` the cap is a no-op; with the
    default cap, it should saturate the shallow cell only when the
    flux is extreme enough.
    """
    temp_ramp = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.3,
        dTdt_max_per_hour=float("inf"),
    )
    temp_legacy = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.0,
        dTdt_max_per_hour=float("inf"),
    )

    flux_w = 600.0  # moderate flux, well below cap on deep cells
    # Patch *both* instances by patching the class once; both share the
    # same underlying method.
    _patch_flux_net(temp_ramp, monkeypatch, flux_w)

    depths = np.array([0.05, 1.0, 1.0, 1.0, 1.0])
    surface_area = np.ones_like(depths)
    volume = depths * surface_area
    water_temp_c = np.full_like(depths, 20.0)

    inputs_ramp = {
        "water_temperature": xr.DataArray(water_temp_c, dims=["nface"]),
        "surface_area": xr.DataArray(surface_area, dims=["nface"]),
        "volume": xr.DataArray(volume, dims=["nface"]),
        "cloudiness": xr.DataArray(np.full(5, 0.1), dims=["nface"]),
        "air_temperature": xr.DataArray(np.full(5, 15.0), dims=["nface"]),
        "solar_flux": xr.DataArray(np.full(5, 0.0), dims=["nface"]),
        "wind_speed": xr.DataArray(np.full(5, 3.0), dims=["nface"]),
        "sediment_temperature": xr.DataArray(np.full(5, 15.0), dims=["nface"]),
        "sediment_thickness": xr.DataArray(np.full(5, 0.1), dims=["nface"]),
        "atmospheric_pressure": xr.DataArray(np.full(5, 1013.0), dims=["nface"]),
        "atmospheric_vapor_pressure": xr.DataArray(np.full(5, 10.0), dims=["nface"]),
    }

    out_ramp = np.asarray(temp_ramp.temperature_change(**inputs_ramp).values)
    out_legacy = np.asarray(temp_legacy.temperature_change(**inputs_ramp).values)

    rho = float(temp_ramp.water_density(20.0))
    cp = float(temp_ramp.water_specific_heat(20.0))
    dt_seconds = temp_ramp.time_step_seconds

    # Deep cells (1-4): ramp factor is 1.0, so ramp and legacy agree exactly.
    expected_deep = flux_w * 1.0 * dt_seconds / (1.0 * rho * cp)
    np.testing.assert_allclose(out_ramp[1:], expected_deep, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(
        out_ramp[1:], out_legacy[1:], rtol=1e-12, atol=0.0,
        err_msg="deep cells must be bit-identical with vs without ramp",
    )

    # Shallow cell (0): ramp factor is 0.05/0.3.
    unramped_shallow = flux_w * 1.0 * dt_seconds / (0.05 * rho * cp)
    expected_shallow = unramped_shallow * (0.05 / 0.3)
    np.testing.assert_allclose(out_ramp[0], expected_shallow, rtol=1e-9, atol=0.0)
    np.testing.assert_allclose(out_legacy[0], unramped_shallow, rtol=1e-12, atol=0.0)

    # Ramp must shrink the shallow-cell delta T meaningfully.
    assert abs(out_ramp[0]) < 0.5 * abs(out_legacy[0])

    # No NaNs or Infs anywhere
    assert np.isfinite(out_ramp).all()
    assert np.isfinite(out_legacy).all()


def test_t6b_mixed_depths_with_cap(monkeypatch):
    """Cap applies per-cell on multi-cell inputs.

    With an extreme flux and the default cap, every cell's |delta_T|
    saturates at the cap regardless of depth or ramp factor.
    """
    temp = Temperature(
        0.3,
        1.5,
        3.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=False,
        q_net_depth_ramp_ref=0.3,
        dTdt_max_per_hour=5.0,
    )

    flux_w = 50_000.0  # extreme: even on deep cells with ramp=1, saturates
    _patch_flux_net(temp, monkeypatch, flux_w)

    depths = np.array([0.05, 1.0, 1.0, 1.0, 1.0])
    surface_area = np.ones_like(depths)
    volume = depths * surface_area
    water_temp_c = np.full_like(depths, 20.0)

    inputs = {
        "water_temperature": xr.DataArray(water_temp_c, dims=["nface"]),
        "surface_area": xr.DataArray(surface_area, dims=["nface"]),
        "volume": xr.DataArray(volume, dims=["nface"]),
        "cloudiness": xr.DataArray(np.full(5, 0.1), dims=["nface"]),
        "air_temperature": xr.DataArray(np.full(5, 15.0), dims=["nface"]),
        "solar_flux": xr.DataArray(np.full(5, 0.0), dims=["nface"]),
        "wind_speed": xr.DataArray(np.full(5, 3.0), dims=["nface"]),
        "sediment_temperature": xr.DataArray(np.full(5, 15.0), dims=["nface"]),
        "sediment_thickness": xr.DataArray(np.full(5, 0.1), dims=["nface"]),
        "atmospheric_pressure": xr.DataArray(np.full(5, 1013.0), dims=["nface"]),
        "atmospheric_vapor_pressure": xr.DataArray(np.full(5, 10.0), dims=["nface"]),
    }

    out = np.asarray(temp.temperature_change(**inputs).values)

    cap_value = 5.0 * (temp.time_step_seconds / 3600.0)
    # Every cell saturates at the cap (positive flux -> +cap).
    np.testing.assert_allclose(out, cap_value, rtol=1e-12, atol=0.0)
