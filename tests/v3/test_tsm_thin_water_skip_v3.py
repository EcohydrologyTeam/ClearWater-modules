"""Tests for v3 TSM's opt-in thin-water heat-flux skip.

Guards the ``q_net_depth_skip_threshold`` parameter on
``clearwater_modules_v3.processes.temperature.Temperature`` (added
2026-05-27 per the Corvallis-Santiam Sept 2008 newly-wet-cell
investigation). The skip is the third, hardest thin-water
regularization, layered on top of the depth ramp
(``q_net_depth_ramp_ref``) and the rate cap (``dTdt_max_per_hour``)
exercised by ``tests/v3/test_tsm_stability_ramp_v3.py``.

Behavior under test (see
``design/clearwater_modules_v3_thin_water_skip.md``):

- ``q_net_depth_skip_threshold = 0.0`` (the default) DISABLES the skip
  and preserves byte-identity with prior runs — even an extremely thin
  cell is not zeroed.
- A positive threshold zeroes the kinetics-side delta AND the ``ramp``
  factor for cells with ``depth < threshold`` (strict). Zeroing ``ramp``
  propagates the skip to the sediment-side delta in ``Temperature.run``
  (``ramp * clip_ratio`` scaling, audit finding F2), so a skipped cell
  freezes BOTH the water and sediment reservoirs on the kinetics side —
  trivially energy-conservative (0 + 0 = 0). The per-component flux
  diagnostics are NOT zeroed (retained for auditing).

The kinetics-side unit tests patch ``flux_components`` to return a
controlled W/m^2 value (same approach as the ramp tests) so the skip
logic is isolated from the full energy budget. The ``run``-driven
sediment-propagation tests use a minimal registry stub.
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.temperature import Temperature


# ---------------------------------------------------------------------------
# Helpers (mirror tests/v3/test_tsm_stability_ramp_v3.py +
# tests/v3/test_tsm_sediment_v3.py)
# ---------------------------------------------------------------------------


def _arr(value, n: int = 1) -> xr.DataArray:
    """Wrap a scalar (or scalar-broadcast) as a 1-D DataArray of length n."""
    return xr.DataArray(np.full(n, value, dtype=float), dims=["nface"])


def _scalar(da) -> float:
    arr = np.asarray(da)
    return float(arr.reshape(-1)[0])


def _stub_temperature_change_inputs(
    *,
    water_temperature=20.0,
    surface_area=1.0,
    volume=1.0,
    n: int = 1,
) -> dict:
    """Build the kwargs ``temperature_change`` /
    ``_temperature_change_with_factors`` expect. Once ``flux_components``
    is patched, only ``water_temperature``, ``surface_area``, and
    ``volume`` (-> depth) matter; the rest just need to be present."""
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


def _patch_flux_components(monkeypatch, flux_w_per_m2):
    """Replace ``flux_components`` (and ``flux_net``) with stubs returning
    a known W/m^2 value broadcast over the input cells."""

    def _stub_net(self, **kwargs):
        wt = kwargs["water_temperature"]
        return xr.zeros_like(wt) + np.asarray(flux_w_per_m2)

    def _stub_components(self, **kwargs):
        wt = kwargs["water_temperature"]
        zeros = xr.zeros_like(wt)
        net = zeros + np.asarray(flux_w_per_m2)
        return {
            "q_sensible": zeros,
            "q_latent": zeros,
            "q_longwave_up": zeros,
            "q_longwave_down": zeros,
            "q_solar": zeros,
            "q_sediment": zeros,
            "q_net": net,
        }

    monkeypatch.setattr(Temperature, "flux_net", _stub_net, raising=True)
    monkeypatch.setattr(
        Temperature, "flux_components", _stub_components, raising=True
    )


def _make(
    *,
    q_net_depth_skip_threshold=0.0,
    q_net_depth_ramp_ref=0.3,
    dTdt_max_per_hour=float("inf"),
    use_sediment_temperature=False,
    evolve_sediment_temperature=True,
) -> Temperature:
    return Temperature(
        wind_a=0.3,
        wind_b=1.5,
        wind_c=2.0,
        time_step=timedelta(minutes=5),
        use_sediment_temperature=use_sediment_temperature,
        evolve_sediment_temperature=evolve_sediment_temperature,
        q_net_depth_ramp_ref=q_net_depth_ramp_ref,
        q_net_depth_skip_threshold=q_net_depth_skip_threshold,
        dTdt_max_per_hour=dTdt_max_per_hour,
    )


def _ramped_delta(temp, *, flux_w, surface_area, volume, water_temp_c, ramp):
    """Analytic ramp(+no-cap) kinetics delta for one cell."""
    rho = float(temp.water_density(water_temp_c))
    cp = float(temp.water_specific_heat(water_temp_c))
    dt_seconds = temp.time_step_seconds
    unramped = flux_w * surface_area * dt_seconds / (volume * rho * cp)
    return unramped * ramp


# ---------------------------------------------------------------------------
# S1 - default off: skip never fires, byte-identity with ramp+cap path
# ---------------------------------------------------------------------------


def test_s1_default_off_does_not_zero_thin_cell(monkeypatch):
    """``q_net_depth_skip_threshold = 0.0`` (default) leaves a very thin
    cell untouched by the skip — output equals the ramp result, nonzero."""
    temp = _make(q_net_depth_skip_threshold=0.0, q_net_depth_ramp_ref=0.3)
    flux_w = 200.0
    _patch_flux_components(monkeypatch, flux_w)

    surface_area, depth = 1.0, 0.02  # depth << any plausible threshold
    volume = depth * surface_area
    inputs = _stub_temperature_change_inputs(
        water_temperature=20.0, surface_area=surface_area, volume=volume
    )

    out = _scalar(temp.temperature_change(**inputs))
    expected = _ramped_delta(
        temp, flux_w=flux_w, surface_area=surface_area, volume=volume,
        water_temp_c=20.0, ramp=depth / 0.3,
    )
    np.testing.assert_allclose(out, expected, rtol=1e-12, atol=0.0)
    assert out != 0.0, "default (skip off) must not zero a thin cell"


# ---------------------------------------------------------------------------
# S2 - skip fires below threshold: delta and ramp both zeroed
# ---------------------------------------------------------------------------


def test_s2_skip_fires_below_threshold(monkeypatch):
    """depth < threshold -> delta == 0 and the returned ramp == 0."""
    temp = _make(q_net_depth_skip_threshold=0.05, q_net_depth_ramp_ref=0.3)
    _patch_flux_components(monkeypatch, 200.0)

    surface_area, depth = 1.0, 0.02  # 0.02 < 0.05 -> skip
    inputs = _stub_temperature_change_inputs(
        water_temperature=20.0, surface_area=surface_area,
        volume=depth * surface_area,
    )

    delta, ramp, _clip, _components = temp._temperature_change_with_factors(
        **inputs
    )
    assert _scalar(delta) == 0.0
    assert _scalar(ramp) == 0.0


# ---------------------------------------------------------------------------
# S3 - strict boundary: a cell exactly at the threshold is NOT skipped
# ---------------------------------------------------------------------------


def test_s3_boundary_is_strict(monkeypatch):
    """``depth < threshold`` is strict: depth == threshold -> not skipped."""
    temp = _make(q_net_depth_skip_threshold=0.05, q_net_depth_ramp_ref=0.3)
    flux_w = 200.0
    _patch_flux_components(monkeypatch, flux_w)

    surface_area, depth = 1.0, 0.05  # exactly at threshold
    volume = depth * surface_area
    inputs = _stub_temperature_change_inputs(
        water_temperature=20.0, surface_area=surface_area, volume=volume
    )

    delta, ramp, _clip, _components = temp._temperature_change_with_factors(
        **inputs
    )
    # Not skipped: ramp is the depth ramp (0.05 / 0.3), delta nonzero.
    np.testing.assert_allclose(_scalar(ramp), 0.05 / 0.3, rtol=1e-12)
    expected = _ramped_delta(
        temp, flux_w=flux_w, surface_area=surface_area, volume=volume,
        water_temp_c=20.0, ramp=0.05 / 0.3,
    )
    np.testing.assert_allclose(_scalar(delta), expected, rtol=1e-12)
    assert _scalar(delta) != 0.0


# ---------------------------------------------------------------------------
# S4 - above threshold: bit-identical with the skip on vs off
# ---------------------------------------------------------------------------


def test_s4_above_threshold_bit_identical_on_vs_off(monkeypatch):
    """A deep cell (depth > threshold) is unaffected by the skip."""
    flux_w = 300.0
    _patch_flux_components(monkeypatch, flux_w)

    surface_area, depth = 1.0, 1.0
    volume = depth * surface_area
    inputs = _stub_temperature_change_inputs(
        water_temperature=20.0, surface_area=surface_area, volume=volume
    )

    temp_off = _make(q_net_depth_skip_threshold=0.0, q_net_depth_ramp_ref=0.3)
    temp_on = _make(q_net_depth_skip_threshold=0.05, q_net_depth_ramp_ref=0.3)

    out_off = _scalar(temp_off.temperature_change(**inputs))
    out_on = _scalar(temp_on.temperature_change(**inputs))
    np.testing.assert_allclose(out_on, out_off, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# S5 - multi-cell mixed depths: per-cell skip / ramp / passthrough
# ---------------------------------------------------------------------------


def test_s5_mixed_depths_per_cell(monkeypatch):
    """4-cell mesh, threshold=0.05, ramp_ref=0.3, cap disabled:

    - cell 0 (0.02 m): skipped       -> delta 0, ramp 0
    - cell 1 (0.05 m): at threshold  -> ramped (0.05/0.3), not skipped
    - cell 2 (0.50 m): deep          -> ramp 1.0
    - cell 3 (1.00 m): deep          -> ramp 1.0
    """
    temp = _make(q_net_depth_skip_threshold=0.05, q_net_depth_ramp_ref=0.3)
    flux_w = 200.0
    _patch_flux_components(monkeypatch, flux_w)

    depths = np.array([0.02, 0.05, 0.50, 1.00])
    surface_area = np.ones_like(depths)
    volume = depths * surface_area
    inputs = {
        "water_temperature": xr.DataArray(np.full(4, 20.0), dims=["nface"]),
        "surface_area": xr.DataArray(surface_area, dims=["nface"]),
        "volume": xr.DataArray(volume, dims=["nface"]),
        "cloudiness": xr.DataArray(np.full(4, 0.1), dims=["nface"]),
        "air_temperature": xr.DataArray(np.full(4, 15.0), dims=["nface"]),
        "solar_flux": xr.DataArray(np.full(4, 0.0), dims=["nface"]),
        "wind_speed": xr.DataArray(np.full(4, 3.0), dims=["nface"]),
        "sediment_temperature": xr.DataArray(np.full(4, 15.0), dims=["nface"]),
        "sediment_thickness": xr.DataArray(np.full(4, 0.1), dims=["nface"]),
        "atmospheric_pressure": xr.DataArray(np.full(4, 1013.0), dims=["nface"]),
        "atmospheric_vapor_pressure": xr.DataArray(np.full(4, 10.0), dims=["nface"]),
    }

    delta, ramp, _clip, _components = temp._temperature_change_with_factors(
        **inputs
    )
    delta = np.asarray(delta.values)
    ramp = np.asarray(ramp.values)

    rho = float(temp.water_density(20.0))
    cp = float(temp.water_specific_heat(20.0))
    dt_seconds = temp.time_step_seconds
    expected_ramp = np.array([0.0, 0.05 / 0.3, 1.0, 1.0])
    np.testing.assert_allclose(ramp, expected_ramp, rtol=1e-12, atol=0.0)

    # cell 0 skipped
    assert delta[0] == 0.0
    # cells 1-3: unramped * ramp
    for i in (1, 2, 3):
        unramped = flux_w * surface_area[i] * dt_seconds / (volume[i] * rho * cp)
        np.testing.assert_allclose(
            delta[i], unramped * expected_ramp[i], rtol=1e-12, atol=0.0
        )
    assert np.isfinite(delta).all()


# ---------------------------------------------------------------------------
# S6 - ramp disabled: scalar ramp promoted to array, zeroed on skip
# ---------------------------------------------------------------------------


def test_s6_ramp_disabled_scalar_promoted_and_zeroed(monkeypatch):
    """With ``q_net_depth_ramp_ref = 0`` the ramp is the scalar 1.0; the
    skip must promote it to a per-cell array and zero the thin cell."""
    temp = _make(q_net_depth_skip_threshold=0.05, q_net_depth_ramp_ref=0.0)
    flux_w = 200.0
    _patch_flux_components(monkeypatch, flux_w)

    depths = np.array([0.02, 1.00])  # cell 0 skipped, cell 1 passthrough
    surface_area = np.ones_like(depths)
    volume = depths * surface_area
    inputs = {
        "water_temperature": xr.DataArray(np.full(2, 20.0), dims=["nface"]),
        "surface_area": xr.DataArray(surface_area, dims=["nface"]),
        "volume": xr.DataArray(volume, dims=["nface"]),
        "cloudiness": xr.DataArray(np.full(2, 0.1), dims=["nface"]),
        "air_temperature": xr.DataArray(np.full(2, 15.0), dims=["nface"]),
        "solar_flux": xr.DataArray(np.full(2, 0.0), dims=["nface"]),
        "wind_speed": xr.DataArray(np.full(2, 3.0), dims=["nface"]),
        "sediment_temperature": xr.DataArray(np.full(2, 15.0), dims=["nface"]),
        "sediment_thickness": xr.DataArray(np.full(2, 0.1), dims=["nface"]),
        "atmospheric_pressure": xr.DataArray(np.full(2, 1013.0), dims=["nface"]),
        "atmospheric_vapor_pressure": xr.DataArray(np.full(2, 10.0), dims=["nface"]),
    }

    delta, ramp, _clip, _components = temp._temperature_change_with_factors(
        **inputs
    )
    # ramp must now be a DataArray (promoted from the scalar 1.0).
    assert isinstance(ramp, xr.DataArray)
    np.testing.assert_allclose(
        np.asarray(ramp.values), np.array([0.0, 1.0]), rtol=1e-12, atol=0.0
    )
    delta = np.asarray(delta.values)
    assert delta[0] == 0.0
    rho = float(temp.water_density(20.0))
    cp = float(temp.water_specific_heat(20.0))
    dt_seconds = temp.time_step_seconds
    unramped = flux_w * 1.0 * dt_seconds / (1.0 * rho * cp)  # ramp disabled -> 1.0
    np.testing.assert_allclose(delta[1], unramped, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# S7 - constructor validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad", [-0.1, float("inf"), float("nan")])
def test_s7_rejects_invalid_threshold(bad):
    """Negative, +inf, and NaN thresholds are rejected at construction."""
    with pytest.raises(ValueError, match="q_net_depth_skip_threshold"):
        _make(q_net_depth_skip_threshold=bad)


@pytest.mark.parametrize("good", [0.0, 0.05, 1.0])
def test_s7b_accepts_valid_threshold(good):
    """0.0 (disable) and positive finite values construct without error."""
    temp = _make(q_net_depth_skip_threshold=good)
    assert temp.q_net_depth_skip_threshold == good


# ---------------------------------------------------------------------------
# S8 - run(): skip propagates to the sediment side (energy pair-cancellation)
# ---------------------------------------------------------------------------


class _StubRegistry:
    """Minimal VariableRegistry stand-in (mirrors test_tsm_sediment_v3.py)."""

    def __init__(self, initial: dict) -> None:
        self._data = dict(initial)

    def get_at_time(self, name, time):
        if name not in self._data:
            raise KeyError(name)
        return self._data[name]

    def set_at_time(self, name, time, value) -> None:
        self._data[name] = value

    def get(self, name):
        return self._data[name]

    def __contains__(self, name) -> bool:
        return name in self._data


def _registry_dict(*, T_water, T_sed, surface_area, volume, h2=0.1, extra=None):
    d = {
        "water_temperature": _arr(T_water),
        "wetted_surface_area": _arr(surface_area),
        "volume": _arr(volume),
        "cloudiness": _arr(0.0),
        "air_temperature": _arr(T_water),
        "solar_radiation": _arr(500.0),  # nonzero so skip-off gives nonzero dT
        "wind_speed": _arr(2.0),
        "atmospheric_pressure": _arr(1013.0),
        "atmospheric_vapor_pressure": _arr(12.0),
        "sediment_temperature": _arr(T_sed),
        "sediment_thickness": _arr(h2),
    }
    if extra:
        d.update(extra)
    return d


def _run_capture_deltas(temp, *, T_water, T_sed, surface_area, volume, extra=None):
    registry = _StubRegistry(
        _registry_dict(
            T_water=T_water, T_sed=T_sed,
            surface_area=surface_area, volume=volume, extra=extra,
        )
    )
    temp._Temperature__skip_first_time_step = False  # type: ignore[attr-defined]
    temp.run(datetime(2026, 1, 1, 0, 0, 0), registry)  # type: ignore[arg-type]
    dT_water = _scalar(registry.get("water_temperature")) - T_water
    dT_sed = _scalar(registry.get("sediment_temperature")) - T_sed
    return dT_water, dT_sed, registry


def test_s8_skip_freezes_both_water_and_sediment():
    """A skipped thin cell freezes BOTH reservoirs (water and sediment
    kinetics deltas exactly 0) -> trivial energy pair-cancellation.

    Sanity: with the skip OFF, the same thin cell exchanges heat on both
    sides (nonzero deltas), so the fixture genuinely exercises a cell the
    skip must suppress.
    """
    surface_area = 100.0
    volume = 2.0  # depth = 0.02 m < 0.05 m threshold
    T_water, T_sed = 10.0, 30.0  # large gradient -> nonzero sediment flux

    # Skip ON: both deltas must be exactly zero.
    temp_on = _make(
        q_net_depth_skip_threshold=0.05,
        q_net_depth_ramp_ref=0.3,
        use_sediment_temperature=True,
        evolve_sediment_temperature=True,
    )
    dT_water_on, dT_sed_on, _ = _run_capture_deltas(
        temp_on, T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume,
    )
    assert dT_water_on == 0.0
    assert dT_sed_on == 0.0

    # Skip OFF: the same thin cell changes on both sides.
    temp_off = _make(
        q_net_depth_skip_threshold=0.0,
        q_net_depth_ramp_ref=0.3,
        use_sediment_temperature=True,
        evolve_sediment_temperature=True,
    )
    dT_water_off, dT_sed_off, _ = _run_capture_deltas(
        temp_off, T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume,
    )
    assert dT_water_off != 0.0, "fixture should drive a nonzero water dT"
    assert dT_sed_off != 0.0, "fixture should drive a nonzero sediment dT"


def test_s8b_deep_cell_unaffected_in_run():
    """A deep cell behaves identically with the skip on vs off in run()."""
    surface_area = 100.0
    volume = 1000.0  # depth = 10 m, well above threshold
    T_water, T_sed = 10.0, 30.0

    common = dict(
        q_net_depth_ramp_ref=0.3,
        use_sediment_temperature=True,
        evolve_sediment_temperature=True,
    )
    temp_on = _make(q_net_depth_skip_threshold=0.05, **common)
    temp_off = _make(q_net_depth_skip_threshold=0.0, **common)

    dTw_on, dTs_on, _ = _run_capture_deltas(
        temp_on, T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume,
    )
    dTw_off, dTs_off, _ = _run_capture_deltas(
        temp_off, T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume,
    )
    np.testing.assert_allclose(dTw_on, dTw_off, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(dTs_on, dTs_off, rtol=1e-12, atol=0.0)
    assert dTw_off != 0.0 and dTs_off != 0.0


# ---------------------------------------------------------------------------
# S9 - diagnostics are NOT zeroed on a skipped cell
# ---------------------------------------------------------------------------


def test_s9_diagnostics_not_zeroed_on_skip(monkeypatch):
    """The per-component flux diagnostics written to the registry are the
    unzeroed fluxes even where the cell is skipped (retained for audit)."""
    flux_w = 250.0
    _patch_flux_components(monkeypatch, flux_w)

    temp = _make(
        q_net_depth_skip_threshold=0.05,
        q_net_depth_ramp_ref=0.3,
        use_sediment_temperature=False,
    )
    surface_area = 100.0
    volume = 2.0  # depth 0.02 m -> skipped
    T_water = 20.0

    # Pre-register q_net so run() writes the diagnostic.
    dT_water, _dT_sed, registry = _run_capture_deltas(
        temp, T_water=T_water, T_sed=15.0,
        surface_area=surface_area, volume=volume,
        extra={"q_net": _arr(-999.0)},
    )
    # Cell was skipped: water delta is zero...
    assert dT_water == 0.0
    # ...but the q_net diagnostic holds the (unzeroed) patched flux.
    np.testing.assert_allclose(
        _scalar(registry.get("q_net")), flux_w, rtol=1e-12, atol=0.0
    )
