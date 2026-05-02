"""Tests for the TSM thin-water stability fix (depth ramp + rate cap).

Covers the changes made to ``dTdt_water_c`` in
``src/clearwater_modules/tsm/processes.py``: a depth-based ramp that
attenuates ``q_net`` as ``depth -> 0`` and a per-substep magnitude cap
expressed as a rate (K/hr).

Design memo: ``design/tsm_stability_thin_water.md``.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.tsm.constants import (
    DEFAULT_METEOROLOGICAL,
    DEFAULT_TEMPERATURE,
)
from clearwater_modules.tsm.model import EnergyBudget
from clearwater_modules.tsm.processes import dTdt_water_c


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _arr(value, n: int = 1) -> xr.DataArray:
    """Wrap a scalar (or scalar-broadcast) as a 1-D DataArray of length n."""
    return xr.DataArray(np.full(n, value, dtype=float), dims=['nface'])


# Standard reference: q_net is the *energy-per-substep-per-m^2*, i.e. it has
# already been multiplied by 86400 * dt by the upstream q_net process. So
# for q_net_W = W/m^2 and dt_days = days, the input to dTdt_water_c is
# q_net = q_net_W * 86400 * dt_days.
def _q_net_input(q_net_w_per_m2: float, dt_days: float) -> float:
    return q_net_w_per_m2 * 86400.0 * dt_days


# Physical constants used in expected-value calculations
RHO = 1000.0     # kg/m^3
CP = 4184.0      # J/kg/K
DT_DAYS = 1.0 / 24.0  # hourly substep
SA = 1.0         # m^2 per cell


# ---------------------------------------------------------------------------
# T1 - parity with prior kernel when both regularisations are disabled
# ---------------------------------------------------------------------------

def test_t1_disabled_matches_legacy_arithmetic():
    """ramp=disabled, cap=disabled -> bit-exact prior formula."""
    depth = 1.0
    volume = depth * SA
    q_net_w = 100.0  # W/m^2

    q_net = _arr(_q_net_input(q_net_w, DT_DAYS))
    surface_area = _arr(SA)
    vol = _arr(volume)
    rho = _arr(RHO)
    cp = _arr(CP)

    # Disable both regularisations
    ramp_ref = _arr(0.0)
    rate_max = _arr(np.inf)
    dt = _arr(DT_DAYS)

    out = dTdt_water_c(
        q_net=q_net,
        surface_area=surface_area,
        volume=vol,
        density_water=rho,
        cp_water=cp,
        q_net_depth_ramp_ref=ramp_ref,
        dTdt_max_per_hour=rate_max,
        dt=dt,
    ).values[0]

    # Legacy formula: dTdt = q_net * SA / (V * rho * cp)
    # where q_net = q_net_W * 86400 * dt
    expected = (q_net_w * 86400.0 * DT_DAYS) * SA / (volume * RHO * CP)
    np.testing.assert_allclose(out, expected, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# T2 - ramp activates on shallow cell
# ---------------------------------------------------------------------------

def test_t2_ramp_activates_on_shallow_cell():
    """Shallow cell: ramp = depth/D_ref damps the per-substep delta T."""
    depth = 0.05
    volume = depth * SA
    q_net_w = 600.0  # W/m^2

    q_net = _arr(_q_net_input(q_net_w, DT_DAYS))
    surface_area = _arr(SA)
    vol = _arr(volume)
    rho = _arr(RHO)
    cp = _arr(CP)

    # Defaults: D_ref = 0.3 m, no rate cap (test ramp in isolation)
    ramp_ref = _arr(0.3)
    rate_max = _arr(np.inf)
    dt = _arr(DT_DAYS)

    out = dTdt_water_c(
        q_net=q_net,
        surface_area=surface_area,
        volume=vol,
        density_water=rho,
        cp_water=cp,
        q_net_depth_ramp_ref=ramp_ref,
        dTdt_max_per_hour=rate_max,
        dt=dt,
    ).values[0]

    # Without ramp: dT_unramped ~= 10.32 K (from memo table). With ramp factor
    # 0.05/0.3 ~= 0.1667 the result should be ~0.1667 * 10.32 = ~1.72 K.
    unramped = (q_net_w * 86400.0 * DT_DAYS) * SA / (volume * RHO * CP)
    expected_ramp = depth / 0.3
    expected = unramped * expected_ramp

    np.testing.assert_allclose(out, expected, rtol=1e-9, atol=0.0)
    # Sanity: well below 10 K per substep
    assert abs(out) < 2.0
    # And specifically: well below the 10 K we'd see without the ramp
    assert abs(out) < abs(unramped) * 0.5


# ---------------------------------------------------------------------------
# T3 - ramp inactive on deep cell (parity with pre-fix value)
# ---------------------------------------------------------------------------

def test_t3_ramp_inactive_on_deep_cell():
    """Deep cell (depth > D_ref): ramp == 1.0, result equals legacy arithmetic."""
    depth = 1.0
    volume = depth * SA
    q_net_w = 300.0

    q_net = _arr(_q_net_input(q_net_w, DT_DAYS))
    surface_area = _arr(SA)
    vol = _arr(volume)
    rho = _arr(RHO)
    cp = _arr(CP)

    ramp_ref = _arr(0.3)        # default
    rate_max = _arr(np.inf)     # disable cap
    dt = _arr(DT_DAYS)

    out = dTdt_water_c(
        q_net=q_net,
        surface_area=surface_area,
        volume=vol,
        density_water=rho,
        cp_water=cp,
        q_net_depth_ramp_ref=ramp_ref,
        dTdt_max_per_hour=rate_max,
        dt=dt,
    ).values[0]

    expected = (q_net_w * 86400.0 * DT_DAYS) * SA / (volume * RHO * CP)
    np.testing.assert_allclose(out, expected, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# T4 - cap activates on extreme flux
# ---------------------------------------------------------------------------

def test_t4_cap_activates_on_extreme_flux():
    """Cap clips |dTdt| to dTdt_max_per_hour * dt_hours."""
    depth = 1.0
    volume = depth * SA
    q_net_w = 50000.0  # unphysical, force the cap

    q_net = _arr(_q_net_input(q_net_w, DT_DAYS))
    surface_area = _arr(SA)
    vol = _arr(volume)
    rho = _arr(RHO)
    cp = _arr(CP)

    ramp_ref = _arr(0.0)        # disable ramp to test cap in isolation
    rate_max = _arr(5.0)        # default 5 K/hr
    dt = _arr(DT_DAYS)

    out = dTdt_water_c(
        q_net=q_net,
        surface_area=surface_area,
        volume=vol,
        density_water=rho,
        cp_water=cp,
        q_net_depth_ramp_ref=ramp_ref,
        dTdt_max_per_hour=rate_max,
        dt=dt,
    ).values[0]

    # cap = 5 K/hr * (DT_DAYS * 24 hr/d) = 5 K/hr * 1 hr = 5.0 K
    cap_value = 5.0 * (DT_DAYS * 24.0)
    np.testing.assert_allclose(cap_value, 5.0, rtol=1e-12)
    assert abs(out) <= cap_value + 1e-9
    # And we expect saturation, not less:
    np.testing.assert_allclose(out, cap_value, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# T5 - cap disabled (np.inf) passes through extreme flux unchanged
# ---------------------------------------------------------------------------

def test_t5_cap_disabled_passes_extreme_flux():
    """dTdt_max_per_hour = +inf -> no clipping, raw value preserved."""
    depth = 1.0
    volume = depth * SA
    q_net_w = 50000.0

    q_net = _arr(_q_net_input(q_net_w, DT_DAYS))
    surface_area = _arr(SA)
    vol = _arr(volume)
    rho = _arr(RHO)
    cp = _arr(CP)

    ramp_ref = _arr(0.0)
    rate_max = _arr(np.inf)
    dt = _arr(DT_DAYS)

    out = dTdt_water_c(
        q_net=q_net,
        surface_area=surface_area,
        volume=vol,
        density_water=rho,
        cp_water=cp,
        q_net_depth_ramp_ref=ramp_ref,
        dTdt_max_per_hour=rate_max,
        dt=dt,
    ).values[0]

    expected_raw = (q_net_w * 86400.0 * DT_DAYS) * SA / (volume * RHO * CP)
    # Should be a very large number (~43 K)
    assert expected_raw > 40.0
    np.testing.assert_allclose(out, expected_raw, rtol=1e-12, atol=0.0)


# ---------------------------------------------------------------------------
# T6 - integration with EnergyBudget kernel: shallow vs deep cells
# ---------------------------------------------------------------------------

def _run_energy_budget(
    depths: np.ndarray,
    n_steps: int,
    *,
    q_solar: float,
    air_temp_c: float,
    init_water_temp_c: float,
    temp_param_overrides: dict,
):
    """Helper: build and run a small EnergyBudget; return water_temp trajectory.

    Uses fully-specified meteo and temperature parameter dicts (not built from
    ``DEFAULT_METEOROLOGICAL`` / ``DEFAULT_TEMPERATURE`` snapshots) because
    other tests in the suite mutate those module-level defaults via the
    ``EnergyBudget.__init__`` merge loop. Building from scratch keeps T6's
    forcing deterministic regardless of test ordering.
    """
    n_cells = len(depths)
    surface_area = np.ones(n_cells)
    volume = depths * surface_area
    water_temp_c0 = np.full(n_cells, init_water_temp_c)

    initial_state = {
        'water_temp_c': xr.DataArray(water_temp_c0, dims=['nface']),
        'surface_area': xr.DataArray(surface_area, dims=['nface']),
        'volume': xr.DataArray(volume, dims=['nface']),
    }

    # Fully-specified meteo params (avoid dependence on possibly-mutated
    # DEFAULT_METEOROLOGICAL).
    meteo_params = {
        'air_temp_c': air_temp_c,
        'q_solar': q_solar,
        'sed_temp_c': 5.0,
        'eair_mb': 1.0,
        'pressure_mb': 1013.0,
        'cloudiness': 0.1,
        'wind_speed': 3.0,
        'wind_a': 0.3,
        'wind_b': 1.5,
        'wind_c': 3.0,
        'wind_kh_kw': 1.0,
    }
    # Fully-specified temperature params (avoid dependence on possibly-mutated
    # DEFAULT_TEMPERATURE).
    temp_params = {
        'stefan_boltzmann': 5.67e-8,
        'cp_air': 1005.0,
        'emissivity_water': 0.97,
        'gravity': -9.806,
        'a0': 6984.505294,
        'a1': -188.903931,
        'a2': 2.133357675,
        'a3': -1.288580973e-2,
        'a4': 4.393587233e-5,
        'a5': -8.023923082e-8,
        'a6': 6.136820929e-11,
        'pb': 1600.0,
        'cps': 1673.0,
        'h2': 0.1,
        'alphas': 0.0432,
        'richardson_option': True,
        'dt': 1.0 / 24.0,  # hourly substep in days
        'q_net_depth_ramp_ref': 0.3,
        'dTdt_max_per_hour': 5.0,
    }
    temp_params.update(temp_param_overrides)

    model = EnergyBudget(
        time_steps=n_steps,
        initial_state_values=initial_state,
        meteo_parameters=meteo_params,
        temp_parameters=temp_params,
        use_sed_temp=False,
    )
    for _ in range(n_steps):
        model.increment_timestep()

    water_temp = np.asarray(model.dataset['water_temp_c'].values)
    if water_temp.shape[0] != n_steps + 1:
        water_temp = water_temp.T
    assert water_temp.shape == (n_steps + 1, n_cells)
    return water_temp


def test_t6_round_trip_with_tsm_kernel_mixed_depths():
    """Integration: shallow cell stays bounded, deep cells unchanged.

    Build a 5-cell EnergyBudget. Cell 0 is shallow (depth = 0.05 m); cells 1-4
    are deep (depth = 1.0 m). Run 10 hourly substeps with constant forcing,
    once with ramp + cap at defaults and once with both disabled (legacy kernel).

    Asserts:
      - Deep cells are bit-identical between the two runs (ramp = 1.0 there).
      - Shallow cell variation is strictly smaller with the ramp than without.
      - Both ramp and legacy runs remain finite at this depth/forcing.
    """
    depths = np.array([0.05, 1.0, 1.0, 1.0, 1.0])
    n_steps = 10

    wt_ramp = _run_energy_budget(
        depths,
        n_steps,
        q_solar=400.0,
        air_temp_c=10.0,
        init_water_temp_c=10.0,
        temp_param_overrides={},  # defaults: ramp_ref=0.3, cap=5 K/hr
    )
    wt_legacy = _run_energy_budget(
        depths,
        n_steps,
        q_solar=400.0,
        air_temp_c=10.0,
        init_water_temp_c=10.0,
        temp_param_overrides={
            'q_net_depth_ramp_ref': 0.0,
            'dTdt_max_per_hour': float('inf'),
        },
    )

    range_ramp = wt_ramp.max(axis=0) - wt_ramp.min(axis=0)
    range_legacy = wt_legacy.max(axis=0) - wt_legacy.min(axis=0)

    # Deep cells (1-4): ramp is identically 1.0 (depth = 1 m > D_ref = 0.3 m)
    # and the cap is far above any per-substep delta T here, so the two runs
    # must be bit-exact on the deep cells.
    np.testing.assert_allclose(
        wt_ramp[:, 1:], wt_legacy[:, 1:], rtol=1e-12, atol=0.0,
        err_msg='deep cells should be bit-identical with vs without ramp',
    )

    # Shallow cell: ramped run shows strictly smaller variation than legacy.
    assert range_ramp[0] < range_legacy[0], (
        f'ramp should reduce shallow-cell variation: '
        f'ramped={range_ramp[0]:.3f} K, legacy={range_legacy[0]:.3f} K'
    )
    # Concrete numbers (q_solar=400, depth=0.05, 10 substeps): the ramp factor
    # is 0.05/0.3 ~= 0.167 so the ramped per-substep flux is ~0.167x the
    # legacy flux. The exact temperature reduction is less because q_net's
    # other components (longwave_up, sensible, latent) depend on water_temp_c
    # and self-regulate on the legacy run. Empirically the legacy range here
    # is ~3.3 K and the ramped range is ~2.2 K. Assert the ramp meaningfully
    # damps the shallow range without locking in fragile exact values.
    assert range_ramp[0] < 3.0
    assert range_legacy[0] > 3.0
    # And the ramp must reduce the range by at least 25%
    assert range_ramp[0] < 0.75 * range_legacy[0]

    # Deep cells: with ramp at defaults, ranges are small (~1.2 K diurnal over
    # this short window with this constant forcing).
    assert (range_ramp[1:] < 5.0).all()

    # No NaNs / Infs anywhere
    assert np.isfinite(wt_ramp).all()
    assert np.isfinite(wt_legacy).all()
