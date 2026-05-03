"""Unit tests for SSM deposition probabilities and rates.

References
----------
- Krone, R. B. (1962). Flume Studies of the Transport of Sediment in
  Estuarial Shoaling Processes.
- Gessler, J. (1965). The beginning of bedload movement of mixtures
  investigated as natural armoring in channels.
- SAND2008-5621 sec 5.5; SSM design spec sec 5.5.

Hand-checked limits used below
------------------------------
Krone (cohesive):

    P_d = max(1 - tau / tau_cs, 0)

      tau = 0          -> P_d = 1
      tau = 0.5*tau_cs -> P_d = 0.5
      tau = tau_cs     -> P_d = 0
      tau = 2*tau_cs   -> P_d = 0  (clamped at 0)

Gessler (non-cohesive):

    P_y = (1/0.57) * (tau_cs/tau - 1)
    P_d = 0.5 * erfc(-P_y / sqrt(2))   (= standard normal CDF at P_y)

      tau = tau_cs        -> P_y = 0     -> P_d = 0.5
      tau -> +inf         -> P_y -> -inf -> P_d -> 0
      tau -> 0+           -> P_y -> +inf -> P_d -> 1   (forced for tau<=0)
      tau = 0.5 * tau_cs  -> P_y = 1.754 -> P_d ~ 0.9603
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import xarray as xr
from scipy.special import erfc

from clearwater_modules_v2.processes.sediment import contracts
from clearwater_modules_v2.processes.sediment.classes import (
    SedimentClass,
    SedimentClassRegistry,
)
from clearwater_modules_v2.processes.sediment.deposition import (
    compute_deposition_flux,
    gessler_probability,
    krone_probability,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _tau_da(values) -> xr.DataArray:
    arr = np.asarray(values, dtype="float64").reshape(-1)
    return xr.DataArray(
        arr,
        dims=(contracts.DIM_NFACE,),
        coords={contracts.DIM_NFACE: np.arange(arr.size)},
    )


# ---------------------------------------------------------------------------
# Krone probability
# ---------------------------------------------------------------------------


def test_krone_at_zero_shear_returns_one():
    tau_cs = 0.5
    pd = krone_probability(_tau_da([0.0]), tau_cs)
    assert pd.values[0] == pytest.approx(1.0, rel=1e-3)


def test_krone_at_critical_shear_returns_zero():
    tau_cs = 0.5
    pd = krone_probability(_tau_da([tau_cs]), tau_cs)
    assert pd.values[0] == pytest.approx(0.0, abs=1e-12)


def test_krone_at_half_critical_returns_half():
    tau_cs = 0.5
    pd = krone_probability(_tau_da([0.5 * tau_cs]), tau_cs)
    assert pd.values[0] == pytest.approx(0.5, rel=1e-3)


def test_krone_above_critical_clamped_to_zero():
    tau_cs = 0.5
    pd = krone_probability(_tau_da([2.0 * tau_cs, 10.0 * tau_cs]), tau_cs)
    np.testing.assert_allclose(pd.values, [0.0, 0.0])


def test_krone_vectorized_over_nface():
    tau_cs = 1.0
    tau = _tau_da([0.0, 0.25, 0.5, 0.75, 1.0, 2.0])
    pd = krone_probability(tau, tau_cs)
    expected = np.array([1.0, 0.75, 0.5, 0.25, 0.0, 0.0])
    np.testing.assert_allclose(pd.values, expected, rtol=1e-3, atol=1e-12)


def test_krone_rejects_nonpositive_tau_cs():
    with pytest.raises(ValueError):
        krone_probability(_tau_da([0.1]), 0.0)
    with pytest.raises(ValueError):
        krone_probability(_tau_da([0.1]), -0.5)


# ---------------------------------------------------------------------------
# Gessler probability
# ---------------------------------------------------------------------------


def test_gessler_at_critical_shear_returns_one_half():
    """tau == tau_cs => P_y = 0 => P_d = 0.5 (standard normal CDF at 0)."""
    tau_cs = 0.4
    pd = gessler_probability(_tau_da([tau_cs]), tau_cs)
    assert pd.values[0] == pytest.approx(0.5, rel=1e-3)


def test_gessler_at_large_shear_approaches_finite_asymptote():
    """tau -> +inf => (tau_cs/tau - 1) -> -1, so P_y -> -1/0.57 = -1.7544
    and P_d -> 0.5 * erfc(1.7544/sqrt(2)) ~ 0.0397.

    This is the true mathematical limit of the Gessler form -- there is no
    'P_d -> 0 at infinite shear' because P_y is bounded below by -1/0.57.
    The function still drops monotonically from 0.5 (at tau = tau_cs) toward
    this small positive asymptote, which is what we verify here.
    """
    tau_cs = 0.4
    pd = gessler_probability(_tau_da([1.0e6 * tau_cs]), tau_cs)
    expected_asymptote = 0.5 * erfc((1.0 / 0.57) / math.sqrt(2.0))
    assert pd.values[0] == pytest.approx(expected_asymptote, rel=1e-3)
    # Sanity: well below 0.5 (the value at tau = tau_cs).
    assert pd.values[0] < 0.05


def test_gessler_at_zero_shear_returns_one():
    """tau -> 0 => P_y -> +inf => P_d -> 1 (forced when tau ~ 0)."""
    tau_cs = 0.4
    pd = gessler_probability(_tau_da([0.0]), tau_cs)
    assert pd.values[0] == pytest.approx(1.0, rel=1e-3)


def test_gessler_at_half_critical_matches_normal_cdf():
    """tau = 0.5 * tau_cs => P_y = 1/0.57 = 1.7544 => P_d ~ 0.9603."""
    tau_cs = 0.4
    p_y_expected = 1.0 / 0.57
    p_d_expected = 0.5 * erfc(-p_y_expected / math.sqrt(2.0))
    pd = gessler_probability(_tau_da([0.5 * tau_cs]), tau_cs)
    assert pd.values[0] == pytest.approx(p_d_expected, rel=1e-3)


def test_gessler_returns_in_unit_interval():
    """P_d must always be in [0, 1] across a realistic shear range."""
    tau_cs = 0.4
    tau = _tau_da(np.geomspace(1e-6, 100.0, 50))
    pd = gessler_probability(tau, tau_cs)
    assert (pd.values >= 0.0).all()
    assert (pd.values <= 1.0).all()


def test_gessler_monotonic_decreasing_with_tau():
    """P_d should decrease as tau increases (more shear, less deposition)."""
    tau_cs = 0.4
    tau = _tau_da(np.geomspace(1e-3, 10.0, 20))
    pd = gessler_probability(tau, tau_cs)
    # Allow for tiny numerical wiggle near the asymptotes.
    diffs = np.diff(pd.values)
    assert (diffs <= 1e-12).all()


def test_gessler_rejects_nonpositive_tau_cs():
    with pytest.raises(ValueError):
        gessler_probability(_tau_da([0.1]), 0.0)


# ---------------------------------------------------------------------------
# compute_deposition_flux — registry dispatch and unit conversions
# ---------------------------------------------------------------------------


def _suspended_da(values) -> xr.DataArray:
    """Build (nface, ssm_class) suspended-concentration DataArray."""
    arr = np.asarray(values, dtype="float64")
    nface, nclass = arr.shape
    return xr.DataArray(
        arr,
        dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
        coords={
            contracts.DIM_NFACE: np.arange(nface),
            contracts.DIM_CLASS: np.arange(nclass),
        },
    )


def _depth_da(values) -> xr.DataArray:
    arr = np.asarray(values, dtype="float64").reshape(-1)
    return xr.DataArray(
        arr,
        dims=(contracts.DIM_NFACE,),
        coords={contracts.DIM_NFACE: np.arange(arr.size)},
    )


def test_compute_deposition_flux_mixed_registry_dispatches_correctly():
    """Mixed cohesive (Krone) + non-cohesive (Gessler) classes use the right
    probability per class.

    Build a registry with one fine class (D50=20 um, cohesive => Krone) and
    one coarse class (D50=200 um, non-cohesive => Gessler). Pick a single
    cell with tau = tau_cs for both classes:
      - Krone:   P_d = 0
      - Gessler: P_d = 0.5
    With matching settling velocities and concentration, the deposition
    ratio class1/class0 should be 0.5/0 = inf -- so we use distinct
    concentrations and verify the absolute values directly.
    """
    tau_cs = 0.5  # Pa
    reg = SedimentClassRegistry.from_iterable(
        [
            SedimentClass(label="silt", d50_um=20.0, tau_cs_pa=tau_cs),  # cohesive
            SedimentClass(label="sand", d50_um=200.0, tau_cs_pa=tau_cs),  # non-cohesive
        ]
    )
    # Sanity: classes labelled correctly.
    assert reg[0].is_cohesive is True
    assert reg[1].is_cohesive is False

    # Single cell, tau = tau_cs => P_d_krone = 0, P_d_gessler = 0.5.
    tau = _tau_da([tau_cs])
    susp = _suspended_da([[100.0, 100.0]])  # mg/L for both classes
    settling = np.array([0.05, 2.5])  # cm/s
    depth = _depth_da([10.0])  # m -- huge so cap is non-binding
    dt = 60.0  # s

    flux = compute_deposition_flux(
        registry=reg,
        suspended_concentration=susp,
        tau_pa=tau,
        settling_velocity_cm_s=settling,
        bottom_water_layer_depth_m=depth,
        dt_seconds=dt,
        max_deposit_fraction=1.0,
    )

    assert flux.dims == (contracts.DIM_NFACE, contracts.DIM_CLASS)
    assert flux.shape == (1, 2)

    # Class 0 (cohesive, Krone): P_d = 0 => flux = 0.
    assert flux.values[0, 0] == pytest.approx(0.0, abs=1e-15)

    # Class 1 (non-cohesive, Gessler): P_d = 0.5.
    # D = P_d * C[g/cm3] * w_s[cm/s] * dt[s]
    #   = 0.5 * (100 mg/L * 1e-6) * 2.5 * 60
    #   = 0.5 * 1e-4 * 150
    #   = 7.5e-3 g/cm^2
    expected = 0.5 * 100.0e-6 * 2.5 * 60.0
    assert flux.values[0, 1] == pytest.approx(expected, rel=1e-3)


def test_compute_deposition_flux_mass_cap_engages():
    """Deposition is capped when the analytic D_s exceeds the available mass
    in the bottom water layer.

    Set up a thin bottom layer (depth = 1 cm = 0.01 m) and a long time step
    with a high settling velocity so the uncapped D_s would exceed
    max_deposit_fraction * C * h.

    For a single cohesive class at tau=0 (P_d=1):
        C = 200 mg/L = 2e-4 g/cm^3
        h = 0.01 m = 1 cm
        cap = max_deposit_fraction * C * h = 0.5 * 2e-4 * 1 = 1.0e-4 g/cm^2
        uncapped = P_d * C * w_s * dt = 1 * 2e-4 * 1.0 * 600 = 0.12 g/cm^2

    The cap (1e-4) is far smaller, so the test verifies it kicks in.
    """
    tau_cs = 0.5
    reg = SedimentClassRegistry.from_iterable(
        [SedimentClass(label="mud", d50_um=10.0, tau_cs_pa=tau_cs)]
    )
    tau = _tau_da([0.0])
    susp = _suspended_da([[200.0]])  # mg/L
    settling = np.array([1.0])  # cm/s
    depth = _depth_da([0.01])  # m
    dt = 600.0
    max_frac = 0.5

    flux = compute_deposition_flux(
        registry=reg,
        suspended_concentration=susp,
        tau_pa=tau,
        settling_velocity_cm_s=settling,
        bottom_water_layer_depth_m=depth,
        dt_seconds=dt,
        max_deposit_fraction=max_frac,
    )

    expected_cap = max_frac * 200.0e-6 * (0.01 * 100.0)  # = 1.0e-4 g/cm^2
    assert flux.values[0, 0] == pytest.approx(expected_cap, rel=1e-3)


def test_compute_deposition_flux_below_cap_uses_analytic_value():
    """When the analytic D_s is smaller than the cap, the analytic value is
    returned unchanged."""
    tau_cs = 0.5
    reg = SedimentClassRegistry.from_iterable(
        [SedimentClass(label="mud", d50_um=10.0, tau_cs_pa=tau_cs)]
    )
    tau = _tau_da([0.0])
    susp = _suspended_da([[10.0]])  # mg/L
    settling = np.array([0.01])  # cm/s
    depth = _depth_da([5.0])  # m
    dt = 60.0

    flux = compute_deposition_flux(
        registry=reg,
        suspended_concentration=susp,
        tau_pa=tau,
        settling_velocity_cm_s=settling,
        bottom_water_layer_depth_m=depth,
        dt_seconds=dt,
        max_deposit_fraction=1.0,
    )

    # Analytic: P_d=1 * C(g/cm3)=10e-6 * w=0.01 * dt=60 = 6.0e-6 g/cm^2.
    # Cap:      C * h = 10e-6 * 500 = 5.0e-3 g/cm^2  (orders of magnitude
    #           larger -- non-binding).
    expected = 1.0 * 10.0e-6 * 0.01 * 60.0
    assert flux.values[0, 0] == pytest.approx(expected, rel=1e-3)


def test_compute_deposition_flux_shape_and_dims_for_multiple_cells():
    """Multi-cell, multi-class output should have shape (nface, ssm_class)
    and broadcast properly."""
    tau_cs = 0.4
    reg = SedimentClassRegistry.from_iterable(
        [
            SedimentClass(label="silt", d50_um=20.0, tau_cs_pa=tau_cs),
            SedimentClass(label="sand", d50_um=200.0, tau_cs_pa=tau_cs),
            SedimentClass(label="coarse_sand", d50_um=500.0, tau_cs_pa=tau_cs),
        ]
    )
    nface = 5
    tau = _tau_da(np.linspace(0.0, 2.0 * tau_cs, nface))
    susp = _suspended_da(np.full((nface, 3), 50.0))
    settling = np.array([0.05, 2.5, 5.0])
    depth = _depth_da(np.full(nface, 1.0))
    flux = compute_deposition_flux(
        registry=reg,
        suspended_concentration=susp,
        tau_pa=tau,
        settling_velocity_cm_s=settling,
        bottom_water_layer_depth_m=depth,
        dt_seconds=30.0,
    )
    assert flux.dims == (contracts.DIM_NFACE, contracts.DIM_CLASS)
    assert flux.shape == (nface, 3)
    # All deposition is non-negative.
    assert (flux.values >= 0.0).all()


def test_compute_deposition_flux_validates_settling_length():
    reg = SedimentClassRegistry.from_iterable(
        [SedimentClass(label="silt", d50_um=20.0, tau_cs_pa=0.3)]
    )
    with pytest.raises(ValueError):
        compute_deposition_flux(
            registry=reg,
            suspended_concentration=_suspended_da([[10.0]]),
            tau_pa=_tau_da([0.1]),
            settling_velocity_cm_s=np.array([0.05, 0.10]),  # wrong length
            bottom_water_layer_depth_m=_depth_da([1.0]),
            dt_seconds=60.0,
        )


def test_compute_deposition_flux_requires_tau_cs_set():
    reg = SedimentClassRegistry.from_iterable(
        [SedimentClass(label="mud", d50_um=10.0, tau_cs_pa=None)]
    )
    with pytest.raises(ValueError):
        compute_deposition_flux(
            registry=reg,
            suspended_concentration=_suspended_da([[10.0]]),
            tau_pa=_tau_da([0.1]),
            settling_velocity_cm_s=np.array([0.01]),
            bottom_water_layer_depth_m=_depth_da([1.0]),
            dt_seconds=60.0,
        )
