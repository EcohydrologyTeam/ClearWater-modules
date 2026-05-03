"""Unit tests for SSM bed shear-stress drivers.

Covers:

* :class:`ExternalShearDriver` — reads τ_b from the mesh; growth limiter
  clamps; missing-input error.
* :class:`CurrentOnlyShearDriver` — Parker (2004) log-law and Manning's
  formulations matched against hand-computed values; edge→face velocity
  reconstruction; ESM composite Manning override.
* :func:`apply_growth_limiter` — vectorized clamping behaviour.
* :class:`WaveCurrentShearDriver` — stub raises NotImplementedError.
"""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment import contracts
from clearwater_modules_v2.processes.sediment.shear import (
    CurrentOnlyShearDriver,
    ExternalShearDriver,
    WaveCurrentShearDriver,
    apply_growth_limiter,
)


# ---------------------------------------------------------------------------
# Synthetic mesh fixtures
# ---------------------------------------------------------------------------

N_FACE = 3
N_EDGE = 2
TIMES = np.array(["2026-01-01T00:00"], dtype="datetime64[ns]")
T0 = TIMES[0]


def _depth(values) -> xr.DataArray:
    """(time, nface) depth field."""
    arr = np.asarray(values, dtype="float64").reshape(1, N_FACE)
    return xr.DataArray(
        arr,
        dims=(contracts.DIM_TIME, contracts.DIM_NFACE),
        coords={contracts.DIM_TIME: TIMES},
    )


def _face_velocity_components(ux_values, uy_values) -> dict[str, xr.DataArray]:
    """Build (time, nface) ux/uy components — the "best path" velocity input."""
    ux = np.asarray(ux_values, dtype="float64").reshape(1, N_FACE)
    uy = np.asarray(uy_values, dtype="float64").reshape(1, N_FACE)
    return {
        "face_velocity_x": xr.DataArray(
            ux,
            dims=(contracts.DIM_TIME, contracts.DIM_NFACE),
            coords={contracts.DIM_TIME: TIMES},
        ),
        "face_velocity_y": xr.DataArray(
            uy,
            dims=(contracts.DIM_TIME, contracts.DIM_NFACE),
            coords={contracts.DIM_TIME: TIMES},
        ),
    }


def _edge_velocity(values) -> xr.DataArray:
    arr = np.asarray(values, dtype="float64").reshape(1, N_EDGE)
    return xr.DataArray(
        arr,
        dims=(contracts.DIM_TIME, contracts.DIM_NEDGE),
        coords={contracts.DIM_TIME: TIMES},
    )


def _make_mesh_with_face_velocity(
    depths,
    ux,
    uy,
    *,
    mannings_n: float | None = None,
    composite_n: float | None = None,
    bed_shear_input: np.ndarray | None = None,
) -> xr.Dataset:
    """3-cell mesh using face-velocity components (skips edge averaging).

    Parameters take per-face vectors (length 3) for depth, ux, uy.
    """
    data: dict[str, xr.DataArray] = {
        contracts.VAR_FACE_HYDRAULIC_DEPTH: _depth(depths),
        **_face_velocity_components(ux, uy),
    }
    if mannings_n is not None:
        data[contracts.VAR_MANNINGS_N] = xr.DataArray(
            np.full(N_FACE, mannings_n, dtype="float64"),
            dims=(contracts.DIM_NFACE,),
        )
    if composite_n is not None:
        data[contracts.VAR_COMPOSITE_MANNINGS_N] = xr.DataArray(
            np.full((1, N_FACE), composite_n, dtype="float64"),
            dims=(contracts.DIM_TIME, contracts.DIM_NFACE),
            coords={contracts.DIM_TIME: TIMES},
        )
    if bed_shear_input is not None:
        data[contracts.VAR_BED_SHEAR_STRESS_INPUT] = xr.DataArray(
            np.asarray(bed_shear_input, dtype="float64").reshape(1, N_FACE),
            dims=(contracts.DIM_TIME, contracts.DIM_NFACE),
            coords={contracts.DIM_TIME: TIMES},
        )
    return xr.Dataset(data)


def _make_mesh_with_edge_velocity(depths, edge_vels) -> xr.Dataset:
    """3-cell mesh wired up so the driver must reconstruct face velocity
    from ``edge_velocity`` via the edges_face1/edges_face2 connectivity.

    Edge layout (2 edges, 3 cells; cells 0–1 share edge 0, cells 1–2 share
    edge 1). With this connectivity, the per-face accumulated speed is:

    * face 0: |U_e0|              (single neighbour)
    * face 1: (|U_e0| + |U_e1|)/2 (two neighbours)
    * face 2: |U_e1|              (single neighbour)
    """
    edges_face1 = xr.DataArray(np.array([0, 1], dtype=np.int64), dims=(contracts.DIM_NEDGE,))
    edges_face2 = xr.DataArray(np.array([1, 2], dtype=np.int64), dims=(contracts.DIM_NEDGE,))
    data = {
        contracts.VAR_FACE_HYDRAULIC_DEPTH: _depth(depths),
        contracts.VAR_EDGE_VELOCITY: _edge_velocity(edge_vels),
        "edges_face1": edges_face1,
        "edges_face2": edges_face2,
    }
    mesh = xr.Dataset(data)
    # The face-velocity averaging in shear.py uses mesh.sizes[DIM_NFACE].
    # Ensure the dim is materialised on the dataset by attaching a coord.
    mesh = mesh.assign_coords({contracts.DIM_NFACE: np.arange(N_FACE)})
    return mesh


def _zeros_face() -> xr.DataArray:
    """(nface,) zero array for previous_tau / d50_surface defaults."""
    return xr.DataArray(np.zeros(N_FACE, dtype="float64"), dims=(contracts.DIM_NFACE,))


# ---------------------------------------------------------------------------
# apply_growth_limiter
# ---------------------------------------------------------------------------


def test_growth_limiter_clamps_when_new_exceeds_threshold():
    """τ_prev=0.5, τ_new=1.0, growth=0.10 → τ_out = 0.5 + 0.10·(1.0−0.5) = 0.55."""
    tau_new = xr.DataArray(np.array([1.0]))
    tau_prev = xr.DataArray(np.array([0.5]))
    out = apply_growth_limiter(tau_new, tau_prev, growth_limit=0.10)
    np.testing.assert_allclose(out.values, [0.55])


def test_growth_limiter_passes_through_below_threshold():
    """If τ_new is within the allowed band, return τ_new unchanged."""
    tau_new = xr.DataArray(np.array([0.52]))      # < 0.55 threshold
    tau_prev = xr.DataArray(np.array([0.50]))
    out = apply_growth_limiter(tau_new, tau_prev, growth_limit=0.10)
    np.testing.assert_allclose(out.values, [0.52])


def test_growth_limiter_passes_through_decreases():
    """Drops in τ are not clamped — only growth is limited."""
    tau_new = xr.DataArray(np.array([0.10]))
    tau_prev = xr.DataArray(np.array([1.00]))
    out = apply_growth_limiter(tau_new, tau_prev, growth_limit=0.10)
    np.testing.assert_allclose(out.values, [0.10])


def test_growth_limiter_disabled_when_zero():
    tau_new = xr.DataArray(np.array([10.0]))
    tau_prev = xr.DataArray(np.array([1.0]))
    out = apply_growth_limiter(tau_new, tau_prev, growth_limit=0.0)
    np.testing.assert_allclose(out.values, [10.0])


def test_growth_limiter_vectorized_per_cell():
    """Different cells: one clamped, one passthrough, one decreasing."""
    tau_new = xr.DataArray(np.array([1.0, 0.52, 0.10]))
    tau_prev = xr.DataArray(np.array([0.5, 0.50, 1.00]))
    out = apply_growth_limiter(tau_new, tau_prev, growth_limit=0.10)
    np.testing.assert_allclose(out.values, [0.55, 0.52, 0.10])


# ---------------------------------------------------------------------------
# ExternalShearDriver
# ---------------------------------------------------------------------------


def test_external_shear_reads_input_field():
    bed_shear = np.array([0.4, 0.8, 1.2])
    mesh = _make_mesh_with_face_velocity(
        depths=[1.0, 1.0, 1.0],
        ux=[0.0, 0.0, 0.0],
        uy=[0.0, 0.0, 0.0],
        bed_shear_input=bed_shear,
    )
    # Disable growth limiter so we just see the raw input.
    driver = ExternalShearDriver(growth_limit=0.0)
    tau = driver.compute(
        mesh,
        T0,
        d50_surface_um=_zeros_face(),
        previous_tau_pa=_zeros_face(),
    )
    np.testing.assert_allclose(tau.values, bed_shear)


def test_external_shear_applies_growth_limiter():
    """Verify the limiter clamps when the imported field jumps above prev."""
    bed_shear = np.array([1.0, 1.0, 1.0])
    mesh = _make_mesh_with_face_velocity(
        depths=[1.0, 1.0, 1.0],
        ux=[0.0, 0.0, 0.0],
        uy=[0.0, 0.0, 0.0],
        bed_shear_input=bed_shear,
    )
    driver = ExternalShearDriver(growth_limit=0.10)
    prev = xr.DataArray(np.array([0.5, 0.5, 0.5]), dims=(contracts.DIM_NFACE,))
    tau = driver.compute(
        mesh,
        T0,
        d50_surface_um=_zeros_face(),
        previous_tau_pa=prev,
    )
    # 0.5 + 0.10·(1.0 − 0.5) = 0.55
    np.testing.assert_allclose(tau.values, [0.55, 0.55, 0.55])


def test_external_shear_missing_input_raises_value_error():
    mesh = _make_mesh_with_face_velocity(
        depths=[1.0, 1.0, 1.0],
        ux=[0.0, 0.0, 0.0],
        uy=[0.0, 0.0, 0.0],
        # bed_shear_input intentionally omitted
    )
    driver = ExternalShearDriver()
    with pytest.raises(ValueError, match="bed_shear_stress_input"):
        driver.compute(
            mesh,
            T0,
            d50_surface_um=_zeros_face(),
            previous_tau_pa=_zeros_face(),
        )


# ---------------------------------------------------------------------------
# CurrentOnlyShearDriver — log-law
# ---------------------------------------------------------------------------


def test_current_only_log_law_matches_hand_computation():
    """Single-cell Parker (2004) log-law hand check:

        U   = 0.5 m/s
        h   = 2.0 m
        k_n = D50 = 1e-3 m  (so ln(11·h / (2·k_n)) = ln(11 000))
        f_c = (0.42 / ln(11·2 / (2·1e-3)))^2
            = (0.42 / 9.3057)^2 ≈ 2.04e-3
        τ   = 1000 · f_c · 0.5^2 ≈ 0.509 Pa

    (The design memo's eyeballed "≈ 1.1 Pa" is off by ~2× because it
    used f_c ≈ 4.4e-3 instead of the actual 2.04e-3; the formula here
    matches Parker 2004 / SAND2008-5621 §S_SHEAR.f90 line 261.)
    """
    # Pick D50 = 1000 μm = 1e-3 m so k_n collapses to D50 (≥ zb_skin).
    mesh = _make_mesh_with_face_velocity(
        depths=[2.0, 2.0, 2.0],
        ux=[0.5, 0.5, 0.5],
        uy=[0.0, 0.0, 0.0],
    )
    d50 = xr.DataArray(np.full(N_FACE, 1000.0, dtype="float64"), dims=(contracts.DIM_NFACE,))
    driver = CurrentOnlyShearDriver(
        formulation="log_law",
        zb_skin_m=0.0,           # rely entirely on D50 = 1 mm
        growth_limit=0.0,        # disable for a clean comparison
    )
    tau = driver.compute(mesh, T0, d50_surface_um=d50, previous_tau_pa=_zeros_face())

    expected_fc = (0.42 / np.log(11.0 * 2.0 / (2.0 * 1.0e-3))) ** 2
    expected_tau = 1000.0 * expected_fc * 0.5 ** 2
    # Sanity: ~0.509 Pa, within 1 % of the hand-computed value.
    assert abs(expected_tau - 0.509) / 0.509 < 0.01
    np.testing.assert_allclose(tau.values, expected_tau, rtol=1e-2)


def test_current_only_log_law_uses_zb_skin_floor():
    """When D50 < zb_skin, k_n should fall back to zb_skin."""
    mesh = _make_mesh_with_face_velocity(
        depths=[2.0, 2.0, 2.0],
        ux=[0.5, 0.5, 0.5],
        uy=[0.0, 0.0, 0.0],
    )
    # D50 = 100 μm = 1e-4 m, smaller than zb_skin = 1.5 mm = 1.5e-3 m.
    d50 = xr.DataArray(np.full(N_FACE, 100.0, dtype="float64"), dims=(contracts.DIM_NFACE,))
    driver = CurrentOnlyShearDriver(
        formulation="log_law",
        zb_skin_m=1.5e-3,
        growth_limit=0.0,
    )
    tau = driver.compute(mesh, T0, d50_surface_um=d50, previous_tau_pa=_zeros_face())

    expected_fc = (0.42 / np.log(11.0 * 2.0 / (2.0 * 1.5e-3))) ** 2
    expected_tau = 1000.0 * expected_fc * 0.5 ** 2
    np.testing.assert_allclose(tau.values, expected_tau, rtol=1e-6)


def test_current_only_log_law_velocity_magnitude_uses_components():
    """|U|² should pick up both x and y components."""
    mesh = _make_mesh_with_face_velocity(
        depths=[2.0, 2.0, 2.0],
        ux=[0.3, 0.3, 0.3],
        uy=[0.4, 0.4, 0.4],   # |U| = 0.5
    )
    d50 = xr.DataArray(np.full(N_FACE, 1000.0, dtype="float64"), dims=(contracts.DIM_NFACE,))
    driver = CurrentOnlyShearDriver(
        formulation="log_law",
        zb_skin_m=0.0,
        growth_limit=0.0,
    )
    tau = driver.compute(mesh, T0, d50_surface_um=d50, previous_tau_pa=_zeros_face())
    expected_fc = (0.42 / np.log(11.0 * 2.0 / (2.0 * 1.0e-3))) ** 2
    expected_tau = 1000.0 * expected_fc * 0.5 ** 2
    np.testing.assert_allclose(tau.values, expected_tau, rtol=1e-6)


def test_current_only_edge_to_face_velocity_reconstruction():
    """Without face_velocity_x/y, the driver must average edge_velocity
    onto faces using edges_face1/edges_face2.

    Using the layout in :func:`_make_mesh_with_edge_velocity`:
      edge speeds = [0.4, 0.6]
      face speeds = [0.4, 0.5, 0.6]   (avg of incident edges)
    """
    mesh = _make_mesh_with_edge_velocity(
        depths=[2.0, 2.0, 2.0],
        edge_vels=[0.4, -0.6],   # absolute values used
    )
    d50 = xr.DataArray(np.full(N_FACE, 1000.0, dtype="float64"), dims=(contracts.DIM_NFACE,))
    driver = CurrentOnlyShearDriver(
        formulation="log_law",
        zb_skin_m=0.0,
        growth_limit=0.0,
    )
    tau = driver.compute(mesh, T0, d50_surface_um=d50, previous_tau_pa=_zeros_face())

    fc = (0.42 / np.log(11.0 * 2.0 / (2.0 * 1.0e-3))) ** 2
    expected = 1000.0 * fc * np.array([0.4, 0.5, 0.6]) ** 2
    np.testing.assert_allclose(tau.values, expected, rtol=1e-6)


def test_current_only_log_law_applies_growth_limiter():
    """End-to-end: a τ jump that exceeds the growth threshold is clamped.

    With h=2 m, U=0.5 m/s, k_n=1e-3 m, the raw τ ≈ 0.509 Pa. Setting
    τ_prev = 0.10 Pa and growth = 0.10 puts the threshold at 0.11 Pa,
    so the new value is clamped to
    τ_out = 0.10 + 0.10·(0.509 − 0.10) ≈ 0.141 Pa.
    """
    mesh = _make_mesh_with_face_velocity(
        depths=[2.0, 2.0, 2.0],
        ux=[0.5, 0.5, 0.5],
        uy=[0.0, 0.0, 0.0],
    )
    d50 = xr.DataArray(np.full(N_FACE, 1000.0, dtype="float64"), dims=(contracts.DIM_NFACE,))
    prev = xr.DataArray(np.full(N_FACE, 0.10, dtype="float64"), dims=(contracts.DIM_NFACE,))
    driver = CurrentOnlyShearDriver(
        formulation="log_law",
        zb_skin_m=0.0,
        growth_limit=0.10,
    )
    tau = driver.compute(mesh, T0, d50_surface_um=d50, previous_tau_pa=prev)
    fc = (0.42 / np.log(11.0 * 2.0 / (2.0 * 1.0e-3))) ** 2
    raw_tau = 1000.0 * fc * 0.5 ** 2
    # Confirm we're truly in the clamping regime.
    assert raw_tau > 0.10 * (1.0 + 0.10)
    expected = 0.10 + 0.10 * (raw_tau - 0.10)
    np.testing.assert_allclose(tau.values, expected, rtol=1e-6)


# ---------------------------------------------------------------------------
# CurrentOnlyShearDriver — Manning
# ---------------------------------------------------------------------------


def test_current_only_manning_matches_hand_computation():
    """Single-cell Manning hand check:

        n = 0.030 s/m^(1/3), h = 2.0 m, U = 0.5 m/s
        f_c = g·n²/h^(1/3) = 9.81 · 0.030² / 2^(1/3)
        τ   = 1000 · f_c · 0.5²
    """
    mesh = _make_mesh_with_face_velocity(
        depths=[2.0, 2.0, 2.0],
        ux=[0.5, 0.5, 0.5],
        uy=[0.0, 0.0, 0.0],
        mannings_n=0.030,
    )
    d50 = _zeros_face()
    driver = CurrentOnlyShearDriver(
        formulation="manning",
        growth_limit=0.0,
        use_composite_manning=False,
    )
    tau = driver.compute(mesh, T0, d50_surface_um=d50, previous_tau_pa=_zeros_face())

    expected_fc = 9.81 * (0.030 ** 2) / (2.0 ** (1.0 / 3.0))
    expected_tau = 1000.0 * expected_fc * 0.5 ** 2
    np.testing.assert_allclose(tau.values, expected_tau, rtol=1e-6)


def test_current_only_manning_uses_composite_when_supplied():
    """When ESM publishes composite_manning_n, it should override the
    static field if ``use_composite_manning=True``."""
    mesh = _make_mesh_with_face_velocity(
        depths=[2.0, 2.0, 2.0],
        ux=[0.5, 0.5, 0.5],
        uy=[0.0, 0.0, 0.0],
        mannings_n=0.025,        # static (ignored when composite is on)
        composite_n=0.080,       # ESM-supplied (vegetation-roughened)
    )
    driver = CurrentOnlyShearDriver(
        formulation="manning",
        growth_limit=0.0,
        use_composite_manning=True,
    )
    tau = driver.compute(mesh, T0, d50_surface_um=_zeros_face(), previous_tau_pa=_zeros_face())
    expected_fc = 9.81 * (0.080 ** 2) / (2.0 ** (1.0 / 3.0))
    expected_tau = 1000.0 * expected_fc * 0.5 ** 2
    np.testing.assert_allclose(tau.values, expected_tau, rtol=1e-6)


def test_current_only_manning_static_when_composite_disabled():
    """``use_composite_manning=False`` keeps the static field even if ESM
    has populated the composite."""
    mesh = _make_mesh_with_face_velocity(
        depths=[2.0, 2.0, 2.0],
        ux=[0.5, 0.5, 0.5],
        uy=[0.0, 0.0, 0.0],
        mannings_n=0.025,
        composite_n=0.080,
    )
    driver = CurrentOnlyShearDriver(
        formulation="manning",
        growth_limit=0.0,
        use_composite_manning=False,
    )
    tau = driver.compute(mesh, T0, d50_surface_um=_zeros_face(), previous_tau_pa=_zeros_face())
    expected_fc = 9.81 * (0.025 ** 2) / (2.0 ** (1.0 / 3.0))
    expected_tau = 1000.0 * expected_fc * 0.5 ** 2
    np.testing.assert_allclose(tau.values, expected_tau, rtol=1e-6)


# ---------------------------------------------------------------------------
# CurrentOnlyShearDriver — error / config paths
# ---------------------------------------------------------------------------


def test_current_only_invalid_formulation_raises():
    with pytest.raises(ValueError, match="formulation"):
        CurrentOnlyShearDriver(formulation="bogus")


def test_current_only_missing_depth_raises():
    mesh = xr.Dataset(
        _face_velocity_components([0.5, 0.5, 0.5], [0.0, 0.0, 0.0])
    )
    driver = CurrentOnlyShearDriver(formulation="log_law", growth_limit=0.0)
    with pytest.raises(ValueError, match="face_hydraulic_depth"):
        driver.compute(
            mesh,
            T0,
            d50_surface_um=_zeros_face(),
            previous_tau_pa=_zeros_face(),
        )


def test_current_only_missing_velocity_raises():
    mesh = xr.Dataset({contracts.VAR_FACE_HYDRAULIC_DEPTH: _depth([1.0, 1.0, 1.0])})
    driver = CurrentOnlyShearDriver(formulation="log_law", growth_limit=0.0)
    with pytest.raises(ValueError, match="velocity"):
        driver.compute(
            mesh,
            T0,
            d50_surface_um=_zeros_face(),
            previous_tau_pa=_zeros_face(),
        )


def test_current_only_manning_missing_n_raises():
    mesh = _make_mesh_with_face_velocity(
        depths=[2.0, 2.0, 2.0],
        ux=[0.5, 0.5, 0.5],
        uy=[0.0, 0.0, 0.0],
        # mannings_n omitted
    )
    driver = CurrentOnlyShearDriver(
        formulation="manning",
        growth_limit=0.0,
        use_composite_manning=False,
    )
    with pytest.raises(ValueError, match="mannings_n"):
        driver.compute(
            mesh,
            T0,
            d50_surface_um=_zeros_face(),
            previous_tau_pa=_zeros_face(),
        )


# ---------------------------------------------------------------------------
# WaveCurrentShearDriver
# ---------------------------------------------------------------------------


def test_wave_current_shear_driver_is_stub():
    """Mode C must be unambiguously deferred to phase 5."""
    driver = WaveCurrentShearDriver()
    mesh = _make_mesh_with_face_velocity(
        depths=[1.0, 1.0, 1.0],
        ux=[0.1, 0.1, 0.1],
        uy=[0.0, 0.0, 0.0],
    )
    with pytest.raises(NotImplementedError, match="phase 5"):
        driver.compute(
            mesh,
            T0,
            d50_surface_um=_zeros_face(),
            previous_tau_pa=_zeros_face(),
        )
