"""Unit tests for SSM armoring helpers.

Tests:
- ``compute_d50_avg`` mass-weighted mean of class D50s.
- ``interpolate_taucrit_from_d50`` piecewise-linear lookup with
  endpoint-clamping (mirrors SEDZLJ ``s_sedzlj.f90:242-263``).
"""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment import contracts
from clearwater_modules_v2.processes.sediment.armoring import (
    compute_d50_avg,
    interpolate_taucrit_from_d50,
)


# ---------------------------------------------------------------------------
# compute_d50_avg
# ---------------------------------------------------------------------------


def test_d50_avg_50_50_two_class_mix():
    """50 % at D50=50 + 50 % at D50=200 should give D50_avg = 125 μm."""
    n_face = 4
    fractions = xr.DataArray(
        np.tile([0.5, 0.5], (n_face, 1)),
        dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
    )
    d50 = np.array([50.0, 200.0])
    out = compute_d50_avg(fractions, d50)
    assert out.dims == (contracts.DIM_NFACE,)
    np.testing.assert_allclose(out.values, np.full(n_face, 125.0))


def test_d50_avg_pure_class_gives_class_d50():
    """A cell that is 100 % class k should report D50_avg = D50_k."""
    fractions = xr.DataArray(
        np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        ),
        dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
    )
    d50 = np.array([10.0, 100.0, 500.0])
    out = compute_d50_avg(fractions, d50)
    np.testing.assert_allclose(out.values, [10.0, 100.0, 500.0])


def test_d50_avg_general_three_class_mix():
    """f = (0.2, 0.3, 0.5), D50 = (50, 100, 300) -> 0.2*50 + 0.3*100 + 0.5*300 = 190."""
    fractions = xr.DataArray(
        np.array([[0.2, 0.3, 0.5]]),
        dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
    )
    d50 = np.array([50.0, 100.0, 300.0])
    out = compute_d50_avg(fractions, d50)
    np.testing.assert_allclose(out.values, [190.0])


def test_d50_avg_preserves_face_dim_and_drops_class_dim():
    fractions = xr.DataArray(
        np.full((5, 3), 1.0 / 3.0),
        dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
    )
    d50 = np.array([50.0, 100.0, 300.0])
    out = compute_d50_avg(fractions, d50)
    assert contracts.DIM_CLASS not in out.dims
    assert contracts.DIM_NFACE in out.dims
    assert out.sizes[contracts.DIM_NFACE] == 5
    np.testing.assert_allclose(out.values, np.full(5, 150.0))


# ---------------------------------------------------------------------------
# interpolate_taucrit_from_d50
# ---------------------------------------------------------------------------


def test_taucrit_midpoint_linear_interp():
    """Lookup [(125, 1.20), (222, 2.27)] at D50=173 (midpoint) should
    give linear interp ≈ 1.74 Pa.

    Hand check: 173 is at fraction (173-125)/(222-125) = 48/97 = 0.4948
    of the interval, so τ_c = 1.20 + 0.4948 * (2.27 - 1.20)
                            = 1.20 + 0.4948 * 1.07 = 1.7295 Pa.
    """
    sizes = np.array([125.0, 222.0])
    taus = np.array([1.20, 2.27])
    d50 = xr.DataArray([173.0], dims=(contracts.DIM_NFACE,))
    out = interpolate_taucrit_from_d50(d50, sizes, taus)
    np.testing.assert_allclose(out.values, [1.7295], atol=1e-3)


def test_taucrit_clamps_below_table_min():
    sizes = np.array([125.0, 222.0])
    taus = np.array([1.20, 2.27])
    d50 = xr.DataArray([10.0, 50.0, 100.0], dims=(contracts.DIM_NFACE,))
    out = interpolate_taucrit_from_d50(d50, sizes, taus)
    np.testing.assert_allclose(out.values, [1.20, 1.20, 1.20])


def test_taucrit_clamps_above_table_max():
    sizes = np.array([125.0, 222.0])
    taus = np.array([1.20, 2.27])
    d50 = xr.DataArray([300.0, 1000.0], dims=(contracts.DIM_NFACE,))
    out = interpolate_taucrit_from_d50(d50, sizes, taus)
    np.testing.assert_allclose(out.values, [2.27, 2.27])


def test_taucrit_exact_endpoints_returned():
    sizes = np.array([125.0, 222.0])
    taus = np.array([1.20, 2.27])
    d50 = xr.DataArray([125.0, 222.0], dims=(contracts.DIM_NFACE,))
    out = interpolate_taucrit_from_d50(d50, sizes, taus)
    np.testing.assert_allclose(out.values, [1.20, 2.27])


def test_taucrit_multi_segment_table():
    """4-point table -- ensure interpolation picks the correct segment."""
    sizes = np.array([50.0, 100.0, 200.0, 400.0])
    taus = np.array([0.10, 0.50, 1.20, 3.00])
    # Midpoint of segment 1 (50-100): 75 -> 0.30
    # Midpoint of segment 2 (100-200): 150 -> 0.85
    # Midpoint of segment 3 (200-400): 300 -> 2.10
    d50 = xr.DataArray([75.0, 150.0, 300.0], dims=(contracts.DIM_NFACE,))
    out = interpolate_taucrit_from_d50(d50, sizes, taus)
    np.testing.assert_allclose(out.values, [0.30, 0.85, 2.10])


def test_taucrit_sorts_unsorted_input_table():
    """Pass an unsorted table; result should match the sorted equivalent."""
    sizes = np.array([222.0, 125.0])
    taus = np.array([2.27, 1.20])
    d50 = xr.DataArray([173.0], dims=(contracts.DIM_NFACE,))
    out = interpolate_taucrit_from_d50(d50, sizes, taus)
    np.testing.assert_allclose(out.values, [1.7295], atol=1e-3)


def test_taucrit_rejects_empty_table():
    with pytest.raises(ValueError):
        interpolate_taucrit_from_d50(
            xr.DataArray([100.0], dims=(contracts.DIM_NFACE,)),
            np.array([]),
            np.array([]),
        )


def test_taucrit_rejects_mismatched_table_lengths():
    with pytest.raises(ValueError):
        interpolate_taucrit_from_d50(
            xr.DataArray([100.0], dims=(contracts.DIM_NFACE,)),
            np.array([100.0, 200.0]),
            np.array([1.0, 2.0, 3.0]),
        )
