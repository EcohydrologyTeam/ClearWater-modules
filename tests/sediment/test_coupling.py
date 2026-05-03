"""Unit tests for SSM/ESM coupling helpers."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment import contracts
from clearwater_modules_v2.processes.sediment.coupling import (
    read_composite_manning_n,
    read_vegetation_feedback,
)


N_FACE = 4
TIMES = np.array(
    ["2026-01-01T00:00", "2026-01-01T01:00", "2026-01-01T02:00"],
    dtype="datetime64[ns]",
)


def _make_time_var(value: float) -> xr.DataArray:
    """(time, nface) DataArray filled with ``value``."""
    return xr.DataArray(
        np.full((TIMES.size, N_FACE), value, dtype="float32"),
        dims=(contracts.DIM_TIME, contracts.DIM_NFACE),
        coords={contracts.DIM_TIME: TIMES},
    )


def _make_static_var(value: float) -> xr.DataArray:
    return xr.DataArray(
        np.full(N_FACE, value, dtype="float32"),
        dims=(contracts.DIM_NFACE,),
    )


# ---------------------------------------------------------------------------
# read_vegetation_feedback
# ---------------------------------------------------------------------------


def test_vegetation_feedback_all_present():
    mesh = xr.Dataset(
        {
            contracts.VAR_VEGETATION_BIOSTABILIZATION: _make_time_var(0.7),
            contracts.VAR_VEGETATION_ROOT_COHESION: _make_time_var(0.3),
            contracts.VAR_VEGETATION_FRONTAL_AREA: _make_time_var(0.2),
        }
    )
    bio, root, frontal = read_vegetation_feedback(mesh, TIMES[1])

    assert bio is not None and root is not None and frontal is not None
    assert bio.dims == (contracts.DIM_NFACE,)
    np.testing.assert_allclose(bio.values, 0.7)
    np.testing.assert_allclose(root.values, 0.3)
    np.testing.assert_allclose(frontal.values, 0.2)


def test_vegetation_feedback_all_missing():
    mesh = xr.Dataset()
    bio, root, frontal = read_vegetation_feedback(mesh, TIMES[0])
    assert bio is None
    assert root is None
    assert frontal is None


def test_vegetation_feedback_partial_missing():
    """Only biostabilization is present; the other two should be None."""
    mesh = xr.Dataset(
        {contracts.VAR_VEGETATION_BIOSTABILIZATION: _make_time_var(0.42)}
    )
    bio, root, frontal = read_vegetation_feedback(mesh, TIMES[2])
    assert bio is not None
    np.testing.assert_allclose(bio.values, 0.42)
    assert root is None
    assert frontal is None


def test_vegetation_feedback_selects_correct_time_slice():
    """Confirm the returned slice is the time we asked for."""
    arr = xr.DataArray(
        np.array(
            [
                np.full(N_FACE, 0.10, dtype="float32"),
                np.full(N_FACE, 0.50, dtype="float32"),
                np.full(N_FACE, 0.90, dtype="float32"),
            ]
        ),
        dims=(contracts.DIM_TIME, contracts.DIM_NFACE),
        coords={contracts.DIM_TIME: TIMES},
    )
    mesh = xr.Dataset({contracts.VAR_VEGETATION_BIOSTABILIZATION: arr})
    bio, _, _ = read_vegetation_feedback(mesh, TIMES[2])
    np.testing.assert_allclose(bio.values, 0.90)


def test_vegetation_feedback_accepts_static_field():
    """If ESM happens to publish a non-time-varying field, return it as-is."""
    mesh = xr.Dataset(
        {contracts.VAR_VEGETATION_ROOT_COHESION: _make_static_var(1.5)}
    )
    _, root, _ = read_vegetation_feedback(mesh, TIMES[0])
    assert root is not None
    assert contracts.DIM_TIME not in root.dims
    np.testing.assert_allclose(root.values, 1.5)


# ---------------------------------------------------------------------------
# read_composite_manning_n
# ---------------------------------------------------------------------------


def test_composite_manning_returns_esm_value_when_present():
    fallback = _make_static_var(0.025)
    mesh = xr.Dataset(
        {contracts.VAR_COMPOSITE_MANNINGS_N: _make_time_var(0.060)}
    )
    out = read_composite_manning_n(mesh, TIMES[0], fallback)
    assert contracts.DIM_TIME not in out.dims
    np.testing.assert_allclose(out.values, 0.060)


def test_composite_manning_falls_back_when_missing():
    fallback = _make_static_var(0.025)
    mesh = xr.Dataset()
    out = read_composite_manning_n(mesh, TIMES[0], fallback)
    # Should be the fallback object, untouched.
    np.testing.assert_allclose(out.values, 0.025)
    assert contracts.DIM_TIME not in out.dims


def test_composite_manning_selects_correct_time_slice():
    arr = xr.DataArray(
        np.array(
            [
                np.full(N_FACE, 0.020, dtype="float32"),
                np.full(N_FACE, 0.040, dtype="float32"),
                np.full(N_FACE, 0.080, dtype="float32"),
            ]
        ),
        dims=(contracts.DIM_TIME, contracts.DIM_NFACE),
        coords={contracts.DIM_TIME: TIMES},
    )
    mesh = xr.Dataset({contracts.VAR_COMPOSITE_MANNINGS_N: arr})
    out = read_composite_manning_n(mesh, TIMES[1], _make_static_var(0.025))
    np.testing.assert_allclose(out.values, 0.040)


def test_composite_manning_static_var_pass_through():
    """If ESM publishes a static composite_manning_n, return it as-is."""
    static_n = _make_static_var(0.045)
    mesh = xr.Dataset({contracts.VAR_COMPOSITE_MANNINGS_N: static_n})
    out = read_composite_manning_n(mesh, TIMES[0], _make_static_var(0.025))
    np.testing.assert_allclose(out.values, 0.045)
