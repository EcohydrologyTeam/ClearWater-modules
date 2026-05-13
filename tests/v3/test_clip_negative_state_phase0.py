"""Phase 0 resolutions for ``clip_negative_state``.

Covers pattern-alignment spec §10 Q1 (step-index source) and §10 Q2
(graceful no-op when diagnostics is None / non-DataArray inputs).

Reference: ``design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md``.
"""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.utils.numerics import Diagnostics, clip_negative_state


# ---------------------------------------------------------------------------
# Q2 — graceful no-op when diagnostics is None
# ---------------------------------------------------------------------------


def test_q2_dataarray_input_with_none_diagnostics_clips_without_counting():
    state = xr.DataArray(np.array([1.0, -0.5, 0.0, -2.0]), dims="cell")
    clipped = clip_negative_state(state, "test", diagnostics=None)
    np.testing.assert_array_equal(clipped.values, np.array([1.0, 0.0, 0.0, 0.0]))
    assert isinstance(clipped, xr.DataArray)


def test_q2_ndarray_input_with_diagnostics_clips_and_counts():
    state = np.array([1.0, -0.5, 0.0, -2.0])
    diagnostics = Diagnostics()
    clipped = clip_negative_state(state, "test", diagnostics)
    np.testing.assert_array_equal(clipped, np.array([1.0, 0.0, 0.0, 0.0]))
    assert isinstance(clipped, np.ndarray)
    assert diagnostics.clip_events == {"test": 2}


def test_q2_ndarray_input_with_none_diagnostics_clips_without_counting():
    state = np.array([1.0, -0.5, 0.0, -2.0])
    clipped = clip_negative_state(state, "test", diagnostics=None)
    np.testing.assert_array_equal(clipped, np.array([1.0, 0.0, 0.0, 0.0]))
    assert isinstance(clipped, np.ndarray)


def test_q2_scalar_input_with_diagnostics_clips_and_counts():
    diagnostics = Diagnostics()
    clipped = clip_negative_state(-1.5, "scalar", diagnostics)
    assert clipped == 0.0
    assert diagnostics.clip_events == {"scalar": 1}


def test_q2_scalar_input_with_none_diagnostics_clips_silently():
    clipped = clip_negative_state(-1.5, "scalar", diagnostics=None)
    assert clipped == 0.0


def test_q2_scalar_positive_input_unchanged():
    assert clip_negative_state(1.5, "scalar", Diagnostics()) == 1.5
    assert clip_negative_state(1.5, "scalar", None) == 1.5


def test_q2_dataarray_all_nonnegative_unchanged_with_none_diagnostics():
    state = xr.DataArray(np.array([1.0, 2.0, 3.0]), dims="cell")
    clipped = clip_negative_state(state, "test", diagnostics=None)
    # No clipping needed; identity preserved.
    np.testing.assert_array_equal(clipped.values, state.values)


def test_q2_dataarray_preserves_coords_and_attrs():
    state = xr.DataArray(
        np.array([1.0, -0.5, 0.0]),
        dims="cell",
        coords={"cell": [0, 1, 2]},
        attrs={"units": "mg/L"},
        name="ammonium",
    )
    clipped = clip_negative_state(state, "ammonium", diagnostics=None)
    assert clipped.dims == state.dims
    np.testing.assert_array_equal(clipped.coords["cell"].values, state.coords["cell"].values)
    assert clipped.attrs == state.attrs
    assert clipped.name == state.name


# ---------------------------------------------------------------------------
# Q1 — step index defaults to diagnostics.current_step
# ---------------------------------------------------------------------------


def test_q1_default_current_step_minus_one_logs_as_none():
    """A fresh Diagnostics() has current_step=-1; clip log records use
    step=None (matching the 1.0.0 convention)."""
    diagnostics = Diagnostics()
    assert diagnostics.current_step == -1
    state = xr.DataArray(np.array([-1.0]), dims="cell")
    clip_negative_state(state, "test", diagnostics)
    assert len(diagnostics.clip_log) == 1
    assert diagnostics.clip_log[0]["step"] is None


def test_q1_current_step_set_by_caller_propagates_to_log():
    diagnostics = Diagnostics()
    diagnostics.current_step = 42
    state = xr.DataArray(np.array([-1.0, -2.0]), dims="cell")
    clip_negative_state(state, "test", diagnostics)
    for record in diagnostics.clip_log:
        assert record["step"] == 42


def test_q1_explicit_step_kwarg_overrides_current_step():
    """Backward-compat: callers can still pass step= explicitly."""
    diagnostics = Diagnostics()
    diagnostics.current_step = 42
    state = xr.DataArray(np.array([-1.0]), dims="cell")
    clip_negative_state(state, "test", diagnostics, step=99)
    assert diagnostics.clip_log[0]["step"] == 99


def test_q1_scalar_clip_propagates_current_step():
    diagnostics = Diagnostics()
    diagnostics.current_step = 7
    clip_negative_state(-0.5, "scalar", diagnostics)
    assert diagnostics.clip_log[0]["step"] == 7


# ---------------------------------------------------------------------------
# Backwards-compatibility — pre-Q1/Q2 call sites still work identically.
# ---------------------------------------------------------------------------


def test_back_compat_dataarray_with_diagnostics_no_step_kwarg():
    """Existing call sites (e.g. Nitrogen._clip) pass diagnostics and
    no step kwarg. Default behaviour is unchanged: clip + count + log
    with step=None when current_step is unset."""
    state = xr.DataArray(np.array([5.0, -1.0, 3.0, -2.0, 0.0]), dims="cell")
    diagnostics = Diagnostics()
    clipped = clip_negative_state(state, "ammonium", diagnostics)

    np.testing.assert_array_equal(clipped.values, [5.0, 0.0, 3.0, 0.0, 0.0])
    assert diagnostics.clip_events == {"ammonium": 2}
    assert len(diagnostics.clip_log) == 2


def test_back_compat_detail_limit_truncation():
    """A clip event larger than detail_limit_per_call still emits one
    aggregate record after the per-cell records."""
    diagnostics = Diagnostics(detail_limit_per_call=3)
    state = xr.DataArray(-np.ones(10), dims="cell")
    clip_negative_state(state, "test", diagnostics)
    assert diagnostics.clip_events == {"test": 10}
    # 3 per-cell records + 1 aggregate
    assert len(diagnostics.clip_log) == 4
    assert diagnostics.clip_log[-1]["cell_index"] == "clipped_aggregate"
    assert diagnostics.clip_log[-1]["n_suppressed"] == 7
