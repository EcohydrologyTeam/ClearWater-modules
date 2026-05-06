"""Regression tests for the shared ``sanitize_rate`` helper.

Three claims from the Gemini code review of the v3 NSM1 Process
implementations are addressed by ``utils.numerics.sanitize_rate``:

1. **inf on dry cells.** Division by zero (``x / depth`` at
   ``depth == 0``) produces ``inf``, not ``NaN``. The previous
   per-Process guards checked only ``NaN``, leaving an ``inf``
   pathway through ``clip_negative_state`` (which clips only
   ``< 0``).
2. **Pathogen lacked any rate sanitization.** The Forward Euler
   integrator went directly from ``rate = self.rate(...)`` to
   ``px_new = px + rate * dt_days`` with no guard.
3. **Scalar evaluation gap in DOX/N2 inline guards.** The
   ``isinstance(rate, xr.DataArray)`` / ``isinstance(rate, np.ndarray)``
   branches did silently fall through for native Python ``float``
   inputs.

These tests pin the behavior of the shared helper and confirm the
five Process classes (DOX, Phosphorus, N2, Carbon, Pathogen,
Alkalinity) all import it.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.utils.numerics import sanitize_rate


# ---------------------------------------------------------------------------
# Container-type coverage
# ---------------------------------------------------------------------------

def test_sanitize_rate_dataarray_replaces_nan_and_inf_with_zero():
    rate = xr.DataArray(np.array([1.0, np.nan, np.inf, -np.inf, 2.0]))
    out = sanitize_rate(rate)
    np.testing.assert_array_equal(
        out.values, np.array([1.0, 0.0, 0.0, 0.0, 2.0])
    )
    # Container type preserved.
    assert isinstance(out, xr.DataArray)


def test_sanitize_rate_ndarray_replaces_nan_and_inf_with_zero():
    rate = np.array([1.0, np.nan, np.inf, -np.inf, 2.0])
    out = sanitize_rate(rate)
    np.testing.assert_array_equal(out, np.array([1.0, 0.0, 0.0, 0.0, 2.0]))
    assert isinstance(out, np.ndarray)


def test_sanitize_rate_scalar_float_handles_nan_inf_and_finite():
    """Claim 3: native Python scalar floats no longer fall through
    silently; they are sanitized just like array-valued rates."""
    assert sanitize_rate(float("nan")) == 0.0
    assert sanitize_rate(float("inf")) == 0.0
    assert sanitize_rate(float("-inf")) == 0.0
    assert sanitize_rate(1.5) == 1.5
    assert sanitize_rate(0.0) == 0.0


def test_sanitize_rate_dataarray_preserves_coords_and_dims():
    rate = xr.DataArray(
        np.array([1.0, np.inf, 2.0]),
        coords={"x": [10.0, 20.0, 30.0]},
        dims=["x"],
    )
    out = sanitize_rate(rate)
    assert out.dims == rate.dims
    np.testing.assert_array_equal(out.x.values, np.array([10.0, 20.0, 30.0]))


# ---------------------------------------------------------------------------
# Claim 1: inf from x/depth at depth == 0 is sanitized
# ---------------------------------------------------------------------------

def test_sanitize_rate_catches_division_by_zero_depth():
    """Claim 1: ``x / depth`` at ``depth == 0`` produces ``inf``,
    which the previous NaN-only guards did not catch."""
    x = xr.DataArray(np.array([1.0, 2.0, 3.0]))
    depth = xr.DataArray(np.array([1.0, 0.0, 1.0]))
    with np.errstate(divide="ignore"):
        rate = x / depth
    # Confirm the unsanitized rate is indeed inf at the dry cell.
    assert np.isinf(rate.values[1])
    out = sanitize_rate(rate)
    np.testing.assert_array_equal(out.values, np.array([1.0, 0.0, 3.0]))


# ---------------------------------------------------------------------------
# Claim 2: Pathogen now has the guard
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# All five v3 Process modules import the shared helper
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("module_name", [
    "clearwater_modules_v3.processes.dox",
    "clearwater_modules_v3.processes.phosphorus",
    "clearwater_modules_v3.processes.n2",
    "clearwater_modules_v3.processes.carbon",
    "clearwater_modules_v3.processes.pathogen",
    "clearwater_modules_v3.processes.alkalinity",
])
def test_process_module_imports_sanitize_rate(module_name):
    """Each affected v3 Process imports the shared ``sanitize_rate``
    helper rather than carrying a per-file inline NaN guard."""
    import importlib
    module = importlib.import_module(module_name)
    assert hasattr(module, "sanitize_rate"), (
        f"{module_name} should import sanitize_rate from utils.numerics "
        f"(post-refactor consolidation)"
    )
