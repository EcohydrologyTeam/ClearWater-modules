"""C5 wet-mask scope test (review-findings 2026-05-04).

The orchestration-layer wet-mask should mask only variables that a
process *writes* (its declared ``output_variables``), not variables it
*reads* as forcings.

Pre-fix behavior: ``Model.__apply_wet_mask`` iterated
``process.variables`` (the full input + output declaration). After
``Temperature.run``, NaN was written into ``wind_speed``,
``air_temperature``, ``solar_radiation``, ``cloudiness``,
``atmospheric_pressure``, ``atmospheric_vapor_pressure``,
``wetted_surface_area``, ``volume``, ``sediment_thickness`` on dry
cells. Subsequent substeps then read NaN forcing data on the dry
margin until the next chunk reload. v1 explicitly avoided this:
NaN-fill state on dry cells, NaN-mask only the *output* slots at
write time.

Post-fix behavior:

- v3 ``Temperature`` declares
  ``output_variables = ["water_temperature", "sediment_temperature"]``.
- v3 ``Model.__apply_wet_mask`` honors ``output_variables`` if defined
  (falls back to ``variables`` for backward compat).
- Result: only ``water_temperature`` and ``sediment_temperature`` are
  NaN-masked on dry cells; all forcing inputs are preserved.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3 import Temperature
from clearwater_modules_v3.model import Model
from clearwater_modules_v3.processes.base import Process


# ---------------------------------------------------------------------------
# Stubs (intentionally local; mirror the patterns in
# test_model_orchestration_v3.py without coupling to that file's stubs).
# ---------------------------------------------------------------------------


class _StubRegistry:
    def __init__(self) -> None:
        self._data: dict[str, xr.DataArray] = {}

    def register(self, name: str, value) -> None:
        self._data[name] = value

    def get(self, name: str):
        return self._data[name]

    def get_at_time(self, name: str, time):
        if name not in self._data:
            raise KeyError(name)
        return self._data[name]

    def set_at_time(self, name: str, time, value) -> None:
        self._data[name] = value

    def get_variable(self, name: str):
        raise KeyError(name)

    def __contains__(self, name: str) -> bool:
        return name in self._data


def _build_model(
    processes,
    *,
    wet_mask_variable: str | None = None,
    wet_mask_threshold: float = 0.0,
) -> Model:
    return Model(
        processes=tuple(processes),
        variable_registry=_StubRegistry(),
        variable_data_sources={},
        start_time=datetime(2026, 1, 1, 0, 0, 0),
        end_time=datetime(2026, 1, 1, 0, 5, 0),
        time_step=timedelta(minutes=5),
        output_variables=[],
        wet_mask_variable=wet_mask_variable,
        wet_mask_threshold=wet_mask_threshold,
    )


# ---------------------------------------------------------------------------
# C5: Temperature declares output_variables containing only writes
# ---------------------------------------------------------------------------


def test_c5_temperature_declares_output_variables():
    """v3 Temperature class must expose an ``output_variables`` list
    so the orchestrator knows what to mask."""
    assert hasattr(Temperature, "output_variables")
    assert Temperature.output_variables == [
        "water_temperature",
        "sediment_temperature",
    ]


def test_c5_temperature_output_variables_is_subset_of_variables():
    """Sanity: every declared output must also be in the full variables
    list, otherwise the registry I/O paths are inconsistent."""
    for name in Temperature.output_variables:
        assert name in Temperature.variables


def test_c5_temperature_output_variables_excludes_forcings():
    """The forcings that Temperature *reads* must NOT be in
    output_variables. Pre-fix, wet-mask corrupted these."""
    forcings = {
        "wind_speed",
        "air_temperature",
        "solar_radiation",
        "cloudiness",
        "atmospheric_pressure",
        "atmospheric_vapor_pressure",
        "wetted_surface_area",
        "volume",
        "sediment_thickness",
    }
    for name in forcings:
        assert name not in Temperature.output_variables, (
            f"{name!r} is a forcing input and should not be in "
            f"Temperature.output_variables"
        )


# ---------------------------------------------------------------------------
# C5: Model.__apply_wet_mask honors output_variables
# ---------------------------------------------------------------------------


def test_c5_wet_mask_respects_output_variables_only():
    """Build a Model with a Temperature, pre-populate the registry with
    forcing data + state, apply the wet-mask, and assert only outputs
    have NaN written into the dry cell."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    model = _build_model(processes=[t])
    registry = model._Model__registry

    # Three cells: cell 0 wet, cell 1 dry, cell 2 wet.
    wet_mask = xr.DataArray(np.array([True, False, True]))
    # Pre-populate the registry with finite values for outputs *and*
    # forcings. The wet-mask should NaN-mask only the dry-cell entries
    # of the outputs and leave forcings unchanged.
    initial_values = {
        # Outputs (should be masked):
        "water_temperature": xr.DataArray(np.array([20.0, 18.0, 22.0])),
        "sediment_temperature": xr.DataArray(np.array([15.0, 16.0, 17.0])),
        # Forcings (should be preserved):
        "wind_speed": xr.DataArray(np.array([3.0, 3.0, 3.0])),
        "air_temperature": xr.DataArray(np.array([25.0, 25.0, 25.0])),
        "solar_radiation": xr.DataArray(np.array([800.0, 800.0, 800.0])),
        "cloudiness": xr.DataArray(np.array([0.0, 0.0, 0.0])),
        "atmospheric_pressure": xr.DataArray(np.array([1013.0, 1013.0, 1013.0])),
        "atmospheric_vapor_pressure": xr.DataArray(np.array([20.0, 20.0, 20.0])),
        "wetted_surface_area": xr.DataArray(np.array([100.0, 0.0, 100.0])),
        "volume": xr.DataArray(np.array([1000.0, 0.0, 1000.0])),
        "sediment_thickness": xr.DataArray(np.array([0.1, 0.1, 0.1])),
    }
    for name, value in initial_values.items():
        registry.register(name, value)

    # Apply the mask.
    model._Model__apply_wet_mask(t, datetime(2026, 1, 1, 0, 0, 0), wet_mask)

    # Outputs: cell 1 must be NaN; cells 0 and 2 must be unchanged.
    for output_name in Temperature.output_variables:
        result = registry.get(output_name).values
        original = initial_values[output_name].values
        assert np.isnan(result[1]), (
            f"{output_name}: dry-cell value should be NaN, got {result[1]}"
        )
        assert result[0] == original[0], f"{output_name}: wet cell 0 changed"
        assert result[2] == original[2], f"{output_name}: wet cell 2 changed"

    # Forcings: ALL cells must be unchanged. This is the C5 fix.
    forcing_names = [
        "wind_speed",
        "air_temperature",
        "solar_radiation",
        "cloudiness",
        "atmospheric_pressure",
        "atmospheric_vapor_pressure",
        "wetted_surface_area",
        "volume",
        "sediment_thickness",
    ]
    for name in forcing_names:
        result = registry.get(name).values
        original = initial_values[name].values
        np.testing.assert_array_equal(
            result, original,
            err_msg=(
                f"{name!r} is a forcing input and must not be modified by "
                f"the wet-mask. C5 fix preserves it; pre-fix code wrote NaN "
                f"into dry cells."
            ),
        )


# ---------------------------------------------------------------------------
# C5: Backward-compatible fallback for processes without output_variables
# ---------------------------------------------------------------------------


class _LegacyProcess(Process):
    """Process that does NOT declare output_variables; should fall back
    to masking process.variables (the prior behavior)."""

    variables: list[str] = ["legacy_state", "legacy_forcing"]
    # No output_variables attribute.

    def __init__(self) -> None:
        super().__init__(timedelta(minutes=5))

    def run(self, time, registry) -> None:
        pass

    def process_name(self) -> str:
        return "Legacy"


def test_c5_fallback_to_variables_when_output_variables_undeclared():
    """Backward-compat: a process without output_variables falls back
    to masking the full variables list. This preserves v2 / pre-C5
    semantics for processes that haven't been migrated yet."""
    p = _LegacyProcess()
    assert not hasattr(p, "output_variables") or getattr(p, "output_variables", None) is None
    model = _build_model(processes=[p])
    registry = model._Model__registry

    wet_mask = xr.DataArray(np.array([True, False, True]))
    registry.register("legacy_state",   xr.DataArray(np.array([1.0, 2.0, 3.0])))
    registry.register("legacy_forcing", xr.DataArray(np.array([4.0, 5.0, 6.0])))

    model._Model__apply_wet_mask(p, datetime(2026, 1, 1, 0, 0, 0), wet_mask)

    # Both variables masked at cell 1 (the legacy / fallback behavior).
    for name in ("legacy_state", "legacy_forcing"):
        result = registry.get(name).values
        assert np.isnan(result[1])
        assert not np.isnan(result[0])
        assert not np.isnan(result[2])


# ---------------------------------------------------------------------------
# C5: No-op when wet_mask is None
# ---------------------------------------------------------------------------


def test_c5_no_op_when_wet_mask_is_none():
    """When no wet_mask is configured, __apply_wet_mask is a no-op
    (preserving v2 backward-compat: every existing v2 config that
    doesn't use the wet-mask runs unchanged on v3)."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    model = _build_model(processes=[t])  # wet_mask_variable defaults to None
    registry = model._Model__registry
    registry.register(
        "water_temperature", xr.DataArray(np.array([20.0, 18.0, 22.0]))
    )
    original = registry.get("water_temperature").values.copy()

    # Compute the mask: should be None.
    mask = model._Model__compute_wet_mask(datetime(2026, 1, 1, 0, 0, 0))
    assert mask is None

    # Apply: should be a no-op.
    model._Model__apply_wet_mask(t, datetime(2026, 1, 1, 0, 0, 0), mask)
    np.testing.assert_array_equal(
        registry.get("water_temperature").values, original
    )
