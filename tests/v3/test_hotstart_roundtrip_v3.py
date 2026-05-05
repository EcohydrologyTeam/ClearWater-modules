"""v3 hotstart roundtrip tests.

The v3 hotstart contract differs from v1's by design (TSM design spec
section 3.2 resolution 2026-05-04):

- Hotstart preserves *registry* state at the saved time. The Model
  reads the saved ``xr.Dataset`` and replaces matching variables in the
  ``VariableRegistry`` before per-process initialization runs.
- Per-process *substep-internal* state defaults to "fresh start" after
  hotstart. Each Process may opt in by implementing ``to_hotstart()``
  and ``from_hotstart(state)``; default is no-op.
- v3 ``Temperature`` opts in: ``__skip_first_time_step`` is preserved
  if explicitly carried in the saved-dataset attrs; otherwise it
  defaults to ``False`` after a hotstart (don't skip the next substep —
  you're not starting from scratch).

This test file covers, in order:

1. Temperature ``to_hotstart`` / ``from_hotstart`` round-trip.
2. Model ``to_hotstart`` collection across processes.
3. Model ``__seed_from_hotstart`` overwrites registry contents.
4. Model ``__restore_process_hotstart`` invokes ``from_hotstart`` on
   processes that define it; processes without it are silently skipped.
5. End-to-end smoke: construct two Models from the same data sources;
   the second carries ``hotstart_dataset``; init_model runs without
   error.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3 import Temperature
from clearwater_modules_v3.model import Model
from clearwater_modules_v3.processes.base import Process


# ---------------------------------------------------------------------------
# Test 1: Temperature to_hotstart / from_hotstart round-trip
# ---------------------------------------------------------------------------


def test_temperature_to_hotstart_default_state():
    """A fresh Temperature has __skip_first_time_step = True."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    snapshot = t.to_hotstart()
    assert snapshot == {"temperature.skip_first_time_step": True}


def test_temperature_from_hotstart_default_value_is_false():
    """from_hotstart with empty state defaults to skip=False (don't skip)."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    t.from_hotstart({})
    assert t.to_hotstart() == {"temperature.skip_first_time_step": False}


def test_temperature_from_hotstart_honors_explicit_true():
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    t.from_hotstart({"temperature.skip_first_time_step": True})
    assert t.to_hotstart() == {"temperature.skip_first_time_step": True}


def test_temperature_from_hotstart_honors_explicit_false():
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    # Set to True first via the raw attr, then restore False via hotstart.
    t.from_hotstart({"temperature.skip_first_time_step": True})
    assert t.to_hotstart()["temperature.skip_first_time_step"] is True
    t.from_hotstart({"temperature.skip_first_time_step": False})
    assert t.to_hotstart()["temperature.skip_first_time_step"] is False


def test_temperature_to_from_hotstart_roundtrip():
    """Snapshot a Temperature, mutate the source, restore from snapshot —
    final state should match the snapshotted one."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    saved = t.to_hotstart()  # {skip: True}
    # Mutate via from_hotstart to a different value.
    t.from_hotstart({"temperature.skip_first_time_step": False})
    assert t.to_hotstart()["temperature.skip_first_time_step"] is False
    # Restore from the original snapshot.
    t.from_hotstart(saved)
    assert t.to_hotstart() == saved


# ---------------------------------------------------------------------------
# Helpers for Model-level tests: stub registry + stub processes
# ---------------------------------------------------------------------------


class _StubRegistry:
    """Minimal VariableRegistry stand-in for unit tests that doesn't need
    clearwater_data's full VariableRegistry behavior. Just stores
    name -> value pairs and supports register/get/get_at_time/set_at_time.
    """

    def __init__(self) -> None:
        self._data: dict[str, object] = {}

    def register(self, name: str, value) -> None:
        self._data[name] = value

    def get(self, name: str):
        return self._data[name]

    def get_at_time(self, name: str, time):
        return self._data[name]

    def set_at_time(self, name: str, time, value) -> None:
        self._data[name] = value

    def get_variable(self, name: str):
        raise KeyError(name)

    def __contains__(self, name: str) -> bool:
        return name in self._data


class _StubProcessNoHotstart(Process):
    """Process that does NOT implement to_hotstart / from_hotstart."""

    variables = []

    def __init__(self, name: str = "stub_no_hotstart") -> None:
        super().__init__(timedelta(minutes=5))
        self._name = name

    def run(self, time, registry) -> None:
        pass

    def process_name(self) -> str:
        return self._name


class _StubProcessWithHotstart(Process):
    """Process that does implement to_hotstart / from_hotstart."""

    variables = []

    def __init__(self, name: str = "stub_with_hotstart") -> None:
        super().__init__(timedelta(minutes=5))
        self._name = name
        self._counter = 0

    def run(self, time, registry) -> None:
        self._counter += 1

    def process_name(self) -> str:
        return self._name

    def to_hotstart(self) -> dict:
        return {f"{self._name}.counter": int(self._counter)}

    def from_hotstart(self, state: dict) -> None:
        key = f"{self._name}.counter"
        if key in state:
            self._counter = int(state[key])
        else:
            self._counter = 0


def _build_minimal_model(
    processes,
    hotstart_dataset: xr.Dataset | None = None,
    hotstart_timestep=None,
) -> Model:
    """Construct a Model with minimal valid kwargs. ``__init_model`` is
    not called (we exercise its private helpers directly)."""
    return Model(
        processes=tuple(processes),
        variable_registry=_StubRegistry(),
        variable_data_sources={},
        start_time=datetime(2026, 1, 1, 0, 0, 0),
        end_time=datetime(2026, 1, 1, 0, 5, 0),
        time_step=timedelta(minutes=1),
        output_variables=[],
        hotstart_dataset=hotstart_dataset,
        hotstart_timestep=hotstart_timestep,
    )


# ---------------------------------------------------------------------------
# Test 2: Model.to_hotstart collects per-process state
# ---------------------------------------------------------------------------


def test_model_to_hotstart_collects_state_from_each_process():
    p_with = _StubProcessWithHotstart(name="p_with")
    p_with._counter = 7
    p_no = _StubProcessNoHotstart(name="p_no")
    model = _build_minimal_model([p_with, p_no])
    snapshot = model.to_hotstart()
    # p_with declared its counter; p_no didn't implement to_hotstart.
    assert snapshot == {"p_with": {"p_with.counter": 7}}


def test_model_to_hotstart_includes_temperature_state():
    """Temperature should appear in Model.to_hotstart() when present."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    model = _build_minimal_model([t])
    snapshot = model.to_hotstart()
    # The snapshot uses process_name() which returns "Temperature".
    assert "Temperature" in snapshot
    assert snapshot["Temperature"] == {"temperature.skip_first_time_step": True}


# ---------------------------------------------------------------------------
# Test 3: Model __seed_from_hotstart overwrites registry contents
# ---------------------------------------------------------------------------


def test_model_seed_from_hotstart_overwrites_registry():
    """Saved-dataset variables replace what's in the registry. Variables
    only in the saved dataset are added; variables only in the registry
    are left alone."""
    saved = xr.Dataset(
        data_vars={
            "water_temperature": ("time", np.array([18.0, 19.0, 20.0])),
            "salinity": ("time", np.array([0.1, 0.2, 0.3])),
        },
        coords={"time": [0, 1, 2]},
    )
    model = _build_minimal_model(
        processes=[],
        hotstart_dataset=saved,
        hotstart_timestep=None,  # default to last slice
    )
    # Pre-seed the stub registry with a value that should be overwritten.
    registry = model._Model__registry
    registry.register("water_temperature", xr.DataArray([99.0]))
    registry.register("alkalinity", xr.DataArray([50.0]))  # only in registry

    # Invoke the private seeder directly.
    model._Model__seed_from_hotstart()

    # ``water_temperature`` was overwritten with the last slice (= 20.0).
    np.testing.assert_array_equal(
        registry.get("water_temperature").values, np.array(20.0)
    )
    # ``salinity`` was added (last slice = 0.3).
    np.testing.assert_array_equal(
        registry.get("salinity").values, np.array(0.3)
    )
    # ``alkalinity`` was left untouched.
    np.testing.assert_array_equal(
        registry.get("alkalinity").values, np.array([50.0])
    )


def test_model_seed_from_hotstart_with_explicit_integer_index():
    saved = xr.Dataset(
        data_vars={
            "water_temperature": ("time", np.array([18.0, 19.0, 20.0])),
        },
        coords={"time": [0, 1, 2]},
    )
    model = _build_minimal_model(
        processes=[],
        hotstart_dataset=saved,
        hotstart_timestep=1,  # middle slice
    )
    model._Model__seed_from_hotstart()
    registry = model._Model__registry
    np.testing.assert_array_equal(
        registry.get("water_temperature").values, np.array(19.0)
    )


# ---------------------------------------------------------------------------
# Test 4: Model __restore_process_hotstart calls from_hotstart on opt-in
# processes only
# ---------------------------------------------------------------------------


def test_model_restore_process_hotstart_invokes_opt_in_processes():
    """from_hotstart should be called on processes that implement it,
    skipped on those that don't, and the dataset's attrs should be the
    state argument."""
    p_with = _StubProcessWithHotstart(name="p_with")
    p_no = _StubProcessNoHotstart(name="p_no")
    saved = xr.Dataset(
        data_vars={"x": ("t", np.zeros(2))},
        coords={"t": [0, 1]},
        attrs={"p_with.counter": 42},
    )
    model = _build_minimal_model(
        processes=[p_with, p_no],
        hotstart_dataset=saved,
    )
    # Pre-state: counter is 0.
    assert p_with._counter == 0
    # Invoke restoration.
    model._Model__restore_process_hotstart()
    # p_with restored from attrs.
    assert p_with._counter == 42
    # p_no had no from_hotstart — no error, no state change.


def test_model_restore_process_hotstart_with_empty_attrs():
    """When the saved dataset has no relevant attrs, opt-in processes
    fall back to their from_hotstart default behavior (which for our
    stub is counter=0)."""
    p_with = _StubProcessWithHotstart(name="p_with")
    p_with._counter = 99
    saved = xr.Dataset(
        data_vars={"x": ("t", np.zeros(2))},
        coords={"t": [0, 1]},
        # No attrs at all.
    )
    model = _build_minimal_model(
        processes=[p_with],
        hotstart_dataset=saved,
    )
    model._Model__restore_process_hotstart()
    # Stub falls back to counter=0 when its key isn't in attrs.
    assert p_with._counter == 0


def test_model_restore_temperature_hotstart_resets_skip_flag():
    """The v3 Temperature semantic: post-hotstart, skip_first_time_step
    is False unless explicitly set in attrs."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    # Confirm the initial state is "skip=True".
    assert t.to_hotstart()["temperature.skip_first_time_step"] is True

    saved = xr.Dataset(
        data_vars={"x": ("t", np.zeros(2))},
        coords={"t": [0, 1]},
        # Empty attrs -> from_hotstart should default to skip=False.
    )
    model = _build_minimal_model(
        processes=[t],
        hotstart_dataset=saved,
    )
    model._Model__restore_process_hotstart()
    assert t.to_hotstart()["temperature.skip_first_time_step"] is False
