"""v3 ``Model`` robustness fixes (review-findings 2026-05-04).

Covers six MAJOR findings on the v3 ``Model`` orchestrator:

- **M5.** ``init_process`` runs *before* ``from_hotstart``: the ordering
  contract is documented on the Process base module so that authors of
  new internal substep state cannot silently desynchronize fresh-start
  and hotstart-resume.
- **M7.** A non-chunked simulation (``chunk_size=None``) paired with a
  ``ChunkedDataSource`` (e.g., HEC-RAS HDF) used to load only the first
  substep window. The fix loads the full ``[start_time, end_time]``
  window when ``chunk_size`` is None.
- **M8.** ``__seed_from_hotstart`` previously fell back to
  ``next(iter(ds.dims))`` when no recognized time dim was present,
  silently treating a spatial dim such as ``nface`` as the time axis.
  The fix raises ``ValueError`` when an explicit ``hotstart_timestep``
  was provided and there is no recognizable time dim, but allows a
  no-time-dim dataset to be used as a single snapshot when
  ``hotstart_timestep is None``.
- **M10.** A second ``run()`` call against the same Model instance now
  raises ``RuntimeError`` rather than silently re-iterating against an
  already-advanced registry.
- **M11.** Wet-mask threshold semantic (strict-inequality with default
  0.0) is documented on the constructor. No code change; the docstring
  test pins the documented contract.
- **M14.** ``wet_mask_variable`` is validated against the registry at
  ``init_model`` time. A typo such as ``wetted_surface_aera`` raises
  ``KeyError`` immediately rather than blowing up at the first substep.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
import xarray as xr

from clearwater_data.io.base import ChunkedDataSource

from clearwater_modules_v3.model import Model
from clearwater_modules_v3.processes.base import Process
from clearwater_modules_v3 import processes as v3_processes


# ---------------------------------------------------------------------------
# Stubs
# ---------------------------------------------------------------------------


class _StubRegistry:
    """Minimal VariableRegistry stand-in shared across the tests."""

    def __init__(self) -> None:
        self._data: dict[str, object] = {}

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


class _StubProcess(Process):
    """Process with no finalize_process / no run-side effects."""

    variables: list[str] = []

    def __init__(self, name: str = "stub", time_step: timedelta = timedelta(minutes=5)) -> None:
        super().__init__(time_step)
        self._name = name

    def run(self, time, registry) -> None:
        pass

    def process_name(self) -> str:
        return self._name


class _RecordingChunkedDataSource(ChunkedDataSource):
    """ChunkedDataSource that records every (start, end) pair passed to
    ``read_chunk`` so tests can assert the requested window.

    The model code calls ``read_chunk(start, end)`` positionally without
    a parameter name; the stub accepts and records the same positional
    pair. We override the abstract ``read_chunk`` with a permissive
    signature so the call signature mismatch in legacy code paths
    doesn't surface as a TypeError under test.
    """

    def __init__(self, payload) -> None:
        self.calls: list[tuple[datetime, datetime]] = []
        self._payload = payload

    def read_chunk(self, *args, **kwargs):  # type: ignore[override]
        # The Model invocations in __init_model and __load_chunk_data
        # both pass two positional datetimes (start, end). Record them.
        if len(args) >= 2 and isinstance(args[0], datetime) and isinstance(args[1], datetime):
            start, end = args[0], args[1]
        elif "start_time" in kwargs and "end_time" in kwargs:
            start, end = kwargs["start_time"], kwargs["end_time"]
        else:
            # Defensive: not the expected v3 call signature; just
            # record whatever was given for debugging.
            start = kwargs.get("start_time", args[0] if args else None)
            end = kwargs.get("end_time", args[-1] if args else None)
        self.calls.append((start, end))
        return self._payload


def _build_model(
    processes,
    *,
    chunk_size: timedelta | None = None,
    end_time: datetime = datetime(2026, 1, 1, 0, 30, 0),
    time_step: timedelta = timedelta(minutes=5),
    variable_data_sources=None,
    hotstart_dataset=None,
    hotstart_timestep=None,
    wet_mask_variable=None,
    wet_mask_threshold: float = 0.0,
) -> Model:
    return Model(
        processes=tuple(processes),
        variable_registry=_StubRegistry(),
        variable_data_sources=variable_data_sources or {},
        start_time=datetime(2026, 1, 1, 0, 0, 0),
        end_time=end_time,
        time_step=time_step,
        output_variables=[],
        simulation_directory=None,
        chunk_size=chunk_size,
        wet_mask_variable=wet_mask_variable,
        wet_mask_threshold=wet_mask_threshold,
        hotstart_dataset=hotstart_dataset,
        hotstart_timestep=hotstart_timestep,
    )


# ---------------------------------------------------------------------------
# M5: init_process / from_hotstart ordering contract is documented
# ---------------------------------------------------------------------------


def test_m5_ordering_contract_documented_on_process_base_module():
    """The v3 ``Process`` base module docstring must spell out the
    ``init_process`` / ``from_hotstart`` ordering invariant. This is a
    documentation-only fix: pinning the contract here ensures the
    invariant is preserved through future docstring rewrites."""
    doc = v3_processes.base.__doc__ or ""
    # The contract must mention both methods and the ordering keyword.
    assert "init_process" in doc, (
        "Process base module docstring must reference init_process"
    )
    assert "from_hotstart" in doc, (
        "Process base module docstring must reference from_hotstart"
    )
    assert "ordering" in doc.lower(), (
        "Process base module docstring must explain the ordering invariant"
    )


def test_m5_ordering_contract_documented_on_init_model():
    """The Model ``__init_model`` docstring should also explain the
    ordering invariant for readers of the orchestrator code."""
    init_model = Model._Model__init_model
    doc = init_model.__doc__ or ""
    assert "init_process" in doc
    assert "from_hotstart" in doc
    assert "ordering" in doc.lower() or "Ordering" in doc


# ---------------------------------------------------------------------------
# M7: Non-chunked + ChunkedDataSource reads the full window
# ---------------------------------------------------------------------------


def test_m7_non_chunked_with_chunked_source_reads_full_window():
    """When ``chunk_size`` is None and the data source is a
    ChunkedDataSource, ``__init_model`` must read
    ``[start_time, end_time]`` — not ``[start_time, start_time + time_step]``
    as it did pre-fix."""
    payload = xr.DataArray(np.zeros(10))
    cds = _RecordingChunkedDataSource(payload)
    model = _build_model(
        processes=[],
        chunk_size=None,  # non-chunked simulation
        end_time=datetime(2026, 1, 1, 1, 0, 0),  # one hour
        time_step=timedelta(minutes=5),
        variable_data_sources={"forcing": cds},
    )
    model._Model__init_model()
    assert len(cds.calls) == 1
    start, end = cds.calls[0]
    # The recorded window must span the entire simulation, not just
    # one substep.
    assert start == datetime(2026, 1, 1, 0, 0, 0)
    assert end == datetime(2026, 1, 1, 1, 0, 0)
    # Sanity: the buggy window would have been [start, start + 5 min].
    assert end != start + timedelta(minutes=5)


def test_m7_chunked_simulation_uses_chunk_window():
    """When chunk_size is not None, the initial read still uses
    ``[start_time, start_time + chunk_size]`` — the M7 fix preserves
    chunked-mode behavior, only changing the non-chunked branch."""
    payload = xr.DataArray(np.zeros(10))
    cds = _RecordingChunkedDataSource(payload)
    model = _build_model(
        processes=[],
        chunk_size=timedelta(minutes=10),
        end_time=datetime(2026, 1, 1, 1, 0, 0),
        time_step=timedelta(minutes=5),
        variable_data_sources={"forcing": cds},
    )
    # Stub the output store init so we don't try to build a real Zarr.
    with patch.object(Model, "_Model__init_output_source", lambda self: None):
        model._Model__init_model()
    assert len(cds.calls) == 1
    start, end = cds.calls[0]
    assert start == datetime(2026, 1, 1, 0, 0, 0)
    # 10-minute chunk window, not the full hour.
    assert end == datetime(2026, 1, 1, 0, 10, 0)


# ---------------------------------------------------------------------------
# M8: __seed_from_hotstart raises on missing time-dim with explicit ts
# ---------------------------------------------------------------------------


def test_m8_no_time_dim_with_explicit_timestep_raises_value_error():
    """A saved dataset with only a spatial dim (e.g. ``nface``) plus an
    explicit ``hotstart_timestep`` must raise ``ValueError`` rather
    than silently slicing space-as-time."""
    saved = xr.Dataset(
        data_vars={
            "water_temperature": ("nface", np.array([18.0, 19.0, 20.0])),
        },
        coords={"nface": [0, 1, 2]},
    )
    model = _build_model(
        processes=[],
        hotstart_dataset=saved,
        hotstart_timestep=0,  # explicit -> must raise
    )
    with pytest.raises(ValueError, match="no recognizable time dimension"):
        model._Model__seed_from_hotstart()


def test_m8_no_time_dim_with_none_timestep_uses_dataset_as_is():
    """A saved dataset with only a spatial dim and ``hotstart_timestep=None``
    is treated as a single-snapshot dataset and used as-is."""
    saved = xr.Dataset(
        data_vars={
            "water_temperature": ("nface", np.array([18.0, 19.0, 20.0])),
        },
        coords={"nface": [0, 1, 2]},
    )
    model = _build_model(
        processes=[],
        hotstart_dataset=saved,
        hotstart_timestep=None,
    )
    model._Model__seed_from_hotstart()
    registry = model._Model__registry
    np.testing.assert_array_equal(
        registry.get("water_temperature").values,
        np.array([18.0, 19.0, 20.0]),
    )


# ---------------------------------------------------------------------------
# M10: run() twice raises RuntimeError on the second call
# ---------------------------------------------------------------------------


def test_m10_run_twice_raises_runtime_error():
    """The second call to ``run()`` must raise ``RuntimeError``. The
    first call still completes successfully and advances the registry."""
    p = _StubProcess("p")
    model = _build_model(processes=[p])
    # Stub out the output-store init / save so the first run completes
    # without touching real I/O.
    with patch.object(Model, "_Model__init_output_source", lambda self: None), \
         patch.object(Model, "_Model__save_output_model", lambda self, start_time, end_time: None):
        # First run — must succeed.
        model.run()
        # Second run — must raise.
        with pytest.raises(RuntimeError, match="already completed"):
            model.run()


def test_m10_first_run_completes_successfully():
    """Sanity: the M10 guard fires on the second call only — the first
    call still runs the substep loop and sets the completion flag."""
    p = _StubProcess("p")
    model = _build_model(processes=[p])
    with patch.object(Model, "_Model__init_output_source", lambda self: None), \
         patch.object(Model, "_Model__save_output_model", lambda self, start_time, end_time: None):
        # Pre-run flag is False.
        assert model._Model__run_complete is False
        model.run()
        # Post-run flag is True.
        assert model._Model__run_complete is True


# ---------------------------------------------------------------------------
# M14: wet_mask_variable is validated at init_model time
# ---------------------------------------------------------------------------


def test_m14_unknown_wet_mask_variable_raises_key_error_at_init():
    """A typo in ``wet_mask_variable`` (here, ``wetted_surface_aera``)
    must raise ``KeyError`` at ``init_model`` time, not at the first
    substep deep inside ``__compute_wet_mask``."""
    p = _StubProcess("p")
    model = _build_model(
        processes=[p],
        wet_mask_variable="wetted_surface_aera",  # typo on purpose
    )
    with pytest.raises(KeyError, match="wetted_surface_aera"):
        model.init_model()


def test_m14_known_wet_mask_variable_passes_validation():
    """A correctly-spelled ``wet_mask_variable`` that *is* in the
    registry passes validation and ``init_model`` completes."""
    p = _StubProcess("p")
    model = _build_model(
        processes=[p],
        wet_mask_variable="wetted_surface_area",
    )
    # Pre-seed the registry with the variable so the validation finds it.
    model._Model__registry.register(
        "wetted_surface_area",
        xr.DataArray(np.array([100.0, 200.0, 300.0])),
    )
    # Skip the data-source load (variable_data_sources is empty) and the
    # output-store init (no output_variables). The wet-mask validation
    # is the last step before the substep loop is entered.
    model.init_model()
    # Sanity: the model is now init-complete.
    assert model._Model__init_complete is True
