"""Integration tests for the v3 ``Riverine`` process MeshView bridge.

These tests build a REAL ``ClearwaterRiverine`` over a bundled HEC-RAS
plan (``plan02_2x1`` in the ClearWater-riverine repo) with the five NSM
constituents (Ap/NH4/NO3/TIP/DOX) and drive the modules bridge against
the resulting transport mesh -- no stub ``riverine_instance``/``MeshView``.

Covers the chunk-safe re-bridge added in the riverine MeshView-compat
change (``design/clearwater_modules_v3_riverine_process_meshview_compat.md``):

- MeshView item access (``mesh["Ap"]``) replaces attribute access.
- The constituent gate is widened: no ``has_process("FloatingAlgae")``
  guard; present constituents bridge unconditionally.
- The inorganic-P name is reconciled to ``tip`` (was
  ``phosphorus_total_inorganic``).
- ``depth`` is bridged from ``mesh["coupling_depth"]`` (the resolved cell
  mean water-column depth, a length), not ``wetted_surface_area`` (area).
  ``init_process`` first calls
  ``riverine_instance.enable_coupling_depth()`` to turn that resolved
  depth on for the coupled run.
- The bridge is re-applied each substep so chunk reloads (which
  re-register FRESH DataArrays) are picked up rather than stranded on the
  previous chunk's buffers.

REQUIRED INVOCATION (the conda ``clearwater`` env has a zarr-3-incompatible
xarray; the riverine pixi ``dev`` env has a working one, and PYTHONPATH
adds the modules source so no install is needed). Run from the modules
repo dir::

    PYTHONPATH=/Users/todd/GitHub/ecohydrology/ClearWater-modules/src \\
      /Users/todd/GitHub/ecohydrology/ClearWater-riverine/.pixi/envs/dev/bin/python \\
      -m pytest tests/v3/test_riverine_process_v3.py -q

Plain ``python``/``conda`` fails on the chunked zarr path with
``'Float64' object has no attribute 'value'``.
"""

from __future__ import annotations

import tempfile
from datetime import timedelta
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import xarray as xr
import yaml

import clearwater_riverine as cwr
from clearwater_data.variables import DataArrayVariable, VariableRegistry
from clearwater_riverine.fork_compat import MeshView
from clearwater_riverine.variables import CHANGE_IN_TIME

from clearwater_modules_v3.processes.riverine import Riverine


# --- Real-fixture plumbing -------------------------------------------------
# The bundled plan lives in the ClearWater-riverine repo, a sibling of the
# modules repo. Resolve it relative to this file so the test is location
# independent.
_RIVERINE_REPO = (
    Path(__file__).resolve().parents[3] / "ClearWater-riverine"
)
PLAN02 = _RIVERINE_REPO / "tests" / "data" / "simple_test_cases" / "plan02_2x1"
PLAN02_HDF = "clearWaterTestCases.p02.hdf"

_RAS_TIME_PATH = (
    "Results/Unsteady/Output/Output Blocks/Base Output/"
    "Unsteady Time Series/Time Date Stamp"
)

# The five NSM constituents the modules bridge maps to canonical names.
_CONSTITUENTS = ["Ap", "NH4", "NO3", "TIP", "DOX"]


pytestmark = pytest.mark.skipif(
    not (PLAN02 / PLAN02_HDF).exists(),
    reason=(
        "ClearWater-riverine plan02 fixture not found at "
        f"{PLAN02 / PLAN02_HDF}; sibling repo checkout required"
    ),
)


def _hdf_time_bounds(hdf_path: Path):
    with h5py.File(hdf_path, "r") as f:
        raw = f[_RAS_TIME_PATH][()]
    stamps = pd.to_datetime(
        pd.Series(raw).str.decode("utf8"), format="%d%b%Y %H:%M:%S"
    )
    return stamps.iloc[0], stamps.iloc[-1]


def _build_real_riverine(
    constituents: list[str] | None = None,
    *,
    chunk_size: str | None = None,
):
    """Build a real ``ClearwaterRiverine`` over plan02 with float IC/BC.

    Returns ``(instance, registry)``. ``water_temperature`` is seeded on
    the shared registry (the modules side requires it) at a constant 15 C.
    The five (or a subset of) NSM constituents land in ``inst.mesh`` and
    the shared registry. When ``chunk_size`` is given (a ``pd.Timedelta``
    string), the instance runs in chunked mode.
    """
    if constituents is None:
        constituents = _CONSTITUENTS
    start, end = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    reg = VariableRegistry()
    consts = {
        c: {
            "initial_conditions": {"provider": "float", "data": {"value": 1.0}},
            "boundary_conditions": {"provider": "float", "data": {"value": 1.0}},
        }
        for c in constituents
    }
    model_cfg = {
        # A scratch sim dir keeps zarr output out of the checked-in fixture.
        "simulation_directory": str(Path(tempfile.mkdtemp())),
        "hydrodynamic_input": str((PLAN02 / PLAN02_HDF).resolve()),
        "start_datetime": str(start),
        "end_datetime": str(end),
        "diffusion_coefficient": 0.01,
        "output_variables": [],
        "mass_flux_calculation": False,
    }
    if chunk_size is not None:
        model_cfg["chunk_size"] = chunk_size
    cfg = {"model": model_cfg, "constituents": consts}
    cfg_path = Path(tempfile.mkdtemp()) / "riv.yml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))

    inst = cwr.ClearwaterRiverine(
        config_filepath=str(cfg_path), variable_registry=reg
    )
    reg.register(
        "water_temperature",
        DataArrayVariable(xr.full_like(reg.get_variable("volume").get_data(), 15.0)),
    )
    return inst, reg


def _even_chunk_size():
    """Derive a chunk_size giving an exact, even >=2-chunk split of plan02.

    Mirrors ``tests/test_coupling_depth.py::_even_chunk_size`` in the
    riverine repo: probe a non-chunked model for its uniform timestep,
    then pick a chunk window that splits the step count into 2 (or 3)
    equal chunks with >=2 slots each. Returns ``(chunk_size_str, dt_s,
    n_steps)`` or ``(None, None, None)``.
    """
    probe, _ = _build_real_riverine()
    dt_s = float(probe.registry.get_variable(CHANGE_IN_TIME).get_data())
    start, end = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    n_steps = round((end - start).total_seconds() / dt_s)
    m = next((k for k in (2, 3) if n_steps % k == 0 and n_steps // k >= 2), None)
    if m is None:
        return None, None, None
    return str(pd.Timedelta(seconds=dt_s) * (n_steps // m)), dt_s, n_steps


class _DummyModel:
    """Stub model. ``init_process`` keeps the ``model`` param for interface
    compatibility but no longer reads it; this just satisfies the call."""

    def has_process(self, process_type) -> bool:  # pragma: no cover - defensive
        return False


def _read(registry: VariableRegistry, name: str) -> np.ndarray:
    return np.asarray(registry.get_variable(name).get_data())


def _init(inst, registry, time_step=timedelta(seconds=30)) -> Riverine:
    process = Riverine(inst, time_step=time_step)
    process.init_process(_DummyModel(), registry)
    return process


# ---------------------------------------------------------------------------
# Full constituent set bridges to canonical names with correct depth.
# ---------------------------------------------------------------------------


def test_init_process_bridges_all_constituents_and_depth():
    inst, registry = _build_real_riverine()
    process = _init(inst, registry)

    # Canonical aliases registered.
    for canonical in (
        "algae_floating",
        "ammonium",
        "nitrate",
        "tip",
        "oxygen_dissolved",
        "depth",
    ):
        assert canonical in registry, f"{canonical!r} not bridged"

    # Inorganic-P reconciled to ``tip``; old name absent.
    assert "tip" in registry
    assert "phosphorus_total_inorganic" not in registry

    # Values track the source mesh constituents (NaN-aware: dry cells are
    # NaN in both, and assert_array_equal treats NaN == NaN as equal).
    np.testing.assert_array_equal(
        _read(registry, "algae_floating"), np.asarray(inst.mesh["Ap"])
    )
    np.testing.assert_array_equal(
        _read(registry, "tip"), np.asarray(inst.mesh["TIP"])
    )

    # depth is sourced from coupling_depth (length), NOT wetted_surface_area
    # (area). plan02 lacks RAS Cell Hydraulic Depth so coupling_depth is
    # resolved on the riverine side as volume / wetted_surface_area.
    np.testing.assert_array_equal(
        _read(registry, "depth"), np.asarray(inst.mesh["coupling_depth"])
    )
    assert not np.allclose(
        np.nan_to_num(_read(registry, "depth"), nan=-1.0),
        np.nan_to_num(_read(registry, "wetted_surface_area"), nan=-2.0),
    )


# ---------------------------------------------------------------------------
# Subset: only constituents present in the mesh bridge.
# ---------------------------------------------------------------------------


def test_init_process_bridges_only_present_constituents():
    # A config declaring only some of the five constituents (Ap omitted).
    inst, registry = _build_real_riverine(
        constituents=["NH4", "NO3", "TIP", "DOX"]
    )
    process = _init(inst, registry)

    assert "Ap" not in inst.mesh  # not declared -> not on the mesh
    assert "algae_floating" not in registry  # Ap absent -> no bridge
    assert "ammonium" in registry
    assert "nitrate" in registry
    assert "tip" in registry
    assert "oxygen_dissolved" in registry
    assert "depth" in registry

    np.testing.assert_array_equal(
        _read(registry, "ammonium"), np.asarray(inst.mesh["NH4"])
    )
    np.testing.assert_array_equal(
        _read(registry, "depth"), np.asarray(inst.mesh["coupling_depth"])
    )


# ---------------------------------------------------------------------------
# Missing depth: fail loud.
#
# Pure error-path check: a real instance cannot easily be made to fail to
# resolve depth, so a MINIMAL stub instance (mesh lacks coupling_depth,
# enable_coupling_depth is a no-op) exercises the guard.
# ---------------------------------------------------------------------------


def _da(values: list[float]) -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=float), dims=["nface"])


class _NoDepthStubInstance:
    """Minimal stand-in whose mesh lacks ``coupling_depth`` and whose
    ``enable_coupling_depth()`` is a no-op, so the bridge must raise."""

    def __init__(self, registry: VariableRegistry) -> None:
        self.mesh = MeshView(registry)
        self.is_chunked = False

    def enable_coupling_depth(self) -> None:  # no-op: depth never appears
        pass


def test_init_process_raises_when_coupling_depth_missing():
    registry = VariableRegistry()
    for name in ("Ap", "NH4", "NO3", "TIP", "DOX"):
        registry.register(name, DataArrayVariable(_da([1.0, 1.1])))
    registry.register("volume", DataArrayVariable(_da([100.0, 110.0])))
    registry.register("water_temperature", DataArrayVariable(_da([20.0, 21.0])))
    registry.register("wetted_surface_area", DataArrayVariable(_da([200.0, 220.0])))

    process = Riverine(_NoDepthStubInstance(registry))

    with pytest.raises(KeyError, match="coupling_depth"):
        process.init_process(_DummyModel(), registry)


# ---------------------------------------------------------------------------
# Chunk-safety: the key regression test for the stale-after-chunk-1 bug,
# driven against a REAL chunked instance across a real chunk boundary.
# ---------------------------------------------------------------------------


def test_rebridge_picks_up_new_chunk_objects():
    chunk_size, dt_s, n_steps = _even_chunk_size()
    if chunk_size is None:
        pytest.skip("plan02 step count has no exact >=2-chunk split")

    inst, registry = _build_real_riverine(chunk_size=chunk_size)
    assert inst.is_chunked
    process = _init(inst, registry, time_step=timedelta(seconds=dt_s))

    # Seed t0 (chunk-1): canonical aliases track the first chunk's mesh.
    np.testing.assert_array_equal(
        _read(registry, "algae_floating"), np.asarray(inst.mesh["Ap"])
    )
    np.testing.assert_array_equal(
        _read(registry, "depth"), np.asarray(inst.mesh["coupling_depth"])
    )
    chunk1_id = id(registry.get_variable("Ap"))

    # Drive across the chunk boundary. update() reloads the next chunk,
    # re-registering FRESH DataArrays; run()'s re-bridge must follow them.
    start, _ = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    t = pd.Timestamp(start)
    crossed = False
    for _ in range(n_steps + 2):
        t = t + pd.Timedelta(seconds=dt_s)
        process.run(t, registry)
        if id(registry.get_variable("Ap")) != chunk1_id:
            crossed = True
            break
    assert crossed, "never crossed a chunk boundary"

    # After the reload the canonical aliases track the CURRENT chunk's
    # mesh, not the stranded chunk-1 buffers.
    np.testing.assert_array_equal(
        _read(registry, "algae_floating"), np.asarray(inst.mesh["Ap"])
    )
    np.testing.assert_array_equal(
        _read(registry, "depth"), np.asarray(inst.mesh["coupling_depth"])
    )
