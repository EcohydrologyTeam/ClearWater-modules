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


# ---------------------------------------------------------------------------
# Full NSM1 state set: the map covers every NSM1 state under its fork name.
# ---------------------------------------------------------------------------

# The full 16-entry fork-name -> canonical mapping the bridge installs. Pinned
# here as an independent copy so the test fails loudly if the production map in
# Riverine._MESH_TO_CANONICAL is changed without updating the contract.
_EXPECTED_MAP = {
    "Ap": "algae_floating",
    "Ab": "benthic_algae",
    "NH4": "ammonium",
    "NO3": "nitrate",
    "OrgN": "organic_nitrogen",
    "N2": "n2",
    "TIP": "tip",
    "OrgP": "organic_phosphorus",
    "POC": "poc",
    "DOC": "doc",
    "DIC": "dic",
    "CBOD": "cbod",
    "POM": "pom",
    "DOX": "oxygen_dissolved",
    "Alk": "alkalinity",
    "PX": "pathogen",
}


def test_mesh_to_canonical_map_contract():
    """The production map matches the pinned full NSM1 contract (no env)."""
    prod = Riverine._MESH_TO_CANONICAL
    assert prod == _EXPECTED_MAP
    # No duplicate fork names or canonical targets.
    assert len(prod) == 16
    assert len(set(prod.values())) == 16
    # Every entry is a non-empty str -> non-empty str.
    for fork, canonical in prod.items():
        assert isinstance(fork, str) and fork
        assert isinstance(canonical, str) and canonical
    # DOX is retained (it was in the original five and must not be dropped).
    assert prod["DOX"] == "oxygen_dissolved"


def test_init_process_bridges_full_nsm1_state_set():
    """A config carrying all 16 fork constituents bridges every canonical
    alias to the registry, each tracking its source mesh array."""
    inst, registry = _build_real_riverine(constituents=list(_EXPECTED_MAP))
    _init(inst, registry)

    for fork, canonical in _EXPECTED_MAP.items():
        assert fork in inst.mesh, f"{fork!r} not on mesh"
        assert canonical in registry, f"{canonical!r} not bridged"
        np.testing.assert_array_equal(
            _read(registry, canonical), np.asarray(inst.mesh[fork])
        )
    assert "depth" in registry


# ---------------------------------------------------------------------------
# Shared-buffer: the bridge points the canonical alias at the SAME buffer as
# the fork-named mesh array (copy(deep=False)), so a write to one side is seen
# on the other. This is what makes the coupling two-way.
# ---------------------------------------------------------------------------


def test_bridge_is_shared_buffer_two_way():
    inst, registry = _build_real_riverine(constituents=list(_EXPECTED_MAP))
    _init(inst, registry)

    # Write via the canonical registry alias -> visible on the fork mesh name.
    canon = registry.get_variable("ammonium").get()
    sentinel_a = np.arange(canon.values.size, dtype=float).reshape(
        canon.values.shape
    )
    canon.values[:] = sentinel_a
    np.testing.assert_array_equal(np.asarray(inst.mesh["NH4"]), sentinel_a)

    # Write via the fork mesh name -> visible on the canonical registry alias.
    mesh_arr = inst.mesh["nitrate"] if "nitrate" in inst.mesh else inst.mesh["NO3"]
    sentinel_b = sentinel_a + 100.0
    mesh_arr.values[:] = sentinel_b
    np.testing.assert_array_equal(_read(registry, "nitrate"), sentinel_b)


# ---------------------------------------------------------------------------
# Per-constituent two-way vs one-way coupling.
#
# A two-way constituent shares the mesh buffer (kinetics writes feed back to
# transport). A one-way constituent (two_way_coupling: false in the riverine
# config) is bridged as an isolated snapshot: transport feeds kinetics, but a
# kinetics write does NOT propagate back to the transport mesh.
# ---------------------------------------------------------------------------


def _build_real_riverine_with_flags(consts_cfg: dict):
    """Build a real ClearwaterRiverine where each constituent's config is
    given explicitly (so two_way_coupling can be set per constituent).

    Mirrors ``_build_real_riverine`` but takes a full constituents-config
    mapping rather than a name list. Returns ``(instance, registry)``.
    """
    start, end = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    reg = VariableRegistry()
    model_cfg = {
        "simulation_directory": str(tempfile.mkdtemp()),
        "hydrodynamic_input": str((PLAN02 / PLAN02_HDF).resolve()),
        "start_datetime": str(start),
        "end_datetime": str(end),
        "diffusion_coefficient": 0.01,
        "output_variables": [],
        "mass_flux_calculation": False,
    }
    cfg = {"model": model_cfg, "constituents": consts_cfg}
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


def _const_block(two_way=None):
    block = {
        "initial_conditions": {"provider": "float", "data": {"value": 1.0}},
        "boundary_conditions": {"provider": "float", "data": {"value": 1.0}},
    }
    if two_way is not None:
        block["two_way_coupling"] = two_way
    return block


def test_one_way_constituent_does_not_write_back_to_mesh():
    # NH4 is one-way; NO3 is two-way (explicit).
    inst, registry = _build_real_riverine_with_flags(
        {
            "NH4": _const_block(two_way=False),
            "NO3": _const_block(two_way=True),
        }
    )
    _init(inst, registry)

    # Transport -> kinetics: both aliases start equal to their mesh source.
    np.testing.assert_array_equal(
        _read(registry, "ammonium"), np.asarray(inst.mesh["NH4"])
    )

    # A kinetics write to the ONE-WAY alias must NOT reach the mesh.
    mesh_before = np.asarray(inst.mesh["NH4"]).copy()
    canon = registry.get_variable("ammonium").get()
    canon.values[:] = np.nan_to_num(canon.values, nan=0.0) + 999.0
    np.testing.assert_array_equal(np.asarray(inst.mesh["NH4"]), mesh_before)

    # The TWO-WAY alias in the same model still writes back (shared buffer).
    canon_no3 = registry.get_variable("nitrate").get()
    sentinel = np.arange(canon_no3.values.size, dtype=float).reshape(
        canon_no3.values.shape
    )
    canon_no3.values[:] = sentinel
    np.testing.assert_array_equal(np.asarray(inst.mesh["NO3"]), sentinel)


def test_one_way_snapshot_reseeds_from_transport_on_rebridge():
    # A one-way constituent's snapshot is refreshed from transport each
    # re-bridge: a kinetics write is discarded, and the alias tracks the
    # mesh again after the next _bridge_mesh_to_registry call.
    inst, registry = _build_real_riverine_with_flags(
        {"NH4": _const_block(two_way=False)}
    )
    process = _init(inst, registry)

    canon = registry.get_variable("ammonium").get()
    canon.values[:] = np.nan_to_num(canon.values, nan=0.0) + 999.0

    # Re-bridge (no transport advance needed): the snapshot is replaced by a
    # fresh deep copy of the current mesh value, discarding the write.
    process._bridge_mesh_to_registry(registry)
    np.testing.assert_array_equal(
        _read(registry, "ammonium"), np.asarray(inst.mesh["NH4"])
    )


def test_coupling_flags_default_two_way_when_accessor_absent():
    # Back-compat: an instance without constituent_coupling() yields {} and
    # the bridge treats every constituent as two-way.
    inst, registry = _build_real_riverine(constituents=["Ap", "NH4"])
    process = _init(inst, registry)

    # A current instance exposes the accessor, so flags reflect the config
    # (every constituent defaults to two-way True).
    assert process._coupling_flags() == {"Ap": True, "NH4": True}

    # Simulate an older instance lacking the accessor.
    class _NoAccessor:
        def __init__(self, wrapped):
            self._w = wrapped

        def __getattr__(self, name):
            if name == "constituent_coupling":
                raise AttributeError(name)
            return getattr(self._w, name)

    process.riverine_instance = _NoAccessor(inst)
    assert process._coupling_flags() == {}
    # And the bridge still works (all two-way / shared buffer).
    process._bridge_mesh_to_registry(registry)
    canon = registry.get_variable("ammonium").get()
    sentinel = np.arange(canon.values.size, dtype=float).reshape(
        canon.values.shape
    )
    canon.values[:] = sentinel
    np.testing.assert_array_equal(np.asarray(inst.mesh["NH4"]), sentinel)
