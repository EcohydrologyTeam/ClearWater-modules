"""Live chunked-coupling integration test for the v3 ``Riverine`` bridge.

Proves the chunk-safe bridge against a REAL multi-chunk transport run:
a real ``ClearwaterRiverine`` is built over the bundled ``plan02_2x1``
HEC-RAS plan with the five NSM constituents and a ``chunk_size`` that
splits the run into two chunks. The modules ``Riverine`` process is
driven across the chunk boundary; each ``run()`` calls ``inst.update()``,
which triggers ``__load_new_chunk`` at the boundary, re-registering fresh
DataArrays. After the boundary the registry's canonical aliases must
track the riverine mesh's CURRENT values -- the regression guard for the
stale-after-chunk-1 bug this whole change fixes.

A non-chunked reference instance (same config, no ``chunk_size``) is
driven in lockstep and must yield the same bridged constituent/depth
values at matching timestamps.

REQUIRED INVOCATION (the conda ``clearwater`` env has a zarr-3-incompatible
xarray; the riverine pixi ``dev`` env has a working one, and PYTHONPATH
adds the modules source). Run from the modules repo dir::

    PYTHONPATH=/Users/todd/GitHub/ecohydrology/ClearWater-modules/src \\
      /Users/todd/GitHub/ecohydrology/ClearWater-riverine/.pixi/envs/dev/bin/python \\
      -m pytest tests/v3/test_riverine_chunked_integration_v3.py -q

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
from clearwater_riverine.variables import CHANGE_IN_TIME, VOLUME

from clearwater_modules_v3.processes.riverine import Riverine


_RIVERINE_REPO = Path(__file__).resolve().parents[3] / "ClearWater-riverine"
PLAN02 = _RIVERINE_REPO / "tests" / "data" / "simple_test_cases" / "plan02_2x1"
PLAN02_HDF = "clearWaterTestCases.p02.hdf"

_RAS_TIME_PATH = (
    "Results/Unsteady/Output/Output Blocks/Base Output/"
    "Unsteady Time Series/Time Date Stamp"
)

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


def _build_real_riverine(*, chunk_size: str | None = None):
    start, end = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    reg = VariableRegistry()
    consts = {
        c: {
            "initial_conditions": {"provider": "float", "data": {"value": 1.0}},
            "boundary_conditions": {"provider": "float", "data": {"value": 1.0}},
        }
        for c in _CONSTITUENTS
    }
    model_cfg = {
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
    """``(chunk_size_str, dt_s, n_steps)`` for an even >=2-chunk split,
    or ``(None, None, None)``. Mirrors the riverine repo helper."""
    probe, _ = _build_real_riverine()
    dt_s = float(probe.registry.get_variable(CHANGE_IN_TIME).get_data())
    start, end = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    n_steps = round((end - start).total_seconds() / dt_s)
    m = next((k for k in (2, 3) if n_steps % k == 0 and n_steps // k >= 2), None)
    if m is None:
        return None, None, None
    return str(pd.Timedelta(seconds=dt_s) * (n_steps // m)), dt_s, n_steps


class _DummyModel:
    def has_process(self, process_type) -> bool:  # pragma: no cover - defensive
        return False


def _init(inst, registry, dt_s) -> Riverine:
    process = Riverine(inst, time_step=timedelta(seconds=dt_s))
    process.init_process(_DummyModel(), registry)
    return process


def _data(registry: VariableRegistry, name: str) -> xr.DataArray:
    return registry.get_variable(name).get_data()


def _vol_window(inst):
    v = inst.registry.get_variable(VOLUME).get_data()
    return (v.time.values[0], v.time.values[-1])


def _equal_nan(a: np.ndarray, b: np.ndarray) -> bool:
    """NaN-aware equality: dry cells are NaN in both meshes."""
    return np.array_equal(np.nan_to_num(a, nan=-1.0), np.nan_to_num(b, nan=-1.0))


def test_chunked_bridge_tracks_current_chunk_and_matches_reference():
    chunk_size, dt_s, n_steps = _even_chunk_size()
    if chunk_size is None:
        pytest.skip("plan02 step count has no exact >=2-chunk split")

    # 1. Build a chunked real instance + the modules process; seed t0.
    inst_c, reg_c = _build_real_riverine(chunk_size=chunk_size)
    assert inst_c.is_chunked
    proc_c = _init(inst_c, reg_c, dt_s)

    # A non-chunked reference, driven in lockstep for the cross-check.
    inst_n, reg_n = _build_real_riverine()
    assert not inst_n.is_chunked
    proc_n = _init(inst_n, reg_n, dt_s)

    # 2. Drive both processes across at least one chunk boundary. Each
    #    run() calls inst.update(); at the boundary the chunked instance
    #    loads a new chunk (its VOLUME time window shifts).
    start, _ = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    t = pd.Timestamp(start)
    boundary_window = _vol_window(inst_c)
    crossed = False
    steps_after_boundary = 0
    for _ in range(n_steps + 2):
        t = t + pd.Timedelta(seconds=dt_s)
        proc_c.run(t, reg_c)
        proc_n.run(t, reg_n)
        if not crossed:
            if _vol_window(inst_c) != boundary_window:
                crossed = True
        else:
            steps_after_boundary += 1
            if steps_after_boundary >= 2:
                break
    assert crossed, "chunked run never crossed a chunk boundary"

    # The chunked instance is now on the SECOND chunk; the reference is
    # at the same simulation time.
    assert pd.Timestamp(inst_c.current_time) == pd.Timestamp(inst_n.current_time)

    # 3. AFTER the boundary, the registry value equals the riverine mesh's
    #    CURRENT (second-chunk) value, for a constituent and for depth.
    #    This fails if the bridge were stale on chunk-1 buffers.
    assert _equal_nan(
        np.asarray(_data(reg_c, "algae_floating")), np.asarray(inst_c.mesh["Ap"])
    ), "algae_floating stale: not tracking the current chunk's mesh['Ap']"
    assert _equal_nan(
        np.asarray(_data(reg_c, "depth")), np.asarray(inst_c.mesh["coupling_depth"])
    ), "depth stale: not tracking the current chunk's mesh['coupling_depth']"

    # The bridged array spans the current chunk only (a strict subset of
    # the full 25-step run): proof the reload actually swapped buffers.
    assert np.asarray(_data(reg_c, "algae_floating")).shape[0] < (n_steps + 1)

    # 4. Cross-check: the non-chunked reference yields the same bridged
    #    constituent/depth values at the matching current timestamp.
    ct = np.datetime64(pd.Timestamp(inst_c.current_time))
    for canonical, mesh_name in (
        ("algae_floating", "Ap"),
        ("depth", "coupling_depth"),
    ):
        slice_c = np.asarray(_data(reg_c, canonical).sel(time=ct))
        slice_n = np.asarray(_data(reg_n, canonical).sel(time=ct))
        np.testing.assert_allclose(
            np.nan_to_num(slice_c, nan=-1.0),
            np.nan_to_num(slice_n, nan=-1.0),
            rtol=1e-9,
            atol=1e-12,
            err_msg=(
                f"chunked vs non-chunked {canonical!r} mismatch at {ct}"
            ),
        )
        # And the chunked slice equals the chunked mesh at the same time.
        mesh_slice = np.asarray(inst_c.mesh[mesh_name].sel(time=ct))
        assert _equal_nan(slice_c, mesh_slice)
