"""Model-output zarr store builds for NON-temperature constituents.

Regression test for the framework limitation where the model-output zarr
store only sized its spatial axis for ``water_temperature`` (a hardcoded
special-case in ``Model.__init_output_source``). Every other
riverine-bridged constituent reported ``space_dimension is None``, so a
config whose ``output_variables`` listed nutrients either produced a
spatially-flat (time-only) template that the ``(time, nface)`` chunk write
could not land in, or -- once a bare ``space_dimension`` was set without
values -- crashed ``_parse_zarr_coordinates`` with ``len(None)``.

The fix wires ``space_dimension="nface"`` onto every bridged constituent
DataArray (``Riverine._bridge_mesh_to_registry``) and generalizes
``Model.__init_output_source`` to resolve the spatial dimension + coordinate
values for ANY constituent, dropping the temperature one-off. This test
drives the full config-driven coupled run with nutrient output variables and
asserts the resulting ``model_outputs.zarr`` carries each nutrient on the
riverine mesh's ``nface`` axis.

REQUIRED INVOCATION (the conda ``clearwater`` env has a zarr-3-incompatible
xarray; the riverine pixi ``dev`` env has a working one, and PYTHONPATH adds
the modules source). Run from the modules repo dir::

    PYTHONPATH=/Users/todd/GitHub/ecohydrology/ClearWater-modules/src \\
      /Users/todd/GitHub/ecohydrology/ClearWater-riverine/.pixi/envs/dev/bin/python \\
      -m pytest tests/v3/test_riverine_output_store_constituents_v3.py -q
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import xarray
import yaml
from packaging.version import Version

from clearwater_data.variables import DataArrayVariable

# Importing the processes package runs the @ProcessFactory.register decorators
# for every NSM1 process, so init_from_file can resolve them by name.
import clearwater_modules_v3.processes  # noqa: F401
from clearwater_modules_v3.config.init import init_from_file
from clearwater_modules_v3.model import Model


_RIVERINE_REPO = Path(__file__).resolve().parents[3] / "ClearWater-riverine"
PLAN02 = _RIVERINE_REPO / "tests" / "data" / "simple_test_cases" / "plan02_2x1"
PLAN02_HDF = "clearWaterTestCases.p02.hdf"

_RAS_TIME_PATH = (
    "Results/Unsteady/Output/Output Blocks/Base Output/"
    "Unsteady Time Series/Time Date Stamp"
)

# The 16 transported constituents the full bridge maps to canonical names.
_FORK_CONSTITUENTS = [
    "Ap", "Ab", "NH4", "NO3", "OrgN", "N2", "TIP", "OrgP",
    "POC", "DOC", "DIC", "CBOD", "POM", "DOX", "Alk", "PX",
]

# The full NSM1 kinetics process set, in dependency order (mirrors the
# init_from_file integration test). The Phosphorus process declares
# FloatingAlgae and BenthicAlgae as upstream, so the full ordered set is
# required to satisfy the Model's process-order DAG check.
_KINETICS_PROCESSES = [
    {"floating_algae": {}},
    {"benthic_algae": {}},
    {"nitrogen": {}},
    {"phosphorus": {}},
    {"carbon": {}},
    {"pom": {}},
    {"cbod": {}},
    {"dox": {"parameters": {"pressure_mb": 1013.25, "kah_20_user": 20.0,
                            "hydraulic_reaeration_option": 1}}},
    {"n2": {"parameters": {"pressure_mb": 1013.25}}},
    {"pathogen": {}},
    {"alkalinity": {}},
]

# NON-temperature constituents requested as model outputs. These are exactly
# the variables that used to break the output store: each is a riverine-bridged
# DataArray on the mesh ``nface`` axis, none is ``water_temperature``.
_NUTRIENT_OUTPUTS = ["ammonium", "nitrate", "tip", "oxygen_dissolved"]


pytestmark = [
    pytest.mark.skipif(
        not (PLAN02 / PLAN02_HDF).exists(),
        reason=(
            "ClearWater-riverine plan02 fixture not found at "
            f"{PLAN02 / PLAN02_HDF}; sibling repo checkout required"
        ),
    ),
    pytest.mark.skipif(
        Version(xarray.__version__) < Version("2025.8"),
        reason=(
            "init_from_file exercises the zarr-backed model-data path; "
            "xarray<2025.8 is incompatible with zarr 3.x (conda clearwater "
            "env). Runs in pixi dev."
        ),
    ),
]


# ---------------------------------------------------------------------------
# Unit coverage for the generic spatial-dimension resolver (no RAS fixture).
#
# These exercise ``Model.__resolve_space_dimensions`` directly: it replaced
# the v2 ``if variable_name == "water_temperature"`` one-off. Both paths must
# yield the ``nface`` axis with its coordinate values:
#   - a constituent that self-reports ``space_dimension="nface"`` (the bridge);
#   - a temperature-like field registered WITHOUT a self-reported spatial dim
#     (the fallback that derives non-time dims from the array itself).
# ---------------------------------------------------------------------------

_resolve = Model._Model__resolve_space_dimensions


def _mesh_like(name: str) -> xarray.DataArray:
    """A (time, nface) field mimicking a bridged mesh array: nface is a
    dimension WITHOUT an explicit coordinate (as the riverine constituents
    are), so the resolver must fall back to the default integer index."""
    return xarray.DataArray(
        np.zeros((2, 3)),
        dims=("time", "nface"),
        coords={"time": pd.date_range("2023-01-01", periods=2, freq="h")},
        name=name,
    )


def test_resolve_space_dimensions_self_reported():
    var = DataArrayVariable(_mesh_like("ammonium"), space_dimension="nface")
    resolved = _resolve(var)
    assert len(resolved) == 1
    dim, values = resolved[0]
    assert dim == "nface"
    assert len(values) == 3


def test_resolve_space_dimensions_fallback_for_unreported_dim():
    # water_temperature-style: registered without space_dimension. The
    # resolver must still discover the non-time 'nface' axis (the v2
    # hardcode's behavior, now generic).
    var = DataArrayVariable(_mesh_like("water_temperature"))
    assert var.space_dimension is None
    resolved = _resolve(var)
    assert len(resolved) == 1
    dim, values = resolved[0]
    assert dim == "nface"
    assert len(values) == 3


def test_resolve_space_dimensions_scalar_contributes_nothing():
    # A float-provider-style scalar yields no spatial dimension.
    class _Scalar:
        space_dimension = None
        space_dimension_values = None

        def get(self):
            return 15.0

    assert _resolve(_Scalar()) == []


def _hdf_stamps(hdf_path: Path):
    with h5py.File(hdf_path, "r") as f:
        raw = f[_RAS_TIME_PATH][()]
    return pd.to_datetime(
        pd.Series(raw).str.decode("utf8"), format="%d%b%Y %H:%M:%S"
    )


def _hdf_time_bounds(hdf_path: Path):
    stamps = _hdf_stamps(hdf_path)
    return stamps.iloc[0], stamps.iloc[-1]


def _hdf_time_step(hdf_path: Path):
    stamps = _hdf_stamps(hdf_path)
    return stamps.iloc[1] - stamps.iloc[0]


def _even_chunk_size(dt: pd.Timedelta, start, end) -> pd.Timedelta:
    """A chunk_size giving an exact, even >=2-chunk split of the run.

    The model-output store is a ChunkedZarrDataStore; it requires a
    chunk_size that is an integer multiple of time_step. A split into 2 or
    3 equal chunks (>=2 steps each) also keeps the per-chunk interior write
    and the trailing partial-chunk write from colliding on the same region.
    Mirrors ``_even_chunk_size`` in the chunked integration test.
    """
    n_steps = round((end - start).total_seconds() / dt.total_seconds())
    m = next((k for k in (2, 3) if n_steps % k == 0 and n_steps // k >= 2), None)
    assert m is not None, f"plan02 step count {n_steps} has no even >=2-chunk split"
    return dt * (n_steps // m)


def _write_configs(simdir: Path, output_variables: list[str]):
    start, end = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    dt = _hdf_time_step(PLAN02 / PLAN02_HDF)
    chunk = _even_chunk_size(dt, start, end)

    riverine_consts = {
        c: {
            "initial_conditions": {"provider": "float", "data": {"value": 1.0}},
            "boundary_conditions": {"provider": "float", "data": {"value": 1.0}},
        }
        for c in _FORK_CONSTITUENTS
    }
    riverine_cfg = {
        "model": {
            "simulation_directory": str(simdir),
            "hydrodynamic_input": str((PLAN02 / PLAN02_HDF).resolve()),
            "start_datetime": str(start),
            "end_datetime": str(end),
            "diffusion_coefficient": 0.01,
            "output_variables": [],
            "mass_flux_calculation": False,
        },
        "constituents": riverine_consts,
    }
    (simdir / "riverine.yml").write_text(
        yaml.safe_dump(riverine_cfg, sort_keys=False)
    )

    processes = [{"riverine": {"configuration_path": "riverine.yml"}}]
    processes.extend(_KINETICS_PROCESSES)

    data_sources = {
        "water_temperature": {"provider": "float", "data": {"value": 15.0}},
        "solar_radiation": {"provider": "float", "data": {"value": 250.0}},
        "atmospheric_pressure": {"provider": "float", "data": {"value": 1013.25}},
    }
    variable_map = {
        "water_temperature": "water_temperature",
        "solar_radiation": "solar_radiation",
        "atmospheric_pressure": "atmospheric_pressure",
    }
    model_cfg = {
        "model": {
            "start_datetime": str(start),
            "end_datetime": str(end),
            "time_step": f"{int(dt.total_seconds())}s",
            "chunk_size": f"{int(chunk.total_seconds())}s",
            "simulation_directory": str(simdir),
            # The crux: NON-temperature constituents as model outputs.
            "output_variables": output_variables,
        },
        "processes": processes,
        "data_sources": data_sources,
        "variable_map": variable_map,
    }
    model_path = simdir / "model.yml"
    model_path.write_text(yaml.safe_dump(model_cfg, sort_keys=False))
    return model_path, start, end


def _build_model_and_init(model_path: str):
    """Assemble the model and run only its initialization phase.

    ``Model.__init_model`` (name-mangled ``_Model__init_model``) loads the
    registry, runs every ``init_process`` (which fires the riverine bridge,
    self-reporting ``space_dimension="nface"`` on each constituent), and then
    builds the model-output zarr template via ``__init_output_source``. That
    template build is exactly the path the old ``water_temperature`` hardcode
    gated and the path that crashed for nutrient outputs pre-fix, so driving
    just the init phase isolates the fix from the unrelated per-substep
    transport/write machinery.
    """
    model = init_from_file(model_path)
    init = getattr(model, "_Model__init_model")
    init()
    return model


def test_output_store_template_builds_for_nutrient_constituents():
    simdir = Path(tempfile.mkdtemp())
    model_path, _, _ = _write_configs(simdir, _NUTRIENT_OUTPUTS)

    # Pre-fix this raised inside __init_output_source: the bridged nutrient
    # DataArrays reported no spatial dimension (and once a bare
    # space_dimension was set without values, _parse_zarr_coordinates hit
    # ``len(None)``), so the (time, nface) output template could not be sized.
    model = _build_model_and_init(str(model_path))

    store_path = simdir / "model_outputs.zarr"
    assert store_path.exists(), "model_outputs.zarr template was not created"

    ds = xarray.open_zarr(store_path, consolidated=False)
    try:
        # Every requested NON-temperature constituent is present in the
        # template and carries the riverine mesh's spatial axis -- via the
        # generic space_dimension path, not a temperature special-case.
        for name in _NUTRIENT_OUTPUTS:
            assert name in ds.data_vars, f"{name!r} missing from output template"
            assert "nface" in ds[name].dims, (
                f"{name!r} not sized on the riverine mesh 'nface' axis; "
                f"got dims {ds[name].dims}"
            )
        assert ds.sizes["nface"] > 0
    finally:
        ds.close()


def test_output_store_nface_matches_mesh_extent():
    """The output template's ``nface`` length equals the transport mesh's
    cell count -- the coordinate values are derived from the bridged
    constituent's own array, not hardcoded for temperature."""
    simdir = Path(tempfile.mkdtemp())
    model_path, _, _ = _write_configs(simdir, ["ammonium"])

    model = _build_model_and_init(str(model_path))

    reg = getattr(model, "registry", None) or getattr(
        model, "_Model__registry", None
    )
    assert reg is not None
    # The bridged array's nface length, read straight from the registry.
    bridged = reg.get_variable("ammonium").get_data()
    expected_nface = bridged.sizes["nface"]

    ds = xarray.open_zarr(simdir / "model_outputs.zarr", consolidated=False)
    try:
        assert ds.sizes["nface"] == expected_nface
    finally:
        ds.close()
