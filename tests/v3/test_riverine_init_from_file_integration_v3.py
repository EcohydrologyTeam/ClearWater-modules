"""End-to-end ``init_from_file`` integration for the Riverine full-state bridge.

Proves the config-driven assembly path: a single top-level YAML wires the
``Riverine`` transport process plus the full NSM1 kinetics process set, and
``init_from_file`` -> ``model.run()`` advances the coupled model to completion
with every canonical state variable populated. This is the path the Tier 2
bridge (``design/clearwater_modules_v3_riverine_full_state_bridge.md``) exists
to enable: a full-NSM1 coupled run without the hand-rolled stepping loop.

(The design spec refers to reproducing a ``build_v3_modules`` trajectory; no
function of that name exists in the tree -- it is shorthand for the manual
``Model(processes=[...])`` construction path. The achievable, durable assertion
is that the config-driven ``init_from_file`` path assembles the same logical
model and runs it to completion, which is what this test checks.)

REQUIRED INVOCATION (the conda ``clearwater`` env has a zarr-3-incompatible
xarray; the riverine pixi ``dev`` env has a working one, and PYTHONPATH adds
the modules source). Run from the modules repo dir::

    PYTHONPATH=/Users/todd/GitHub/ecohydrology/ClearWater-modules/src \\
      /Users/todd/GitHub/ecohydrology/ClearWater-riverine/.pixi/envs/dev/bin/python \\
      -m pytest tests/v3/test_riverine_init_from_file_integration_v3.py -q
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import yaml

# Importing the processes package runs the @ProcessFactory.register decorators
# for every NSM1 process, so init_from_file can resolve them by name. Without
# this import only transitively-imported processes are registered.
import clearwater_modules_v3.processes  # noqa: F401
from clearwater_modules_v3.config.init import init_from_file


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

# The full NSM1 kinetics process set (riverine prepended at build time).
_KINETICS_PROCESSES = [
    {"floating_algae": {}},
    {"benthic_algae": {}},
    {"nitrogen": {}},
    {"phosphorus": {}},
    {"carbon": {}},
    {"pom": {}},
    {"cbod": {}},
    # Process parameter overrides nest under ``parameters:`` (the form
    # ``ProcessFactory.from_config`` -> ``DOX(**config)`` accepts; the
    # constructor takes ``parameters`` and ``time_step`` only).
    {"dox": {"parameters": {"pressure_mb": 1013.25, "kah_20_user": 20.0,
                            "hydraulic_reaeration_option": 1}}},
    {"n2": {"parameters": {"pressure_mb": 1013.25}}},
    {"pathogen": {}},
    {"alkalinity": {}},
]

# Canonical state variables the bridge feeds the kinetics, spot-checked post-run.
_CANONICAL_SPOT_CHECK = [
    "oxygen_dissolved", "ammonium", "nitrate", "tip",
    "algae_floating", "benthic_algae", "alkalinity", "pathogen",
]


pytestmark = pytest.mark.skipif(
    not (PLAN02 / PLAN02_HDF).exists(),
    reason=(
        "ClearWater-riverine plan02 fixture not found at "
        f"{PLAN02 / PLAN02_HDF}; sibling repo checkout required"
    ),
)


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
    """The transport grid's native step. The bridged constituent arrays
    carry this time index, and the coupled kinetics substep must land on
    it exactly (``get_at_time`` does an exact ``.sel``), so the model
    ``time_step`` must equal the transport step."""
    stamps = _hdf_stamps(hdf_path)
    return stamps.iloc[1] - stamps.iloc[0]


def _write_configs(simdir: Path):
    """Write the riverine sub-config and the top-level model config.

    Returns the path to the top-level model YAML for ``init_from_file``.
    """
    start, end = _hdf_time_bounds(PLAN02 / PLAN02_HDF)
    dt = _hdf_time_step(PLAN02 / PLAN02_HDF)

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

    # Forcings the kinetics read that the riverine bridge does not supply.
    # Float providers are keyed by name through variable_map.
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
            # Must equal the transport grid's native step: the bridged
            # constituent arrays are time-indexed on that grid and the
            # coupled kinetics read them with an exact get_at_time.
            "time_step": f"{int(dt.total_seconds())}s",
            "simulation_directory": str(simdir),
            "output_variables": [],
        },
        "processes": processes,
        "data_sources": data_sources,
        "variable_map": variable_map,
    }
    model_path = simdir / "model.yml"
    model_path.write_text(yaml.safe_dump(model_cfg, sort_keys=False))
    return model_path


def _registry(model):
    reg = getattr(model, "registry", None)
    if reg is None:
        reg = getattr(model, "_Model__registry", None)
    assert reg is not None, "could not locate the model registry"
    return reg


def test_init_from_file_full_nsm1_runs_to_completion():
    simdir = Path(tempfile.mkdtemp())
    model_path = _write_configs(simdir)

    # Config-driven assembly: riverine transport + full NSM1 kinetics.
    model = init_from_file(str(model_path))

    # Run the coupled model to completion (no exception => the bridge fed
    # every kinetics input under its canonical name at every substep).
    model.run()

    # Every spot-checked canonical state is present and carries finite values.
    reg = _registry(model)
    for name in _CANONICAL_SPOT_CHECK:
        assert name in reg, f"{name!r} not registered after run"
        arr = np.asarray(reg.get_variable(name).get_data())
        assert arr.size > 0, f"{name!r} is empty"
        assert np.isfinite(arr).any(), f"{name!r} has no finite values"
