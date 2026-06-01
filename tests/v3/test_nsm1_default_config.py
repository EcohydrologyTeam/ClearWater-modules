"""The shipped ``nsm1_default.yml`` is config-driven-loadable.

Guards two things:

1. Each process's parameter block is nested under ``parameters:``. A
   top-level param (e.g. ``dox: {pressure_mb: …}``) makes the factory call
   ``DOX(pressure_mb=…)`` — but the constructor only accepts ``parameters``
   and ``time_step``, so it raises ``TypeError`` during process
   construction. This regression guard parses the YAML directly and runs
   in any environment.
2. The config assembles into a ``Model`` through ``init_from_file`` and the
   nested params reach the constructor. ``init_from_file`` exercises the
   zarr-backed model-data setup, which is broken under xarray < 2025.8
   (zarr 3.x incompatibility — the conda ``clearwater`` env); that test is
   skipped there and runs in the pixi ``dev`` env.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml
import xarray
from packaging.version import Version

import clearwater_modules_v3

_CONFIG = (
    Path(clearwater_modules_v3.__file__).parent / "config" / "nsm1_default.yml"
)
# Processes in nsm1_default.yml that carry an explicit parameter block.
_PARAMETERIZED = {"dox", "n2"}


def _process_blocks() -> dict:
    """Return ``{process_name: config_dict}`` from the shipped YAML."""
    doc = yaml.safe_load(_CONFIG.read_text())
    blocks = {}
    for entry in doc["processes"]:
        (name, cfg), = entry.items()
        blocks[name] = cfg or {}
    return blocks


# --- always-on YAML-structure guard ---------------------------------------


def test_process_params_nested_under_parameters():
    """Any process that carries parameters nests them under ``parameters:``
    (so the factory's ``Process(**config)`` gets ``parameters=…``, not loose
    kwargs the constructor would reject)."""
    blocks = _process_blocks()
    for name in _PARAMETERIZED:
        cfg = blocks[name]
        assert "parameters" in cfg, (
            f"{name!r} block must nest its params under 'parameters:'; got "
            f"top-level keys {sorted(cfg)}"
        )
        assert isinstance(cfg["parameters"], dict) and cfg["parameters"]
        # No constructor-incompatible loose keys remain alongside.
        assert set(cfg) <= {"parameters", "time_step"}, (
            f"{name!r} has loose top-level keys: {sorted(set(cfg) - {'parameters', 'time_step'})}"
        )
    # dox's nested params carry the documented non-default reaeration.
    assert blocks["dox"]["parameters"]["kah_20_user"] == 20.0


# --- end-to-end build via init_from_file (pixi dev env) --------------------


@pytest.mark.skipif(
    Version(xarray.__version__) < Version("2025.8"),
    reason=(
        "init_from_file exercises the zarr-backed model-data path; xarray<2025.8 "
        "is incompatible with zarr 3.x (conda clearwater env). Runs in pixi dev."
    ),
)
def test_shipped_config_builds_via_init_from_file(tmp_path, monkeypatch):
    import clearwater_modules_v3.processes  # noqa: F401 — register factories
    from clearwater_modules_v3.config.init import init_from_file

    # The shipped config uses ``simulation_directory: '.'``, so init_from_file
    # writes its model_inputs.zarr store into CWD. Run in a tmp dir so the test
    # leaves no artifact in the repo; _CONFIG is an absolute path, so chdir is
    # safe for config resolution.
    monkeypatch.chdir(tmp_path)
    model = init_from_file(str(_CONFIG))
    assert model is not None
    # All 11 NSM1 processes from the YAML assembled (get_process/has_process
    # key by class name).
    expected = [
        "FloatingAlgae", "BenthicAlgae", "Nitrogen", "Phosphorus", "Carbon",
        "POM", "CBOD", "DOX", "N2", "Pathogen", "Alkalinity",
    ]
    for name in expected:
        assert model.has_process(name), f"config did not assemble process {name!r}"
    # The nested dox parameters reached the DOX constructor (a top-level
    # param would have raised TypeError in from_config before init returned).
    dox = model.get_process("DOX")
    assert dox.kah_20_user == 20.0
    assert dox.hydraulic_reaeration_option == 1
