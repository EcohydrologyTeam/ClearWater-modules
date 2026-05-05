"""v3 config-loader robustness tests.

Covers three error-reporting improvements to
``clearwater_modules_v3.config.init`` documented as findings M12, M13,
and M15 in ``design/clearwater_modules_v3_review_findings.md``:

* M12: ``_resolve_hotstart`` rejects unsupported file suffixes upfront
  and wraps xarray open failures with a message naming the YAML key.
* M13: ``_resolve_hotstart`` validates a non-integer ``timestep`` via
  ``pandas.to_datetime`` so YAML typos surface before the value reaches
  ``Dataset.sel``.
* M15: ``init_from_config`` reports missing required keys with the full
  YAML path (top-level block, list index, or source name) rather than
  bubbling up a bare ``KeyError``.

The M12/M13 tests exercise ``_resolve_hotstart`` directly. The full
``init_from_config`` reaches ``_resolve_hotstart`` only after process
construction and model-data wiring, both of which require process
plug-ins and a populated ``data_sources`` block. Testing the helper
directly satisfies the spec's "use a minimal config" guidance without
pulling in those concerns.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.config.init import (
    _resolve_hotstart,
    init_from_config,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def valid_nc_path(tmp_path: Path) -> Path:
    """Write a minimal valid netCDF dataset and return its absolute path."""
    ds = xr.Dataset(
        data_vars={
            "water_temperature": ("time", np.array([18.0, 19.0, 20.0])),
        },
        coords={"time": np.array([0, 1, 2])},
    )
    out_path = tmp_path / "hotstart.nc"
    ds.to_netcdf(out_path)
    return out_path


def _minimal_model_block(tmp_path: Path) -> dict:
    """Build a minimally valid ``config['model']`` block.

    The block is sufficient for ``init_from_config`` to clear the
    ``model``-section validation. Tests that aim to fail later sections
    can copy this and add or omit the relevant keys.
    """
    return {
        "start_datetime": "2026-01-01 00:00:00",
        "end_datetime": "2026-01-01 01:00:00",
        "time_step": "5min",
        "simulation_directory": str(tmp_path),
    }


# ---------------------------------------------------------------------------
# M12: unsupported suffix and open-failure error wrapping
# ---------------------------------------------------------------------------


def test_resolve_hotstart_unsupported_suffix_names_yaml_key(tmp_path: Path):
    """A ``.xyz`` suffix raises ValueError that mentions the YAML key
    and the offending suffix, instead of a deep xarray failure."""
    cfg = {"dataset_path": str(tmp_path / "no_such_file.xyz")}
    with pytest.raises(ValueError, match="hotstart"):
        _resolve_hotstart(cfg, tmp_path)


def test_resolve_hotstart_unsupported_suffix_lists_supported_options(
    tmp_path: Path,
):
    cfg = {"dataset_path": str(tmp_path / "x.xyz")}
    with pytest.raises(ValueError) as excinfo:
        _resolve_hotstart(cfg, tmp_path)
    msg = str(excinfo.value)
    assert ".xyz" in msg
    assert ".nc" in msg
    assert ".zarr" in msg


def test_resolve_hotstart_supported_suffix_missing_file_names_yaml_key(
    tmp_path: Path,
):
    """A supported suffix (``.nc``) but a file that does not exist still
    produces a ValueError naming ``hotstart.dataset_path`` rather than a
    bare ``FileNotFoundError`` from xarray."""
    missing_nc = tmp_path / "definitely_does_not_exist.nc"
    cfg = {"dataset_path": str(missing_nc)}
    with pytest.raises(ValueError, match="hotstart.dataset_path"):
        _resolve_hotstart(cfg, tmp_path)


def test_resolve_hotstart_happy_path_valid_nc(valid_nc_path: Path):
    """The happy path still works: a real ``.nc`` file opens and the
    helper returns ``(Dataset, timestep)`` unchanged."""
    cfg = {"dataset_path": str(valid_nc_path)}
    ds, timestep = _resolve_hotstart(cfg, valid_nc_path.parent)
    assert isinstance(ds, xr.Dataset)
    assert timestep is None
    assert "water_temperature" in ds.data_vars
    ds.close()


# ---------------------------------------------------------------------------
# M13: timestep parseability
# ---------------------------------------------------------------------------


def test_resolve_hotstart_invalid_timestep_string_names_yaml_key(
    valid_nc_path: Path,
):
    """A YAML typo like ``2022-13-01`` (invalid month) raises a
    ValueError that names ``hotstart.timestep``."""
    cfg = {"dataset_path": str(valid_nc_path), "timestep": "2022-13-01"}
    with pytest.raises(ValueError, match="hotstart.timestep"):
        _resolve_hotstart(cfg, valid_nc_path.parent)


def test_resolve_hotstart_integer_timestep_passes_through(
    valid_nc_path: Path,
):
    """An integer ``timestep`` is treated as a positional index for
    downstream ``isel`` and must not be coerced through
    ``pd.to_datetime``."""
    cfg = {"dataset_path": str(valid_nc_path), "timestep": 0}
    ds, timestep = _resolve_hotstart(cfg, valid_nc_path.parent)
    assert isinstance(ds, xr.Dataset)
    assert timestep == 0
    assert isinstance(timestep, int)
    ds.close()


# ---------------------------------------------------------------------------
# M15: per-section error messages with full YAML path
# ---------------------------------------------------------------------------


def test_init_from_config_missing_model_start_datetime_names_path(
    tmp_path: Path,
):
    """Missing ``model.start_datetime`` raises ValueError mentioning both
    the top-level block name and the missing key."""
    config = {
        "model": {
            # start_datetime intentionally omitted
            "end_datetime": "2026-01-01 01:00:00",
            "time_step": "5min",
            "simulation_directory": str(tmp_path),
        },
        "processes": [],
        "data_sources": {},
        "variable_map": {},
    }
    with pytest.raises(ValueError) as excinfo:
        init_from_config(config)
    msg = str(excinfo.value)
    assert "model" in msg
    assert "start_datetime" in msg


def test_init_from_config_missing_process_spec_mentions_index(
    tmp_path: Path,
):
    """An entry in ``processes`` that is not a single-key dict naming
    the process must produce an error message identifying the offending
    list index."""
    config = {
        "model": _minimal_model_block(tmp_path),
        # Index 0 is a malformed entry. v3 used to surface this as a
        # bare KeyError or a confusing AttributeError downstream.
        "processes": [
            {"some_process": "not-a-mapping"},
        ],
        "data_sources": {},
        "variable_map": {},
    }
    with pytest.raises(ValueError) as excinfo:
        init_from_config(config)
    msg = str(excinfo.value)
    # The error must point at the offending list element so users can
    # navigate straight to it in the YAML.
    assert "processes" in msg
    assert "0" in msg
    assert "some_process" in msg


def test_init_from_config_missing_data_sources_data_block_names_source(
    tmp_path: Path,
):
    """A data_sources entry missing the required ``data`` block must
    produce an error message mentioning ``data_sources`` and the
    offending source name. (The spec's M15 example targets ``provider``;
    in the v3 implementation, the earliest required key inside a source
    spec is ``data``, and that is what the per-section validator
    reports.)"""
    source_name = "my_source"
    config = {
        "model": _minimal_model_block(tmp_path),
        "processes": [],
        "data_sources": {
            source_name: {
                # ``data`` block intentionally omitted; ``provider`` is
                # ignored at this stage because v3 checks ``data`` first.
                "provider": "csv",
            },
        },
        "variable_map": {},
    }
    with pytest.raises(ValueError) as excinfo:
        init_from_config(config)
    msg = str(excinfo.value)
    assert "data_sources" in msg
    assert source_name in msg
