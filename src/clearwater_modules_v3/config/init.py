"""v3 ``init_from_file`` entry point.

Backward-compatible extension of v2's ``init_from_file``. Two new optional
top-level YAML keys are honored:

```yaml
hotstart:                     # optional
  dataset_path: hotstart.nc   # path resolved relative to simulation_directory
  timestep: '2022-05-13 12:00:00'   # optional; default is the last timestep

wet_mask:                     # optional
  variable: wetted_surface_area     # registry variable to threshold
  threshold: 1.0              # cells with value > threshold are considered wet
```

When neither key is present, v3 behavior matches v2 exactly: every
existing v2 configuration runs against v3 unchanged.

Implementation reuses v2's helper functions (process construction, data
source construction, variable map parsing) so v3 inherits any changes
LimnoTech makes to those helpers. v3 only owns the entry-point wrapper
and the construction of v3 ``Model``.
"""

from __future__ import annotations

import warnings
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import xarray as xr

from clearwater_data.io.base import ChunkedDataSource, DataSource
from clearwater_data.io.pathing import resolve_path
from clearwater_data.io.zarr import ZarrDataSource
from clearwater_data.variables import VariableRegistry

# v2 helpers are still authoritative for process construction, data-source
# wiring, and variable-map parsing. v3 reuses them so it inherits any
# upstream improvements without forking the implementation.
from clearwater_modules_v2.config import init as _v2_init
from clearwater_modules_v2.config.read import read_config

from clearwater_modules_v3.model import Model


def init_from_file(file_path: Path | str) -> Model:
    """Build a v3 ``Model`` from a YAML configuration file."""
    config = read_config(file_path)
    return init_from_config(config)


def init_from_config(config: dict) -> Model:
    """Build a v3 ``Model`` from an already-parsed configuration dict."""
    # ----- v2-compatible model-level fields -----
    # Per-section validation produces error messages that name the
    # offending YAML path (top-level block, list index, or key). See
    # finding M15 in design/clearwater_modules_v3_review_findings.md:
    # the previous monolithic ``try/except KeyError`` lost the deep key
    # path and made sponsor-support triage harder than necessary.

    # config['model'] required keys
    start_time = pd.to_datetime(_required(config, "model", "start_datetime"))
    end_time = pd.to_datetime(_required(config, "model", "end_datetime"))
    time_step = pd.Timedelta(_required(config, "model", "time_step"))

    simulation_directory_path = resolve_path(
        Path(_required(config, "model", "simulation_directory"))
    )
    config["model"]["simulation_directory"] = simulation_directory_path

    # config['processes'] required: list of single-key dicts
    processes_section = _required(config, "processes")
    for index, process_dict in enumerate(processes_section):
        if not isinstance(process_dict, dict) or not process_dict:
            raise ValueError(
                f"Invalid entry at config['processes'][{index}]: expected a "
                f"single-key dict naming the process; got {process_dict!r}"
            )
        process_name = next(iter(process_dict.keys()))
        process_spec = process_dict[process_name]
        if process_spec is None:
            # v2 tolerates a None spec; nothing more to resolve here.
            continue
        if not isinstance(process_spec, dict):
            raise ValueError(
                f"Invalid spec at config['processes'][{index}][{process_name!r}]: "
                f"expected a mapping; got {process_spec!r}"
            )
        file_path_value = process_spec.get("configuration_path")
        if file_path_value is not None:
            config["processes"][index][process_name]["configuration_path"] = (
                resolve_path(simulation_directory_path / file_path_value)
            )

    # config['data_sources'] required: name -> spec mapping
    data_sources_section = _required(config, "data_sources")
    if not isinstance(data_sources_section, dict):
        raise ValueError(
            f"Invalid value at config['data_sources']: expected a mapping; "
            f"got {data_sources_section!r}"
        )
    for source_name in data_sources_section.keys():
        # Each source must declare a ``data`` block; the ``file_path`` key
        # within ``data`` is optional (e.g. float sources have no file).
        source_data = _required(
            config, "data_sources", source_name, "data"
        )
        if not isinstance(source_data, dict):
            raise ValueError(
                f"Invalid value at config['data_sources'][{source_name!r}]['data']: "
                f"expected a mapping; got {source_data!r}"
            )
        file_path_value = source_data.get("file_path")
        if file_path_value is not None:
            config["data_sources"][source_name]["data"]["file_path"] = resolve_path(
                simulation_directory_path / file_path_value
            )

    # Chunking (v2-compatible).
    chunk_size = config["model"].get("chunk_size", None)
    if chunk_size is None:
        chunk_size = config["model"].get("chunk_time_step", None)
        if chunk_size is not None:
            warnings.warn(
                "The `chunk_time_step` configuration option is deprecated and "
                "will be removed in a future release. Please use `chunk_size` "
                "instead.",
                DeprecationWarning,
                stacklevel=2,
            )
    if chunk_size is not None:
        chunk_size = pd.Timedelta(chunk_size)

    variable_registry = VariableRegistry()

    # Process and data-source wiring through v2's helpers. v2's helpers are
    # double-underscore-prefixed at module level (private by convention,
    # but module-level dunder names are NOT class-mangled, so they are
    # accessible from here via direct attribute access). v3 reuses them
    # rather than forking so any LimnoTech changes to those helpers flow
    # through to v3 unchanged. See finding C9 in
    # design/clearwater_modules_v3_review_findings.md for context.
    init_processes = _resolve_v2_helper("__init_processes")
    init_model_data = _resolve_v2_helper("__init_model_data")

    processes = init_processes(
        config, variable_registry, default_time_step=time_step
    )
    variables = {v for p in processes for v in p.variables}
    variable_data_sources = init_model_data(
        config=config,
        variables=variables,
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
    )

    # ----- v3 additions: hotstart -----
    hotstart_dataset, hotstart_timestep = _resolve_hotstart(
        config.get("hotstart"), simulation_directory_path
    )

    # ----- v3 additions: wet-mask -----
    wet_mask_variable, wet_mask_threshold = _resolve_wet_mask(config.get("wet_mask"))

    return Model(
        processes=processes,
        variable_registry=variable_registry,
        variable_data_sources=variable_data_sources,
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
        output_variables=config["model"].get("output_variables", []),
        simulation_directory=simulation_directory_path,
        chunk_size=chunk_size,
        hotstart_dataset=hotstart_dataset,
        hotstart_timestep=hotstart_timestep,
        wet_mask_variable=wet_mask_variable,
        wet_mask_threshold=wet_mask_threshold,
    )


# ----- v3-specific helpers -----


def _required(d: dict, *path):
    """Look up ``d[path[0]][path[1]]...`` and raise ``ValueError`` with the
    full key path on ``KeyError`` or on a non-mapping intermediate value.

    Used by ``init_from_config`` to produce error messages that name the
    offending YAML location (top-level block, list index, or key) instead
    of bubbling up a bare ``KeyError`` whose message lacks the parent
    path. See finding M15 in
    ``design/clearwater_modules_v3_review_findings.md``.
    """
    cur = d
    for i, key in enumerate(path):
        try:
            cur = cur[key]
        except (KeyError, TypeError, IndexError):
            location = ".".join(repr(p) for p in path[: i + 1])
            raise ValueError(
                f"Missing required key {key!r} at config path: {location}"
            ) from None
    return cur


_HOTSTART_NETCDF_SUFFIXES = (".nc", ".cdf", ".netcdf4")
_HOTSTART_ZARR_SUFFIXES = (".zarr", "")


def _resolve_hotstart(
    hotstart_cfg: dict | None,
    simulation_directory_path: Path,
) -> tuple[xr.Dataset | None, datetime | str | int | None]:
    """Resolve the optional ``hotstart`` YAML block to constructor kwargs.

    Returns ``(hotstart_dataset, hotstart_timestep)``. Both ``None`` when
    the YAML block is absent. ``dataset_path`` is resolved relative to
    ``simulation_directory_path`` if not absolute. Supports netCDF
    (``.nc``, ``.cdf``, ``.netcdf4``) via ``xarray.open_dataset`` and
    Zarr (``.zarr`` or directory) via ``xarray.open_zarr``.

    Raises a descriptive ``ValueError`` naming the YAML key when the
    suffix is unsupported, when the dataset cannot be opened, or when
    a non-integer ``timestep`` cannot be parsed by ``pandas``. See
    findings M12 and M13 in
    ``design/clearwater_modules_v3_review_findings.md``.
    """
    if hotstart_cfg is None:
        return (None, None)
    dataset_path = hotstart_cfg.get("dataset_path")
    if dataset_path is None:
        raise ValueError(
            "hotstart YAML block must include `dataset_path` when present"
        )

    # Resolve relative to simulation_directory_path without invoking
    # ``resolve_path``'s existence validation: we want our own error
    # messages naming the YAML key to take precedence over a bare
    # ``FileNotFoundError``. The suffix check below also runs before any
    # filesystem access so unsupported suffixes are caught even when the
    # file is missing (finding M12).
    raw_path = Path(dataset_path)
    if raw_path.is_absolute():
        dataset_path = raw_path
    else:
        dataset_path = simulation_directory_path / raw_path
    suffix = dataset_path.suffix.lower()

    # Enumerate supported suffixes upfront so users get a clear breadcrumb
    # naming the YAML key rather than a deep xarray traceback (finding M12).
    if suffix in _HOTSTART_ZARR_SUFFIXES:
        opener = xr.open_zarr
    elif suffix in _HOTSTART_NETCDF_SUFFIXES:
        opener = xr.open_dataset
    else:
        supported = ", ".join(
            repr(s) for s in (*_HOTSTART_NETCDF_SUFFIXES, *_HOTSTART_ZARR_SUFFIXES)
        )
        raise ValueError(
            f"Unsupported file suffix {suffix!r} for "
            f"hotstart.dataset_path={str(dataset_path)!r}. "
            f"Supported suffixes: {supported}."
        )

    try:
        ds = opener(dataset_path)
    except (FileNotFoundError, OSError, ValueError) as exc:
        raise ValueError(
            f"Failed to open hotstart.dataset_path={str(dataset_path)!r}: {exc}"
        ) from exc

    timestep = hotstart_cfg.get("timestep")
    # Eagerly validate non-integer timesteps so a typo surfaces here with a
    # breadcrumb to ``hotstart.timestep`` rather than deep inside
    # ``Dataset.sel(method='nearest')`` (finding M13). Integer indexing is
    # routed through ``isel`` downstream and is left untouched. ``bool`` is
    # a subclass of ``int`` in Python; reject it explicitly because
    # ``timestep: true`` in YAML is never meaningful.
    is_int_index = isinstance(timestep, int) and not isinstance(timestep, bool)
    if timestep is not None and not is_int_index:
        try:
            pd.to_datetime(timestep)
        except (ValueError, TypeError, OverflowError) as exc:
            raise ValueError(
                f"Could not parse hotstart.timestep={timestep!r} as a "
                f"datetime: {exc}"
            ) from exc
    return (ds, timestep)


def _resolve_wet_mask(
    wet_mask_cfg: dict | None,
) -> tuple[str | None, float]:
    """Resolve the optional ``wet_mask`` YAML block to constructor kwargs.

    Returns ``(wet_mask_variable, wet_mask_threshold)``. ``(None, 0.0)``
    when the YAML block is absent.
    """
    if wet_mask_cfg is None:
        return (None, 0.0)
    variable = wet_mask_cfg.get("variable")
    if variable is None:
        raise ValueError(
            "wet_mask YAML block must include `variable` when present"
        )
    threshold = float(wet_mask_cfg.get("threshold", 0.0))
    return (variable, threshold)


def _resolve_v2_helper(name: str):
    """Resolve a helper function on the v2 init module by exact name.

    v2's helpers use leading double-underscore names at module scope.
    Python name-mangling applies only inside class bodies, so module-level
    ``__init_processes`` is exposed under that exact attribute name. v3
    looks the helper up directly so any rename or removal upstream
    surfaces as a loud, descriptive ``AttributeError`` at startup rather
    than as a silent failure mid-simulation. See finding C9 in
    design/clearwater_modules_v3_review_findings.md.
    """
    try:
        return getattr(_v2_init, name)
    except AttributeError as exc:
        raise AttributeError(
            f"v3 expected `{name}` on clearwater_modules_v2.config.init; "
            f"if v2 has been refactored, update v3's reuse contract in "
            f"clearwater_modules_v3/config/init.py and "
            f"tests/v3/test_v2_helper_contract.py."
        ) from exc


__all__ = ["init_from_file", "init_from_config"]
