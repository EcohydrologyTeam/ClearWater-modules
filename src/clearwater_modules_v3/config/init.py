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

Implementation. v3 owns the process and data-source wiring directly via
the module-private helpers ``_init_processes``, ``_init_model_data``,
``_init_data_sources``, ``_parse_variable_map``, and
``_rename_parameter``. These were originally lifted from v2 and v3
reused them by attribute lookup; the v3-self-sufficient refactor moved
them in-tree so the v3 config layer is independent of
``clearwater_modules_v2.config.init``. Process instantiation goes
through ``clearwater_modules_v3.processes.base.ProcessFactory``.
"""

from __future__ import annotations

import warnings
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import xarray as xr

from clearwater_data.custom_types import ArrayLike
from clearwater_data.io.base import ChunkedDataSource, DataSource
from clearwater_data.io.csv import CSVDataSource
from clearwater_data.io.float import FloatDataSource
from clearwater_data.io.pathing import resolve_path
from clearwater_data.io.zarr import ZarrDataSource, ZarrDataStore
from clearwater_data.variables import Variable, VariableRegistry

# read_config is still reused from v2 (YAML parsing helper that has not
# been forked into v3). The process-construction and data-source-wiring
# helpers were ported in-tree by the v3-self-sufficient refactor; see
# the module-private functions further down in this file.
from clearwater_modules_v3.config.read import read_config

from clearwater_modules_v3.model import Model
from clearwater_modules_v3.processes.base import Process, ProcessFactory


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

    # Process and data-source wiring through v3-native helpers. These were
    # originally reused from v2 via attribute lookup; the v3-self-sufficient
    # refactor ported them in-tree (see the module-private helpers below).
    processes = _init_processes(
        config, variable_registry, default_time_step=time_step
    )
    variables = {v for p in processes for v in p.variables}
    variable_data_sources = _init_model_data(
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


def _init_processes(
    config: dict,
    variable_registry: VariableRegistry,
    default_time_step: timedelta,
) -> list[Process]:
    """Construct Process instances from the YAML ``processes:`` block.

    Routes through ``ProcessFactory.from_config`` keyed on the process
    name. ``riverine`` is special-cased to receive the model-level
    ``start_datetime`` / ``end_datetime`` because the underlying
    ``ClearwaterRiverine`` is a standalone solver that needs the time
    range up front; all other Process classes derive time information
    from the registry and the per-substep ``time`` argument to ``run``.

    Originally lifted verbatim from
    ``clearwater_modules_v2.config.init.__init_processes``; the
    v3-self-sufficient refactor moved it in-tree and changed
    ``ProcessFactory`` to come from
    ``clearwater_modules_v3.processes.base`` so v3 owns its own factory
    registry.
    """
    process_instances: list[Process] = []
    for process_spec in config["processes"]:
        process_name, process_config = (
            *process_spec.keys(),
            *process_spec.values(),
        )

        process_config = process_config if process_config is not None else {}
        if "time_step" not in process_config:
            process_config["time_step"] = default_time_step
        if not isinstance(process_config["time_step"], timedelta):
            process_config["time_step"] = pd.to_timedelta(
                process_config["time_step"]
            )

        # riverine is a standalone solver and needs the datetime range
        # from the model config.
        if process_name.lower() == "riverine":
            process_config["start_datetime"] = pd.to_datetime(
                config["model"]["start_datetime"]
            )
            process_config["end_datetime"] = pd.to_datetime(
                config["model"]["end_datetime"]
            )

        process_instance = ProcessFactory.from_config(
            process_name, process_config, variable_registry
        )
        process_instances.append(process_instance)

    return process_instances


def _init_model_data(
    config: dict,
    variables: set[str],
    start_time: datetime,
    end_time: datetime,
    time_step: timedelta,
) -> dict[str, DataSource | ChunkedDataSource]:
    """Wire raw data sources into a model-input zarr store and per-variable
    DataSource handles.

    For each variable consumed by any Process, locate its source in the
    YAML ``data_sources:`` block, read the raw values, resample to the
    model time step, and write into a per-simulation zarr store. The
    result is a dict mapping variable names to DataSource instances that
    the Model uses for time-indexed reads at run time.

    Float data sources bypass the zarr store and are passed through as
    direct DataSource objects (constants don't need resampling).

    Originally lifted verbatim from
    ``clearwater_modules_v2.config.init.__init_model_data``.
    """
    sources = _init_data_sources(config)
    source_variable_map = _parse_variable_map(config["variable_map"])

    data_store = ZarrDataStore(
        store_path=config["model"]["simulation_directory"] / "model_inputs.zarr",
        start_date=start_time,
        end_date=end_time,
        time_step=time_step,
        variables=variables,
    )

    model_input_data_source = ZarrDataSource(store_path=data_store.store_path)

    variable_data_sources: dict[str, DataSource | ChunkedDataSource] = {}

    for source_name, variable_parameter_map in source_variable_map.items():
        if source_name not in sources:
            raise KeyError(f"Source {source_name} not found in configuration")

        source = sources[source_name]
        if isinstance(source, FloatDataSource):
            variable_data_sources[source_name] = source
            continue

        for variable_name, parameter_name in variable_parameter_map.items():
            if variable_name not in variables:
                warnings.warn(
                    f"Variable not required for any processes: {variable_name} "
                    f"will not be written to the data store"
                )
                continue

            data = source.read(parameter_name)
            data.subset_time(start_time, end_time)
            data.resample(new_time_frequency=time_step)
            data = _rename_parameter(
                variable=data,
                parameter_name=parameter_name,
                variable_name=variable_name,
            )
            data_store.write(data, variable_name)
            variable_data_sources[variable_name] = model_input_data_source

    return variable_data_sources


def _parse_variable_map(
    variable_map: dict[str, str],
) -> dict[str, dict[str, str | None]]:
    """Convert the user-facing ``variable_map`` block to a source-keyed dict.

    The YAML's ``variable_map`` is keyed by variable name and values
    are ``"source_name|parameter_name"`` strings (or just ``"source_name"``
    if the source uses the same name for the parameter). This helper
    inverts the mapping so the caller can iterate by source.
    """
    parsed_map: dict[str, dict[str, str | None]] = {}
    for variable_name, source_specification in variable_map.items():
        if len(source_specification.split("|")) == 2:
            source_name, parameter_name = source_specification.split("|")
        else:
            source_name, parameter_name = source_specification, None

        if parsed_map.get(source_name) is None:
            parsed_map[source_name] = {}
        parsed_map[source_name][variable_name] = parameter_name

    return parsed_map


def _init_data_sources(config: dict) -> dict[str, DataSource]:
    """Instantiate the DataSource subclasses listed in ``data_sources:``.

    Currently supports ``csv`` and ``float`` providers. Source names
    cannot contain ``|`` (that character is reserved for the
    ``variable_map`` parser).
    """
    data_source: dict[str, DataSource] = {}
    for source_name, source_config in config["data_sources"].items():
        provider_name = source_config["provider"]
        if "|" in source_name:
            raise ValueError(
                f"Invalid source name: {source_name}. Source names cannot "
                f"contain the '|' character."
            )
        if provider_name.lower() == "csv":
            data_source[source_name] = CSVDataSource(**source_config["data"])
        elif provider_name.lower() == "float":
            data_source[source_name] = FloatDataSource(**source_config["data"])
        else:
            raise ValueError(
                f"Unknown data or unsupported data provider type: "
                f"`{provider_name}` for data_source {source_name}"
            )
    return data_source


def _rename_parameter(
    variable: Variable,
    parameter_name: str,
    variable_name: str,
) -> ArrayLike:
    """Rename a raw variable's data to the canonical model-side name."""
    if isinstance(variable.get(), xr.Dataset):
        return variable.get().rename({parameter_name: variable_name})
    elif isinstance(variable.get(), xr.DataArray):
        return variable.get().rename(variable_name)


__all__ = ["init_from_file", "init_from_config"]
