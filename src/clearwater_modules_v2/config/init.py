from model import Model
from processes.base import Process, ProcessFactory
from pathlib import Path
from config.read import read_config
from datetime import timedelta, datetime

from clearwater_data.variables import DataArrayVariable, Variable, VariableRegistry
from clearwater_data.io.zarr import ZarrDataStore, ZarrDataSource
from clearwater_data.io.csv import CSVDataSource
from clearwater_data.io.base import DataSource, ChunkedDataSource
from clearwater_data.io.float import FloatDataSource

import pandas as pd
import xarray as xr

from clearwater_data.custom_types import ArrayLike

import warnings


def init_from_file(file_path: Path | str) -> Model:
    config = read_config(file_path)
    return init_from_config(config)


def init_from_config(config: dict) -> Model:
    # pull out model level configuration information
    try:
        start_time = pd.to_datetime(config["model"]["start_datetime"])
        end_time = pd.to_datetime(config["model"]["end_datetime"])
        time_step = pd.Timedelta(config["model"]["time_step"])
        root_directory = Path(config["model"]["root_directory"])
    except KeyError as e:
        raise ValueError(f"Missing key in config: {e}")

    # initialize the model process instances based on the configuration
    processes = __init_processes(config, default_time_step=time_step)

    # initialize the data store from data
    variables = {v for p in processes for v in p.variables}

    # TODO: this needs to be replaced by ZarrDataSource
    # store_path = data.init_data_store(
    #    root_directory=root_directory,
    #    start_time=start_time,
    #    end_time=end_time,
    #    time_step=time_step,
    #    variables=variables,
    # )

    # read data from original sources and map to an input zarr store
    variable_data_sources = __init_model_data(
        config=config,
        variables=variables,
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
    )

    model_data_source = ZarrDataSource(store_path=root_directory / "model_inputs.zarr")

    output_data_store = ZarrDataStore(
        store_path=root_directory / "model_outputs.zarr",
        start_date=start_time,
        end_date=end_time,
        time_step=time_step,
        variables=config["model"].get("output_variables", []),
    )

    # TODO: read data sources from conf
    return Model(
        processes=processes,
        variable_registry=VariableRegistry(),
        variable_data_sources=variable_data_sources,
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
        output_variables=config["model"].get("output_variables", []),
        output_store=output_data_store,
    )


def __init_model_data(
    config: dict,
    variables: set[str],
    start_time: datetime,
    end_time: datetime,
    time_step: timedelta,
) -> dict[str, DataSource | ChunkedDataSource]:
    # read and load all original sources
    # init data source providers
    sources = __init_data_sources(config)

    # crosswalk the variables associated with each source to their model naming equivalent
    source_variable_map = __parse_variable_map(config["variable_map"])

    # init model data store
    # this is an intermediate data storage solution for model inputs
    data_store = ZarrDataStore(
        store_path=Path(config["model"]["root_directory"]) / "model_inputs.zarr",
        start_date=start_time,
        end_date=end_time,
        time_step=time_step,
        variables=variables,
    )

    # init model input data source
    # this is an instance of a ZarrDataSource that points to the model inputs data store
    model_input_data_source = ZarrDataSource(store_path=data_store.store_path)

    # init model variablle map
    variable_data_sources = {}

    for source_name, variable_parameter_map in source_variable_map.items():
        # check if the user failed to provide a source definition
        if source_name not in sources:
            raise KeyError(f"Source {source_name} not found in configuration")

        # float data can be provided directly to the model
        source = sources[source_name]
        # TODO we need to come back for float data sources
        if isinstance(source, FloatDataSource):
            variable_data_sources[source_name] = source
            continue

        # other data sources need to be read from the original data source
        # and written to the intermediate data store
        for variable_name, parameter_name in variable_parameter_map.items():
            # handle case where a varibale is not required for any of the defined processes
            if variable_name not in variables:
                warnings.warn(
                    f"Variable not required for any processes: {variable_name} will not be written to the data store"
                )
                continue

            # read the data from the original data source
            data = source.read(parameter_name)

            # To store all input data in a common zarr store, we need to align the raw data
            # to the time dimensions of the model input data store
            data.subset_time(start_time, end_time)
            data.resample(
                new_time_frequency=time_step,
            )  # TODO: consider support for different interpolation methods

            # validate and transform the data
            data = __rename_parameter(
                variable=data,
                parameter_name=parameter_name,
                variable_name=variable_name,
            )

            # write the data to the intermediate data store
            data_store.write(data, variable_name)

            # record the datasource in the variable_data_sources dictionary
            variable_data_sources[variable_name] = model_input_data_source

    return variable_data_sources


def __parse_variable_map(
    variable_map: dict[str, str],
) -> dict[str, dict[str, str | None]]:
    # the variable map as specified by the user will map variables to their sources and potential parameter names
    # for data load, we'll want to loop data sources and then save them out as their respective variables
    # this method converts the user provided variable mapping to {source : {variable_name : parameter_name|None}}
    parsed_map = {}
    for variable_name, source_specification in variable_map.items():
        # split to source definition to put out source name and parameter name
        if len(source_specification.split("|")) == 2:
            source_name, parameter_name = source_specification.split("|")
        else:
            source_name, parameter_name = source_specification, None

        # create the dictionary for the source variables if needed
        if parsed_map.get(source_name) is None:
            parsed_map[source_name] = {}
        # add the variable to the source dictionary
        parsed_map[source_name][variable_name] = parameter_name

    return parsed_map


def __init_data_sources(
    config: dict,
) -> dict[str, DataSource]:
    data_source: dict[str, DataSource] = {}
    for source_name, source_config in config["data_sources"].items():
        provider_name = source_config["provider"]
        if "|" in source_name:
            raise ValueError(
                f"Invalid source name: {source_name}. Source names cannot contain the '|' character."
            )
        if provider_name.lower() == "csv":
            data_source[source_name] = CSVDataSource(**source_config["data"])
        elif provider_name.lower() == "float":
            data_source[source_name] = FloatDataSource(**source_config["data"])
        else:
            raise ValueError(
                f"Unknown data or unsupported data provider type: `{provider_name}` for data_source {source_name}"
            )
    return data_source


def __init_processes(config: dict, default_time_step: timedelta) -> list[Process]:
    process_instances = []
    for process_spec in config["processes"]:
        # TODO: this does not seem scalable to me
        # Maybe we should set up factory pattern or something
        process_name, process_config = *process_spec.keys(), *process_spec.values()

        process_config = process_config if process_config is not None else {}
        if "time_step" not in process_config:
            process_config["time_step"] = default_time_step
        if not isinstance(process_config["time_step"], timedelta):
            process_config["time_step"] = pd.to_timedelta(process_config["time_step"])

        # riverine is a standalone model and needs the datetime range from the model config
        if process_name.lower() == "riverine":
            process_config["start_datetime"] = pd.to_datetime(
                config["model"]["start_datetime"]
            )
            process_config["end_datetime"] = pd.to_datetime(
                config["model"]["end_datetime"]
            )

        process_instance = ProcessFactory.from_config(process_name, process_config)
        process_instances.append(process_instance)

    return process_instances


def __rename_parameter(
    variable: Variable,
    parameter_name: str,
    variable_name: str,
) -> ArrayLike:
    if isinstance(variable.get(), xr.Dataset):
        return variable.get().rename({parameter_name: variable_name})
    elif isinstance(variable.get(), xr.DataArray):
        return variable.get().rename(variable_name)
