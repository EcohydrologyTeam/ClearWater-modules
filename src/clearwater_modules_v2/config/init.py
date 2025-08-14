from model import Model
import processes
from processes.base import Process
from pathlib import Path
from .read import read_config
from datetime import timedelta, datetime
from variables import VariableRegistry
import data

import pandas as pd

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
    store_path = data.init_data_store(
        root_directory=root_directory,
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
        variables=variables,
    )
    __populate_data_store(
        config=config,
        store_path=store_path,
        variables=variables,
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
    )

    # TODO: read data sources from conf
    return Model(
        processes=processes,
        variable_registry=VariableRegistry(),
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
        output_variables=config["model"].get("output_variables", []),
    )


def __populate_data_store(
    config: dict,
    store_path: Path,
    variables: set[str],
    start_time: datetime,
    end_time: datetime,
    time_step: timedelta,
) -> None:
    # read and load all sources
    # init data source providers
    sources = __init_data_sources(config)
    # crosswalk the variabels assoicated with each source
    source_variable_map = __parse_variable_map(config["variable_map"])
    for source_name, variable_field_map in source_variable_map.items():
        # check if the user failed to provide a source definition
        if source_name not in sources:
            raise KeyError(f"Source {source_name} not found in configuration")

        # load in the data from the source
        # temporarily skip float data providers
        source = sources[source_name]
        if isinstance(source, data.FloatDataProvider):
            continue

        # for each variable in the source
        for variable_name, field_name in variable_field_map.items():
            if variable_name not in variables:
                warnings.warn(
                    f"Variable not required for any processes: {variable_name} will not be written to the data store"
                )
                continue
            source.write_to_store(
                store_path=store_path,
                start_time=start_time,
                end_time=end_time,
                time_step=time_step,
                variable_name=variable_name,
                field_name=field_name,
            )


def __parse_variable_map(
    variable_map: dict[str, str],
) -> dict[str, dict[str, str | None]]:
    # the variable map as specified by the user will map variables to their sources and potential field names
    # for data load, we'll want to loop data sources and then save them out as their respective variables
    # this method converts the user provided variable mapping to {source : {variable_name : field_name|None}}
    parsed_map = {}
    for variable_name, source_specification in variable_map.items():
        # split to source definition to put out source name and field name
        if len(source_specification.split("|")) == 2:
            source_name, field_name = source_specification.split("|")
        else:
            source_name, field_name = source_specification, None

        # create the dictionary for the source variables if needed
        if parsed_map.get(source_name) is None:
            parsed_map[source_name] = {}
        # add the variable to the source dictionary
        parsed_map[source_name][variable_name] = field_name

    return parsed_map


def __init_data_sources(config: dict) -> dict[str, data.DataProvider]:
    data_source: dict[str, data.DataProvider] = {}
    for source_name, source_config in config["data_sources"].items():
        provider_name = source_config["provider"]
        if "|" in source_name:
            raise ValueError(
                f"Invalid source name: {source_name}. Source names cannot contain the '|' character."
            )
        if provider_name.lower() == "csv":
            data_source[source_name] = data.CSVDataProvider(**source_config["data"])
        elif provider_name.lower() == "float":
            data_source[source_name] = data.FloatDataProvider(**source_config["data"])
        else:
            raise ValueError(
                f"Unknown data or unsupported data provider type: `{provider_name}` for data_source {source_name}"
            )
    return data_source


def __init_processes(config: dict, default_time_step: timedelta) -> list[Process]:
    process_instances = []
    for process in config["processes"]:
        process_name, process_config = *process.keys(), *process.values()
        if process_name.lower() == "riverine":
            process_instances.append(
                __init_riverine(
                    process_config,
                    default_time_step,
                    config,
                )
            )
        elif process_name.lower() == "temperature":
            process_instances.append(
                __init_temperature(process_config, default_time_step)
            )
        else:
            raise ValueError(f"Unknown process type: {process_name}")
    return process_instances


def __init_riverine(
    process_config: dict,
    default_time_step: timedelta,
    config: dict,
) -> processes.Riverine:
    configuration_path = process_config["configuration_path"]
    if "time_step" in process_config:
        time_step_frequency = pd.Timedelta(process_config["time_step"])
    else:
        time_step_frequency = default_time_step
    return processes.Riverine.from_file_path(
        configuration_path,
        start_datetime=config["model"]["start_datetime"],
        end_datetime=config["model"]["end_datetime"],
        time_step_frequency=time_step_frequency,
    )


def __init_temperature(
    process_config: dict,
    default_time_step: timedelta,
) -> processes.Temperature:
    if "time_step" in process_config:
        time_step_frequency = pd.Timedelta(process_config["time_step"])
        process_config.pop("time_step")
    else:
        time_step_frequency = default_time_step
    return processes.Temperature(
        **process_config, time_step_frequency=time_step_frequency
    )


"""
    init_model_data -> generate a zarr store which I can write data to
    loop for provider in sources:
        provider -> data from source -> source_dataset
        write source_dataset to zarr store at mapped_variable_name

    we also need a method to read data from the zarr into the registry 
    but figure that out next

"""
