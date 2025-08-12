from model import Model
import processes
from processes.base import Process
from pathlib import Path
from .read import read_config
from datetime import timedelta, datetime
from variables import VariableRegistry

import pandas as pd


def init_from_file(file_path: Path | str) -> Model:
    config = read_config(file_path)
    return init_from_config(config)


def init_from_config(config: dict) -> Model:
    try:
        start_time = pd.to_datetime(config["model"]["start_datetime"])
        end_time = pd.to_datetime(config["model"]["end_datetime"])
        time_step = pd.Timedelta(config["model"]["time_step"])
    except KeyError as e:
        raise ValueError(f"Missing key in config: {e}")
    processes = __get_init_processes(config, default_time_step=time_step)
    # TODO: read data sources from conf
    return Model(
        processes=processes,
        variable_registry=VariableRegistry(),
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
        output_variables=config["model"].get("output_variables", []),
    )


def __get_init_processes(config: dict, default_time_step: timedelta) -> list[Process]:
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
