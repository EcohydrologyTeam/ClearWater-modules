import logging
import os

from clearwater_modules_v2.processes.base import Process
from clearwater_data.variables import VariableRegistry
from datetime import datetime, timedelta
from typing import Iterable
from clearwater_data.io.base import (
    DataSource,
    ChunkedDataSource,
)
from clearwater_data.io.zarr import ChunkedZarrDataStore

from logging import getLogger, basicConfig

LOGGER = getLogger(__name__)


def set_logging_config(log_level: str = "DEBUG") -> None:
    basicConfig(
        level=log_level,
        format="[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


class Model:
    def __init__(
        self,
        processes: tuple[Process],
        variable_registry: VariableRegistry,
        variable_data_sources: dict[str, DataSource | ChunkedDataSource],
        start_time: datetime,
        end_time: datetime,
        time_step: timedelta,
        output_variables: Iterable[str],
        simulation_directory: os.PathLike | None = None,
        # output_store: DataStore | ChunkedDataStore,
        chunk_size: timedelta | None = None,
    ) -> None:
        self.__processes: tuple[Process] = processes
        self.__registry: VariableRegistry = variable_registry
        self.__variable_data_sources: dict[str, DataSource | ChunkedDataSource] = (
            variable_data_sources
        )
        self.__start_time: datetime = start_time
        self.__end_time: datetime = end_time
        self.__time_step: timedelta = time_step
        self.__output_variables: Iterable[str] = output_variables
        self.__simulation_directory: os.Pathlike = simulation_directory if simulation_directory else "."
        # self.__output_store: DataStore | ChunkedDataStore = output_store

        # check if we are running chunked
        self.__chunked_mode: bool = chunk_size is not None
        self.__chunk_size: timedelta | None = chunk_size

        # TODO: if no output warning and don't set up __output_file

        self.__init_complete: bool = False

        # TODO have the configuration provided from the configuration file
        set_logging_config()

    def validate(self) -> None:
        if self.__start_time >= self.__end_time:
            raise ValueError("Start time must be before end time.")

    def init_model(self) -> None:
        self.__init_model()
        self.__init_complete = True

    def run(self) -> None:
        if not self.__init_complete:
            self.__init_model()
        if self.__chunked_mode:
            self.__process_loop_chunked()
        else:
            self.__process_loop_full()

    def has_process(self, process_type: type[Process] | str) -> bool:
        if isinstance(process_type, str):
            return any(
                p.process_name().lower() == process_type.lower()
                for p in self.__processes
            )
        return any(isinstance(p, process_type) for p in self.__processes)

    def get_process(self, process_type: type[Process] | str) -> Process:
        if not self.has_process(process_type):
            raise ValueError(f"Process {process_type} not found in model.")
        if isinstance(process_type, str):
            return next(
                p
                for p in self.__processes
                if p.process_name().lower() == process_type.lower()
            )
        return next(p for p in self.__processes if isinstance(p, process_type))

    def __init_model(self) -> None:
        # load model or first chunk
        for variable_name, data_source in self.__variable_data_sources.items():
            if isinstance(data_source, ChunkedDataSource):
                data = data_source.read_chunk(
                    self.__start_time, self.__start_time + self.__chunk_size
                )
            else:
                data = data_source.read(variable_name)

            self.__registry.register(
                variable_name,
                data,
            )

        # initialize any processes that require initialization
        for process in self.__processes:
            process.init_process(self, self.__registry)

        # set up the output data store
        self.__init_output_source()

    def __finalize_model(self) -> None:
        for process in self.__processes:
            process.finalize_process(self, self.__registry)

    def __init_output_source(self) -> None:
        # we need to have the output variables match the signature of the variable themselves,
        # which means we need to have the dimensions of the variable and provide that
        # to the output store initialization method.

        if self.__output_variables is None or len(self.__output_variables) == 0:
            # TODO: warning
            return

        space_dimensions = {}
        for variable_name in self.__output_variables:
            variable = self.__registry.get_variable(variable_name)

            # We need the dimensions of the output store to match the dimensions of the variables we are writing out.
            # If there are any spatial dimensions associated with the variable, we need to pass that information to the output store.
            if (
                variable.space_dimension is not None
                and variable.space_dimension not in space_dimensions
            ):
                space_dimensions[variable.space_dimension] = (
                    variable.space_dimension_values
                )

            # TODO: This is a manual override to set spatial dimension until we have
            # fully implemented variable space dimension in riverine
            if variable_name == "water_temperature":
                data = variable.get()
                space_dimensions["nface"] = data["nface"].values

        self.__output_data_store = ChunkedZarrDataStore(
            store_path=self.__simulation_directory / "model_outputs.zarr",
            start_date=self.__start_time,
            end_date=self.__end_time,
            time_step=self.__time_step,
            variables=self.__output_variables,
            spatial_field=(
                list(space_dimensions.keys()) if len(space_dimensions) > 0 else None
            ),
            spatial_field_values=(
                list(space_dimensions.values()) if len(space_dimensions) > 0 else None
            ),
            chunk_size=self.__chunk_size,
        )

    def __load_chunk_data(self, chunk_start: datetime, chunk_end: datetime) -> None:
        for variable_name, data_source in self.__variable_data_sources.items():
            # Non-chunked data sources should already be loaded in the registry as part of model initialization.
            # we should only need to loaded chunked data sources here.
            if not isinstance(data_source, ChunkedDataSource):
                continue

            data = data_source.read_chunk(chunk_start, chunk_end)
            self.__registry.register(
                variable_name,
                data,
            )

    def __process_loop_chunked(self) -> None:
        # TODO: this need actual chunking logic
        current_time = self.__start_time
        chunk_end_time = self.__start_time + self.__chunk_size

        while current_time <= self.__end_time and chunk_end_time <= self.__end_time:
            # TODO: look at riverine's code and mirror where applicable.
            # TODO: align with riverine

            # load next chunk's data
            # the minus timestep ensures that previous timestep data is available for process logic
            self.__load_chunk_data(current_time - self.__time_step, chunk_end_time)

            while current_time <= chunk_end_time and current_time < self.__end_time:
                current_time_seconds = current_time.timestamp()
                LOGGER.info(f"Running timestep: {current_time}")
                for process in self.__processes:
                    # check if this process should be updated at this timestamp
                    if current_time_seconds % process.time_step_seconds == 0:
                        process.run(current_time, self.__registry)
                current_time += self.__time_step

            # write out chunk's data
            self.__save_output_model(
                start_time=chunk_end_time - self.__chunk_size,  # start of the chunk
                end_time=chunk_end_time,
            )

            self.__finalize_model()

            # iterate to next chunk
            chunk_end_time += self.__chunk_size

        # TODO: confirm if this is necessary to write out the last chunk or if it will be handled in the loop above.
        # output last chunk
        self.__save_output_model(
            start_time=self.__end_time - self.__chunk_size, end_time=self.__end_time
        )

    def __process_loop_full(self) -> None:
        current_time = self.__start_time

        while current_time < self.__end_time:
            current_time_seconds = current_time.timestamp()
            LOGGER.info(f"Running timestep: {current_time}")
            for process in self.__processes:
                # check if this process should be updated at this timestamp
                # Process should be calculated is current_time + process.time_step_frequency
                if current_time_seconds % process.time_step_seconds == 0:
                    process.run(current_time, self.__registry)
            current_time += self.__time_step
        self.__save_output_model(self.__start_time, self.__end_time)

    def __save_output_model(self, start_time: datetime, end_time: datetime) -> None:
        for var in self.__output_variables:
            var = self.__registry.get(var)

            self.__output_data_store.write_chunk(
                data=var,
                parameter_name=var.name,
                start_time=start_time,
                end_time=end_time,
            )
