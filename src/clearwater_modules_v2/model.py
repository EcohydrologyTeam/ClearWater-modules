from processes.base import Process
from clearwater_data.variables import VariableRegistry
from datetime import datetime, timedelta
from typing import Iterable
from clearwater_data.io.base import (
    DataSource,
    ChunkedDataSource,
    DataStore,
    ChunkedDataStore,
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
        output_store: DataStore | ChunkedDataStore,
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
        self.__output_store: DataStore | ChunkedDataStore = output_store

        # check if we are running chunked
        self.__chunked_mode: bool = chunk_size is not None
        self.__chunk_size: timedelta | None = chunk_size
        # if chunked mode, we also want to make sure the output is chunked
        if self.__chunked_mode:
            if not isinstance(self.__output_store, ChunkedDataStore):
                raise ValueError(
                    "Output store must be chunked when running the model in chunked mode. "
                    + "Either set chunk_size to None or provide ChunkedDataStore for model output."
                )

    def validate(self) -> None:
        if self.__start_time >= self.__end_time:
            raise ValueError("Start time must be before end time.")

    def run(self) -> None:
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

        for process in self.__processes:
            process.init_process(self, self.__registry)

    def __process_loop_chunked(self) -> None:
        # TODO: this need actual chunking logic
        current_time = self.__start_time
        while current_time < self.__end_time:
            current_time_seconds = current_time.timestamp()
            print(f"Running timestep: {current_time}")
            for process in self.__processes:
                # check if this process should be updated at this timestamp
                if current_time_seconds % process.time_step_seconds == 0:
                    process.run(current_time, self.__registry)
            current_time += self.__time_step
        self.__save_output_model()

    def __process_loop_full(self) -> None:
        current_time = self.__start_time

        while current_time < self.__end_time:
            current_time_seconds = current_time.timestamp()
            for process in self.__processes:
                # check if this process should be updated at this timestamp
                # Process should be calculated is current_time + process.time_step_frequency
                if current_time_seconds % process.time_step_seconds == 0:
                    process.run(current_time, self.__registry)
            current_time += self.__time_step
        self.__save_output_model()

    def __save_output_model(self) -> None:
        for var in self.__output_variables:
            var = self.__registry.get(var)
            self.__output_store.write(
                data=var,
                parameter_name=var.name,
            )
