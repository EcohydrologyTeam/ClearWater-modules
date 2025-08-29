from processes.base import Process
from clearwater_data.variables import VariableRegistry
from datetime import datetime, timedelta
from typing import Iterable
from clearwater_data.io.base import DataSource, ChunkedDataSource


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
        chunk_size: timedelta | None = None,
    ) -> None:
        self.__processes = processes
        self.__registry = variable_registry
        self.__variable_data_sources = variable_data_sources
        self.__start_time = start_time
        self.__end_time = end_time
        self.__time_step = time_step
        self.__output_variables = output_variables
        self.__chunk_size = chunk_size

    def validate(self) -> None:
        if self.__start_time >= self.__end_time:
            raise ValueError("Start time must be before end time.")

    def run(self) -> None:
        self.__init_model()
        self.__process_loop()
        self.__finalize_model()

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
            process.init_process(self.__registry)

    def __process_loop(self) -> None:
        current_time = self.__start_time
        while current_time < self.__end_time:
            current_time_seconds = current_time.timestamp()
            for process in self.__processes:
                # check if this process should be updated at this timestamp
                if current_time_seconds % process.time_step_seconds == 0:
                    process.run(current_time, self.__registry)
            current_time += self.__time_step

    def __finalize_model(self) -> None:
        for var in self.__output_variables:
            var = self.__registry.get(var)
            var.save()
