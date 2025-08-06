from processes.base import Process
from variables import VariableRegistry, Variable
from datetime import datetime, timedelta
from typing import Iterable


class Model:
    def __init__(
        self,
        processes: tuple[Process],
        variable_registry: VariableRegistry,
        start_time: datetime,
        end_time: datetime,
        time_step: timedelta,
        output_variables: Iterable[str | Variable],
    ) -> None:
        self.__processes = processes
        self.__variables = variable_registry
        self.__start_time = start_time
        self.__end_time = end_time
        self.__time_step = time_step
        self.__output_variables = output_variables

    def validate(self) -> None:
        if self.__start_time >= self.__end_time:
            raise ValueError("Start time must be before end time.")

    def run(self) -> None:
        self.__init_model()
        self.__process_loop()
        self.__finalize_model()

    def __init_model(self) -> None:
        for process in self.__processes:
            process.init_process(self.__variables)

    def __process_loop(self) -> None:
        current_time = self.__start_time
        while current_time < self.__end_time:
            current_time_seconds = current_time.timestamp()
            for process in self.__processes:
                # check if this process should be updated at this timestamp
                if current_time_seconds % process.time_step_seconds == 0:
                    process.run(current_time, self.__variables)
            current_time += self.__time_step

    def __finalize_model(self) -> None:
        for var in self.__output_variables:
            if isinstance(var, str):
                var = self.__variables.get(var)
            var.save()
