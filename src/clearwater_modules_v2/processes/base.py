from abc import ABC, abstractmethod
from datetime import datetime, timedelta
import functools
from typing import Callable

from clearwater_data.variables import VariableRegistry

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from model import Model


class Process(ABC):
    """
    Base class for processes. Defines a class-level variable registry and a run method to be implemented by subclasses.
    """

    # Class-level definition of variables associated with this process
    variables = []
    time_step_seconds: int

    def __init__(self, time_step: timedelta) -> None:
        self.time_step = time_step
        self.time_step_seconds = self.time_step.total_seconds()

    @classmethod
    def from_config(
        cls, config: dict, variable_registry: VariableRegistry
    ) -> "Process":
        return cls(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """
        Initialize of the process.
        """
        # base method assumes no initialization is needed
        pass

    def validate(self, registry: VariableRegistry) -> None:
        """
        Validate the process.
        """
        for variable in self.variables:
            if variable not in registry:
                raise ValueError(
                    f"Variable {variable} not found. Are you sure you provided a valid configuration for {variable}?"
                )

    @abstractmethod
    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """
        Run the process. To be implemented by subclasses.
        """
        raise NotImplementedError

    def process_name(self) -> str:
        return self.__class__.__name__


class ProcessFactory:
    processes: dict[str, callable] = {}

    @classmethod
    def from_config(
        cls, process_name: str, config: dict, variable_registry: VariableRegistry
    ) -> "Process":
        if process_name not in cls.processes:
            raise ValueError(
                f"Process type {process_name} not registered. Did you register the process at the {__name__}"
            )
        return cls.processes[process_name](config, variable_registry)

    @classmethod
    def register(cls, process_name: str):
        def init_method(from_config_method: Callable) -> Callable:
            cls.processes[process_name] = from_config_method

            @functools.wraps(from_config_method)
            def from_config_method(config: dict) -> "Process":
                return from_config_method(config)

            return from_config_method

        return init_method
