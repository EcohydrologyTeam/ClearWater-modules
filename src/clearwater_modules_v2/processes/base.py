from abc import ABC, abstractmethod
from datetime import datetime, timedelta
from variables import VariableRegistry


class Process(ABC):
    """
    Base class for processes. Defines a class-level variable registry and a run method to be implemented by subclasses.
    """

    # Class-level definition of variables associated with this process
    variables = []
    time_step_seconds: int

    def __init__(self, time_step_frequency: timedelta) -> None:
        self.time_step_frequency = time_step_frequency
        self.time_step_seconds = time_step_frequency.total_seconds()

    def init_process(self, registry: VariableRegistry) -> None:
        """
        Initialize of the process.
        """

    @abstractmethod
    def run(self, time_step: datetime, registry: VariableRegistry) -> None:
        """
        Run the process. To be implemented by subclasses.
        """
        raise NotImplementedError
