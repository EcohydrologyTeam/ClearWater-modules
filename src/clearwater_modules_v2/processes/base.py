from abc import ABC, abstractmethod
from datetime import datetime
from variables import VariableRegistry

class Process(ABC):
    """
    Base class for processes. Defines a class-level variable registry and a run method to be implemented by subclasses.
    """
    # Class-level definition of variables associated with this process
    variables = []
    time_step_frequency = "5min"

    def init_process(self, variables:VariableRegistry) -> None:
        """
        Initialize of the process.
        """
        pass    

    @abstractmethod
    def run(self, time_step:datetime, variables:VariableRegistry) -> None:
        """
        Run the process. To be implemented by subclasses.
        """
        raise NotImplementedError
