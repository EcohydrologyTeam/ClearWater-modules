from variables.base import Variable
from datetime import datetime


class FloatVariable(Variable):
    """
    A variable that stores a single float value.
    """

    def __init__(self, value: float):
        self.value = value

    def get(self) -> float:
        """
        Get a reference to the variable's value
        """
        return self.value

    def get_at_time(self, time: datetime) -> float:
        """
        Get a reference to the variable's value at a specific time
        """
        # single floating value is time independent
        return self.get()
