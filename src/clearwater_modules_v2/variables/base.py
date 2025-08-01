from abc import ABC, abstractmethod

from custom_types import ArrayLike
from datetime import datetime


class Variable(ABC):
    """
    Base class for variables.
    """

    def get(self) -> object:
        """
        Get a reference to the variable's value
        """
        raise NotImplementedError

    def get_at_time(self, time: datetime) -> object:
        """
        Get a reference to the variable's value at a specific time
        """
        raise NotImplementedError

    # TODO: Consider the notion of units
    # look at metpy
    # look at pint
