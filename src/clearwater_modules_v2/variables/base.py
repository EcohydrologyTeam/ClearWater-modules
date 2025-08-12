from abc import ABC, abstractmethod

from custom_types import ArrayLike
from datetime import datetime


class Variable(ABC):
    """
    Base class for variables.
    """

    @property
    @abstractmethod
    def time_dimension(self) -> str | None:
        """
        Get the time dimension of the variable.
        """
        raise NotImplementedError

    def get(self) -> ArrayLike:
        """
        Get a reference to the variable's value
        """
        raise NotImplementedError

    def get_at_time(self, time: datetime) -> ArrayLike:
        """
        Get a reference to the variable's value at a specific time
        """
        raise NotImplementedError

    # TODO: Consider the notion of units
    # look at metpy
    # look at pint
