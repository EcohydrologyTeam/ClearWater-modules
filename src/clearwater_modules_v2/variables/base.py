from abc import ABC, abstractmethod
from enum import Enum


class Variable(ABC):
    """
    Base class for variables. Defines a class-level variable registry and a run method to be implemented by subclasses.
    """

    pass


# TODO: Consider the notion of units
# look at metpy
# look at pint
