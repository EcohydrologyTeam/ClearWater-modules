from abc import ABC, abstractmethod
from enum import Enum


class Variable(ABC):
    """
    Base class for variables. Defines a class-level variable registry and a run method to be implemented by subclasses.
    """
    pass

#look at metpy
#look at pint
class Units(Enum):
    CM = ("cm", "100")
    M = ("m", "1")
    MM = ("mm", "1000")
    
Units.CM.values[1]

def convert_units(value: float, from_unit: Units, to_unit: Units):
    return value * (from_unit.value[1] / to_unit.value[1])