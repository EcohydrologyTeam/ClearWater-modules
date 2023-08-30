"""Stored base types shared by all sub-modules."""

import abc
import xarray as xr
from typing import (
    TypedDict,
    Callable,
    Optional,
)

Equation = Callable[[], float]

class Process(abc.ABC):
    """Abstract base class for all processes."""
    
    __required_inputs: TypedDict
    __default_constants: TypedDict
    __default_equations: list[Equation]
    
    def __init__(
        self,
        required_inputs: TypedDict | xr.Dataset,
        update_constants: Optional[dict[str, float | int | bool]],
        update_equations: Optional[list[Equation]],
        **kwargs,
    ) -> None:
        ...

    @abc.abstractmethod
    def validate_inputs(self) -> None:
        """Validate inputs."""
        ...

    @abc.abstractproperty
    def constants(self) -> TypedDict | xr.Dataset:
        """Return the constants."""
        ... 

    @abc.abstractmethod
    def run(self) -> TypedDict | xr.Dataset:
        """Run the process."""
        ...

