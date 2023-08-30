"""Stored base types shared by all sub-modules."""

import abc
import xarray as xr
from typing import (
    TypedDict,
    Callable,
    Optional,
)

Equation = Callable[..., float]
ConstantsDict = dict[str, float | int | bool]


class Model(abc.ABC):
    """Abstract base class for all processes."""
    
    __required_inputs: TypedDict
    __default_constants: TypedDict
    __default_equations: list[Equation]
    
    def __init__(
        self,
        required_inputs: ConstantsDict | xr.Dataset,
        update_constants: Optional[ConstantsDict],
        update_equations: Optional[list[Equation]],
        **kwargs,
    ) -> None:
        ...

    @abc.abstractmethod
    def validate_inputs(self) -> None:
        """Validate inputs."""
        ...

    @abc.abstractproperty
    def constants(self) -> ConstantsDict| xr.Dataset:
        """Return the constants."""
        ... 

    @abc.abstractmethod
    def run(self) -> ConstantsDict | xr.Dataset:
        """Run the process."""
        ...

