"""Stored base types shared by all sub-modules."""

import abc
import xarray as xr
from dataclasses import dataclass
from typing import (
    TypedDict,
    Callable,
    Optional,
    Literal,
)

Process = Callable[..., float]
ParametersDict = dict[str, float | int | bool]
VariableTypes = Literal['static', 'dynamic', 'state']

@dataclass(slots=True, frozen=True)
class Variable:
    """Variable type."""
    name: str
    long_name: str
    units: str
    description: str
    use: VariableTypes
    process: Optional[Process] = None


class Model(abc.ABC):
    """Abstract base class for all processes."""

    __required_inputs: TypedDict
    __default_parameters: TypedDict
    __default_variables: list[Variable]

    def __init__(
        self,
        required_inputs: ParametersDict | xr.Dataset,
        update_parameters: Optional[ParametersDict] = None,
        track_non_state_variables: bool = False,
        ignore_variables: Optional[list[Variable | str]] = None,
        **kwargs,
    ) -> None:
        """Initialize the model.

        Args:
            required_inputs: Required inputs for the model.
            update_parameters: Optional parameters to update.
            track_non_state_variables: Track non-state variables.
            ignore_variables: Variables to ignore.
        """
        ...

    @abc.abstractmethod
    def validate_inputs(self) -> None:
        """Validate inputs."""
        ...

    @abc.abstractproperty
    def parameters(self) -> ParametersDict | xr.Dataset:
        """Return the parameters."""
        ...

    @abc.abstractmethod
    def run(self) -> xr.Dataset:
        """Run the process."""
        ...
