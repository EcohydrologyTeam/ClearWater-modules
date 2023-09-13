"""Stored base types shared by all sub-modules."""

import xarray as xr
from dataclasses import dataclass
from typing import (
    Iterable,
    Protocol,
    TypedDict,
    Callable,
    Optional,
    Literal,
    runtime_checkable,
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


class SplitVariablesDict(TypedDict):
    """A dict containing all variables split by type.

    Attributes:
        static: A list of static variables (out).
        dynamic: A list of dynamic variables (in).
        state: A list of state variables (in/out).
    """
    static: list[Variable]
    dynamic: list[Variable]
    state: list[Variable]


@runtime_checkable
class CanRegisterVariable(Protocol):

    def register_variable(self, variable: Variable) -> None:
        """Register a variable with the model."""
        ...

class Model(CanRegisterVariable):
    _variables: list[Variable] = []
    
    @classmethod
    def get_variable_names(cls) -> list[str]:
        """Return a list of default variable names."""
        return [var.name for var in cls._variables]

    @classmethod
    def register_variable(cls, variable: Variable) -> None:
        """Register a variable with the model."""
        if variable.name not in cls.get_variable_names():
            cls._variables.append(variable)
    
    @classmethod
    def unregister_variables(cls, variables: str | list[str]) -> None:
        """Unregister a variable with the model."""
        if isinstance(variables, str):
            variables = [variables]
        cls._variables = [
            var for var in cls._variables if var.name not in variables
        ]
    
    @property
    def all_variables(self) -> list[Variable]:
        """Return a list of variables."""
        return self._variables
    
    @property
    def static_variables(self) -> list[Variable]:
        """Return a list of static variables."""
        return [var for var in self._variables if var.use == 'static']

    @property
    def dynamic_variables(self) -> list[Variable]:
        """Return a list of dynamic variables."""
        return [var for var in self._variables if var.use == 'dynamic']

    @property
    def state_variables(self) -> list[Variable]:
        """Return a list of state variables."""
        return [var for var in self._variables if var.use == 'state']

    def validate_inputs(self) -> None:
        """Validate inputs."""
        ...

    def run(self) -> xr.Dataset:
        """Run the process."""
        ...
    
def register_variable(models: CanRegisterVariable | Iterable[CanRegisterVariable]):
    """A decorator to register a variable with a model."""
    if not isinstance(models, Iterable):
        models = [models]
    def decorator(cls):
        def wrapper(*args, **kwargs):
            variable = cls(*args, **kwargs)
            if not issubclass(cls, Variable):
                raise TypeError(
                    f'Expected a Variable, got {type(cls)} instead.'
                )
            for model in models:
                model.register_variable(variable)
            return variable
        return wrapper
    return decorator


