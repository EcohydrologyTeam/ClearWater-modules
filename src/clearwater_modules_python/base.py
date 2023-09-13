"""Stored base types shared by all sub-modules."""
import warnings
from re import S
import xarray as xr
from dataclasses import dataclass
import clearwater_modules_python.utils as utils
import clearwater_modules_python.sorter as sorter
from clearwater_modules_python.shared.types import (
    InitialVariablesDict,
    Variable,
)
from typing import (
    runtime_checkable,
    Protocol,
    Optional,
    Iterable,
)

@runtime_checkable
class CanRegisterVariable(Protocol):

    def register_variable(self, variable: Variable) -> None:
        """Register a variable with the model."""
        ...

class Model(CanRegisterVariable):
    _variables: list[Variable] = []

    def __init__(
        self, 
        initial_state_values: InitialVariablesDict,
        static_variable_values: InitialVariablesDict,
        track_dynamic_variables: bool = False,
        hotstart_dataset: xr.Dataset | None = None,
        time_dim: Optional[str] = None,
    ) -> None:
        """Initialize the model, should be accessed by subclasses.

        Args:
            initial_values: A dict with variable names as keys, and initial 
                values for static and (optionally) state variables as values.
        """
        self.initial_state_values = initial_state_values
        for state_var in self.state_variables:
            if state_var.name not in self.initial_state_values.keys():
                 raise ValueError(
                     f'No initial value found for state variable: {state_var.name}.'
                 )


        self.static_variable_values = static_variable_values
        self.track_dynamic_variables = track_dynamic_variables
        if not time_dim:
            time_dim = 'time_step'
        self.time_dim = time_dim
        
        if hotstart_dataset is not None:
            if self.time_dim not in hotstart_dataset.dims:
                raise ValueError(
                    f'Hotstart dataset must have a {self.time_dim} dimension.'
                )
            self.dataset = hotstart_dataset
        else:
            # initialize the main model dataset
            self.dataset: xr.Dataset = self.init_state_arrays()

            # make 2-dimensional static variable arrays
            self.init_static_arrays()
        
        self._sorted_variables: list[Variable] = []


    @classmethod
    def get_variable_names(cls) -> list[str]:
        """Return a list of default variable names."""
        return [var.name for var in cls._variables]

    @classmethod
    def register_variable(cls, variable: Variable) -> None:
        """Register a variable with the model."""
        if variable.name not in cls.get_variable_names():
            cls._variables.append(variable)
            cls._sorted_variables = []

    
    @classmethod
    def unregister_variables(cls, variables: str | list[str]) -> None:
        """Unregister a variable with the model."""
        if isinstance(variables, str):
            variables = [variables]
        cls._variables = [
            var for var in cls._variables if var.name not in variables
        ]
    
    @classmethod
    def get_variable(cls, name: str) -> Variable:
        """Returns a variable dataclass by name"""
        for var in cls._variables:
            if var.name == name:
                return var
        raise ValueError(f'No variable found with name: {name}')

    @property
    def all_variables(self) -> list[Variable]:
        """Return a list of variables."""
        return self._variables
    
    @property
    def static_variables(self) -> list[Variable]:
        """Return a list of static variables."""
        return [var for var in self.all_variables if var.use == 'static']

    @property
    def dynamic_variables(self) -> list[Variable]:
        """Return a list of dynamic variables."""
        return [var for var in self.all_variables if var.use == 'dynamic']

    @property
    def state_variables(self) -> list[Variable]:
        """Return a list of state variables."""
        return [var for var in self.all_variables if var.use == 'state']
    
    @property
    def static_variables_names(self) -> list[str]:
        """Return a list of static variable names."""
        return [var.name for var in self.static_variables]
    
    @property
    def dynamic_variables_names(self) -> list[str]:
        """Return a list of dynamic variable names."""
        return [var.name for var in self.dynamic_variables]

    @property
    def state_variables_names(self) -> list[str]:
        """Return a list of state variable names."""
        return [var.name for var in self.state_variables]

    @property
    def sorted_dynamic_variables(self) -> list[Variable]:
        """Return a list of variables in order."""
        if len(self._sorted_variables) == 0:
            self._sorted_variables = sorter.sort_dynamic_variables(
                sorter.split_variables(self.all_variables),
            )
        return self.all_variables


    def run(self) -> xr.Dataset:
        """Run the process."""
        ...
    
    def init_state_arrays(self) -> xr.Dataset:
        """Initializes the state arrays."""
        match_dims: list[str] = []
        data_arrays: dict[str, xr.DataArray] = {}
        
        for k, v in self.initial_state_values.items():
            if k not in self.state_variables_names:
                warnings.warn(
                    f'Variable {k} is not a state variable, skipping.',
                )
                continue
            if not isinstance(v, xr.DataArray):
                match_dims.append(k)
            else:
                utils.validate_arrays(v, *list(data_arrays.values()))
                data_arrays[k] = v
        if len(data_arrays) > 0:
            array_i = list(data_arrays.values())[0]
        else:
            array_i = xr.DataArray(
                [[1.0]],
                dims=['x', 'y'],
                coords=[[1.0], [1.0]],
            )
        for var_name in match_dims:
            variable = self.get_variable(var_name)
            attrs = {
                'long_name': variable.long_name,
                'units': variable.units,
                'description': variable.description,
            }
            data_arrays[var_name] = xr.full_like(
                array_i, 
                self.initial_state_values[var_name],
                dtype=type(self.initial_state_values[var_name]),
                attrs=attrs,
            )
        ds = xr.Dataset(
            data_vars=data_arrays,
            coords=array_i.coords,
        )
        return ds.expand_dims({self.time_dim: [0]}) 
   

    def init_static_arrays(self) -> None:
        """Return a static dataset."""
        for var in self.static_variables:
            if var.name not in self.static_variable_values.keys():
                raise ValueError(
                    f'No initial value found for static variable: {var.name}.'
                )
            if var.name in self.dataset.coords.keys():
                raise ValueError(
                    f'Variable name {var.name} already exists in coords.'
                )
            attrs = {
                'long_name': var.long_name,
                'units': var.units,
                'description': var.description,
            }
            self.dataset[var.name] = xr.full_like(
                self.dataset[list(self.dataset.data_vars)[0]],
                self.static_variable_values[var.name],
                dtype=type(self.static_variable_values[var.name]),
            ).sel({self.time_dim: 0})
            self.dataset[var.name].attrs = attrs


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


