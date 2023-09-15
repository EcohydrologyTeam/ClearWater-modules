"""Stored base types shared by all sub-modules."""
import warnings
import xarray as xr
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
        initial_state_values: Optional[InitialVariablesDict] = None,
        static_variable_values: Optional[InitialVariablesDict] = None,
        track_dynamic_variables: bool = True,
        hotstart_dataset: Optional[xr.Dataset] = None,
        time_dim: Optional[str] = None,
    ) -> None:
        """Initialize the model, should be accessed by subclasses.

        Args:
            initial_state_values: A dict with variable names as keys, and initial 
                state variables as values.
            static_variable_values: A dict with variable names as keys, and static
                values as values. All Model.static_variables must be present.
                NOTE: This may be triggered in a subclass __init__ method.
            track_dynamic_variables: If True, dynamic variables will be tracked
                in the model dataset. If False, they will not be tracked.
            hotstart_dataset: An optional dataset to use as a hotstart. 
                This skips the initialization of the model dataset, and uses the 
                provided dataset instead after validating the presence of all static 
                and state variables.
            time_dim: The name of the time dimension. If not provided, defaults
                to 'time_step'.
        """
        self.initial_state_values = initial_state_values
        self.static_variable_values = static_variable_values
        self.hotstart_dataset = hotstart_dataset
        self.track_dynamic_variables = track_dynamic_variables

        if not time_dim:
            time_dim = 'time_step'
        self.time_dim = time_dim
       
        if isinstance(self.initial_state_values, dict) and isinstance(self.static_variable_values, dict):
            print('Initializing from dicts...')
            self.dataset: xr.Dataset = self._init_dataset_from_dicts(
                initial_state_values=self.initial_state_values,
                static_variable_values=self.static_variable_values,
            )

        elif isinstance(hotstart_dataset, xr.Dataset):
            print('Initializing from hotstart dataset...')
            self.dataset: xr.Dataset = self._init_from_dataset(hotstart_dataset)
            self.hotstart_dataset = None

        else:
            raise ValueError(
                'Must provide either initial state and static values, or a hotstart dataset.'
            )

        self._sorted_variables: list[Variable] = []
    
    def _init_dataset_from_dicts(
        self,
        initial_state_values: InitialVariablesDict,
        static_variable_values: InitialVariablesDict,
    ) -> None:
        """Initialize Model.dataset from dicts."""
        if not isinstance(initial_state_values, dict):
            raise TypeError(
                f'Expected initial state value dict, got {type(self.initial_state_values)} instead.'
            )
        if not isinstance(static_variable_values, dict):
            raise TypeError(
                f'Expected static variable value dict, got {type(self.static_variable_values)} instead.'
            )
        for state_var in self.state_variables:
            if state_var.name not in initial_state_values.keys():
                raise ValueError(
                    f'No initial value found for state variable: {state_var.name}.'
                )

        # initialize the main model dataset
        dataset: xr.Dataset = self._init_state_arrays(initial_state_values)
        dataset: xr.Dataset = self._init_static_arrays(
            dataset, 
            static_variable_values,
        )

        print('Model initialized from input dicts successfully!.')
        return dataset

    def _init_from_dataset(self, hotstart_dataset: xr.Dataset) -> xr.Dataset:
        """Initialize the model from a hotstart dataset."""
        if self.time_dim not in hotstart_dataset.dims:
            raise ValueError(
                f'Hotstart dataset must have a {self.time_dim} dimension.'
            )
        return hotstart_dataset

    def _init_state_arrays(
        self,
        initial_state_values: InitialVariablesDict,
    ) -> xr.Dataset:
        """Initializes the state arrays."""
        match_dims: list[str] = []
        data_arrays: dict[str, xr.DataArray] = {}

        for k, v in initial_state_values.items():
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
                initial_state_values[var_name],
                dtype=type(initial_state_values[var_name]),
            )
            data_arrays[var_name].attrs = attrs
        ds = xr.Dataset(
            data_vars=data_arrays,
            coords=array_i.coords,
        )
        return ds.expand_dims({self.time_dim: [0]})

    def _init_static_arrays(
        self, 
        dataset: xr.Dataset,
        static_variable_values: InitialVariablesDict,
    ) -> xr.Dataset:
        """Broadcasts static variables to an existing dataset.
        
        Args:
            dataset: The dataset to broadcast to.
            static_variable_values: A dictionary of static variable names and 
                values (either float/bool/int or a xarray.DataArray).

        Returns:
            The dataset with the static variables added.
        """
        for var in self.static_variables:
            if var.name not in static_variable_values.keys():
                raise ValueError(
                    f'No initial value found for static variable: {var.name}.'
                )
            if var.name in dataset.coords.keys():
                raise ValueError(
                    f'Variable name {var.name} already exists in coords.'
                )
            attrs = {
                'long_name': var.long_name,
                'units': var.units,
                'description': var.description,
            }
            dataset[var.name] = xr.full_like(
                dataset[
                    list(dataset.data_vars)[0]
                ].isel({self.time_dim: 0}),
                static_variable_values[var.name],
                dtype=type(static_variable_values[var.name]),
            )
            dataset[var.name].attrs = attrs
        return dataset


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
        cls._sorted_variables = []

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
    def computation_order(self) -> list[Variable]:
        """Return a list of variables to compute in order (dynamic + state)."""
        if len(self._sorted_variables) == 0:
            self._sorted_variables = sorter.sort_variables_for_computation(
                sorter.split_variables(self.all_variables),
            )
        return self._sorted_variables
    
    @property
    def _update_vars(self)-> list[str]:
        """Return a list of variables to update."""
        if self.track_dynamic_variables:
            return self.dynamic_variables_names + self.state_variables_names 
        else:
            return self.state_variables_names

    def increment_timestep(
        self,
        update_state_values: Optional[dict[str, xr.DataArray]] = None,
    ) -> xr.Dataset:
        """Run the process."""
        if update_state_values is None:
            update_state_values = {}
        last_timestep: int = self.dataset[self.time_dim].values[-1]
        timestep_ds: xr.Dataset = self.dataset.isel(
            {self.time_dim: -1}).copy(deep=True)

        # update the state variables as necessary (i.e. interacting w/ other models)
        for var_name, value in update_state_values.items():
            utils.validate_arrays(value, timestep_ds[var_name])
            timestep_ds[var_name] = value

        # compute the dynamic variables in order
        timestep_ds = utils.iter_computations(
            timestep_ds,
            self.computation_order,
        )
        if not self.track_dynamic_variables:
            timestep_ds = timestep_ds.drop_vars(self.dynamic_variables_names)

        timestep_ds = timestep_ds.drop_vars(self.static_variables_names)
        timestep_ds= timestep_ds.expand_dims(
            {self.time_dim: [last_timestep + 1]},
        )
        
        self.dataset = xr.concat(
            [
                self.dataset,
                timestep_ds,
            ],
            dim=self.time_dim,
            data_vars='minimal',
        )
        
        # add dynamic variable attributes
        if self.track_dynamic_variables:
            for var in self.dynamic_variables:
                if var.name in self.dataset.data_vars.keys():
                    self.dataset[var.name].attrs = {
                        'long_name': var.long_name,
                        'units': var.units,
                        'description': var.description,
                    }

        return self.dataset

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
