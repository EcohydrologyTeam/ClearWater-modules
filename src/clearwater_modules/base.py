"""Stored base types shared by all sub-modules."""
import warnings
import xarray as xr
import numpy as np
import clearwater_modules.utils as utils
import clearwater_modules.sorter as sorter
from clearwater_modules.shared.types import (
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
        time_steps: int,
        initial_state_values: Optional[InitialVariablesDict] = None,
        static_variable_values: Optional[InitialVariablesDict] = None,
        updateable_static_variables: Optional[list[str]] = None,
        track_dynamic_variables: bool = False,
        hotstart_dataset: Optional[xr.Dataset] = None,
        time_dim: Optional[str] = None,
        timestep: Optional[int] = 0,
    ) -> None:
        """Initialize the model, should be accessed by subclasses.

        Args:
            time_steps: An integer to indicate the number of timesteps to run.
            initial_state_values: A dict with variable names as keys, and initial
                state variables as values.
            static_variable_values: A dict with variable names as keys, and static
                values as values. All Model.static_variables must be present.
                NOTE: This may be triggered in a subclass __init__ method.
            updateable_static_variables: A list of static variable names that should
                be converted to state variables. This allows them to be updated
                between timesteps. Note that these still don't have process functions.
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
        self._track_dynamic_variables = track_dynamic_variables
        self.timestep = timestep
        self.time_steps = time_steps + 1  # xarray indexing
        self.temporal_variables: list = []

        if not time_dim:
            time_dim = 'time_step'
        self.time_dim = time_dim

        if not isinstance(updateable_static_variables, list):
            updateable_static_variables = []
        self.updateable_static_variables = updateable_static_variables
        self.__non_updateable_static_variables: list[str] | None = None

        # create list of temporal variables
        if self._track_dynamic_variables:
            self.temporal_variables = self.state_variables_names + \
                    self.updateable_static_variables + self.dynamic_variables_names
        else:
            self.temporal_variables = self.state_variables_names + \
                    self.updateable_static_variables

        if isinstance(self.initial_state_values, dict) and isinstance(self.static_variable_values, dict):
            print('Initializing from dicts...')
            self.dataset: xr.Dataset = self._init_dataset_from_dicts(
                initial_state_values=self.initial_state_values,
                static_variable_values=self.static_variable_values,
                updateable_static_variables=self.updateable_static_variables,
                time_steps=self.time_steps,
            )

        elif isinstance(hotstart_dataset, xr.Dataset):
            print('Initializing from hotstart dataset...')
            self.dataset: xr.Dataset = self._init_from_dataset(
                hotstart_dataset,
                time_steps,
            )
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
        updateable_static_variables: list[str],
        time_steps: int,
    ) -> xr.Dataset:
        """Initialize Model.dataset from dicts."""
        if not isinstance(initial_state_values, dict):
            raise TypeError(
                f'Expected initial state value dict, got {type(initial_state_values)} instead.'
            )
        if not isinstance(static_variable_values, dict):
            raise TypeError(
                f'Expected static variable value dict, got {type(static_variable_values)} instead.'
            )
        static_variable_values = static_variable_values.copy()

        for state_var in self.state_variables:
            if state_var.name not in initial_state_values.keys():
                raise ValueError(
                    f'No initial value found for state variable: {state_var.name}.'
                )

        # reassign updateable_static_variables to state variables
        for static in updateable_static_variables:
            if static not in static_variable_values.keys():
                warnings.warn(
                    f'Variable name = {static} is not a static variable, skipping.'
                )
                continue
            else:
                initial_state_values[static] = static_variable_values.pop(static)

        # initialize the main model dataset
        dataset: xr.Dataset = self._init_state_arrays(
            initial_state_values,
            time_steps,
        )
        dataset: xr.Dataset = self._init_static_arrays(
            dataset,
            static_variable_values,
            time_steps,
        )
        if self._track_dynamic_variables:
            dataset: xr.Dataset = self._init_dynamic_arrays(
                dataset,
            )

        print('Model initialized from input dicts successfully!.')
        return dataset

    def _init_from_dataset(
        self,
        hotstart_dataset: xr.Dataset,
        time_steps: int
    ) -> xr.Dataset:
        """Initialize the model from a hotstart dataset."""
        if self.time_dim not in hotstart_dataset.dims:
            raise ValueError(
                f'Hotstart dataset must have a {self.time_dim} dimension.'
            )
        else:
            coords = {
                key: value if key != self.time_dim
                else np.arange(time_steps)
                for key, value in hotstart_dataset.coords.items()
            }

            new_hotstart_dataset = xr.Dataset(
                data_vars={
                    var_name: (
                        hotstart_dataset[var_name].dims,
                        np.full(
                            tuple(
                                hotstart_dataset[var_name].sizes[dim]
                                for dim in hotstart_dataset[var_name].dims
                            ),
                            np.nan
                        )
                    )
                    for var_name, _ in hotstart_dataset.data_vars.items()
                },
                coords={
                    **coords
                }
            )

            # set temporal variables to the last timestep of the hotstart dataset
            new_hotstart_dataset[self.temporal_variables].loc[
                {self.time_dim: 0}
            ] = hotstart_dataset[self.temporal_variables].isel(
                {self.time_dim: -1}
            )

            new_hotstart_dataset[self._non_updateable_static_variables] = hotstart_dataset[self._non_updateable_static_variables]

        return new_hotstart_dataset

    def _init_state_arrays(
        self,
        initial_state_values: InitialVariablesDict,
        time_steps: int,
    ) -> xr.Dataset:
        """Initializes the state arrays."""
        match_dims: list[str] = []
        data_arrays: dict[str, xr.DataArray] = {}
        coords: dict = {}
        add_data: list[str] = []

        for k, v in initial_state_values.items():
            if k not in (self.state_variables_names + self.updateable_static_variables):
                warnings.warn(
                    f'Variable {k} is not a state variable, skipping.',
                )
                continue
            if not isinstance(v, xr.DataArray):
                match_dims.append(k)
            else:
                utils.validate_arrays(v, *list(data_arrays.values()))
                data_arrays[k] = v
                coords = coords | dict(data_arrays[k].coords.items())
                add_data.append(k)
        if len(data_arrays) > 0:
            ds = xr.Dataset(
                data_vars={
                    k: (
                        (self.time_dim,) + data_arrays[k].dims,
                        np.full(
                            (time_steps,) + tuple(data_arrays[k].sizes[dim] for dim in data_arrays[k].dims),
                            np.nan
                        )
                    )
                    for k in data_arrays.keys()
                },
                coords={
                    self.time_dim: np.arange(time_steps),
                    **coords,
                }
            )
        else:
            ds = xr.Dataset(
                data_vars={
                    k: (
                        (self.time_dim, 'x', 'y'),
                        np.full((time_steps, 1, 1), np.nan)
                    )
                    for k in match_dims
                },
                coords={
                    self.time_dim: np.arange(time_steps),
                    'x': [1.0],
                    'y': [1.0],
                }
            )

        for var_name in match_dims + add_data:
            variable = self.get_variable(var_name)
            attrs = {
                'long_name': variable.long_name,
                'units': variable.units,
                'description': variable.description,
            }

            if var_name not in data_arrays.keys():
                ds[var_name] = xr.DataArray(
                    np.full(
                        tuple(ds.sizes[dim] for dim in ds.dims),
                        np.nan
                    ),
                    dims=ds.dims
                )

                ds[var_name].loc[{self.time_dim: 0}] = xr.full_like(
                    ds[var_name].isel({self.time_dim: 0}),
                    initial_state_values[var_name],
                    dtype=type(initial_state_values[var_name]),
                )

            else:
                ds[var_name].loc[{self.time_dim: 0}] = initial_state_values[var_name]

            ds[var_name].attrs = attrs

        return ds  # ds.expand_dims({self.time_dim: np.arange(time_steps)})

    def _init_static_arrays(
        self,
        dataset: xr.Dataset,
        static_variable_values: InitialVariablesDict,
        time_steps: int,
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
            if var.name in self.updateable_static_variables:
                continue
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

    def _init_dynamic_arrays(
            self,
            dataset: xr.Dataset,
    ) -> xr.Dataset:
        """Initialize dynamic variables."""
        k = self.state_variables_names[0]
        for dynamic_variable in self.dynamic_variables_names:
            dataset[dynamic_variable] = xr.DataArray(
                np.full(
                        tuple(
                            dataset[k].sizes[dim]
                            for dim in dataset[k].dims),
                        np.nan
                    ),
                dims=dataset[k].dims
            )

        for var in self.dynamic_variables:
            if var.name in dataset.data_vars.keys():
                dataset[var.name].attrs = {
                    'long_name': var.long_name,
                    'units': var.units,
                    'description': var.description,
                }

        return dataset

    @classmethod
    def get_variable_names(cls) -> list[str]:
        """Return a list of default variable names."""
        return [var.name for var in cls._variables]

    @classmethod
    def get_state_variables(cls) -> list[Variable]:
        """Returns a list of state variable names and types.
        This can be used to inform the 'initial_state_values' argument
        pre-initialization.
        """
        return [var for var in cls._variables if var.use == 'state']

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
        return self.get_state_variables()

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
    def _update_vars(self) -> list[str]:
        """Return a list of variables to update."""
        if self._track_dynamic_variables:
            return self.dynamic_variables_names + self.state_variables_names
        else:
            return self.state_variables_names

    @property
    def _non_updateable_static_variables(self) -> list[str]:
        """Return a list of static variable names that are non-updateable (2-D)."""
        if self.__non_updateable_static_variables is None:
            self.__non_updateable_static_variables = [
                var.name for var in self.static_variables if var.name not in self.updateable_static_variables
            ]
        return self.__non_updateable_static_variables

    @property
    def track_dynamic_variables(self) -> bool:
        """Track dynamic variables property."""
        return self._track_dynamic_variables

    @track_dynamic_variables.setter
    def track_dynamic_variables(self, value: bool) -> bool:
        if self._track_dynamic_variables == value:
            pass
        elif value:
            self._track_dynamic_variables = value
            self.dataset = self._init_dynamic_arrays(
                self.dataset,
            )
            self.temporal_variables = self.temporal_variables + self.dynamic_variables_names

    def _iter_computations(self):
        """Iterate over the computation order."""
        inputs = map(
            lambda x: utils._prep_inputs(
                self.timestep_ds,
                x),
            self.computation_order
        )
        dims = self.timestep_ds.dims

        for name, func, arrays in inputs:
            array: np.ndarray = func(*arrays)
            self.timestep_ds[name] = (dims, array)

    def increment_timestep(
        self,
        update_state_values: Optional[dict[str, xr.DataArray]] = None,
    ) -> xr.Dataset:
        """Run the process."""
        self.timestep += 1

        if update_state_values is None:
            update_state_values = {}

        # by default, set current timestep equal to last timestep
        self.timestep_ds: xr.Dataset = self.dataset.isel(
            {self.time_dim: self.timestep - 1}
        )

        # update the state variables as necessary (i.e. interacting w/ other models)
        for var_name, value in update_state_values.items():
            if var_name not in (self.state_variables_names + self.updateable_static_variables):
                raise ValueError(
                    f'Variable {var_name} cannot be updated between timesteps, skipping.',
                )
            utils.validate_arrays(value, self.timestep_ds[var_name])
            self.timestep_ds[var_name] = value

        # compute the dynamic variables in order
        self._iter_computations()

        if not self._track_dynamic_variables:
            self.timestep_ds = self.timestep_ds.drop_vars(
                self.dynamic_variables_names
            )
        self.timestep_ds = self.timestep_ds.drop_vars(
            self._non_updateable_static_variables
        )

        self.dataset[self.temporal_variables].loc[
            {self.time_dim: self.timestep}
        ] = self.timestep_ds
    
        return self.dataset


def register_variable(
    models: CanRegisterVariable | Iterable[CanRegisterVariable]
):
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
