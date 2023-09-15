"""This tests makes sure that shared model functionality works as expected."""
import pytest
import numpy as np
from tests.conftest import initial_array
import xarray as xr
from clearwater_modules_python.base import (
    Model,
    Variable,
    InitialVariablesDict,
)


class MockModel(Model):
    ...


def state_process(dynamic_2: float) -> float:
    return dynamic_2 * 2


@pytest.fixture(scope='module')
def initial_static_values() -> InitialVariablesDict:
    return {'a': 1.0, 'b': 2.0}


@pytest.fixture(scope='module')
def initial_state_values(
    initial_array: xr.DataArray,
) -> InitialVariablesDict:
    return {'state_variable': initial_array}


@pytest.fixture(scope='module')
def state_variable() -> Variable:
    return Variable(
        name='state_variable',
        long_name='State Variable',
        units='m',
        description='A state variable.',
        use='state',
        process=state_process,
    )


@pytest.fixture(scope='function')
def model(
    static_variables: list[Variable],
    dynamic_variables: list[Variable],
    state_variable: Variable,
    initial_static_values: InitialVariablesDict,
    initial_state_values: InitialVariablesDict,
) -> Model:
    for var in static_variables + dynamic_variables:
        MockModel.register_variable(var)
    MockModel.register_variable(state_variable)

    model_instance = MockModel(
        initial_state_values=initial_state_values,
        static_variable_values=initial_static_values,
    )
    assert isinstance(model_instance, Model)
    return model_instance


def test_model_static_variables(
    model: Model,
    static_variables: list[Variable],
) -> None:
    """Test the model."""
    assert isinstance(model.static_variables, list)
    print(f'Static variables: {model.static_variables}\n')
    assert len(model.static_variables) == len(static_variables)
    assert isinstance(model.static_variables[0], Variable)
    assert isinstance(model.static_variables_names, list)
    assert len(model.static_variables_names) == len(static_variables)
    assert isinstance(model.static_variables_names[0], str)
    assert model.static_variables_names[0] == model.static_variables[0].name


def test_model_dynamic_variables(
    model: Model,
    dynamic_variables: list[Variable],
) -> None:
    assert isinstance(model.dynamic_variables, list)
    print(f'Dynamic variables {model.dynamic_variables}\n')
    assert len(model.dynamic_variables) == len(dynamic_variables)
    assert isinstance(model.dynamic_variables[0], Variable)
    assert isinstance(model.dynamic_variables_names, list)
    assert len(model.dynamic_variables_names) == len(dynamic_variables)
    assert isinstance(model.dynamic_variables_names[0], str)
    assert model.dynamic_variables_names[0] == model.dynamic_variables[0].name


def test_model_state_variables(
    model: Model,
) -> None:
    assert isinstance(model.state_variables, list)
    print(f'State variables: {model.state_variables}\n')
    assert len(model.state_variables) == 1
    assert isinstance(model.state_variables[0], Variable)
    assert isinstance(model.state_variables_names, list)
    assert len(model.state_variables_names) == 1
    assert isinstance(model.state_variables_names[0], str)
    assert model.state_variables_names[0] == model.state_variables[0].name


def test_initial_state(model: Model, state_variable: Variable) -> None:
    """Test the initial state."""
    assert state_variable.name in model.dataset.data_vars
    assert len(model.dataset[state_variable.name].dims) == 3
    assert 'time_step' in model.dataset[state_variable.name].dims
    assert model.dataset[state_variable.name].shape == (1, 10, 10)


def test_static_array(model: Model) -> None:
    """Test the static array."""
    for var in model.static_variables:
        assert var.name in model.dataset.data_vars
        assert len(model.dataset[var.name].dims) == 2
        assert model.dataset[var.name].shape == (10, 10)
        assert 'long_name' in model.dataset[var.name].attrs
        assert 'units' in model.dataset[var.name].attrs
        assert 'description' in model.dataset[var.name].attrs


def test_state_array(model: Model) -> None:
    """Test the state array."""
    for var in model.state_variables:
        assert var.name in model.dataset.data_vars
        assert len(model.dataset[var.name].dims) == 3
        assert model.dataset[var.name].shape == (1, 10, 10)
        assert 'long_name' in model.dataset[var.name].attrs
        assert 'units' in model.dataset[var.name].attrs
        assert 'description' in model.dataset[var.name].attrs


def test_variable_compution_order(model: Model) -> None:
    """Test that dynamic variables are sorted correctly."""
    sorted = model.computation_order
    print(sorted)
    assert len(sorted) == (
        len(model.dynamic_variables) +
        len(model.state_variables)
    )
    assert sorted[0].name == 'dynamic_0'
    assert sorted[-1].name == 'state_variable'


def test_computation_with_dynamics(model: Model) -> None:
    """Test that dynamic variables are saved to the dataset.""" 
    model.track_dynamic_variables = True
    ds: xr.Dataset = model.increment_timestep()
    assert isinstance(ds, xr.Dataset)
    assert len(model.dataset[model.time_dim]) == 2
    assert 'dynamic_0' in model.dataset.data_vars


def test_computation_no_dynamics(model: Model) -> None:
    """Test that dynamic variables not are saved to the dataset.""" 
    model.track_dynamic_variables = False
    ds = model.increment_timestep()
    assert isinstance(ds, xr.Dataset)
    assert len(model.dataset[model.time_dim]) == 2
    assert 'dynamic_0' not in model.dataset.data_vars


def test_static_variable_dims(model: Model) -> None:
    """Test that static variables remain 2-dimensional."""
    ds = model.increment_timestep()
    for var_name in model.static_variables_names:
        assert var_name in ds.data_vars
        assert len(ds[var_name].dims) == 2


def test_variable_attributes(model: Model) -> None:
    """Tests that all variables in the Model.dataset have the correct attributes."""
    model.increment_timestep()
    for var in model.dataset.data_vars:
        assert 'long_name' in model.dataset[var].attrs
        assert 'units' in model.dataset[var].attrs
        assert 'description' in model.dataset[var].attrs


def test_model_hotstart(model: Model) -> None:
    """Test if the hotstart works."""
    ds = model.increment_timestep()
    ds.attrs['hotstart'] = True

    hotstart_model = MockModel(
        hotstart_dataset=ds,
    )
    assert isinstance(hotstart_model, Model)
    assert len(hotstart_model.dataset[model.time_dim]) == 2
    assert model.dataset.attrs.get('hotstart') == True
