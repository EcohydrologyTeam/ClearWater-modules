"""This tests makes sure that shared model functionality works as expected."""
import pytest
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


@pytest.fixture(scope='module')
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
    return model_instance

def test_model_instance(
    model: Model,
    static_variables: list[Variable],
    dynamic_variables: list[Variable],
) -> None:
    """Test the model."""
    assert isinstance(model, Model)
    assert isinstance(model.static_variables, list)
    assert isinstance(model.dynamic_variables, list)
    assert isinstance(model.state_variables, list)
    print(f'Static variables: {model.static_variables}\n')
    assert len(model.static_variables) == len(static_variables)
    print(f'Dynamic variables {model.dynamic_variables}\n')
    assert len(model.dynamic_variables) == len(dynamic_variables)
    print(f'State variables: {model.state_variables}\n')
    assert len(model.state_variables) == 1

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

def test_state_array(model: Model) -> None:
    """Test the state array."""
    for var in model.state_variables:
        assert var.name in model.dataset.data_vars
        assert len(model.dataset[var.name].dims) == 3
        assert model.dataset[var.name].shape == (1, 10, 10)
