"""This tests makes sure that shared model functionality works as expected."""
import pytest
from tests.conftest import initial_array
import xarray as xr
from clearwater_modules_python.base import (
    Model,
    Variable,
    CanRegisterVariable,
)

class MockModel(Model):
    def __init__(self, initial_state: xr.DataArray):
        self.initial_state = initial_state    

def state_process(dynamic_2: float) -> float:
    return dynamic_2 * 2

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
    initial_array: xr.DataArray,
) -> CanRegisterVariable:


    for var in static_variables + dynamic_variables:
        MockModel.register_variable(var)
    MockModel.register_variable(state_variable)
    
    model_instance = MockModel(initial_state=initial_array)
    return model_instance

def test_model_instance(
    model: CanRegisterVariable,
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
