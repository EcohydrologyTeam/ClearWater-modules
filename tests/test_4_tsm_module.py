import pytest
import numpy as np

from clearwater_modules.tsm.model import (
    EnergyBudget
)

from clearwater_modules.tsm.constants import (
    DEFAULT_METEOROLOGICAL,
    DEFAULT_TEMPERATURE,
)


@pytest.fixture(scope='module')
def initial_tsm_state(initial_array) -> dict[str, float]:
    """Return initial state values for the model."""
    return {
        'water_temp_c': initial_array,
        'surface_area': initial_array ** 2,
        'volume': initial_array ** 3,
    }


@pytest.fixture(scope='module')
def tsm_state_variable_names(initial_state_values) -> list[str]:
    """Return the names of the state variables."""
    return list(initial_state_values.keys())


@pytest.fixture(scope='module')
def energy_budget_instance(initial_tsm_state) -> EnergyBudget:
    """Return an instance of the TSM class."""
    return EnergyBudget(
        initial_state_values=initial_tsm_state,
        updateable_static_variables=['a0'],
        time_dim='tsm_time_step',
    )


def test_tsm_specific_attributes(energy_budget_instance) -> None:
    """Checks that all TSM variables are present."""
    assert energy_budget_instance.time_dim == 'tsm_time_step'
    assert isinstance(energy_budget_instance.met_parameters, dict)
    assert isinstance(energy_budget_instance.temp_parameters, dict)
    assert energy_budget_instance.updateable_static_variables == ['a0']
    assert energy_budget_instance.met_parameters == DEFAULT_METEOROLOGICAL
    assert energy_budget_instance.temp_parameters == DEFAULT_TEMPERATURE


def test_tsm_variable_sorting(energy_budget_instance) -> None:
    """Checks that we can auto-sort our TSM variable"""
    assert isinstance(energy_budget_instance.computation_order, list)
    assert energy_budget_instance.computation_order[-1].name == 'water_temp_c'


def test_tsm_timestep(energy_budget_instance) -> None:
    """Checks that we can auto-sort our TSM variable"""
    energy_budget_instance.increment_timestep()
    assert len(energy_budget_instance.dataset.tsm_time_step) == 2


def test_use_sed_temp(initial_tsm_state) -> None:
    """Tests that when we set use_sed_temp to False we get a False boolean array."""
    if 'a0' in initial_tsm_state:
        del initial_tsm_state['a0']
    no_sed_temp = EnergyBudget(
        initial_state_values=initial_tsm_state,
        use_sed_temp=False,
    )
    assert no_sed_temp.dataset.use_sed_temp.dtype == bool
    assert bool(np.all(no_sed_temp.dataset.use_sed_temp == False))
