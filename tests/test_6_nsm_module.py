import numpy as np
import pytest

from clearwater_modules.nsm1.model import (
    NutrientBudget
)

from clearwater_modules.nsm1.constants import (
    DEFAULT_ALGAE,
    DEFAULT_ALKALINITY,
    DEFAULT_BALGAE,
    DEFAULT_NITROGEN,
    DEFAULT_CARBON,
    DEFAULT_CBOD,
    DEFAULT_DOX,
    DEFAULT_N2,
    DEFAULT_POM,
    DEFAULT_PATHOGEN,
    DEFAULT_PHOSPHORUS,
    DEFAULT_GLOBALPARAMETERS,
    DEFAULT_GLOBALVARS

)


@pytest.fixture(scope='module')
def initial_nsm1_state(initial_array) -> dict[str, float]:
    """Return initial state values for the model."""
    return {
        'Ap': initial_array,
        'Ab': initial_array,
        'NH4': initial_array,
        'NO3': initial_array,
        'OrgN': initial_array,
        'N2': initial_array, 
        'TIP': initial_array,
        'OrgP': initial_array,
        'POC': initial_array,
        'DOC': initial_array,
        'DIC': initial_array,
        'POM': initial_array, 
        'CBOD': initial_array,
        'DOX': initial_array,
        'PX': initial_array,
        'Alk': initial_array,      
    }


@pytest.fixture(scope='module')
def time_steps() -> int:
    return 1


@pytest.fixture(scope='module')
def nsm1_state_variable_names(initial_state_values) -> list[str]:
    """Return the names of the state variables."""
    return list(initial_state_values.keys())


@pytest.fixture(scope='module')
def nutrient_budget_instance(
    time_steps,
    initial_nsm1_state
) -> NutrientBudget:
    """Return an instance of the NSM1 class."""
    return NutrientBudget(
        time_steps=time_steps,
        initial_state_values=initial_nsm1_state,
        updateable_static_variables=['a0'],
        time_dim='tsm_time_step',
    )


def test_nsm1_specific_attributes(nutrient_budget_instance) -> None:
    """Checks that all NSM1  variables are present."""
    assert nutrient_budget_instance.time_dim == 'tsm_time_step'
    assert isinstance(nutrient_budget_instance.algae_parameters, dict)
    assert isinstance(nutrient_budget_instance.alkalinity_parameters, dict)
    assert isinstance(nutrient_budget_instance.Balgae_parameters, dict)
    assert isinstance(nutrient_budget_instance.nitrogen_parameters, dict)
    assert isinstance(nutrient_budget_instance.carbon_parameters, dict)
    assert isinstance(nutrient_budget_instance.CBOD_parameters, dict)
    assert isinstance(nutrient_budget_instance.DOX_parameters, dict)
    assert isinstance(nutrient_budget_instance.N2_parameters, dict)
    assert isinstance(nutrient_budget_instance.POM_parameters, dict)
    assert isinstance(nutrient_budget_instance.pathogen_parameters, dict)
    assert isinstance(nutrient_budget_instance.phosphorus_parameters, dict)
    assert isinstance(nutrient_budget_instance.gp_parameters, dict)
    assert isinstance(nutrient_budget_instance.gvars_parameters, dict)

    assert nutrient_budget_instance.updateable_static_variables == ['a0']

    assert nutrient_budget_instance.algae_parameters == DEFAULT_ALGAE
    assert nutrient_budget_instance.alkalinity_parameters == DEFAULT_ALKALINITY
    assert nutrient_budget_instance.Balgae_parameters == DEFAULT_BALGAE
    assert nutrient_budget_instance.nitrogen_parameters == DEFAULT_NITROGEN
    assert nutrient_budget_instance.carbon_parameters == DEFAULT_CARBON
    assert nutrient_budget_instance.CBOD_parameters == DEFAULT_CBOD
    assert nutrient_budget_instance.DOX_parameters == DEFAULT_DOX
    assert nutrient_budget_instance.N2_parameters == DEFAULT_N2
    assert nutrient_budget_instance.POM_parameters == DEFAULT_POM
    assert nutrient_budget_instance.pathogen_parameters == DEFAULT_PATHOGEN
    assert nutrient_budget_instance.phosphorus_parameters == DEFAULT_PHOSPHORUS
    assert nutrient_budget_instance.gp_parameters == DEFAULT_GLOBALPARAMETERS
    assert nutrient_budget_instance.gvars_parameters == DEFAULT_GLOBALVARS

def test_nsm1_variable_sorting(nutrient_budget_instance) -> None:
    """Checks that we can auto-sort our NSM1 variable"""
    assert isinstance(nutrient_budget_instance.computation_order, list)
    assert nutrient_budget_instance.computation_order[-1].name == 'water_temp_c' #TODO what is first?


def test_nsm1_timestep(nutrient_budget_instance) -> None:
    """Checks that we can auto-sort our NSM1 variable"""
    nutrient_budget_instance.increment_timestep()
    assert len(nutrient_budget_instance.dataset.tsm_time_step) == 2
    assert nutrient_budget_instance.dataset.sel(tsm_time_step=1).isnull().any() == False


#TODO do we need this for NSM?
def test_use_sed_temp(
    initial_tsm_state,
    time_steps,
) -> None:
    """Tests that when we set use_sed_temp to False we get a False boolean array."""
    if 'a0' in initial_tsm_state:
        del initial_tsm_state['a0']
    no_sed_temp = NutrientBudget(
        time_steps=time_steps,
        initial_state_values=initial_tsm_state,
        use_sed_temp=False,
    )
    assert no_sed_temp.dataset.use_sed_temp.dtype == bool
    assert bool(np.all(no_sed_temp.dataset.use_sed_temp == False))
