from numba import (
    types,
    typed,
)
import pytest

from clearwater_modules.tsm import EnergyBudget
from clearwater_modules.tsm.constants import (
    Meteorological,
    Temperature,
)


@pytest.fixture(scope='function')
def initial_tsm_state() -> dict[str, float]:
    """Return initial state values for the model."""
    return {
        'water_temp_c': 20.0,
        'surface_area': 1.0,
        'volume': 1.0,
    }


@pytest.fixture(scope='function')
def default_meteo_params() -> Meteorological: 
    """Returns default meteorological static variable values for the model.

    NOTE: As of now (11/17/2023) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float values.
    """
    return Meteorological(
        air_temp_c=20,
        q_solar=400,
        sed_temp_c=5.0,
        eair_mb=1.0,
        pressure_mb=1013.0,
        cloudiness=0.1,
        wind_speed=3.0,
        wind_a=0.3,
        wind_b=1.5,
        wind_c=1.0,
        wind_kh_kw=1.0,
    )


@pytest.fixture(scope='function')
def default_temp_params() -> Temperature: 
    """Returns default temperature static variable values for the model.

    NOTE: As of now (11/17/2023) these match the built in defaults, but are 
    copied here to allow for easy modification of the defaults in the future.

    Returns a typed dictionary, with string keys and float or bool values.
    """
    return Temperature(
        stefan_boltzmann=5.67e-8,
        cp_air=1005,
        emissivity_water=0.97,
        gravity=-9.806,
        a0=6984.505294,
        a1=-188.903931,
        a2=2.133357675,
        a3=-1.288580973E-2,
        a4=4.393587233E-5,
        a5=-8.023923082E-8,
        a6=6.136820929E-11,
        pb=1600.0,
        cps=1673.0,
        h2=0.1,
        alphas=0.0432,
        richardson_option=True,
    )

def get_energy_budget_instance(
    initial_tsm_state,
    default_meteo_params,
    default_temp_params,
) -> EnergyBudget:
    """Return an instance of the TSM class."""
    return EnergyBudget(
        initial_state_values=initial_tsm_state,
        meteo_parameters=default_meteo_params,
        temp_parameters=default_temp_params,
        time_dim='tsm_time_step',
    )


@pytest.fixture(scope='module')
def tolerance() -> float:
    """Controls the precision of the pytest.approx() function."""
    return 100 # TODO: flip back to 0.0001 after we finish debugging

def test_defaults(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary


    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # Run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 15.3389 


def test_changed_water_temp_c(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_tsm_state
    initial_state_dict['water_temp_c'] = 40.0 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # Run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == -12.1871 


def test_changed_surface_area(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_tsm_state
    initial_state_dict['surface_area'] = 2.0 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_state_dict,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # Run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 10.6779 



def test_changed_volume(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    initial_state_dict = initial_tsm_state
    initial_state_dict['volume'] = 2.0

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_state_dict,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # Run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 17.66947 


def test_changes_air_temp_c(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['air_temp_c'] = 30.0

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # Run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.48758 


def test_changed_sed_temp_c(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['sed_temp_c'] = 10.0

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # Run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 18.10903 


def test_changed_q_solar(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """Test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['q_solar'] = 450.0

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # Run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 16.37379 


def test_changed_wind_kh_kw(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['wind_kh_kw'] = 0.5 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 15.33893


def test_changed_eair_mb(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['eair_mb'] = 2.0 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 15.48259


def test_changed_pressure_mb(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['pressure_mb'] = 970 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 15.16238149

def test_changed_cloudiness(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['cloudiness'] = 0.0 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 15.32707


def test_changed_wind_a(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['wind_a'] = 1.0e-7 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 15.4728 


def test_changed_wind_b(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['wind_b'] = 1.0e-6 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 16.3432


def test_changed_wind_c(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """test the model with default parameters."""
    # alter parameters as necessary
    default_meteo_params['wind_c'] = 0.5 

    # instantiate the model
    tsm: EnergyBudget = get_energy_budget_instance(
        initial_tsm_state=initial_tsm_state,
        default_meteo_params=default_meteo_params,
        default_temp_params=default_temp_params,
    )

    # run the model
    tsm.increment_timestep()
    water_temp_c = tsm.dataset.isel(tsm_time_step=-1).water_temp_c.values.item()
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 16.6123


def test_use_sed_temp(
    initial_tsm_state, 
    default_meteo_params,
    default_temp_params,
    tolerance,
) -> None:
    """test the model with default parameters."""
    assert True == True
    #TODO: implement this test
    ...

