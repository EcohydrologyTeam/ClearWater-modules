import copy
from datetime import datetime, timedelta
import pytest
from clearwater_modules_v2.processes.temperature import Temperature
from clearwater_data.variables.registry import VariableRegistry
from clearwater_data.variables.float import FloatVariable


@pytest.fixture(scope='module')
def default_temperature_dict() -> dict[str, float]:
    """Return default values for the pytest model simulation."""
    return {
        'water_temperature': 20.0,
        'wetted_surface_area': 1.0,
        'volume': 1.0,
        'cloudiness': 0.1,
        'air_temperature': 20.0,
        'solar_radiation': 400.0,
        'wind_speed': 3.0,
        'atmospheric_pressure': 1013.0,
        'atmospheric_vapor_pressure': 1.0,
        'sediment_temperature': 5.0,
        'sediment_thickness': 0.1,
        'wind_a': 0.3,
        'wind_b': 1.5,
        'wind_c': 1.0,
        'sediment_density': 1600.0,
        'sediment_specific_heat': 1673.0,
        'air_diffusivity_ratio': 1.0,
        'sediment_diffusivity': 0.0432,
        'time_step': 1.0,
        }


CASES_DICT = {
    "test_defaults": (
        {"water_temperature": 20.0},
        19.9999461,
    ),
    "test_changed_water_temp_c": (
        {"water_temperature": 40.0},
        39.99939598,
    ),
    "def test_changed_surface_area": (
        {"wetted_surface_area": 2.0},
        19.9998921,
    ),
    "test_changed_volume": (
        {"volume": 2.0},
        19.99997303,
    ),
    "test_changed_air_temp_c": (
        {"air_temperature": 30.0},
        19.99999407,
    ),
    "test_changed_sed_temp_c": (
        {"sediment_temperature": 10.0},
        19.99997811,
    ),
    "test_changed_q_solar": (
        {"solar_radiation": 450.0},
        19.99995803,
    ),
    "test_changed_wind_kh_kw": (
        {"air_diffusivity_ratio": 0.5},
        19.99994605,
    ),
    "test_changed_eair_mb": (
        {"atmospheric_vapor_pressure": 2.0},
        19.99994772,
    ),
    "test_changed_pressure_mb": (
        {"atmospheric_pressure": 970.0},
        19.99994401,
    ),
    "test_changed_cloudiness": (
        {"cloudiness": 0.0},
        19.99994592,
    ),
    "test_changed_wind_a": (
        {"wind_a": 1.0e-7},
        19.9999476,
    ),
    "test_changed_wind_b": (
        {"wind_b": 1.0},
        19.99995768,
    ),
    "test_changed_wind_c": (
        {"wind_c": 0.5},
        19.99996079,
    ),
    #"test_use_sed_temp": ( #TODO: need to implement use_sed_temp parameter in Temperature process
    #    {"use_sed_temp": False},
    #    20.0000422364348,
    #),
}

CASES = list(CASES_DICT.values())
CASE_IDS = list(CASES_DICT.keys())


@pytest.mark.parametrize(
    "overrides, expected_temperature",
    CASES, #CASES is a list of tuples (overrides dict, expected temperature)
    ids=CASE_IDS,
)
def test_temperature_process(
    overrides,
    expected_temperature,
    default_temperature_dict,
):
    data = copy.deepcopy(default_temperature_dict)
    data.update(overrides)

    #slice the updated dict into registry variables and process parameters
    items = list(data.items())
    split_index = 11 #note: adjusted to match number of registry variables in the process
    data_registry = dict(items[:split_index])
    data_process = dict(items[split_index:])

    variable_registry = VariableRegistry()
    for name, value in data_registry.items():
        variable_registry.register(name, FloatVariable(value))

    process = Temperature(
        wind_a=data_process["wind_a"],
        wind_b=data_process["wind_b"],
        wind_c=data_process["wind_c"],
        sediment_density=data_process["sediment_density"],
        sediment_specific_heat=data_process["sediment_specific_heat"],
        sediment_diffusivity=data_process["sediment_diffusivity"],
        time_step=timedelta(seconds=data_process["time_step"])
    )

    date_time = datetime(2026, 1, 1, 0, 0, 0)
    process.run(date_time, variable_registry)
    result = variable_registry.get("water_temperature")
    assert pytest.approx(result, rel=1e-7) == expected_temperature