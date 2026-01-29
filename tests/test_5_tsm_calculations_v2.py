




@pytest.fixture(scope='function')
def default_temperature_dict() -> dict[str, float]:
    """Return default values for the pytest model."""
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
        'water_temperature_compare': 19.9999461,
        }


@pytest.fixture(scope='function')
def parametrize_changes_dict() -> dict[str, float]:
    """Return default values for the pytest model."""
    return {
        'test_defaults': {
            'water_temperature': 20.0,
            'water_temperature_compare': 19.9999461,
        },

        'test_changed_water_temp': {
            'water_temperature': 40.0,
            'water_temperature_compare': 39.99939598,
        },

    }