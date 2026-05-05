"""v3 port of ``tests/test_5_tsm_calculations.py``.

This is a faithful port of the v1 TSM calculation regression tests against
the v3-native ``Temperature`` class. The expected post-substep water
temperatures are reused verbatim from the v1 file (rebaselined 2026-05-01
under the corrected latent-heat-of-vaporization formula). v3's
``latent_heat_vaporization`` already uses the corrected (Celsius) formula,
so it produces the same Lv as v1.

v3-specific notes
-----------------

1. **API shape.** v1 exercises module-level math via the ``EnergyBudget``
   driver (``increment_timestep``); v3 exercises the same math through
   ``Temperature.temperature_change(...)``. The post-substep water
   temperature is reconstructed as
   ``water_temperature + delta_temperature``.

2. **Guard convention.** v3 adds a depth ramp + per-hour rate cap on the
   per-substep delta T. Both are disabled here
   (``q_net_depth_ramp_ref=0.0``, ``dTdt_max_per_hour=float('inf')``) so
   the comparison reduces to v1's pre-guard math. The guards are tested
   separately in the stability-ramp test.

3. **Sediment parameters.** v1's ``Temperature`` constants used
   ``cps=1673.0`` and ``alphas=0.0432``. v3's class defaults are
   ``sediment_specific_heat=1000.0`` and ``sediment_diffusivity=0.0061``.
   The fixture overrides v3 defaults to match v1 so the expected values
   line up.

4. **Time step.** v1 ran with ``dt=1/86400`` (1 second). v3's
   ``time_step`` is a ``timedelta``; we use ``timedelta(seconds=1)``.
"""

from datetime import timedelta

import pytest

from clearwater_modules_v3.processes.temperature import Temperature


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='function')
def initial_state() -> dict[str, float]:
    """Return initial water-column state values (mirrors v1 ``initial_tsm_state``)."""
    return {
        'water_temperature': 20.0,
        'surface_area': 1.0,
        'volume': 1.0,
    }


@pytest.fixture(scope='function')
def meteo_inputs() -> dict[str, float]:
    """Return default meteorological inputs (mirrors v1 ``default_meteo_params``).

    v3 names: ``solar_radiation`` (was ``q_solar``),
    ``atmospheric_vapor_pressure`` (was ``eair_mb``), ``atmospheric_pressure``
    (was ``pressure_mb``), ``air_temperature`` (was ``air_temp_c``),
    ``sediment_temperature`` (was ``sed_temp_c``), ``sediment_thickness``
    (was ``h2``).
    """
    return {
        'air_temperature': 20.0,
        'solar_flux': 400.0,
        'sediment_temperature': 5.0,
        'atmospheric_vapor_pressure': 1.0,
        'atmospheric_pressure': 1013.0,
        'cloudiness': 0.1,
        'wind_speed': 3.0,
        'sediment_thickness': 0.1,
    }


@pytest.fixture(scope='function')
def temperature_kwargs() -> dict:
    """Return ``Temperature`` constructor kwargs.

    Defaults reproduce v1's ``Meteorological``/``Temperature`` constants:
        wind_a=0.3, wind_b=1.5, wind_c=1.0
        sediment_density=1600.0, sediment_specific_heat=1673.0
        sediment_diffusivity=0.0432
        air_diffusivity_ratio=1.0  (v1 ``wind_kh_kw``)

    Guards disabled so the math matches v1's pre-guard form.
    """
    return {
        'wind_a': 0.3,
        'wind_b': 1.5,
        'wind_c': 1.0,
        'sediment_density': 1600.0,
        'sediment_specific_heat': 1673.0,
        'sediment_diffusivity': 0.0432,
        'air_diffusivity_ratio': 1.0,
        'time_step': timedelta(seconds=1),
        'use_sediment_temperature': True,
        'q_net_depth_ramp_ref': 0.0,
        'dTdt_max_per_hour': float('inf'),
    }


@pytest.fixture(scope='module')
def tolerance() -> float:
    """Controls the precision of the pytest.approx() function (matches v1)."""
    return 0.0000001


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _run_one_substep(
    temperature_kwargs: dict,
    initial_state: dict,
    meteo_inputs: dict,
) -> float:
    """Build a ``Temperature`` and compute one substep's updated water temp.

    Equivalent to v1's ``tsm.increment_timestep()`` followed by reading
    ``water_temp_c`` for the latest step.
    """
    temp = Temperature(**temperature_kwargs)
    delta_t = temp.temperature_change(
        water_temperature=initial_state['water_temperature'],
        surface_area=initial_state['surface_area'],
        volume=initial_state['volume'],
        cloudiness=meteo_inputs['cloudiness'],
        air_temperature=meteo_inputs['air_temperature'],
        solar_flux=meteo_inputs['solar_flux'],
        wind_speed=meteo_inputs['wind_speed'],
        sediment_temperature=meteo_inputs['sediment_temperature'],
        sediment_thickness=meteo_inputs['sediment_thickness'],
        atmospheric_pressure=meteo_inputs['atmospheric_pressure'],
        atmospheric_vapor_pressure=meteo_inputs['atmospheric_vapor_pressure'],
    )
    # delta_t may be a 0-d numpy array; coerce to float.
    return float(initial_state['water_temperature'] + float(delta_t))


# ---------------------------------------------------------------------------
# Tests (one per v1 test, same names where reasonable)
# ---------------------------------------------------------------------------

def test_defaults(temperature_kwargs, initial_state, meteo_inputs, tolerance) -> None:
    """Default parameters: water temp after one 1-second substep."""
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999932590607706  # rebaselined 2026-05-01


def test_changed_water_temp_c(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Increase initial water temperature."""
    initial_state['water_temperature'] = 40.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 39.99926924100962  # rebaselined 2026-05-01


def test_changed_surface_area(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Increase wetted surface area."""
    initial_state['surface_area'] = 2.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999865181215412  # rebaselined 2026-05-01


def test_changed_volume(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Increase volume."""
    initial_state['volume'] = 2.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999966295303853  # rebaselined 2026-05-01


def test_changed_air_temp_c(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Increase air temperature."""
    meteo_inputs['air_temperature'] = 30.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999989667533338  # rebaselined 2026-05-01


def test_changed_sed_temp_c(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Increase sediment temperature."""
    meteo_inputs['sediment_temperature'] = 10.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999964651929712  # rebaselined 2026-05-01


def test_changed_q_solar(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Increase solar radiation."""
    meteo_inputs['solar_flux'] = 450.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.99994456808903  # rebaselined 2026-05-01


def test_changed_wind_kh_kw(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Decrease the diffusivity ratio (v1 ``wind_kh_kw`` -> v3 ``air_diffusivity_ratio``)."""
    temperature_kwargs['air_diffusivity_ratio'] = 0.5
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999932590607706  # rebaselined 2026-05-01


def test_changed_eair_mb(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Increase atmospheric vapor pressure (v1 ``eair_mb``)."""
    meteo_inputs['atmospheric_vapor_pressure'] = 2.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999935347790892  # rebaselined 2026-05-01


def test_changed_pressure_mb(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Decrease atmospheric pressure (v1 ``pressure_mb``)."""
    meteo_inputs['atmospheric_pressure'] = 970
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.99992980763556  # rebaselined 2026-05-01


def test_changed_cloudiness(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Zero cloudiness."""
    meteo_inputs['cloudiness'] = 0.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999932453268244  # rebaselined 2026-05-01


def test_changed_wind_a(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Decrease wind-function coefficient ``a``."""
    temperature_kwargs['wind_a'] = 1.0e-7
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.99993575672672  # rebaselined 2026-05-01


def test_changed_wind_b(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Decrease wind-function coefficient ``b``."""
    temperature_kwargs['wind_b'] = 1.0
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999948421208064  # rebaselined 2026-05-01


def test_changed_wind_c(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Decrease wind-function exponent ``c``."""
    temperature_kwargs['wind_c'] = 0.5
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 19.999952663004645  # rebaselined 2026-05-01


def test_use_sed_temp(
    temperature_kwargs, initial_state, meteo_inputs, tolerance
) -> None:
    """Disable sediment temperature coupling."""
    temperature_kwargs['use_sediment_temperature'] = False
    water_temp_c = _run_one_substep(temperature_kwargs, initial_state, meteo_inputs)
    assert isinstance(water_temp_c, float)
    assert pytest.approx(water_temp_c, tolerance) == 20.000028774573728  # rebaselined 2026-05-01
