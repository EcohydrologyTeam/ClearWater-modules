"""Regression tests for v3 TSM sign-convention refactor (audit 2026-05-05).

The audit identified that ``flux_upwelling_longwave`` and
``flux_latent_heat`` were pre-negated inside the methods, while
``flux_atmospheric_longwave`` returned a magnitude and
``flux_sensible`` and ``flux_sediment`` carried sign through the
temperature gradient. ``flux_net`` then summed all six. The audit
recommended the v1 / Fortran convention: every method returns a
magnitude (or signed-by-gradient value) and signs are applied at
composition time::

    q_net = (q_sensible + q_solar + q_sediment + q_atmospheric_LW
             - q_upwelling_LW - q_latent)

These tests pin:

1. The new sign convention (no internal pre-negation).
2. ``flux_net`` numerical equivalence to the pre-refactor form.

Refs:
    design/clearwater_modules_v3_tsm_audit_2026-05-05.md
    section 1, finding F-sign-convention.
"""
from __future__ import annotations

import numpy as np
import pytest

from clearwater_modules_v3.processes.temperature import Temperature
from clearwater_modules_v3.utils import constants, conversions


@pytest.fixture(scope="module")
def temperature_module() -> Temperature:
    return Temperature()  # v1 default wind params (post-C1).


# ---------- Sign convention: magnitudes only ----------


def test_upwelling_longwave_returns_positive_magnitude(
    temperature_module: Temperature,
) -> None:
    """``flux_upwelling_longwave`` must return a positive magnitude.

    Stefan-Boltzmann emission is always nonnegative; this test pins
    that the method returns the magnitude, not the negated form. The
    sign is applied in ``flux_net`` via ``- upwelling``.
    """
    water_temperature = 20.0
    result = temperature_module.flux_upwelling_longwave(water_temperature)
    expected = (
        constants.EMISSIVITY_WATER
        * constants.STEFAN_BOLTZMANN
        * conversions.celsius_to_kelvin(water_temperature) ** 4
    )
    np.testing.assert_allclose(float(result), float(expected), rtol=1e-12)
    assert float(result) > 0.0, "flux_upwelling_longwave must be a positive magnitude"


def test_latent_heat_returns_positive_when_evaporative(
    temperature_module: Temperature,
) -> None:
    """``flux_latent_heat`` must be positive in the evaporative regime.

    With ``e_sat > e_air`` (water warmer / drier-air case), the
    polynomial term ``(e_sat - e_air)`` is positive, all other factors
    are positive, so the magnitude is positive. The sign is applied at
    composition via ``- latent``.
    """
    water_temperature = 25.0  # warm water -> high e_sat
    atmospheric_pressure = 1013.0
    wind_speed = 3.0
    atmospheric_vapor_pressure = 5.0  # dry air -> low e_air
    richardson_function = 1.0

    result = temperature_module.flux_latent_heat(
        atmospheric_pressure=atmospheric_pressure,
        water_temperature=water_temperature,
        wind_speed=wind_speed,
        atmospheric_vapor_pressure=atmospheric_vapor_pressure,
        richardson_function=richardson_function,
    )
    assert float(result) > 0.0, (
        "flux_latent_heat must be positive when e_sat > e_air "
        "(evaporative regime); the negative sign is applied in flux_net"
    )


def test_latent_heat_can_be_negative_for_condensation(
    temperature_module: Temperature,
) -> None:
    """``flux_latent_heat`` must be negative when ``e_sat < e_air``.

    Condensation regime: humid air over cold water. The product
    ``(e_sat - e_air)`` is negative; the magnitude propagates that
    sign. Subtracting a negative in ``flux_net`` correctly adds energy
    to the water column (latent heat of condensation).
    """
    water_temperature = 5.0  # cold water -> low e_sat
    atmospheric_pressure = 1013.0
    wind_speed = 3.0
    atmospheric_vapor_pressure = 50.0  # very humid air -> high e_air
    richardson_function = 1.0

    result = temperature_module.flux_latent_heat(
        atmospheric_pressure=atmospheric_pressure,
        water_temperature=water_temperature,
        wind_speed=wind_speed,
        atmospheric_vapor_pressure=atmospheric_vapor_pressure,
        richardson_function=richardson_function,
    )
    assert float(result) < 0.0, (
        "flux_latent_heat must be negative when e_sat < e_air "
        "(condensation onto a cold surface)"
    )


def test_atmospheric_longwave_returns_positive_magnitude(
    temperature_module: Temperature,
) -> None:
    """``flux_atmospheric_longwave`` must be positive (always incoming)."""
    air_temperature = 15.0
    cloudiness = 0.3
    result = temperature_module.flux_atmospheric_longwave(
        air_temperature, cloudiness
    )
    assert float(result) > 0.0


def test_sensible_signed_by_air_minus_water(
    temperature_module: Temperature,
) -> None:
    """``flux_sensible`` carries sign from ``T_air - T_water``.

    Hot air over cold water: positive (heats water). Cold air over hot
    water: negative (cools water). This is the v1 / Fortran convention
    and is preserved unchanged by the audit refactor.
    """
    common = dict(wind_speed=3.0, richardson_function=1.0)
    hot_air = float(
        temperature_module.flux_sensible(
            water_temperature=5.0, air_temperature=25.0, **common
        )
    )
    cold_air = float(
        temperature_module.flux_sensible(
            water_temperature=25.0, air_temperature=5.0, **common
        )
    )
    assert hot_air > 0.0
    assert cold_air < 0.0
    # Symmetry: equal magnitude in opposite directions (uses
    # density_water at the water-temperature argument, so not exact
    # symmetry, but should be opposite signs and same order of
    # magnitude).
    assert hot_air * cold_air < 0.0


def test_sediment_signed_by_sed_minus_water(
    temperature_module: Temperature,
) -> None:
    """``flux_sediment`` carries sign from ``T_sed - T_water``."""
    hot_sed = float(
        temperature_module.flux_sediment(
            water_temperature=5.0,
            sediment_temperature=20.0,
            sediment_thickness=0.1,
        )
    )
    cold_sed = float(
        temperature_module.flux_sediment(
            water_temperature=20.0,
            sediment_temperature=5.0,
            sediment_thickness=0.1,
        )
    )
    assert hot_sed > 0.0
    assert cold_sed < 0.0


# ---------- flux_net composition: signs applied here ----------


def test_flux_net_uses_subtraction_for_upwelling_and_latent(
    temperature_module: Temperature,
) -> None:
    """``flux_net`` must subtract upwelling and latent magnitudes.

    Refactor invariant. With the magnitudes-only convention, the
    composition is::

        q_net = sensible + solar + sediment + atmospheric
                - upwelling - latent

    Verify that ``flux_net`` matches a manual reconstruction that
    explicitly subtracts the two outgoing-energy magnitudes.
    """
    inputs = dict(
        water_temperature=20.0,
        cloudiness=0.3,
        air_temperature=22.0,
        solar_flux=400.0,
        wind_speed=3.0,
        atmospheric_pressure=1013.0,
        atmospheric_vapor_pressure=15.0,
        sediment_temperature=18.0,
        sediment_thickness=0.1,
    )
    actual = float(temperature_module.flux_net(**inputs))

    # Manual reconstruction. Mirror the internal call sequence in
    # flux_net so we exercise the same math path.
    mixing_ratio = temperature_module.mixing_ratio_air(
        inputs["atmospheric_vapor_pressure"], inputs["atmospheric_pressure"]
    )
    density_air = temperature_module.density_air(
        inputs["atmospheric_pressure"], inputs["air_temperature"], mixing_ratio
    )
    density_air_sat = temperature_module.density_air_sat(
        inputs["water_temperature"], inputs["atmospheric_pressure"]
    )
    _, ri_function = temperature_module.richardson_number(
        inputs["wind_speed"],
        density_air_sat=density_air_sat,
        density_air=density_air,
    )
    sensible = float(
        temperature_module.flux_sensible(
            water_temperature=inputs["water_temperature"],
            air_temperature=inputs["air_temperature"],
            wind_speed=inputs["wind_speed"],
            richardson_function=ri_function,
        )
    )
    latent = float(
        temperature_module.flux_latent_heat(
            atmospheric_pressure=inputs["atmospheric_pressure"],
            water_temperature=inputs["water_temperature"],
            wind_speed=inputs["wind_speed"],
            atmospheric_vapor_pressure=inputs["atmospheric_vapor_pressure"],
            richardson_function=ri_function,
        )
    )
    sediment = float(
        temperature_module.flux_sediment(
            water_temperature=inputs["water_temperature"],
            sediment_temperature=inputs["sediment_temperature"],
            sediment_thickness=inputs["sediment_thickness"],
        )
    )
    atmospheric = float(
        temperature_module.flux_atmospheric_longwave(
            inputs["air_temperature"], inputs["cloudiness"]
        )
    )
    upwelling = float(
        temperature_module.flux_upwelling_longwave(inputs["water_temperature"])
    )
    expected = (
        sensible
        + inputs["solar_flux"]
        + sediment
        + atmospheric
        - upwelling
        - latent
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-12)
