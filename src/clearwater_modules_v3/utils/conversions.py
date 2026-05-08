"""v3 unit and kinetics conversions.

``celsius_to_kelvin`` is overridden locally in v3 to use the SI
convention ``T_K = T_C + 273.15`` rather than v2's ``+273.16``
(triple-point) value. See ``utils/constants.py:KELVIN_OFFSET`` and
audit 2026-05-05 (open question 5) for the rationale.

``arrhenius_correction`` is defined in-place here (previously re-exported
from v2). The function body matches v2's verbatim; severing the
re-export is part of the v2 retirement plan.
"""

from clearwater_data.custom_types import ArrayLike

from clearwater_modules_v3.utils.constants import KELVIN_OFFSET


def arrhenius_correction(
    water_temperature: ArrayLike,
    reaction_kinetics: ArrayLike,
    theta: ArrayLike,
) -> ArrayLike:
    """Adjusted kinetics rate at the specified water temperature.

    Uses the van't Hoff form of the Arrhenius equation:
    ``k(T) = k_20 * theta ** (T - 20)``.

    Parameters
    ----------
    water_temperature : ArrayLike
        Water temperature in degrees Celsius.
    reaction_kinetics : ArrayLike
        Kinetics reaction (decay) coefficient at 20 degrees Celsius.
    theta : ArrayLike
        Temperature correction factor.

    Returns
    -------
    ArrayLike
        Adjusted kinetics rate for the specified water temperature.
    """
    return reaction_kinetics * theta ** (water_temperature - 20.0)


def celsius_to_kelvin(celsius: ArrayLike) -> ArrayLike:
    """Convert Celsius to Kelvin using the SI offset (273.15 K = 0 deg C)."""
    return celsius + KELVIN_OFFSET


def kelvin_to_celsius(kelvin: ArrayLike) -> ArrayLike:
    """Convert Kelvin to Celsius using the SI offset (0 deg C = 273.15 K)."""
    return kelvin - KELVIN_OFFSET


__all__ = ["arrhenius_correction", "celsius_to_kelvin", "kelvin_to_celsius"]
