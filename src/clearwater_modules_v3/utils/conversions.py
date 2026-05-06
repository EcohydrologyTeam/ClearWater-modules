"""v3 unit and kinetics conversions.

``celsius_to_kelvin`` is overridden locally in v3 to use the SI
convention ``T_K = T_C + 273.15`` rather than v2's ``+273.16``
(triple-point) value. See ``utils/constants.py:KELVIN_OFFSET`` and
audit 2026-05-05 (open question 5) for the rationale.
``arrhenius_correction`` is re-exported from v2 unchanged.
"""

from clearwater_data.custom_types import ArrayLike
from clearwater_modules_v2.utils.conversions import arrhenius_correction

from clearwater_modules_v3.utils.constants import KELVIN_OFFSET


def celsius_to_kelvin(celsius: ArrayLike) -> ArrayLike:
    """Convert Celsius to Kelvin using the SI offset (273.15 K = 0 deg C)."""
    return celsius + KELVIN_OFFSET


def kelvin_to_celsius(kelvin: ArrayLike) -> ArrayLike:
    """Convert Kelvin to Celsius using the SI offset (0 deg C = 273.15 K)."""
    return kelvin - KELVIN_OFFSET


__all__ = ["arrhenius_correction", "celsius_to_kelvin", "kelvin_to_celsius"]
