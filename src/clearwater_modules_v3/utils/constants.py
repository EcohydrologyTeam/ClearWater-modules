"""v3 physical constants.

Phase 1: re-export the four module-level constants from v2.
"""

from clearwater_modules_v2.utils.constants import (
    AIR_SPECIFIC_HEAT,
    EMISSIVITY_WATER,
    GRAVITY,
    STEFAN_BOLTZMANN,
)

__all__ = [
    "AIR_SPECIFIC_HEAT",
    "EMISSIVITY_WATER",
    "GRAVITY",
    "STEFAN_BOLTZMANN",
]
