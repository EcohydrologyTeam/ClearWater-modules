"""v3 physical constants.

Constants are defined in-place. ``KELVIN_OFFSET`` follows the SI
convention (273.15 K = 0 deg C); historical v1/v2 utilities used
273.16 (the triple point of water) for v1-utility-parity. v3 uses the
SI value canonically. See ``utils/conversions.py:celsius_to_kelvin``
for the C-to-K conversion function and
``design/clearwater_modules_v3_tsm_audit_2026-05-05.md`` section 5 for
the rationale.

The other physical constants (GRAVITY, STEFAN_BOLTZMANN,
EMISSIVITY_WATER, AIR_SPECIFIC_HEAT) were originally inherited from v2
by re-export; the v2-retirement refactor moved them in-tree.
"""

GRAVITY: float = -9.806
"""m/s^2; sign-convention is downward-positive in the TSM/NSM kinetics."""

STEFAN_BOLTZMANN: float = 5.67e-8
"""W/m^2/K^4. Note this is the truncated v1-parity value; the SI value
is 5.670374e-8."""

EMISSIVITY_WATER: float = 0.97
"""Dimensionless. Used in the upwelling-longwave term."""

AIR_SPECIFIC_HEAT: float = 1005.0
"""J/(kg*K). Used in the sensible-heat flux term."""

KELVIN_OFFSET: float = 273.15
"""SI definition of the absolute-temperature offset for 0 deg C.

Use ``T_K = T_C + KELVIN_OFFSET`` and ``T_C = T_K - KELVIN_OFFSET``.
The previous v3 form (re-exported from v2) used ``+273.16``, the
triple point of water; that is non-canonical as a 0-degC offset and
introduces a ``+0.01 K`` systematic bias in every Kelvin-evaluated
quantity (Stefan-Boltzmann T^4, Brutsaert e_sat, latent-heat
polynomial inputs, air-density polynomial). The bias is tiny in
absolute terms (~3.4e-5 relative at 293 K) but is the wrong unit
choice for SI temperature physics.

Audit 2026-05-05 (open question 5) resolved to switch v3 to 273.15.
"""

__all__ = [
    "AIR_SPECIFIC_HEAT",
    "EMISSIVITY_WATER",
    "GRAVITY",
    "KELVIN_OFFSET",
    "STEFAN_BOLTZMANN",
]
