"""v3 NSM1 sorption partitioning utility.

Stateless partitioning primitive ported from v1 NSM1's shared library
(``clearwater_modules/shared/processes.py``). Provides the dissolved
fraction of total inorganic phosphorus (TIP) given a partition coefficient
and particulate solids concentration.
"""

import xarray as xr


def fdp(
    use_TIP: xr.DataArray,
    Solid: xr.DataArray,
    kdpo4: xr.DataArray,
) -> xr.DataArray:
    """Dissolved fraction of total inorganic phosphorus.

    Computes the dissolved fraction via the linear-equilibrium sorption
    isotherm ``fdp = 1 / (1 + kdpo4 * Solid / 1e-6)``. When ``use_TIP`` is
    false, returns 0 to suppress all TIP-related kinetics downstream.

    Args:
        use_TIP | bool | total inorganic phosphorus module switch.
        Solid | mg/L | inorganic suspended solids concentration.
        kdpo4 | L/kg | solid-water partition coefficient for orthophosphate.

    Returns:
        dimensionless | dissolved fraction of TIP, in [0, 1].
    """
    return xr.where(use_TIP, 1.0 / (1.0 + kdpo4 * Solid / 0.000001), 0.0)
