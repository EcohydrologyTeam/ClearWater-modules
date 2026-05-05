"""v3 NSM1 sediment oxygen demand utility.

Stateless SOD primitive ported from v1 NSM1's shared library
(``clearwater_modules/shared/processes.py``). Returns the
temperature-corrected sediment oxygen demand. The optional Monod limitation
on dissolved oxygen present in the v1 ``SOD_tc`` is NOT applied here; v3
NSM1 keeps that limitation as a separate step inside the DOX Process so
this primitive remains a pure Arrhenius correction.
"""

import xarray as xr

from clearwater_modules_v3.utils.conversions import arrhenius_correction


def SOD_tc(
    SOD_20: xr.DataArray,
    SOD_theta: xr.DataArray,
    T_water_C: xr.DataArray,
) -> xr.DataArray:
    """Temperature-corrected sediment oxygen demand.

    Args:
        SOD_20 | g-O2/m^2/d | sediment oxygen demand at 20 deg C.
        SOD_theta | dimensionless | Arrhenius temperature-correction factor.
        T_water_C | deg C | water temperature.

    Returns:
        g-O2/m^2/d | temperature-corrected sediment oxygen demand.
    """
    return arrhenius_correction(T_water_C, SOD_20, SOD_theta)
