import numpy as np
import xarray as xr

def pwv(
    TwaterK: xr.DataArray
) -> xr.DataArray:
    """Calculate partial pressure of water vapor

    Args:
        TwaterK: Water temperature kelvin
    """
    return np.exp(11.8571 - 3840.70 / TwaterK - 216961 / TwaterK ** 2)


#TwaterC??
def DOs_atm_alpha(
    TwaterK: xr.DataArray
) -> xr.DataArray:
    """Calculate DO saturation atmospheric correction coefficient

    Args:
        TwaterK: Water temperature kelvin
    """
    return .000975 - 1.426 * 10 ** -5 * TwaterK + 6.436 * 10 ** -8 * TwaterK ** 2

def DOX_sat(
    TwaterK: xr.DataArray,
    pressure_mb: xr.DataArray,
    pwv: xr.DataArray,
    DOs_atm_alpha: xr.DataArray
) -> xr.DataArray:
    """Calculate DO saturation value

    Args:
        TwaterK: Water temperature kelvin
        pressure_mb: Atmospheric pressure (mb)
        pwv: Patrial pressure of water vapor (atm)
        DOs_atm_alpha: DO saturation atmospheric correction coefficient
    """
    no_exp = (-139.34410 + 1.575701 * 10 ** 5 / TwaterK - 6.642308 * 10 ** 7 / TwaterK ** 2
                                 + 1.243800 * 10 ** 10 / TwaterK - 8.621949 * 10 ** 11 / TwaterK)
    DOX_sat_uncorrected = np.exp(-139.34410 + 1.575701 * 10 ** 5 / TwaterK - 6.642308 * 10 ** 7 / TwaterK ** 2
                                 + 1.243800 * 10 ** 10 / TwaterK - 8.621949 * 10 ** 11 / TwaterK)

    DOX_sat_uncorrected = -139.34410 + ( 1.575701E05 / TwaterK ) - ( 6.642308E07 / (TwaterK**2.0) ) + ( 1.243800E10 / (TwaterK**3.0) ) - ( 8.621949E11 / (TwaterK**4.0) )

    DOX_sat_corrected = DOX_sat_uncorrected * pressure_mb * \
        (1 - pwv / pressure_mb) * (1 - DOs_atm_alpha * pressure_mb) / \
        ((1 - pwv) * (1 - DOs_atm_alpha))
    
    print(no_exp)
    print(fortran_copy)
    print(DOX_sat_uncorrected)
    print(DOX_sat_uncorrected)
    
    return DOX_sat_corrected

pwv_1 = pwv(298.15)
correction = DOs_atm_alpha(298.15)

DOX_sat(298.15,1013.25,pwv_1,correction)