import math
from clearwater_modules.shared.processes import arrhenius_correction, celsius_to_kelvin
import numba
import xarray as xr


@numba.njit
def KHN2_tc(
    TwaterK : xr.DataArray,
) -> xr.DataArray :
    
    """Calculate Henry's law constant (mol/L/atm)
    
    Constant values found on NIST

    Args:
        TwaterK: water temperature kelvin (K)
        Henry's law constant for solubility in water at 298.15K: 0.00065 (mol/(kg*bar))
        Temperature dependence constant: 1300 (K) 
        Reference temperature: 298.15 (K) 
    """

    return 0.00065 * math.exp(1300.0 * (1.0 / TwaterK - 1 / 298.15))   
        
@numba.njit
def P_wv(
    TwaterK : xr.DataArray,
) -> xr.DataArray :
        
    """Calculate partial pressure water vapor (atm)

    Constant values found in documentation

    Args:
        TwaterK: water temperature kelvin (K)

    """
    return math.exp(11.8571  - (3840.70 / TwaterK) - (216961.0 / (TwaterK**2)))

@numba.njit     
#N2 saturation
def N2sat(
    KHN2_tc : xr.DataArray,
    pressure_atm: xr.DataArray,
    P_wv: xr.DataArray
) -> xr.DataArray:
    
    """Calculate N2 at saturation f(Twater and atm pressure) (mg-N/L)

    Args:
        KHN2_tc: Henry's law constant (mol/L/atm)
        pressure_atm: atmosphric pressure in atm (atm)
        P_wv: Partial pressure of water vapor (atm)
    """
        
    N2sat = 2.8E+4 * KHN2_tc * 0.79 * (pressure_atm - P_wv)  
    if (N2sat < 0.0) :  #Trap saturation concentration to ensure never negative
        N2sat = 0.0 

    return N2sat

@numba.njit    
def dN2dt(
    ka_tc : xr.DataArray, #TODO this should be calculated in Carbon based on kah_tc and kaw_tc 
    N2sat : xr.DataArray,
    N2: xr.DataArray,
) -> xr.DataArray: 
    
    """Calculate change in N2 air concentration (mg-N/L/d)

    Args:
        ka_tc: Oxygen re-aeration rate (1/d)
        N2sat: N2 at saturation f(Twater and atm pressure) (mg-N/L)
        N2: Nitrogen concentration air (mg-N/L)
    """
        
    return 1.034 * ka_tc * (N2sat - N2)

@numba.njit    
def N2_new(
    N2: xr.DataArray,
    dN2dt : xr.DataArray,
) -> xr.DataArray: 
    
    """Calculate change in N2 air concentration (mg-N/L/d)

    Args:
        N2: Nitrogen concentration air (mg-N/L)
        dN2dt: Change in nitrogen concentration air
    """
        
    return N2 + dN2dt

@numba.njit    
def TDG(
    N2: xr.DataArray,
    N2sat : xr.DataArray,
    DOX: xr.DataArray,
    DOX_sat: xr.DataArray,
    use_DOX: bool,
) -> xr.DataArray: 
    
    """Calculate total dissolved gas (%)

    Args:
        N2: Nitrogen concentration air (mg-N/L)
        N2sat: N2 at saturation f(Twater and atm pressure) (mg-N/L)
        DOX: Dissolved oxygen concentration (mg-O2/L)
        DOX_sat: O2 at saturation f(Twater and atm pressure) (mg-O2/L)
        use_DOX: true/false use dissolved oxygen module (true/false)
    """
    if use_DOX :
        TDG = (79.0 * N2 / N2sat) + (21.0 * DOX / DOX_sat)
    else:
        TDG = N2/N2sat

    return TDG