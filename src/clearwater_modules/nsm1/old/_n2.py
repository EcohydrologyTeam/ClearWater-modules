"""
=======================================================================================
Nutrient Simulation Module 1 (NSM1): N2 Kinetics
=======================================================================================

Developed by:
* Dr. Todd E. Steissberg (ERDC-EL)
* Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

This module computes the water quality of a single computational cell. The algorithms 
and structure of this program were adapted from the Fortran 95 version of this module, 
developed by:
* Dr. Billy E. Johnson (ERDC-EL)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

Version 1.0

Initial Version: June 12, 2021
Last Revision Date: June 13, 2021
"""

'''
variables in:

depth
TwaterC
pressure_atm
N2
dN2dt
DOX
TDG
use_N2
use_DOX

ka_tc
O2sat

TwaterK 
'''




import math
from clearwater_modules.shared.processes import arrhenius_correction, celsius_to_kelvin
import numba
@numba.njit
def TwaterK(
    TwaterC: float,
) -> float:

    return celsius_to_kelvin(TwaterK)


@numba.njit
def KHN2_tc(
    TwaterK : float,
) -> float :
    
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
    TwaterK : float,
) -> float :
        
    """Calculate partial pressure water vapor (atm)

    Constant values found in documentation

    Args:
        TwaterK: water temperature kelvin (K)

    """
    return math.exp(11.8571  - (3840.70 / TwaterK) - (216961.0 / (TwaterK**2)))

    return 0.00065 * math.exp(1300.0 * (1.0 / TwaterK - 1 / 298.15))


@numba.njit
# Correct N2 saturation for atmopsheric pressure
def P_wv(
    TwaterK: float,
) -> float:

    return math.exp(11.8571 - (3840.70 / TwaterK) - (216961.0 / (TwaterK**2)))


@numba.njit
# N2 saturation
def N2sat(
    KHN2_tc: float,
    pressure_atm: float,
    P_wv: float
) -> float:
    N2sat = 2.8E+4 * KHN2_tc * 0.79 * (pressure_atm - P_wv)
    if (N2sat < 0.0):  # Trap saturation concentration to ensure never negative
        N2sat = 0.0

    return N2sat


@numba.njit
def dN2dt(
    ka_tc: float,
    N2sat: float,
    N2: float,
) -> float:

    Args:
        ka_tc: Oxygen re-aeration rate (1/d)
        N2sat: N2 at saturation f(Twater and atm pressure) (mg-N/L)
        N2: Nitrogen concentration air (mg-N/L)
    """
        
    return 1.034 * ka_tc * (N2sat - N2)

@numba.njit    
def N2_new(
    N2: float,
    dN2dt : float,
) -> float: 
    
    """Calculate change in N2 air concentration (mg-N/L/d)

    Args:
        N2: Nitrogen concentration air (mg-N/L)
        dN2dt: Change in nitrogen concentration air
    """
        
    return N2 + dN2dt

@numba.njit    
def TDG(
    N2: float,
    N2sat : float,
    DOX: float,
    O2sat: float,
    use_DOX: bool,
) -> float: 
    
    """Calculate total dissolved gas (%)

    Args:
        N2: Nitrogen concentration air (mg-N/L)
        N2sat: N2 at saturation f(Twater and atm pressure) (mg-N/L)
        DOX: Dissolved oxygen concentration (mg-O2/L)
        O2sat: O2 at saturation f(Twater and atm pressure) (mg-O2/L)
        use_DOX: true/false use dissolved oxygen module (true/false)
    """
    if use_DOX :
        TDG = (79.0 * N2 / N2sat) + (21.0 * DOX / O2sat)
    else:
        TDG = N2/N2sat

    return TDG