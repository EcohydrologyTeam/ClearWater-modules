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
    TwaterK: float,
) -> float:

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

    return 1.034 * ka_tc * (N2sat - N2)
