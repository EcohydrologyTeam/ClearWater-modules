"""
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Temperature Correction
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

Initial Version: June 5, 2021
"""

from numba import njit


@njit(cache=True, fastmath=True)
def TempCorrection(rc20, theta, TwaterC):
    """
    Computes an adjusted kinetics reaction rate coefficient for the specified water 
    temperature using the van't Hoff form of the Arrhenius equation

    Parameters:
        TwaterC (float): Water temperature in degrees Celsius
    """

    return rc20 * theta**(TwaterC - 20.0)
