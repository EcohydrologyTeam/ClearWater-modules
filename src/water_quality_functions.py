'''
=======================================================================================
ClearWater Modules: Water quality equations
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

Initial Version: April 10, 2021
Last Revision Date: April 11, 2021
'''


def ArrheniusCorrection(TwaterC: float, rc20: float, theta: float):
    '''
    Computes an adjusted kinetics reaction rate coefficient for the specified water 
    temperature using the van't Hoff form of the Arrhenius equation

    Parameters
    ----------
    TwaterC : float
        Water temperature in degrees Celsius
    rc20 : float
        Kinetics reaction (decay) coefficient at 20 degrees Celsius
    theta : float
        Temperature correction factor

    Returns
    ----------
    float
        Adjusted kinetics rate for the specified water temperature
    '''
    return rc20 * theta**(TwaterC - 20.0)


class TempCorrection:
    '''
    Temperature correction class
    '''

    def __init__(self, rc20: float, theta: float):
        self.rc20 = rc20
        self.theta = theta

    def arrhenius_correction(self, TwaterC: float):
        return ArrheniusCorrection(TwaterC, self.rc20, self.theta)
