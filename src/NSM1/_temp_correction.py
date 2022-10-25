'''
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
'''
from numba import njit, jitclass, types, typed

@jitclass([('rc20', types.float64),('theta', types.float64)])
class TempCorrection:
    '''
    Temperature correction class
    '''

    def __init__(self, rc20: types.float64, theta: types.float64):
        '''
        Initialize temperature correction

        rc20 : float
            Kinetics reaction (decay) coefficient at 20 degrees Celsius
        theta : float
            Temperature correction factor
        '''
        self.rc20 = rc20
        self.theta = theta

    def arrhenius_correction(self, TwaterC: types.float64) -> types.float64:
        '''
        Computes an adjusted kinetics reaction rate coefficient for the specified water 
        temperature using the van't Hoff form of the Arrhenius equation

        Parameters:
            TwaterC (float): Water temperature in degrees Celsius
        '''
        return_value : types.float64 = self.rc20 * self.theta**(TwaterC - 20.0)

        return return_value
