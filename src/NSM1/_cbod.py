'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Carbonaceous Biochemical Oxygen Demand (CBOD)
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

Initial Version: June 13, 2021
'''

import math
from ._temp_correction import TempCorrection
from ._globals import Globals


class CBOD:

    def __init__(self, CBOD: float, TwaterC: float, DOX: float):
        '''
        Compute CBOD kinetics

        Parameters
        ----------
        CBOD : float
            Current CBOD concentration
        TwaterC : float
            Water temperature in degrees Celsius
        DOX: float
            Current dissolved oxygen concentration?

        Returns
        ----------
        dCBODdt: float
            Change in CBOD concentration

        '''

        gv = Globals()

        # CBOD parameters
        # Oxidation rate (1/day), Range {0.02-3.4}
        kbod = TempCorrection(0.0, 0.0)  # TODO Need to define

        # Sedimentation rate (1/day), Range {-0.36-0.36}
        ksbod = TempCorrection(0.0, 0.0)  # TODO Need to define

        # Half-saturation oxygen attenuation constant for CBOD oxidation (mg-O/L)
        KsOxbod = 0.0

        # CBOD pathway
        # CBOD oxidation into DIC (mg-O2/L/day)
        CBOD_Oxidation = 0.0
        CBOD_Sediment = 0.0  # CBOD lost to sediment
        CBOD_Oxidation_index = 0
        CBOD_Sediment_index = 0

        # Initialize parameters and pathways
        kbod = TempCorrection(0.12, 1.047)
        ksbod = TempCorrection(0.0, 1.024)

        if gv.globals['use_DOX']:
            KsOxbod = 0.5
        CBOD_Oxidation_index = 0
        CBOD_Sediment_index = 0

        # Note: I am going to code this to compute a single CBOD value
        # to be consistent with the other functions in this module. Then
        # this function can be called multiple times for multiple CBOD values
        kbod_tc = kbod.arrhenius_correction(TwaterC)
        ksbod_tc = ksbod.arrhenius_correction(TwaterC)

        # CBOD Kinetics: compute kinetic rate
        # dBOD/dt = - CBOD decay due to oxidation - CBOD Sedimentation

        if gv.globals['use_DOX']:
            CBOD_Oxidation = DOX / (DOX + KsOxbod) * kbod_tc * CBOD
        if math.isnan(CBOD_Oxidation):
            CBOD_Oxidation = 0.0
        else:
            CBOD_Oxidation = kbod_tc * CBOD

        CBOD_Sediment = ksbod_tc * CBOD
        dCBODdt = - CBOD_Oxidation - CBOD_Sediment

        # TODO: Add code to handle output pathways

        self.dCBODdt = dCBODdt
