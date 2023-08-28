"""
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
"""

import math
from ._temp_correction import TempCorrection
from collections import OrderedDict


class CBOD:

    def __init__(self, CBOD: float, TwaterC: float, DOX: float, use_DOX: bool, CBOD_constant_changes: dict):
        self.CBOD = CBOD
        self.DOX = DOX
        self.TwaterC = TwaterC
        self.use_DOX = use_DOX
        self.CBOD_constant_changes = CBOD_constant_changes

        self.CBOD_constants = OrderedDict()
        self.CBOD_constants = {
            'kbod': 0.12,
            'ksbod': 0,
            'KsOxbod': 0.5,
        }

        for key in self.CBOD_constant_changes.keys():
            if key in self.CBOD_constants:
                self.CBOD_constants[key] = self.CBOD_constant_changes[key]

    def Calculation(self):
        # CBOD parameters
        # Oxidation rate (1/day), Range {0.02-3.4}
        kbod_tc = TempCorrection(
            self.CBOD_constants['kbod'], 1.047).arrhenius_correction(self.TwaterC)

        # Sedimentation rate (1/day), Range {-0.36-0.36}
        ksbod_tc = TempCorrection(
            self.CBOD_constants['ksbod'], 1.047).arrhenius_correction(self.TwaterC)

        KsOxbod = self.CBOD_constants['KsOxbod']

        # CBOD pathway
        # CBOD oxidation into DIC (mg-O2/L/day)
        CBOD_Oxidation = 0.0
        CBOD_Sediment = 0.0  # CBOD lost to sediment
        CBOD_Oxidation_index = 0
        CBOD_Sediment_index = 0

        # Note: I am going to code this to compute a single CBOD value
        # to be consistent with the other functions in this module. Then
        # this function can be called multiple times for multiple CBOD values

        # CBOD Kinetics: compute kinetic rate
        # dBOD/dt = - CBOD decay due to oxidation - CBOD Sedimentation

        if self.use_DOX:
            CBOD_Oxidation = self.DOX / \
                (self.DOX + KsOxbod) * kbod_tc * self.CBOD
        if math.isnan(CBOD_Oxidation):
            CBOD_Oxidation = 0.0
        else:
            CBOD_Oxidation = kbod_tc * self.CBOD

        CBOD_Sediment = ksbod_tc * self.CBOD
        dCBODdt = - CBOD_Oxidation - CBOD_Sediment

        # TODO: is this all the output pathways needed?
        return dCBODdt
