# Changes: lambda is not calcualted
# Some options such as light calculation method was not changeable in FORTRAN but is now changeable

'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Algae Kinetics
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
'''

import math
from collections import OrderedDict
from _temp_correction import TempCorrection

from numba import types
from numba.typed import Dict


class Algae:
    def __init__(self, Ap, depth, vsap, ApGrowth,ApDeath,ApRespiration):
        self.Ap=Ap
        self.depth=depth
        self.vsap=vsap
        self.ApGrowth=ApGrowth
        self.ApDeath=ApDeath
        self.ApRespiration=ApRespiration

    def Calculations (self):
        '''
        Compute algae kinetics (Main function)

        The growth rate of phytplankton algae is limited by
        (1) Light (FL)
        (2) Nitrogen (FN)
        (3) Phosphorous (FP)

        Each is computed individually and the applied to the maximum growth rate to obtain the local growth rate

        Parameters
        ----------
        (Global) module_choices T/F
        'use_Algae'
        'use_NH4'
        'use_NO3'
        'use_TIP'
        'use_OrgN'

        (Global) global_vars
        'Ap'       Algae concentration                      [ug/Chla/L]
        'NH4'      Ammonium concentration                   [mg-N/L]
        'NO3'      Nitrate concentration                    [mg-N/L]
        'TIP'      Total inorganic phosphorus               [mg-P/L]
        'TwaterC'  Water temperature                        [C]
        'depth'    Depth from water surface                 [m]

        (Global) global_par
        'lambda'   Light attenuation coefficient            [unitless]
        'fdp'      Fraction P dissolved                     [unitless]
        'PAR'      Surface light intensity                  [W/m^2]

        algae_const
        'Awd'      Algal dry                                [mg]
        'Awc'      Algal carbon                             [mg]
        'Awn'      Algal nitrogen                           [mg]
        'Awp'      Algal phosphorus                         [mg]
        'Awa'      Algal chlorophyll                        [mg]

        'KL'      light limiting consant for growth         [W/m^2]
        'KSN'     Half-Sat N limiting constant for growth   [mg-N/L]
        'KSP'     Half-Sat P limiting constant for growth   [mg-P/L]
        'mu_max'  Max algae growth                          [1/d]
        'kdp'     Algae mortality rate                      [1/d]
        'krp'     Algal respiration rate                    [1/d]
        'vsap'    Settling velocity                         [m/d]

        algae_options
        'growth_rate_option' Algal growth rate cal option            [no units]
        'light_limitation_option' Algal light limitation cal option  [no units]

        Returns
        ----------
        dApdt: float
            Change in algae concentration

        '''

        print("Calculating change in algae concentration")

        # Algal settling
        self.ApSettling = self.vsap / self.depth * self.Ap          # [ug-Chla/L/d]

        # Algal Biomass Concentration
        # dA/dt = A*(AlgalGrowthRate - AlgalRespirationRate - AlgalDeathRate - AlgalSettlingRate)(mg/L/day)
        dApdt = self.ApGrowth - self.ApRespiration - self.ApDeath - self.ApSettling     # [ug-Chla/L/d]
        

        print (dApdt)

        return dApdt


        