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
from ._globals import Globals
from ._temp_correction import TempCorrection


class Algae:

    def __init__(self, TwaterC: float, depth: float, PAR: float, fdp: float, TIP: float):
        '''
        Compute algae kinetics (Main function)

        The growth rate of phytplankton algae is limited by
        (1) Light (FL)
        (2) Nitrogen (FN)
        (3) Phosphorous (FP)

        Each is computed individually and the applied to the maximum growth rate to obtain the local growth rate

        Parameters
        ----------
        TwaterC : float
            Water temperature in degrees Celsius
        depth: float
            Depth
        PAR: float
            PAR
        fdp: float
            fdp
        TIP: float
            TIP

        Returns
        ----------
        dApdt: float
            Change in algae concentration

        '''

        gv = Globals()

        # *** Note: depth was a global variable ***

        # Initialize algae
        # Algae stoichiometric ratio
        AWd = 100.0
        AWc = 40.0
        AWn = 7.2
        AWp = 1.0
        AWa = 1000.0
        rna = 7.2 / 1000.0
        rpa = 1.0 / 1000.0
        rca = 40.0 / 1000.0
        rda = 100.0 / 1000.0

        # Parameters related to algae growth
        mu_max = TempCorrection(1.0, 1.047)
        KL = 10.0
        growth_rate_option = 1
        light_limitation_option = 1

        # Parameters related to algae death
        kdp = TempCorrection(0.15, 1.047)

        # Parameters related to algae respiration
        krp = TempCorrection(0.2, 1.047)
        vsap = 0.15

        if gv.initial_values['use_NH4'] or gv.initial_values['use_NO3']:
            KsN = 0.04
            PN = 0.5

        if gv.initial_values['use_TIP']:
            KsP = 0.0012

        if gv.initial_values['use_POC'] or gv.initial_values['use_DOC']:
            Fpocp = 0.9

        # Temperature correction
        mu_max_tc = mu_max.arrhenius_correction(TwaterC)
        krp_tc = krp.arrhenius_correction(TwaterC)
        kdp_tc = kdp.arrhenius_correction(TwaterC)

        sqrt1 = 0.0
        sqrt2 = 0.0

        Ap = gv.globals['Ap']

        # Depth averaged light function
        KEXT = (gv.globals['lambda'] * gv.globals['depth'])

        # (1) Algal light limitation (FL)

        if (Ap <= 0.0 or KEXT <= 0.0 or PAR <= 0.0):
            # After sunset or if there is no algae present
            FL = 0.0
        elif gv.globals['light_limitation_option'] == 1:
            # Half-saturation formulation
            FL = (1.0 / KEXT) * math.log((KL + PAR) / (KL + PAR * math.exp(-KEXT)))
        elif light_limitation_option == 2:
            # Smith's model
            if abs(KL) < 1.0E-10:
                FL = 1.0
            else:
                sqrt1 = (1.0 + (PAR / KL)**2.0)**0.5
                sqrt2 = (1.0 + (PAR * math.exp(-KEXT) / KL)**2.0)**0.5
                FL = (1.0 / KEXT) * math.log((PAR / KL + sqrt1) / (PAR * math.exp(-KEXT) / KL + sqrt2))
        elif light_limitation_option == 3:
            # Steele's model
            if abs(KL) < 1.0E-10:
                FL = 0.0
            else:
                FL = (2.718/KEXT) * (math.exp(-PAR/KL * math.exp(-KEXT)) - math.exp(-PAR/KL))

        # Limit factor to between 0.0 and 1.0.
        # This should never happen, but it would be a mess if it did.
        if FL > 1.0:
            FL = 1.0
        if FL < 0.0:
            FL = 0.0
        # (2) Algal nitrogen limitation (FN)
        # KsN = Michaelis-Menton half-saturation constant (mg N/L) relating inorganic N to algal growth
        if gv.globals['use_NH4'] or gv.globals['use_NO3']:
            NH4 = gv.globals['NH4']
            NO3 = gv.globals['NO3']
            FN = (NH4 + NO3) / (KsN + NH4 + NO3)
            if math.isnan(FN):
                FN = 0.0
            if FN > 1.0:
                FN = 1.0
        else:
            FN = 1.0
        # (3) Algal phosphorous limitation (FP)
        # PO4 = Dissolved (inorganic) phosphorous (mg-P/L)
        # KsP = Michaelis-Menton half-saturation constant (mg-P/L) relating inorganic P to algal growth
        if gv.globals['use_TIP']:
            FP = fdp * TIP / (KsP + fdp * TIP)
            if math.isnan(FP):
                FP = 0.0
            if FP > 1.0:
                FP = 1.0
        else:
            FP = 1.0

        # Algal growth rate with three options
        # (a) Multiplicative (b) Limiting nutrient (c) Harmonic Mean
        if growth_rate_option == 1:
            # (a) Multiplicative (day-1)
            mu = mu_max_tc * FL * FP * FN
        elif growth_rate_option == 2:
            # (b) Limiting nutrient (day-1)
            mu = mu_max_tc * FL * min(FP, FN)
        elif growth_rate_option == 3:
            # (c) Harmonic Mean Option (day-1)
            if FN == 0.0 or FP == 0.0:
                mu = 0.0
            else:
                mu = mu_max_tc * FL * 2.0 / (1.0 / FN + 1.0 / FP)

        # Algal growth
        ApGrowth = mu * Ap

        # Algal respiration
        ApRespiration = krp_tc * Ap

        # Algal mortality
        ApDeath = kdp_tc * Ap

        # Algal settling
        ApSettling = vsap / depth * Ap

        # Algal Biomass Concentration
        # dA/dt = A*(AlgalGrowthRate - AlgalRespirationRate - AlgalDeathRate - AlgalSettlingRate)(mg/L/day)
        dApdt = ApGrowth - ApRespiration - ApDeath - ApSettling

        self.dApdt = dApdt
