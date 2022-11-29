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

from _temp_correction_numba import TempCorrection
from numba import njit
from numba.types import float64, unicode_type
from numba.typed import Dict


@njit(cache = True, fastmath = True)
def Calculations (global_module_choices, global_vars, algae_constant_changes):
        
    global_module_choices=global_module_choices
    global_vars = global_vars
    algae_constant=algae_constant_changes

    '''
    # Initialize genrally constant parameters
    algae_constant=Dict.empty(key_type=unicode_type, value_type=float64)
    algae_constant['AWd'] = 100
    algae_constant['AWc'] = 40
    algae_constant['AWn'] = 7.2
    algae_constant['AWp'] = 1
    algae_constant['AWa'] = 1000

    algae_constant['KL'] = 10
    algae_constant['KsN'] = 0.04
    algae_constant['KsP'] = 0.0012
    algae_constant['mu_max'] = 1
    algae_constant['kdp'] = 0.15
    algae_constant['krp'] = 0.2
    algae_constant['vsap'] = 0.15
    algae_constant['growth_rate_option'] = 1
    algae_constant['light_limitation_option'] = 1
  


    for key in algae_constant_changes.keys() :
        if key in algae_constant:
            algae_constant[key] = algae_constant_changes[key]
    '''
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

    #print("Calculating change in algae concentration")
    
    #Updating generally constant parameters

    rna : float64 = algae_constant['AWn'] / algae_constant['AWa']             # Algal N : Chla ratio [mg-N/ugChla]
    rpa : float64 = algae_constant['AWp'] / algae_constant['AWa']             # Algal P : Chla ratio [mg-P/ugChla]
    rca : float64 = algae_constant['AWc'] / algae_constant['AWa']             # Algal C : Chla ratio [mg-C/ugChla]
    rda : float64 = algae_constant['AWd'] / algae_constant['AWa']             # Algal D : Chla ratio [mg-D/ugChla]

    # Parameters related to algae growth and settling
    KL = algae_constant['KL']                                     # Light limiting constant for algal growth [W/m^2]
    growth_rate_option = algae_constant['growth_rate_option']                    
    light_limitation_option = algae_constant['light_limitation_option']           

    # Temperature correction
    mu_max_tc = TempCorrection(algae_constant['mu_max'], 1.047,global_vars['TwaterC'])           # Maximum algal growth rate [1/d]
    krp_tc = TempCorrection(algae_constant['krp'], 1.047,global_vars['TwaterC'])                 # algae respiration rate [1/d]
    kdp_tc = TempCorrection(algae_constant['kdp'], 1.047,global_vars['TwaterC'])                 # algae mortality rate [1/d]

    sqrt1 = 0.0
    sqrt2 = 0.0

    Ap = global_vars['Ap']           # algae [ug-Chla/L]
    PAR = global_vars['PAR']
    # Depth averaged light function
    KEXT = (global_vars['lambda'] * global_vars['depth'])    # lambda is light attenuation coefficient (1/m). depth is depth from water surface (m) [unitless] TODO: other depth not initalized method

    # (1) Algal light limitation (FL)

    if (Ap <= 0.0 or KEXT <= 0.0 or PAR <= 0.0):
        # After sunset or if there is no algae present
        FL = 0.0                                                                        # light limiting factor for algal growth [unitless]
    elif light_limitation_option == 1:
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
    if global_module_choices['use_NH4'] or global_module_choices['use_NO3']:
        NH4 = global_vars['NH4']                         # Ammonium [mg-N/L]
        NO3 = global_vars['NO3']                         # Nitrate [mg-N/L]
        FN = (NH4 + NO3) / (algae_constant['KsN'] + NH4 + NO3)                  # [unitless]
        if math.isnan(FN):
            FN = 0.0
        if FN > 1.0:
            FN = 1.0
    else:
        FN = 1.0
    # (3) Algal phosphorous limitation (FP)
    # PO4 = Dissolved (inorganic) phosphorous (mg-P/L)
    # KsP = Michaelis-Menton half-saturation constant (mg-P/L) relating inorganic P to algal growth
    fdp = global_vars['fdp']
    TIP = global_vars['TIP']

    if global_module_choices['use_TIP']:
        FP = fdp * TIP / (algae_constant['KsP'] + fdp * TIP)              # [unitless]
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
        mu = mu_max_tc * FL * FP * FN                   # [1/d]
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
    ApGrowth = mu * Ap                      # [ug-Chla/L/d]

    # Algal respiration
    ApRespiration = krp_tc * Ap             # [ug-Chla/L/d]

    # Algal mortality
    ApDeath = kdp_tc * Ap                   # [ug-Chla/L/d]

    # Algal settling
    ApSettling = algae_constant['vsap'] / global_vars['depth'] * Ap          # [ug-Chla/L/d]

    # Algal Biomass Concentration
    # dA/dt = A*(AlgalGrowthRate - AlgalRespirationRate - AlgalDeathRate - AlgalSettlingRate)(mg/L/day)
    dApdt = ApGrowth - ApRespiration - ApDeath - ApSettling     # [ug-Chla/L/d]
    
    algae_pathway = Dict.empty(key_type=unicode_type, value_type=float64)
    algae_pathway = {
        'ApGrowth': ApGrowth,
        'ApRespiration' : ApRespiration,
        'ApDeath': ApDeath,
        'rna': rna,
        'rda' : rda,         
        'rca': rca,
        'rpa': rpa
    }

    print ("print dApdt", dApdt)

    return dApdt, algae_pathway


