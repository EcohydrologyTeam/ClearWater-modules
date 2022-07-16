'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Benthic Algae Kinetics
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


class BenthicAlgae:

    def __init__(self, TwaterC: float, depth: float, PAR: float, 
        fdp: float, TIP: float, Ab: float, KLbr: float):
        '''
        Compute benthic algae kinetics (Main function)

        The growth rate of phytplankton algae is limited by
        (1) Light (FL)
        (2) Nitrogen (FN)
        (3) Phosphorous (FP)

        The growth rate of benthic algae is limited by
        (1b)  Light (FLb)
        (2b)  Nitrogen (FNb)
        (3b)  Phosphorous (FPb)
        (4b)  Bottom Area Density (FSb)

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
        # *** Note: depth was a global variable ***

        gv = Globals()

        # User-supplied maximum benthic algal growth rate (1/day), Range {1.0 - 2.25}
        mub_max = TempCorrection(0.0, 0.0) # TODO Need to define

        # Respiration rate (1/day), Range {0.1 - 0.8} 
        krb = TempCorrection(0.0, 0.0) # TODO Need to define

        # Mortality rate (1/day), Range {0.0 - 0.8} 
        kdb = TempCorrection(0.0, 0.0) # TODO Need to define

        rnb = 0.0     # Ratio of nitrogen to bottom algal biomass (mg-N/mg-D)
        rpb = 0.0     # Ratio of phosphorus to bottom algal biomass (mg-P/mg-D)
        rcb = 0.0     # Ratio of carbon to bottom algal biomass        (mg-C/mg-D) 
        rab = 0.0     # Ratio of chlorophyll-a to bottom algal biomass (µg-Chl-ab/mg-D) 

        BWd = 0.0     # Algal biomass stoichiometry                    (unitless)
        BWc = 0.0     # Algal carbon stoichiometry                     (unitless)
        BWn = 0.0     # Algal nitrogen stoichiometry                   (unitless)
        BWp = 0.0     # Algal phosphorus stoichiometry                 (unitless)
        BWa = 0.0     # Algal chlorophyll-a stoichiometry              (unitless)
        KLb = 0.0     # Light constant for benthic algae growth                        (W m-2)     Range {2-20}
        KsNb = 0.0    # Half-saturation N constant forbenthic algae growth             (mg-N/L)    Range {0.01-0.766}
        KsPb = 0.0    # Half-saturation P constant for benthic algae growth            (mg-P/L)    Range {0.001-0.08} 
        KSb = 0.0     # Half-saturation density constant for benthic algae growth      (g-D/m2)    Range {10 - 20} 
        PNb = 0.0     # NH4 preference factor for benthic algae growth                 (unitless)  Range {0.0 - 1.0} 
        Fw = 0.0      # Fraction of benthic algae mortality into water column          (unitless)  Range {0.0-1.0}                
        Fb = 0.0      # Fraction of bottom area available for bottom algae growth      (unitless)  Range {0.0-1.0}                    
        Fpocb = 0.0   # Fraction of benthic algal mortality that is converted to POC   (unitless)  Range {0.0-1.0}

        b_growth_rate_option = 0.0      # Local specific growth rate of (benthic) algae
        b_light_limitation_option = 0.0 # Local light limitation formulation

        # Pathway      
        AbGrowth = 0.0        # Benthic algae growth             (g-Ab/m2/day)
        AbRespiration = 0.0   # Benthic algae respiration        (g-Ab/m2/day)
        AbDeath = 0.0         # Benthic algae death              (g-Ab/m2/day)
        FLb = 0.0             # Benthic algae light limitation   (unitless)
        FNb = 0.0             # Benthic algae N limitation       (unitless)
        FPb = 0.0             # Benthic algae P limitation       (unitless)
        FSb = 0.0             # Benthic algae density limitation (unitless)

        AbGrowth_index = 0
        AbRespiration_index = 0
        AbDeath_index = 0        
        FLb_index = 0
        FNb_index = 0
        FPb_index = 0
        FSb_index = 0

        # Local variables
        mub = 0.0

        # Initialize all input parameters with default values

        # Stoichiometric ratio
        BWd = 100.0
        BWc = 40.0
        BWn = 7.2
        BWp = 1.0
        BWa = 5000.0

        # Derived paramters
        rnb = 7.2 / 100.0
        rpb = 1.0 / 100.0
        rcb = 40.0 / 100.0
        rab = 5000.0 / 100.0

        # Parameters related to benthic algae growth
        mub_max_rc20 = 0.4
        mub_max_theta = 1.047
        KLb = 10.0
        KSb = 10.0

        # Parameters related to benthic algae respiration
        krb_rc20 = 0.2
        krb_theta = 1.06

        # Parameters related to benthic algae death
        kdb_rc20 = 0.3
        kdb_theta = 1.047
        Fb = 0.9    

        if gv.globals['use_OrgN'] or gv.globals['use_OrgP'] or gv.globals['use_POC'] or gv.globals['use_DOC']:
            Fw = 0.9

        if gv.globals['use_NH4'] or gv.globals['use_NO3']:
            KsNb = 0.25

        if gv.globals['use_TIP']:
            KsPb = 0.125

        if gv.globals['use_NH4'] and gv.globals['use_NO3']:
            PNb = 0.5

        if gv.globals['use_POC'] or gv.globals['use_DOC']:
            Fpocb = 0.9

        # Integer parameters
        b_growth_rate_option = 1
        b_light_limitation_option = 1

        # Arrhenius temperature correction
        mub_max_tc = mub_max.arrhenius_correction(TwaterC)
        krb_tc = krb.arrhenius_correction(TwaterC)
        kdb_tc = kdb.arrhenius_correction(TwaterC)

        # Compute kinetic rate

        KEXT = 0.0

        # Note that KENT is defined differently here than it was for the algal equations.
        # The equations are different, this expression is more convenient here.
        KEXT = math.exp(-gv.globals['lambda'].value * depth)

        # Benthic algae growth
        # The growth rate of benthic algae is limited by
        # (1b)  Light (FLb)
        # (2b)  Nitrogen (FNb)
        # (3b)  Phosphorous (FPb)
        # (4b)  Bottom Area Density (FSb)

        # Each is computed individually and then applied to the maximum growth rate to obtain the local specific growth rate

        if Ab <= 0.0 or KEXT <= 0.0 or PAR <= 0.0:
        # After sunset, no growth
            FLb = 0.0                                         
        elif b_light_limitation_option == 1:
            # Use half-saturation formulation
            FLb = PAR * KEXT / (KLb + PAR * KEXT)
        elif b_light_limitation_option == 2:
            # Use Smith's equation
            FLb = PAR * KEXT / ( (KLb**2.0 + (PAR * KEXT)**2.0)**0.5 )
        elif b_light_limitation_option == 3:
            # Use Steele's equation
            if abs(KLbr) < 1.0E-10:
                FLb = 0.0 
        else:
            FLb = PAR * KEXT / KLb * math.exp(1.0 - PAR * KEXT / KLb)

        # Limit the benthic light limitation factor to between 0.0 and 1.0
        if FLb > 1.0: FLb = 1.0
        if FLb < 0.0: FLb = 0.0

        # (2b) Benthic Nitrogen Limitation (FNb)
        if gv.globals['use_NH4'] or gv.globals['use_NO3']:
            NH4 = gv.globals['NH4']
            NO3 = gv.globals['NO3']
            FNb = (NH4 + NO3) / (KsNb + NH4 + NO3)
            if math.isnan(FNb): FNb = 0.0
            if FNb > 1.0: FNb = 1.0
        else:
            FNb = 1.0

        # (3b) Benthic Phosphorous Limitation (FPb)
        if gv.globals['use_TIP']:
            FPb = fdp * TIP / (KsPb + fdp * TIP)
            if math.isnan(FPb): FPb = 0.0
            if FPb > 1.0: FPb = 1.0
        else:
            FPb = 1.0

        # (4b) Benthic Density Attenuation (FSb)
        FSb   = 1.0 - (Ab / (Ab + KSb))
        if math.isnan(FSb): FSb = 1.0 
        if FSb > 1.0: FSb = 1.0

        # Benthic Local Specific Growth Rate
        if b_growth_rate_option == 1:
            # (a) Multiplicative (day-1)
            mub = mub_max_tc * FLb * FPb * FNb * FSb
        elif b_growth_rate_option == 2:
            # (b) Limiting nutrient (day-1)
            mub = mub_max_tc * FLb * FSb * min(FPb, FNb)

        # Benthic growth
        AbGrowth = mub * Ab

        # Benthic respiration
        AbRespiration = krb_tc * Ab

        # Benthic death
        AbDeath = kdb_tc * Ab

        # Benthic Algal Biomass Concentration
        # dAb/dt = Ab*(BenthicGrowthRate - BenthicRespirationRate - BenthicDateRate)
        # (g/m2/day)
        dAbdt = AbGrowth - AbRespiration - AbDeath

        # Compute derived variables
        # Chlorophyll-a
        Chlb = rab * Ab

        # TODO: Add output pathways ***
        self.dAbdt = dAbdt