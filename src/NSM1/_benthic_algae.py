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
from _temp_correction import TempCorrection
from collections import OrderedDict


class BenthicAlgae:

    def __init__(self, global_module_choices: OrderedDict, global_vars: OrderedDict, Balgae_constant_changes: OrderedDict ):

        self.Balgae_constant_changes=Balgae_constant_changes
        self.global_module_choices=global_module_choices
        self.global_vars = global_vars
        
        # Initialize genrally constant parameters
        self.Balgae_constant = OrderedDict()
        self.Balgae_constant = {
            'BWd': 100,       
            'BWc': 40,      
            'BWn' : 7.2,      
            'BWp' : 1,      
            'BWa' : 3500,       

            'KLb': 10,       
            'KsNb' : 0.25,       
            'KsPb' : 0.125,    
            'Ksb' : 10,       

            'mub_max' : 0.4,   
            'krb' : 0.2,     
            'kdb': 0.3,   
            
            'b_growth_rate_option' : 1,     
            'b_light_limitation_option' : 1,

            'Fw' : 0.9,
            'Fb' : 0.9 
        }

        for key in self.Balgae_constant_changes.keys() :
            if key in self.Balgae_constant:
                self.Balgae_constant[key] = self.Balgae_constant_changes[key]
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
        'BWd'       Benthic algae dry weight        [unitless]
        'BWc'       Benthic algae carbon            [unitless]
        'Bwn'       Benthic algae nitrogen          [unitless]
        'BWp'       Benthic algae phosphorus        [unitless]
        'BWa'       Benthic algae Chla              [unitless]

        'KLb'       Light limiting constant for benthic algae growth            [W/m^2]
        'KsNb       Half-Saturation N limiting constant for Benthic algae       [mg-N/L]
        'KsPb'      Half-Saturation P limiting constant for Benthic algae       [mg-P/L]
        'Ksb'       Half-Saturation density constant for benthic algae growth   [g-D/m^2]

        'mub_max'   maximum benthic algal growth rate                           [1/d]
        'krb'       Respiration rate                                            [1/d]
        'kdb'       Mortality rate                                              [1/d]   

        (Global) global_variables
        'Ab'        Benthic algae concentration         [g/m^2]
        'NH4'       Ammonium concentration              [mg-N/L]
        'NO3'       Nitrate                             [mg-N/L]
        'TIP'       Total inorganic P                   [mg-P/L]
        'TwaterC'   Water Temperature                   [C]
        'depth'     Water depth                         [depth]

        'lambda'    Light attenuation coefficient       [1/m]
        'fdp'       fraction dissolved P                [unitless]
        'PAR'       Surface light intensity             [W/m^2]
        
        (Global) module_choices T/F
        'use_POC'
        'use_DOC'
        'use_NH4'
        'use_NO3'
        'use_TIP'
        'use_OrgN'
        'use_OrgP'

        'b_growth_rate_option'      Benthic algal growth rate options
        'b_light_limitation_option' Benthic algal light limitation

        ----------
        dApdt: float
            Change in algae concentration

        '''
        # *** Note: depth was a global variable ***
    def Calculations (self) : 
        print("Calculating change in benthic algae concentration")
        
        self.Fw = self.Balgae_constant['Fw']
        self.Fb = self.Balgae_constant['Fb']

        # User-supplied maximum benthic algal growth rate (1/day), Range {1.0 - 2.25}
        mub_max_tc = TempCorrection(self.Balgae_constant['mub_max'], 1.047).arrhenius_correction(self.global_vars['TwaterC']) # TODO Need to define

        # Respiration rate (1/day), Range {0.1 - 0.8} 
        krb_tc = TempCorrection(self.Balgae_constant['krb'], 1.06).arrhenius_correction(self.global_vars['TwaterC']) # TODO Need to define

        # Mortality rate (1/day), Range {0.0 - 0.8} 
        kdb_tc = TempCorrection(self.Balgae_constant['kdb'], 1.047).arrhenius_correction(self.global_vars['TwaterC']) # TODO Need to define

        self.rnb = self.Balgae_constant['BWn']/self.Balgae_constant['BWd']     # Ratio of nitrogen to bottom algal biomass (mg-N/mg-D)
        self.rpb = self.Balgae_constant['BWp']/self.Balgae_constant['BWd']     # Ratio of phosphorus to bottom algal biomass (mg-P/mg-D)
        self.rcb = self.Balgae_constant['BWc']/self.Balgae_constant['BWd']     # Ratio of carbon to bottom algal biomass        (mg-C/mg-D) 
        self.rab = self.Balgae_constant['BWa']/self.Balgae_constant['BWd']     # Ratio of chlorophyll-a to bottom algal biomass (µg-Chl-ab/mg-D) 

        # Note that KENT is defined differently here than it was for the algal equations.
        # The equations are different, this expression is more convenient here.
        KEXT = math.exp(-self.global_vars['lambda'] * self.global_vars['depth'])

        # Benthic algae growth
        # The growth rate of benthic algae is limited by
        # (1b)  Light (FLb)
        # (2b)  Nitrogen (FNb)
        # (3b)  Phosphorous (FPb)
        # (4b)  Bottom Area Density (FSb)

        # Each is computed individually and then applied to the maximum growth rate to obtain the local specific growth rate

        if self.global_vars['Ab'] <= 0.0 or KEXT <= 0.0 or self.global_vars['PAR'] <= 0.0:
        # After sunset, no growth
            FLb = 0.0                                         
        elif self.Balgae_constant['b_light_limitation_option'] == 1:
            # Use half-saturation formulation
            FLb = self.global_vars['PAR'] * KEXT / (self.Balgae_constant['KLb'] + self.global_vars['PAR'] * KEXT)
        elif self.Balgae_constant['b_light_limitation_option'] == 2:
            # Use Smith's equation
            FLb = self.global_vars['PAR'] * KEXT / ( (self.Balgae_constant['KLb']**2.0 + (self.global_vars['PAR'] * KEXT)**2.0)**0.5 )
        elif self.Balgae_constant['b_light_limitation_option'] == 3:
            # Use Steele's equation
            if abs(self.Balgae_constant['KLbr']) < 1.0E-10:
                FLb = 0.0 
        else:
            FLb = self.global_vars['PAR'] * KEXT / self.Balgae_constant['KLb'] * math.exp(1.0 - self.global_vars['PAR'] * KEXT / self.Balgae_constant['KLb'])

        # Limit the benthic light limitation factor to between 0.0 and 1.0
        if FLb > 1.0: FLb = 1.0
        if FLb < 0.0: FLb = 0.0

        # (2b) Benthic Nitrogen Limitation (FNb)
        if self.global_module_choices['use_NH4'] or self.global_module_choices['use_NO3']:
            FNb = (self.global_vars['NH4'] + self.global_vars['NO3']) / (self.Balgae_constant['KsNb'] + self.global_vars['NH4'] + self.global_vars['NO3'])
            if math.isnan(FNb): FNb = 0.0
            if FNb > 1.0: FNb = 1.0
        else:
            FNb = 1.0

        # (3b) Benthic Phosphorous Limitation (FPb)
        if self.global_module_choices['use_TIP']:
            FPb = self.global_vars['fdp'] * self.global_vars['TIP'] / (self.Balgae_constant['KsPb'] + self.global_vars['fdp'] * self.global_vars['TIP'])
            if math.isnan(FPb): FPb = 0.0
            if FPb > 1.0: FPb = 1.0
        else:
            FPb = 1.0

        # (4b) Benthic Density Attenuation (FSb)
        FSb   = 1.0 - (self.global_vars['Ab'] / (self.global_vars['Ab'] + self.Balgae_constant['Ksb']))
        if math.isnan(FSb): FSb = 1.0 
        if FSb > 1.0: FSb = 1.0

        # Benthic Local Specific Growth Rate
        if self.Balgae_constant['b_growth_rate_option'] == 1:
            # (a) Multiplicative (day-1)
            mub = mub_max_tc * FLb * FPb * FNb * FSb
        elif self.Balgae_constant['b_growth_rate_option'] == 2:
            # (b) Limiting nutrient (day-1)
            mub = mub_max_tc * FLb * FSb * min(FPb, FNb)

        # Benthic growth
        AbGrowth = mub * self.global_vars['Ab']

        # Benthic respiration
        AbRespiration = krb_tc * self.global_vars['Ab']

        # Benthic death
        AbDeath = kdb_tc * self.global_vars['Ab']

        # Benthic Algal Biomass Concentration
        # dAb/dt = Ab*(BenthicGrowthRate - BenthicRespirationRate - BenthicDateRate)
        # (g/m2/day)
        dAbdt = AbGrowth - AbRespiration - AbDeath

        # Compute derived variables
        # Chlorophyll-a
        Chlb = self.rab * self.global_vars['Ab']

     
        self.dAbdt = dAbdt

        Balgae_pathways = {
        'rnb' : self.rnb,
        'rpb' : self.rpb,
        'AbGrowth' : AbGrowth,
        'AbDeath' : AbDeath,
        'AbRespiration' : AbRespiration,
        'Fw' : self.Fw,
        'Fb' : self.Fb
        }

        print(dAbdt)
        return dAbdt, Balgae_pathways