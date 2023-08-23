"""
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Alkalinity Kinetics
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
"""


class Alkalinity:

    def __init__(self, global_module_choices, global_vars, alkalinity_constant_changes):
        
        Alk - state variable
        Must also know - pull in from globals 
        DO - state variable 
        NO3 - state variable
        NH4 - state variable
        Ap - state variable
        Ab - state variable
        
        


        Alk_AlgalGrowth = (ralkaa * ApUptakeFr_NH4 - ralkan * (1 - ApUptakeFr_NH4)) * mu * Ap * rca 
        #mu is calculated in the algae module
        #ApUptakeFr_NH4 is calculated in the nitrogen module
        #ralkaa -- constant
        #ralkan -- constant
        #Ap -- state variable
        #rca -- constant - regional

        Alk_AlgalRespiration = ralkaa * Ap * krp_tc * rca 
        #krp_tc temperature adjusted respiration rate calculated in algae module
        #ralkaa -- constant
        #Ap -- state variable
        #rca -- constant - regional

        Alk_BenthicGrowth = Fb * (ralkba * AbUptakeFr_NH4 - ralkbn * (1 - AbUptakeFr_NH4)) * mub * Ab * rcb / depth  
        #Fb -- fraction of area available for benthic algae - regional 
        #ralkba -- constant
        #AbUptakeFr_NH4 calculated in nitrogen module
        #ralkbn -- constant
        #mub calculated in benthic algae module
        #Ab -- state variable
        #rcb -- constant - regional
        #depth - from model output - regional

        Alk_BenthicRespiration = Fb * ralkba * Ab * krb_tc * rcb / depth
        #Fb -- fraction of area available for benthic algae - regional 
        #ralkba -- constant
        #krb_tc temperature adjust respiration rate calculated in benthic algae module
        #Ab -- state variable
        #rcb -- constant - regional
        #depth - from model output - regional
       
        Alk_Nitrification = ralkn * NH4_Nitrification
        #ralkn -- constant
        #NH4_Nitrification calculated in the nitrogen module

        Alk_Denitrification = ralkden * NO3_Denit
        #ralkden -- constant
        #NO3_Denit calculated in nitrogen module


        dAlkdt = Alk_AlgalRespiration - Alk_AlgalGrowth + Alk_BenthicRespiration - Alk_BenthicGrowth + Alk_Denitrification - Alk_Nitrification

        
       # RATIOS - separate rca and rcb from other ratios
        ralkaa # ratio translating algal growth into Alk if NH4 is the N source (eq ug-Chla -1)
        ralkan # ratio translating algal growth into Alk if NO3 is the N source (eq ug-Chla -1)
        ralkba # ratio translating benthic algae growth into Alk if NH4 is the N source (eq mg-D -1)
        ralkbn # ratio translating benthic algae growth into Alk if NO3 is the N source (eq mg-D -1)
        ralkn # ratio translating NH4 nitrification into Alk (eq mg-N -1)
        ralkden # ratio translating NO3 denitrification into Alk (eq mg-N -1)
        rca # ratio of carbon to algae 
        rcb # ratio of carbon to benthic algae 

        RATES
        kdnit
        knit

        
        pass