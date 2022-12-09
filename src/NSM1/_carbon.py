'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Carbon Cycling
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
from collections import OrderedDict

class Carbon:
    def __init__(self, global_module_choices, global_vars, carbon_constant_changes):
        self.global_vars = global_vars
        self.global_module_choices = global_module_choices

        self.Carbon_constants = OrderedDict()
        self.Carbon_constants = {
            'F_pocp' : .9
            'vsoc' : .01
            'kpoc' : .005
            'kdoc' : .01
            'KsOxmc' : 1
            'Fco2' : .2
            'Pco2' : 383 
        }
    def Calculations():
        """
        Compute carbon kinetics (Main function)

        Parameters
        ----------
        (Global) module_choices
        'use_Algae'
        'use_BAlgae'
        'use_POC'
        'use_DOC'
        'use_DIC'
        'use_DOX'
        'use_SedFlux'

        (Global) global_vars
        'AP'    Algae concentration             [ug/Chla/L]
        'Ab'    Benthic algae concentration
        'POC'   Particulate organic carbon concentration    
        'DOC'
        'DOX'
        'CBOD'




        """

        #State variables
        POC #non-living particulate detrital carbon
        DOC #dissolved organic carbon
        DIC #dissolved inorganic carbon, does not include solid-phase calcium carbonate

        #POC#################################################################
        Algal_mortality = F_pocp * kdp * rca * Ap
        #F_pocp -- constant regional
        #kdp -- algae_constant
        #rca -- global_constant 
        #Ap -- state variable

        POC_hydrolysis = kpoc * POC
        #kpoc -- carbon_constant
        #POC -- state variable

        POC_settling = vsoc * POC / depth
        #vsoc -- carbon_constant
        #POC -- state variable
        #depth -- from hydraulic model 

        Benthic_algae_mortality = F_pocb * kdp * rcb * Ab * Fw * Fb / depth
        #F_pocb -- carbon_constant regional
        #kdp -- algae_constant
        #rcb -- global_constant
        #Ab -- state variable
        #Fw -- constant regional
        #Fb -- constant regional
        #depth -- from hydraulic model
        #######################################################################


        #DOC###################################################################
        Algal_mortality = (1 - F_pocp) * kdp * rca * Ap
        #F_pocp -- carbon_constant regional
        #kdp -- algae_constant
        #rca -- global constant
        #Ap -- state variable
        POC_hydrolysis = kpoc * POC
        #kpoc -- carbon_constant
        #POC -- state variable
        DOC_oxidation = DOX / (DOX + KsOxmc) * kdoc * DOC 
        #DOX -- state variable
        #KsOxmc -- carbon_constant
        #kdoc -- CBOD_constant
        #DOC -- state variable
        Benthic_algae_mortality = (1 - F_pocb) * rcb * Fb * Fw * Ab * kdb / depth 
        #F_pocb -- carbon_constant regional
        #rcb -- global constant
        #Fb -- benthic_algae_constant regional
        #Fw -- benthic_algae_constant regional
        #Ab -- state variable
        #kdp -- benthic_algae_constant
        #depth -- from hydraulics

        # NSM II -- POM_dissolution = kpom * fcom * POM 
        # NSM II -- DOC_denitrification = 
        #######################################################################

    
        #DIC####################################################################
        Atmospheric_CO2_Reaeration = 12 * kac * (KH * pco2 10**-3 - Fco2 * DIC) 
        #kac -- carbon_constant
        #KH -- carbon_constant
        #pco2 -- carbon_constant -- met input... 
        #Fco2 -- carbon_constant
        #DIC -- state variable
        Algal_respiration = krp * rca * Ap 
        #krp -- algae_constant
        #rca -- global_constant
        #Ap -- state variable
        Algal_photosynthesis = mu * rca * Ap 
        #mu -- algae_pathway
        #rca -- global constant
        #Ap -- state variable
        Benthis_algae_respiration = krb * rcb * Ab * Fb / depth 
        #krb -- benthic_algae_constant
        #rcb -- global_constant
        #Ab -- state variable
        #Fb -- benthic_algae_constant
        #depth -- from hydraulics
        Benthic_algae_photosynthesis = mub * rcb * Ab * Fb / depth 
        #mub -- benthic_algae_pathway
        #rcb -- global_constant
        #Ab -- state variable
        #Fb -- benthic_algae_constant
        #depth -- from hydraulics
        CBOD_oxidation = 0
        for i in CBOD:
            CBOD_oxidation = CBOD_oxidation + CBOD_pathways[CBOD_Ox] 
        CBOD_oxidation = CBOD_oxidation / roc
        #CBOD_Ox -- CBOD_pathway
        #roc -- global_constant 
        Sediment_release = SOD_tc / roc / depth
        #SOD_tc -- global_constant
        #roc -- global_constant
        #depth -- from hydraulics   
        # NSM II -- DOC_mineralization = 
        ##########################################################################


        pass
