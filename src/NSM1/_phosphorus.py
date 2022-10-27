'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Phosphorus Kinetics
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
from _temp_correction import TempCorrection

class Phosphorus:

    def __init__(self, global_module_choices, global_vars, algae_pathways, Balgae_pathways, sedFlux_pathways, P_constant_changes):
        self.global_module_choices = global_module_choices
        self.global_vars = global_vars
        self.P_constant_changes=P_constant_changes
        self.algae_pathways = algae_pathways
        self.Balgae_pathways = Balgae_pathways
        self.sedFlux_pathways = sedFlux_pathways

        self.P_constants = OrderedDict()
        self.P_constants = {

        }
        
    def Calculation(self) :

    # use modGlobal,       only: r, nGS, depth, TwaterC, Solid, TIP, dTIPdt, OrgP, dOrgPdt, DIP, TOP, TP, Ap,use_Algae, use_BAlgae, use_OrgP, use_TIP, use_SedFlux
    # use modGlobalParam,  only: vs, vsop, fdp, kdpo4
    # use modAlgae,        only: ApGrowth, ApRespiration, ApDeath, rpa
    # use modBenthicAlgae, only: Fb, Fw, AbGrowth, AbRespiration, AbDeath, rpb
    # use modSedFlux,      only: JDIP
    # use modDLL,          only: R8, nRegion, list_class, SetOptionalIndex, TempCorrectionStruct, Arrhenius_TempCorrection  
        kop=0.1
        rop4=0

        if self.global_module_choices['use_OrgP'] :
            kop_tc=TempCorrection(self.P_constants['kop'], 1.047).arrhenius_correction(self.global_vars['TwaterC'])

        if self.global_module_choices['use_TIP'] and not self.global_module_choices['use_SedFlux'] :
            rpo4_tc=TempCorrection(self.P_constants['rpo4'], 1.074).arrhenius_correction(self.global_vars['TwaterC'])

        #Total Organic Phosphorus   
        '''
            Organic Phosphorus  (mgP/day)
        
            dOrgP/dt =    Algae_OrgP              (Algae -> OrgP)
                        - OrgP Decay              (OrgP -> DIP)	
                        - OrgP Settling           (OrgP -> bed) 
                        + BenthicAlgalDeath       (Benthic Algae -> OrgP)    	
            
        '''

        if self.global_module_choices['use_OrgP'] :
            OrgP_DIP_decay  = kop_tc  * self.global_vars['OrgP']
            OrgP_Settling   = self.global_vars['vsop']/ self.global_vars['depth'] * self.global_vars['OrgP']
        
        if self.global_module_choices['use_Algae'] :
            ApDeath_OrgP = self.algae_pathways['rpa'] * self.algae_pathways['ApDeath']
        else :
            ApDeath_OrgP    = 0
      
        if self.global_module_choices['use_BAlgae'] :
            AbDeath_OrgP = self.Balgae_pathways['rpb'] * self.P_constants['Fw'] * self.P_constants['Fb'] * self.Balgae_pathways['AbDeath'] / self.global_vars['depth'] 
        else :
            AbDeath_OrgP = 0.0

        dOrgPdt = ApDeath_OrgP + AbDeath_OrgP - OrgP_DIP_decay - OrgP_Settling
     
        # Total Inorganic Phosphorus  (mg/P/day)
        '''
         dTIP/dt =     OrgP Decay                (OrgP -> DIP)
                     - DIP AlgalUptake           (DIP -> Floating Algae)
                     - DIP BenthicAlgae Uptake   (DIP -> Floating Algae)	
                     - TIP Settling              (TIP -> bed)
                     + DIP From Benthos          (Benthos -> DIP) 
        '''
        if self.global_module_choices['use_TIP'] :
            if self.global_module_choices['use_SedFlux'] :
                DIPfromBed = self.sedFlux_pathways['JDIP'] / self.global_vars['depth']
        else :
            DIPfromBed = rpo4_tc / self.global_vars['depth']

        TIP_Settling = self.global_vars['vs'] / self.global_vars['depth'] * (1.0 - self.global_vars['fdp']) * self.global_vars['TIP']  
        
        if not self.global_module_choices['use_OrgP'] :
            OrgP_DIP_decay = 0.0

        if self.global_module_choices['use_Algae'] :    
            DIP_ApRespiration = self.algae_pathways['rpa'] * self.algae_pathways['ApRespiration']      
            DIP_ApGrowth= self.algae_pathways['rpa'] * self.algae_pathways['ApGrowth']
        else :
            DIP_ApRespiration = 0.0
            DIP_ApGrowth = 0.0

        if self.global_module_choices['use_BAlgae'] : 
            DIP_AbRespiration = self.Balgae_pathways['rpb'] * self.P_constants['Fb'] * self.Balgae_pathways['AbRespiration'] / self.global_vars['depth']
            DIP_AbGrowth = self.Balgae_pathways['rpb'] * self.P_constants['Fb'] * self.Balgae_pathways['AbGrowth'] / self.global_vars['depth']
        else :
            DIP_AbRespiration = 0.0
            DIP_AbGrowth = 0.0

        dTIPdt = OrgP_DIP_decay - TIP_Settling + DIPfromBed + DIP_ApRespiration - DIP_ApGrowth + DIP_AbRespiration - DIP_AbGrowth 

        TOP = 0.0
        if self.global_module_choices['use_OrgP'] :
            TOP = TOP + self.global_vars['OrgP'] 
        if self.global_module_choices['use_Algae'] :
            TOP = TOP + self.algae_pathways['rpa'] * self.global_vars['Ap'] 
        
        TP = TOP
        if self.global_module_choices['use_TIP'] :
            TP  = TP + self.global_vars['TIP'] 

            '''
            fdp = 1.0
            do i = 1, nGS
                fdp = fdp + kdpo4(i,r) * Solid(i) / 1.0E6
            end do
            '''
            DIP = (self.global_vars['TIP']) * self.global_vars['fdp'] #TODO Check formula, check all with the change of d/dt
