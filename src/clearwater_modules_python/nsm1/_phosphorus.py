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
            'kop' : 0.1,
            'rpo4' : 0
        }
        
        for key in self.P_constant_changes.keys() :
            if key in self.P_constants:
                self.P_constants[key] = self.P_constant_changes[key]

    def Calculation(self) :

        '''
        (Global) module_choices T/F
        'use_Algae'
        'use_BAlgae'
        'use_OrgP'
        'use_TIP'
        'use_SedFlux'

        (Global) global_vars
        'Ap'       Algae concentration                              [ug/Chla/L]
        'TwaterC'  Water temperature                                [C]
        'depth'    Depth from water surface                         [m]
        'OrgP'     Organic phosphorus                               [mg-P/L]
        'TIP'      Total inorganic phosphorus                       [mg-P/L]
        'vs'       Sediment settling velocity                       [m/d]
        'fdp'      fraction P dissolved                             [unitless]

        phosphorus_constant_changes
        'kop'       Decay rate or orgnaic P to DIP                  [1/d]
        'rpo4'      Benthic sediment release rate of DIP            [g-P/m2*d]

        from Algae
        'rPa'           AlgalP : Chla ratio                              [mg-P/ugChla]
        'ApGrowth'      Algal growth rate                                [ug-chla/L/d]
        'ApDeath'       Algal death rate                                 [ug-chla/L/d]
        'ApRespiration' AlgalRespiration rate                            [ug-chla/L/d]

        from Benthic Algae
        'rpb'           Benthic Algal P: Benthic Algal Dry Weight        [mg-P/mg-D]
        'AbGrowth'      Benthic Algal growth rate                        [g/m^2*d]
        'AbDeath'       Benthic Algal death rate                         [g/m^2*d]
        'AbRespiration' Benthic Algal respiration rate                   [g/m^2*d]
        
        from SedFlux
        'JDIP'     Sediment water flux of phosphate                   [g-P/m^2*d]      

        '''
        print("Calculating change in phosphorus concentration")

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
            OrgP_Settling   = (self.global_vars['vsop']/ self.global_vars['depth']) * self.global_vars['OrgP']
        
            if self.global_module_choices['use_Algae'] :
                ApDeath_OrgP = self.algae_pathways['rpa'] * self.algae_pathways['ApDeath']
            else :
                ApDeath_OrgP    = 0
        
            if self.global_module_choices['use_BAlgae'] :
                AbDeath_OrgP = (self.Balgae_pathways['rpb'] * self.Balgae_pathways['Fw'] * self.Balgae_pathways['Fb'] * self.Balgae_pathways['AbDeath']) / self.global_vars['depth'] 
            else :
                AbDeath_OrgP = 0.0

            dOrgPdt = ApDeath_OrgP + AbDeath_OrgP - OrgP_DIP_decay - OrgP_Settling
        else: 
            dOrgPdt = 0
        
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
            
            if self.global_module_choices['use_OrgP'] :
                OrgP_DIP_decay = kop_tc * self.global_vars['OrgP']
            else :
                OrgP_DIP_decay = 0.0

            if self.global_module_choices['use_Algae'] :    
                DIP_ApRespiration = self.algae_pathways['rpa'] * self.algae_pathways['ApRespiration']      
                DIP_ApGrowth= self.algae_pathways['rpa'] * self.algae_pathways['ApGrowth']
            else :
                DIP_ApRespiration = 0.0
                DIP_ApGrowth = 0.0

            if self.global_module_choices['use_BAlgae'] : 
                DIP_AbRespiration = self.Balgae_pathways['rpb'] * self.Balgae_pathways['AbRespiration']
                DIP_AbGrowth = self.Balgae_pathways['rpb'] * self.Balgae_pathways['Fb'] * self.Balgae_pathways['AbGrowth'] / self.global_vars['depth']
            else :
                DIP_AbRespiration = 0.0
                DIP_AbGrowth = 0.0

            dTIPdt = OrgP_DIP_decay - TIP_Settling + DIPfromBed + DIP_ApRespiration - DIP_ApGrowth + DIP_AbRespiration - DIP_AbGrowth 
        else :
            dTIPdt = 0

        #Derived variable calculations
        TOP = 0.0
        if self.global_module_choices['use_OrgP'] :
            TOP = TOP + self.global_vars['OrgP'] 
        if self.global_module_choices['use_Algae'] :
            TOP = TOP + self.algae_pathways['rpa'] * self.global_vars['Ap'] 
        
        TP = TOP
        if self.global_module_choices['use_TIP'] :
            TP  = TP + self.global_vars['TIP'] 

            '''
            Residual FORTRAN code I did not know what it ment
            fdp = 1.0
            do i = 1, nGS
                fdp = fdp + kdpo4(i,r) * Solid(i) / 1.0E6
            end do
            '''

            DIP = (self.global_vars['TIP']) * self.global_vars['fdp'] #TODO Check formula, check all with the change of d/dt
        
        print("dOrgPdt", dOrgPdt)
        print("dTIPdt", dTIPdt)
        print("TOP", TOP)
        print("TP", TP)

        return dOrgPdt, dTIPdt, TOP, TP
