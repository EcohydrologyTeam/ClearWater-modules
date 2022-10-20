'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Nitrogen Kinetics
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


class Nitrogen:

    def __init__(self,global_module_choices: OrderedDict, global_vars: OrderedDict, algae_pathways: OrderedDict, Balgae_pathways: OrderedDict, sedFlux_pathways: OrderedDict, nitrogen_constant_changes: OrderedDict):
        self.nitrogen_constant_changes=nitrogen_constant_changes
        self.global_module_choices=global_module_choices
        self.global_vars = global_vars
        self.algae_pathways = algae_pathways
        self.Balgae_pathways = Balgae_pathways
        self.sedFlux_pathways = sedFlux_pathways

    def Calculations (self):
        #use modGlobal,       only: r, depth, TwaterC, NH4, dNH4dt, NO3, dNO3dt, OrgN, dOrgNdt, DIN, TON, TKN, TN, DOX, Ap,  & use_OrgN, use_NH4, use_NO3, use_DOX, use_Algae, use_BAlgae, use_SedFlux
        #use modGlobalParam,  only: vson
        #use modAlgae,        only: PN, ApGrowth, ApRespiration, ApDeath, rna  
        #use modBenthicAlgae, only: Fb, Fw, PNb, AbGrowth, AbRespiration, AbDeath, rnb
        #use modSedFlux,      only: JNH4, JNO3
        '''
        (Global) module_choices T/F
        'use_Algae'
        'use_BAlgae'
        'use_SedFlux'
        'use_NH4'
        'use_NO3'
        'use_DOC
        'use_DOX'
        'use_OrgN'

        (Global) global_vars
        'Ap'       Algae concentration                              [ug/Chla/L]
        'NH4'      Ammonium concentration                           [mg-N/L]
        'NO3'      Nitrate concentration                            [mg-N/L]
        'TwaterC'  Water temperature                                [C]
        'depth'    Depth from water surface                         [m]
        'DOX'      Dissolved oxygen                                 [mg-O2/L]

        nitrogen_constant_changes
        'vson'     Organic N settling velocity                      [m/d]
        'KNR'      Oxygen inhibitation factor for nitrification     [mg-O2/L]
        'knit'     Nitrification rate ammonia decay NH4-->NO2       [1/d]
        'kon'      Decay rate of OrgN --> NH4                       [1/d]
        'kdnit'    Denitrification rate                             [1/d]
        'rnh4'     Sediment release rate of NH4                     [g-N/m^2 * d]
        'KsOxdn'   Half-saturation oxygen inhibition constant denitrification       [mg-O2/L]
        'PN'        NH4 preference factor (1=full NH4)              [unitless]
        'PNb'      NH4 preference factor (1=full NH4)               [unitless]
        'Fw'       Fraction Benthic algae mortality into water      [unitless]
        'Fb'       Fraciton of bottom arae available                [untiless]

        from Algae
        'rna'           AlgalN:Chla ratio                                [mg-N/ugChla]
        'ApGrowth'      Algral growth rate                               [ug-chla/L/d]
        'ApDeath'       Algal death rate                                 [ug-chla/L/d]
        'ApRespiration' AlgalRespiration rate                            [ug-chla/L/d]



        from Benthic Algae
        'rnb'           Benthic Algal N: Benthic Algal Dry Weight        [mg-N/mg-D]
        'AbGrowth'      Benthic Algal growth rate                        [g/m^2*d]
        'AbDeath'       Benthic Algal death rate                         [g/m^2*d]
        'AbRespiration' Benthic Algal respiration rate                   [g/m^2*d]
        
        from SedFlux
        'JNH4'     Sediment water flux of ammonia                   [g-N/m^2*d]
        'JNO3'     Sediment water flux of NO3                       [g-N/m^2*d]
      

        '''
        print("Calculating change in Nitrogen concentration")

        # local variables                
        NitrificationInhibition =0                       # Nitrification Inhibitation (limits nitrification under low DO conditions)	
        ApUptakeFr_NH4  =0                               # fraction of actual floating algal uptake that is from ammonia pool
        ApUptakeFr_NO3 =0                                # fraction of actual floating algal uptake that is from nitrate pool
        AbUptakeFr_NO3 =0                                # fraction of actual benthic algal uptake that is from nitrate pool
        AbUptakeFr_NH4 =0                                # fraction of actual benthic algal uptake that is from ammonia pool

       #Parameters
        self.nitrogen_constant = OrderedDict()
        self.nitrogen_constant = {
            'knit':0.1,              
            'kon' : 0.1,              
            'kdnit': 0.002,	          
            'rnh4': 0,                
            'vno3': 1,	              
            'KsOxdn': 0.1,              
            'KNR' : 0.6,
            'PN' : 0.5,
            'PNb' : 0.5,
            'Fw' : 0.9,
            'Fb': 0.9
        }

        for key in self.nitrogen_constant_changes.keys() :
            if key in self.nitrogen_constant:
                self.nitrogen_constant[key] = self.nitrogen_constant_changes[key]

        if self.global_module_choices['use_NH4'] :
            knit_tc=TempCorrection(self.nitrogen_constant['knit'], 1.083).arrhenius_correction(self.global_vars['TwaterC'])

        if not self.global_module_choices['use_SedFlux'] :
            rnh4_tc=TempCorrection(self.nitrogen_constant['rnh4'], 1.074).arrhenius_correction(self.global_vars['TwaterC'])
            vno3_tc = TempCorrection(self.nitrogen_constant['vno3'], 1.08).arrhenius_correction(self.global_vars['TwaterC'])

        if self.global_module_choices['use_OrgN'] :
            kon_tc=TempCorrection(self.nitrogen_constant['kon'], 1.074).arrhenius_correction(self.global_vars['TwaterC'])

        if self.global_module_choices['use_NO3'] :
            kdnit_tc = TempCorrection(self.nitrogen_constant['kdnit'], 1.045).arrhenius_correction(self.global_vars['TwaterC'])


        # set value of UptakeFr_NH4/NO3 for special conditions
        if self.global_module_choices['use_NH4'] and not self.global_module_choices['use_NO3'] :
            ApUptakeFr_NH4      = 1.0
            ApUptakeFr_NO3      = 0.0
            AbUptakeFr_NH4      = 1.0
            AbUptakeFr_NO3      = 0.0
        if not self.global_module_choices['use_NH4'] and self.global_module_choices['use_NO3'] :
            ApUptakeFr_NH4      = 0.0
            ApUptakeFr_NO3      = 1.0
            AbUptakeFr_NH4      = 0.0
            AbUptakeFr_NO3      = 1.0
        if not self.global_module_choices['use_NH4'] and not self.global_module_choices['use_NO3'] :
            ApUptakeFr_NH4      = 0.5
            ApUptakeFr_NO3      = 0.5
            AbUptakeFr_NH4      = 0.5
            AbUptakeFr_NO3      = 0.5

    #  Calculating Nitrogen Kinetics
        if self.global_module_choices['use_Algae'] and self.global_module_choices['use_NH4'] and self.global_module_choices['use_NO3'] : 
            ApUptakeFr_NH4 = self.nitrogen_constant['PN'] * self.global_vars['NH4'] / (self.nitrogen_constant['PN'] * self.global_vars['NH4']  + (1.0 - self.nitrogen_constant['PN']) * self.global_vars['NO3'])          # [unitless]
            ApUptakeFr_NO3 = 1 - ApUptakeFr_NH4
    # Check for case when NH4 and NO3 are very small.  If so, force uptake_fractions appropriately. 
            if math.isnan(ApUptakeFr_NH4) :
                ApUptakeFr_NH4 = self.nitrogen_constant['PN']
                ApUptakeFr_NO3 = 1.0 - ApUptakeFr_NH4

    # Check for benthic and recompute if necessary
        if self.global_module_choices['use_BAlgae'] and self.global_module_choices['use_NH4'] and self.global_module_choices['use_NO3'] :
            AbUptakeFr_NH4 = (self.nitrogen_constant['PNb'] * self.global_vars['NH4']) / (self.nitrogen_constant['PNb'] * self.global_vars['NH4']  + (1.0 - self.nitrogen_constant['PNb']) * self.global_vars['NO3'])
            AbUptakeFr_NO3 = 1 - AbUptakeFr_NH4

    # Check if NH4 and NO3 are very small.  If so, force uptake_fractions appropriately. 
            if math.isnan(AbUptakeFr_NH4) :
                AbUptakeFr_NH4 = self.nitrogen_constant['PNb']
                AbUptakeFr_NO3 = 1.0 - AbUptakeFr_NH4

        '''
     Organic Nitrogen                     (mg-N/d*L)
     dOrgN/dt =   Algae_OrgN              (Floating Algae -> OrgN)
                  - OrgN_NH4_Decay        (OrgN -> NH4)
                  - OrgN_Settling         (OrgN -> bed)
                  + Benthic Death         (Benthic Algae -> OrgN)														 
        ''' 

        if self.global_module_choices['use_OrgN'] :
            OrgN_NH4_Decay = kon_tc * self.global_vars['OrgN']
            OrgN_Settling = self.global_vars['vson'] / self.global_vars['depth'] * self.global_vars['OrgN']

            if self.global_module_choices['use_Algae'] :
                ApDeath_OrgN  = self.algae_pathways['rna'] * self.algae_pathways['ApDeath']
            else :
                ApDeath_OrgN  = 0.0

            if self.global_module_choices['use_BAlgae'] :  
                AbDeath_OrgN = self.Balgae_pathways['rnb'] * self.nitrogen_constant['Fw'] * self.nitrogen_constant['Fb'] * self.Balgae_pathways['AbDeath'] / self.global_vars['depth']
            else : 
                AbDeath_OrgN = 0.0

            dOrgNdt = ApDeath_OrgN + AbDeath_OrgN - OrgN_NH4_Decay - OrgN_Settling
        else:
            dOrgNdt = 0 
        '''
        Ammonia Nitrogen (NH4)                 (mg-N/day*L)
        dNH4/dt   =    OrgN_NH4_Decay          (OrgN -> NH4)  
                       - NH4 Oxidation         (NH4 -> NO3)
                       - NH4AlgalUptake        (NH4 -> Floating Algae)
                       + Benthos NH4           (Benthos -> NH4)
                       - Benthic Algae Uptake  (NH4 -> Benthic Algae)														 
        '''
    # Compute nitrification inhibition coefficient used to retard oxidation 
    # rate in case of low dissolved oxygen.  The coefficient should range between zero and one.
     	
    # Modify Nitrogren Inhibition Factor. If the function is disabled, make the function linear at DO concentrations 
    # greater than zero, and shut off nitrogen oxidation completely once DO is depleted
        if self.global_module_choices['use_NH4'] : #TODO this looks different
            if self.global_module_choices['use_DOX'] :
                NitrificationInhibition = 1.0 - math.exp(-self.nitrogen_constant['KNR'] * self.global_vars['DOX'])
            else :
                NitrificationInhibition = 1.0

            NH4_Nitrification = NitrificationInhibition * knit_tc * self.global_vars['NH4']
    
            if self.global_module_choices['use_SedFlux'] :
                NH4fromBed = self.sedFlux_pathways['JNH4'] / self.global_vars['depth']
            else :
                NH4fromBed = rnh4_tc / self.global_vars['depth']
    
            if self.global_module_choices['use_Algae'] :
                NH4_ApRespiration = self.algae_pathways['rna'] * self.algae_pathways['ApRespiration']
                NH4_ApGrowth= ApUptakeFr_NH4 * self.algae_pathways['rna'] * self.algae_pathways['ApGrowth']
            else : 
                NH4_ApRespiration = 0.0
                NH4_ApGrowth      = 0.0

            if self.global_module_choices['use_BAlgae'] : 
                #TODO changed the calculation for respiration from the inital FORTRAN due to conflict with the reference guide
                NH4_AbRespiration = self.Balgae_pathways['rnb'] * self.Balgae_pathways['AbRespiration']
                NH4_AbGrowth      = (AbUptakeFr_NH4 * self.Balgae_pathways['rnb'] *  self.nitrogen_constant['Fb'] * self.Balgae_pathways['AbGrowth']) / self.global_vars['depth'] 
            else :
                NH4_AbRespiration = 0.0
                NH4_AbGrowth      = 0.0

            if not self.global_module_choices['use_OrgN'] :
                OrgN_NH4_Decay = 0.0 

            dNH4dt = OrgN_NH4_Decay - NH4_Nitrification + NH4fromBed + NH4_ApRespiration - NH4_ApGrowth + NH4_AbRespiration - NH4_AbGrowth
        else:
            dNH4dt = 0
 
        '''
        Nitrite Nitrogen  (NO3)                       (mg-N/day*L)
        dNO3/dt  =      NH4 Oxidation                 (NH4 -> NO3) 
                        - NO3 Sediment Denitrification
                        - NO3 Algal Uptake            (NO3-> Floating Algae) 
                        - NO3 Benthic Algal Uptake    (NO3-> Benthic  Algae) 
        '''
        if self.global_module_choices['use_NO3'] :
            if self.global_module_choices['use_DOX'] :
                NO3_Denit = (1.0 - (self.global_vars['DOX'] / (self.global_vars['DOX'] + self.nitrogen_constant['KsOxdn']))) * kdnit_tc * self.global_vars['NO3']
                if math.isnan(NO3_Denit) : 
                    NO3_Denit = kdnit_tc * self.global_vars['NO3']
            else :
                NO3_Denit    = 0.0

            if self.global_module_choices['use_SedFlux'] :
                NO3_BedDenit = self.sedFlux_pathways['JNO3'] / self.global_vars['depth']
            else :
                NO3_BedDenit = vno3_tc * self.global_vars['NO3'] / self.global_vars['depth']
   
            if self.global_module_choices['use_Algae'] :
                NO3_ApGrowth  = ApUptakeFr_NO3 * self.algae_pathways['rna'] * self.algae_pathways['ApGrowth']
            else :
                NO3_ApGrowth  = 0.0
 
            if self.global_module_choices['use_BAlgae'] : 
                NO3_AbGrowth  = (AbUptakeFr_NO3 * self.Balgae_pathways['rnb'] *  self.nitrogen_constant['Fb'] * self.Balgae_pathways['AbGrowth'] )/ self.global_vars['depth']  
            else :
                NO3_AbGrowth  = 0.0

            if not self.global_module_choices['use_NH4'] :
                 NH4_Nitrification = 0.0

            dNO3dt = NH4_Nitrification - NO3_Denit - NO3_BedDenit - NO3_ApGrowth - NO3_AbGrowth 
        else:
            dNO3dt = 0

        # output pathways
        '''
    ApDeath_OrgN
    OrgN_Settling
    OrgN_NH4_Decay
    H4_Nitrification
    NH4_ApRespiration
    NH4_ApGrowth
    NH4fromBed
    NO3_ApGrowth
    NO3_Denit
    NO3_BedDenit
    AbDeath_OrgN
    NH4_AbRespiration
    NH4_AbGrowth
    NO3_AbGrowth
    '''
  
        # compute derived variables   
        DIN = 0.0
        TON = 0.0
        TKN = 0.0
        if self.global_module_choices['use_NH4'] :
            DIN = DIN + self.global_vars['NH4']
            TKN = TKN + self.global_vars['NH4']

        if self.global_module_choices['use_NO3'] :
            DIN = DIN + self.global_vars['NO3']

        if self.global_module_choices['use_OrgN'] :
            TON = TON + self.global_vars['OrgN']

        if self.global_module_choices['use_Algae']:
            TON = TON + self.algae_pathways['rna'] * self.global_vars['Ap']
        
        TKN = TKN + TON
        TN  = DIN + TON

        print("dNH4dt", dNH4dt)
        print("dNO3dt", dNO3dt)
        print("dOrgNdt", dOrgNdt)

        print("DIN", DIN)
        print("TON", TON)
        print("TKN", TKN)
        print("TN", TN)

        return DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt