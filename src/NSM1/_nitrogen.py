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

from cmath import exp
import math
from collections import OrderedDict

from pandas import isna
from _temp_correction import TempCorrection


class Nitrogen:

    def __init__(self):
        pass
    def Calculations (self):
        # use modGlobal,       only: r, depth, TwaterC, NH4, dNH4dt, NO3, dNO3dt, OrgN, dOrgNdt, DIN, TON, TKN, TN, DOX, Ap,  & use_OrgN, use_NH4, use_NO3, use_DOX, use_Algae, use_BAlgae, use_SedFlux
        #use modGlobalParam,  only: vson
        #use modAlgae,        only: PN, ApGrowth, ApRespiration, ApDeath, rna  
        #use modBenthicAlgae, only: Fb, Fw, PNb, AbGrowth, AbRespiration, AbDeath, rnb
        #use modSedFlux,      only: JNH4, JNO3
        #use modDLL,          only: R8, nRegion, list_class, SetOptionalIndex, TempCorrectionStruct, Arrhenius_TempCorrection 
 
        TwaterC = 0
        PN =0   #NH4 preference factor for algal growth
        PNb =0 #No3 preference factor for benthic algal growth
        NH4 = 0 #Concentration NH4
        NO3 = 0 #Concentration NO3
        OrgN = 0 #Concentration Orgnaic Nitrogen
        vson = 0 # Organic N settlign velocity (m/d)
        depth = 0 #depth
        rna = 0 #algal N: Chla ratio (mg-N/ug-Chla)
        ApDeath = 0 # algal mortality rate (1/d)
        rnb = 0 # benthic algae N: D ration (mg-N/mg-D)
        Fw = 0 # Fraction of benthic algae mortality into the water column (0-1)
        Fb = 0 #fraction of bottom area available for benthic algae growth (0-1)
        AbDeath = 0 # benthic algal mortality rate (1/d)
        KNR = 0 # oxygen inhibitation factor for nitrification (0.6-0.7) [mg-O2/L]
        DOX = 0 # dissolved oxygen
        JNH4 = 0 # sediment release rate for NH4 (g-N/m^2 * d)
        ApRespiration = 0 # Ap Resipration rate
        ApGrowth = 0 # Algae growth rate
        AbRespiration = 0 # Ab resipration rate
        AbGrowth = 0 # Benthic algae growth rate
        JNO3 = 0 # sediment release rate for NO3 (g-N/<^2 * d)
        Ap = 0 

        #pathways
        ApDeath_OrgN  =0                # Floating Algae -> OrgN     (mg-N/L/day)
        OrgN_Settling =0                # OrgN -> bed (settling)     (mg-N/L/day) 
        OrgN_NH4_Decay  =0              # OrgN -> NH4 Hydrolysis     (mg-N/L/day) 
        NH4_Nitrification =0            # NH4 -> NO3  Nitrification  (mg-N/L/day)
        NH4_ApRespiration =0            # Floating algae -> NH4      (mg-N/L/day)   
        NH4_ApGrowth =0                 # NH4 -> Floating algae      (mg-N/L/day) 
        NH4fromBed =0                   # bed ->  NH4 (diffusion)    (mg-N/L/day)
        NO3_ApGrowth =0                 # NO3 -> Floating algae      (mg-N/L/day)
        NO3_Denit =0                    # NO3 -> Loss                (mg-N/L/day) 
        NO3_BedDenit =0                 # Sediment denitrification   (mg-N/L/day) 
        AbDeath_OrgN =0                 # Benthic Algae -> OrgN      (mg-N/L/day)
        NH4_AbRespiration =0            # Benthic algae -> NH4       (mg-N/L/day) 
        NH4_AbGrowth =0                 # NH4 -> Benthic Algae       (g-N/L/day)
        NO3_AbGrowth =0                 # NO3 -> Benthic Algae       (g-N/L/day)

        ApDeath_OrgN_index = 0
        OrgN_Settling_index = 0
        OrgN_NH4_Decay_index = 0
        NH4_Nitrification_index = 0
        NH4_ApRespiration_index = 0
        NH4_ApGrowth_index = 0
        NH4fromBed_index = 0
        NO3_ApGrowth_index = 0
        NO3_Denit_index = 0
        NO3_BedDenit_index = 0
        AbDeath_OrgN_index = 0
        NH4_AbRespiration_index = 0
        NH4_AbGrowth_index = 0
        NO3_AbGrowth_index = 0

        # local variables        
        KNR = 0.6                                        # QCLL Hardwire (mg-O2/L)  Range {0.6-0.7}            
        NitrificationInhibition =0                       # Nitrification Inhibitation (limits nitrification under low DO conditions)	
        ApUptakeFr_NH4  =0                               # fraction of actual floating algal uptake that is from ammonia pool
        ApUptakeFr_NO3 =0                                # fraction of actual floating algal uptake that is from nitrate pool
        AbUptakeFr_NO3 =0                                # fraction of actual benthic algal uptake that is from nitrate pool
        AbUptakeFr_NH4 =0                                # fraction of actual benthic algal uptake that is from ammonia pool

        initial_values: dict = {
            'use_Algae': True,
            'use_BAlgae': False,
            'use_OrgN': True,
            'use_NH4': True,
            'use_NO3': True,
            'use_TIP': True,
            'use_OrgP': True,
            'use_POC': False,
            'use_DOC': False,
            'use_DIC': False,
            'use_DOX': True,
            'use_N2': False,
            'use_Pathogen': False,
            'use_Alk': False,
            'use_POM2': False,
            'use_SedFlux': False
        }

       #Parameters
        knit=0.1              # ammonia decay  NH4 -> NO3 (1/day)  Range {0.1-1.0}
        knit_tc=0           # TODO temp correction
        kon=0.1              # OrgN -> NH4 (1/day)  Range {0.02-0.4}
        kon_tc=0            # TODO temp correction
        kdnit=0.002	            # denitrification rate (1/day)
        kdnit_tc=0          # TODO temp correction
        rnh4=0              # Benthos source rate  Benthos -> NH4 (mg/m2day)  Range {variable}	
        rnh4_tc =0          # TODO temp correction
        vno3=0	            # Sediment denitrification reaction velocity (m/d)	
        vno3_tc	=0          # TODO temp correction    
        KsOxdn=0            # Half-saturation oxygen attenuation constant for denitrification	  (mg-O2/L)

        if initial_values['use_NH4'] :
            knit_tc=TempCorrection(knit, 1.083).arrhenius_correction(TwaterC)

        if not initial_values['use_SedFlux'] :
            rnh4_tc=TempCorrection(rnh4, 1.074).arrhenius_correction(TwaterC)
            vno3_tc = TempCorrection(vno3, 1.08).arrhenius_correction(TwaterC)

        if initial_values['use_OrgN'] :
            kon_tc=TempCorrection(kon, 1.074).arrhenius_correction(TwaterC)

        if initial_values['use_NO3'] :
            kdnit_tc = TempCorrection(kdnit, 1.045).arrhenius_correction(TwaterC)

        if initial_values['use_DOX'] :
            KsOxdn=0.1


        # set value of UptakeFr_NH4/NO3 for special conditions
        if initial_values['use_NH4'] and not initial_values['use_NO3'] :
            ApUptakeFr_NH4      = 1.0
            ApUptakeFr_NO3      = 0.0
            AbUptakeFr_NH4      = 1.0
            AbUptakeFr_NO3      = 0.0
        if not initial_values['use_NH4'] and initial_values['use_NO3'] :
            ApUptakeFr_NH4      = 0.0
            ApUptakeFr_NO3      = 1.0
            AbUptakeFr_NH4      = 0.0
            AbUptakeFr_NO3      = 1.0
        if not initial_values['use_NH4'] and not initial_values['use_NO3'] :
            ApUptakeFr_NH4      = 0.5
            ApUptakeFr_NO3      = 0.5
            AbUptakeFr_NH4      = 0.5
            AbUptakeFr_NO3      = 0.5

        '''
  ! set a real parameter
  logical function SetNitrogenRealParameter(name, paramValue)
    character(len=*), intent(in) :: name
    real(R8),         intent(in) :: paramValue(nRegion) 
    !

    SetNitrogenRealParameter = .true.                   
    select case (name)
      case ('knit_rc20')
        knit%rc20 = paramValue
      case ('knit_theta')
        knit%theta = paramValue
      case ('kon_rc20')
        kon%rc20 = paramValue
      case ('kon_theta')
        kon%theta = paramValue
      case ('rnh4_rc20')
        rnh4%rc20 = paramValue
      case ('rnh4_theta')
        rnh4%theta = paramValue
      case ('kdnit_rc20')
        kdnit%rc20 = paramValue
      case ('kdnit_theta')  
        kdnit%theta = paramValue
      case ('vno3_rc20')
        vno3%rc20 = paramValue
      case ('vno3_theta')
        vno3%theta = paramValue
      case ('KsOxdn')
        KsOxdn = paramValue
      case default
        ! did not find the parameter, return false
        SetNitrogenRealParameter = .false.
    end select
  end function
  !

  !===========================================================================================================================
  ! set pathway indexes
  subroutine SetPathwayIndexes(list)
    type(list_class), intent(inout) :: list
    logical success 
    !

    success = SetOptionalIndex(list, 'ApDeath_OrgN',      ApDeath_OrgN_index) 
    success = SetOptionalIndex(list, 'OrgN_Settling',     OrgN_Settling_index) 
    success = SetOptionalIndex(list, 'OrgN_NH4_Decay',    OrgN_NH4_Decay_index) 
    success = SetOptionalIndex(list, 'NH4_Nitrification', NH4_Nitrification_index) 
    success = SetOptionalIndex(list, 'NH4_ApRespiration', NH4_ApRespiration_index) 
    success = SetOptionalIndex(list, 'NH4_ApGrowth',      NH4_ApGrowth_index) 
    success = SetOptionalIndex(list, 'NH4fromBed',        NH4fromBed_index) 
    !
    success = SetOptionalIndex(list, 'NO3_ApGrowth',      NO3_ApGrowth_index) 
    success = SetOptionalIndex(list, 'NO3_Denit',         NO3_Denit_index) 
    success = SetOptionalIndex(list, 'NO3_BedDenit',      NO3_BedDenit_index) 
    success = SetOptionalIndex(list, 'AbDeath_OrgN',      AbDeath_OrgN_index)     
    success = SetOptionalIndex(list, 'NH4_AbRespiration', NH4_AbRespiration_index) 
    success = SetOptionalIndex(list, 'NH4_AbGrowth',      NH4_AbGrowth_index)   
    success = SetOptionalIndex(list, 'NO3_AbGrowth',      NO3_AbGrowth_index)     
  end subroutine 
  !
'''

    #  subroutine NitrogenKinetics()
    # TODO make sure the ApUptakeFr_NO3 is used as well
        if initial_values['use_Algae'] and initial_values['use_NH4'] and initial_values['use_NO3'] : 
            ApUptakeFr_NH4 = PN * NH4 / (PN * NH4  + (1.0 - PN) * NO3)
    # Check for case when NH4 and NO3 are very small.  If so, force uptake_fractions appropriately. 
        if math.isnan(ApUptakeFr_NH4) :
            ApUptakeFr_NH4 = PN
            ApUptakeFr_NO3 = 1.0 - ApUptakeFr_NH4

    # Check for benthic and recompute if necessary
        if initial_values['use_BAlgae'] and initial_values['use_NH4'] and initial_values['use_NO3'] :
            AbUptakeFr_NH4 = (PNb * NH4) / (PNb * NH4  + (1.0 - PNb) * NO3)
    
    # Check if NH4 and NO3 are very small.  If so, force uptake_fractions appropriately. 
        if math.isnan(AbUptakeFr_NH4) :
            AbUptakeFr_NH4 = PNb
            AbUptakeFr_NO3 = 1.0 - AbUptakeFr_NH4

        '''
     Organic Nitrogen                     (mgN/day)
     dOrgN/dt =   Algae_OrgN              (Floating Algae -> OrgN)
                  - OrgN_NH4_Decay        (OrgN -> NH4)
                  - OrgN_Settling         (OrgN -> bed)
                  + Benthic Death         (Benthic Algae -> OrgN)														 
    ''' 
        if initial_values['use_OrgN'] :
            OrgN_NH4_Decay = kon_tc * OrgN
            OrgN_Settling = vson / depth * OrgN #TODO check the negative

            if initial_values['use_Algae'] :
                ApDeath_OrgN  = rna * ApDeath
            else :
                ApDeath_OrgN  = 0.0

            if initial_values['use_BAlgae'] :  
                AbDeath_OrgN = rnb * Fw * Fb * AbDeath / depth
            else : 
                AbDeath_OrgN = 0.0

            dOrgNdt = ApDeath_OrgN + AbDeath_OrgN - OrgN_NH4_Decay - OrgN_Settling 
        '''
        Ammonia Nitrogen (NH4)                 (mgN/day)
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
        if initial_values['use_NH4'] : #TODO this looks different
            if initial_values['use_DOX'] :
                NitrificationInhibition = 1.0 - math.exp(-KNR * DOX)
            else :
                NitrificationInhibition = 1.0

            NH4_Nitrification = NitrificationInhibition * knit_tc * NH4
    
            if initial_values['use_SedFlux'] :
                NH4fromBed = JNH4 / depth
            else :
                NH4fromBed = rnh4_tc / depth
    
            if initial_values['use_Algae'] :
                NH4_ApRespiration = rna * ApRespiration
                NH4_ApGrowth= ApUptakeFr_NH4 * rna * ApGrowth
            else : 
                NH4_ApRespiration = 0.0
                NH4_ApGrowth      = 0.0

            if initial_values['use_BAlgae'] : 
                NH4_AbRespiration = rnb *  Fb * AbRespiration / depth 
                NH4_AbGrowth      = AbUptakeFr_NH4 * rnb *  Fb * AbGrowth / depth 
            else :
                NH4_AbRespiration = 0.0
                NH4_AbGrowth      = 0.0

            if not initial_values['use_OrgN'] :
                OrgN_NH4_Decay = 0.0 
    
        dNH4dt = OrgN_NH4_Decay - NH4_Nitrification + NH4fromBed + NH4_ApRespiration - NH4_ApGrowth + NH4_AbRespiration - NH4_AbGrowth
 
        '''
        Nitrite Nitrogen  (NO3)                       (mgN/day)
        dNO3/dt  =      NH4 Oxidation                 (NH4 -> NO3) 
                        - NO3 Sediment Denitrification
                        - NO3 Algal Uptake            (NO3-> Floating Algae) 
                        - NO3 Benthic Algal Uptake    (NO3-> Benthic  Algae) 
        '''
        if initial_values['use_NO3'] :
            if initial_values['use_DOX'] :
                NO3_Denit = (1.0 - DOX / (DOX + KsOxdn)) * kdnit_tc * NO3
                if math.isnan(NO3_Denit) : 
                    NO3_Denit = kdnit_tc * NO3
            else :
                NO3_Denit    = 0.0

            if initial_values['use_SedFlux'] :
                NO3_BedDenit = JNO3 / depth
            else :
                    NO3_BedDenit = vno3_tc * NO3 / depth
   
            if initial_values['use_Algae'] :
                NO3_ApGrowth  = ApUptakeFr_NO3 * rna * ApGrowth
            else :
                NO3_ApGrowth  = 0.0
 
      # Benthic contribution
            if initial_values['use_BAlgae'] : 
                NO3_AbGrowth  = AbUptakeFr_NO3 * rnb *  Fb * AbGrowth / depth  
            else :
                NO3_AbGrowth  = 0.0
            if not initial_values['use_NH4'] :
                 NH4_Nitrification = 0.0
      
        dNO3dt = NH4_Nitrification - NO3_Denit - NO3_BedDenit - NO3_ApGrowth	- NO3_AbGrowth 


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
        if initial_values['use_NH4'] :
            DIN = DIN + NH4
            TKN = TKN + NH4

        if initial_values['use_NO3'] :
            DIN = DIN + NO3

        if initial_values['use_OrgN'] :
            TON = TON + OrgN

        if initial_values['use_Algae']:
            TON = TON + rna * Ap
            TKN = TKN + TON
            TN  = DIN + TON