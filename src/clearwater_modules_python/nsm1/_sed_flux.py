'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Sediment Flux
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
Last Revision Date: June 13, 2021
'''

import math
from collections import OrderedDict
from _temp_correction import TempCorrection

# Note: OrgN, NH4, NO3, OrgP, TIP, POC and DOX must be on for SedFlux.
class SedFlux:

  def __init__(self, global_module_choices, global_vars, global_vars_sedflux, sedFlux_constant_changes, CBOD_pathways, algae_pathways, Balgae_pathways):
      self.global_module_choices=global_module_choices
      self.global_vars=global_vars
      self.global_vars_sedflux = global_vars_sedflux

      self.algae_pathways=algae_pathways
      self.Balgae_pathways=Balgae_pathways
      self.CBOD_pathways = CBOD_pathways 

      self.sedFlux_constant_changes = sedFlux_constant_changes

      self.sedFlux_constants = OrderedDict()
      self.sedFlux_constants = {
        'Dd' : 0.0025,
        'Dp' : 0.00006,
        'vnh41' : 0.1313,
        'vno31' : 0.1,
        'vch41' : 0.7,
        'vh2sd': 0.2,
        'vh2sp': 0.4,
        'KPONG1' : 0.035,
        'KPONG2' : 0.0018,
        'KPOPG1' : 0.035,
        'KPOPG2' :0.0018,
        'KPOCG1' : 0.035,
        'KPOCG2' : 0.0018,
        'vno32' : 0.25,

        'KsDp' : 4.0, 
        'POCr' : 0.1,
        'res' : 0.001,
        'maxit' : 500,
        'kst' : 0.03,
        'So4_fresh' : 2,
        'focm2' : 0.4,
        'd_kpo41' : 20.0,
        'DOcr' : 2,
        'KsOxch' : 0.1,
        'KSh2s' : 4.0,
        'KsOxna1' : 0.37,
        'KsNh4' : 0.728,
        'Css1' : 0.5,
        'FPON1' : 0.15,
        'FPON2' : 0.35,
        'FPOP1' : 0.15,
        'FPOP2' : 0.35,
        'FPOC1' : 0.15, 
        'FPOC2' : 0.35,
        'FAP1' : 0.6,
        'FAP2' : 0.2,
        'FAB1' : 0.6,
        'FAB2' : 0.2,
        'FCBOD1' : 0.6,
        'FCBOD2' : 0.2,
        'kdnh42' : 1.0,
        'kdh2s2' : 100.0,
        'kdpo42' :20.0,
        'KsSO4' : 1.08,
        'Css2' : 0.5,
        'SedFlux_solution_option' : 1,
        'Methane_solution_option' : 1,
        'POCdiagenesis_part_option' : 1,
        'BFORmax' : 0.0,
      }

      for key in self.sedFlux_constant_changes.keys() :
        if key in self.sedFlux_constants:
          self.sedFlux_constants[key] = self.sedFlux_constant_changes[key]

  def Calculation (self) :
    '''
      use modGlobal,       only: r, nCBOD, use_Algae, use_BAlgae, use_OrgN, use_NH4, use_NO3, use_OrgP, use_TIP, use_POC,  & 
                              use_DIC, use_DOX, CBOD, Ap, OrgN, NH4, NO3, OrgP, TIP, POC, DIC, DOX, NH41, NO31, TIP1,   &
                              CH41, SO41, TH2S1, DIC1, POC2, PON2, POP2, NH42, NO32, TIP2, CH42, SO42, TH2S2, DIC2, ST, & 
                              HSO4, H2S1, DIP1, H2S2, DIP2, TPOC2, TPON2, TPOP2, POM2, depth, TsedC, Salinity, t, dt
      use modGlobalParam,  only: vs, vson, vsop, vsoc, fdp, h2, vb
      use modAlgae,        only: vsap, rca, rna, rpa
      use modBenthicAlgae, only: AbDeath, rnb, rpb, rcb, Fw
      use modCBOD,         only: ksbod_tc
      use modDLL,          only: R8, nRegion, AddIndex, list_class, SetOptionalIndex, TempCorrectionStruct, Arrhenius_TempCorrection
      '''

    '''
      #correction coefficients (14)
      Dd = 0         #pore-water diffusion coefficient between layer 1 and 2                 (m2/d)  
      Dd_tc = 0   
      Dp = 0          # particle mixing diffusion coefficient between layer 1 and 2           (m2/d)  
      Dp_tc = 0
      vnh41 = 0       #nitrification reaction velocity in sediment layer 1                    (m/d)  
      vnh41_tc = 0  
      vno31  = 0      # denitrification reaction velocity in sediment layer 1                 (m/d)  
      vno31_tc= 0   
      vch41   = 0     # methane oxidation reaction velocity in sediment layer 1               (m/d)  
      vch41_tc= 0   
      vh2sd  = 0      # dissolved sulfide oxidation reaction velocity in sediment layer 1     (m/d)  
      vh2sd_tc = 0  
      vh2sp  = 0       # particulate sulfide oxidation reaction velocity in sediment layer 1   (m/d)  
      vh2sp_tc = 0  

      KPONG1  = 0     # diagenesis rate of PON G1 in sediment layer 2                         (1/d)  
      KPONG2 = 0      # diagenesis rate of PON G2 in sediment layer 2                         (1/d)  
      KPON_tc = 0
      KPOPG1 = 0      # diagenesis rate of POP G1 in sediment layer 2                         (1/d)  
      KPOPG2 = 0      # diagenesis rate of POP G2 in sediment layer 2                         (1/d)  
      KPOP_tc= 0   
      KPOCG1 = 0      # diagenesis rate of POC G1 in sediment layer 2                         (1/d)  
      KPOCG2 = 0      # diagenesis rate of POC G2 in sediment layer 2                         (1/d)  
      KPOC_tc= 0   
      vno32  = 0      # denitrification reaction velocity in sediment layer 2                 (m/d)
      vno32_tc = 0  

      res = 0            # relative error of SOD solution
      maxit = 0          # max allowable iteration number for SOD solution
      kst= 0             # decay rate of benthic stress                                      (1/d)
      SO4_fresh= 0       # SO4 concentration of overlaying water column in freshwater        (mg-O2/L)
      focm2= 0           # ratio of carbon to organic matter in bed sediment                 (unitless)

      d_kpo41= 0         # factor that increases the aerobic layer phosphate partition coefficient  (unitless)
      DOcr= 0            # critical oxygen concentration for incremental phosphate sorption  (mg-O2/L) [avoid to repeat with DOC]
      KsOxna1= 0         # half-saturation oxygen constant for sediment nitrification        (mg-O2/L)
      KsNh4= 0           # half-saturation ammonia constant for sediment nitrification       (mg-N/L)
      KsOxch= 0          # half-saturation coefficient for oxygen in oxidation of methane    (mg-O2/L)
      KSh2s= 0           # sulfide oxidation normalization constant                          (mg-O2/L)
      KsDp= 0            # half-saturation constant for oxygen in particle mixing            (mg-O2/L)
      POCr= 0            # reference POC G1 concentration for bioturbation                   (mg-C/L)
      Css1= 0            # solids concentration in sediment layer 1                          (kg/L)
      
      FPON1= 0           # fraction of settled RPON to sediment PON G1                 (unitless)
      FPON2= 0           # fraction of settled RPON to sediment PON G2                 (unitless)
      FPOP1= 0           # fraction of settled RPOP to sediment POP G1                 (unitless)
      FPOP2= 0           # fraction of settled RPOP to sediment POP G2                 (unitless)
      FPOC1= 0           # fraction of settled RPOC to sediment POC G1                 (unitless)
      FPOC2= 0           # fraction of settled RPOC to sediment POC G2                 (unitless)
      FAP1= 0            # fraction of settled algae to G1                             (unitless) 
      FAP2= 0            # fraction of settled algae to G2                             (unitless)
      FAB1= 0            # fraction of benthic algae death to G1                       (unitless)
      FAB2= 0            # fraction of benthic algae death to G2                       (unitless)
      FCBOD1= 0          # fraction of CBOD sedimentation to G1                        (unitless)
      FCBOD2= 0          # fraction of CBOD sedimentation to G2                        (unitless)
      kdnh42= 0          # partition coefficient for ammonium in sediment layer 2      (L/kg)
      kdh2s2= 0          # partition coefficient for sulfide in sediment layer 2       (L/kg)
      kdpo42= 0          # partition coefficient for inorganic P in sediment layer 2   (L/kg)
      KsSO4= 0           # half-saturation constant for sulfate in sulfate reduction   (mg-O2/L)
      Css2= 0            # solids concentration in sediment layer 2                    (kg/L)

  # integer parameter (3)
      SedFlux_solution_option=1        # numerical method (1 steady and 2 unsteady)
      Methane_solution_option=1        # method for solving methane concentration (1 analytical and 2 numerical)
      POCdiagenesis_part_option=1      # method for partitioing carbon diagenesis flux into methane and sulfide (1 half-saturation and 2 sulfate reduction depth)

  # pathway (72)
      JPOC= 0                    # depositional flux of POC G1/G2/G3 from overlying water column to bed sediments (g-C/m2/d)
      POC2_Diagenesis[3]= 0         # sediment disgenesis of POC G1/G2/G3 in sediment layer 2                        (g-C/m2/d)
      POC2_Burial= 0             # burial of POC G1/G2/G3 in sediment layer 2                                     (g-C/m2/d)
      JC= 0                      # total sediment diagenesis flux of POC                                          (g-C/m2/d)
      JC_dn= 0                   # carbon (oxygen equivalents) consumed by denitrification  (g-O2/m2/d)

      JPON= 0                    # depositional flux of PON G1/G2/G3 from overlying water column to bed sediments (g-N/m2/d)
      PON2_Diagenesis= 0         # sediment disgenesis of PON G1/G2/G3 in sediment layer 2                        (g-N/m2/d)
      PON2_Burial= 0             # burial of PON G1/G2/G3 in sediment layer 2                                     (g-N/m2/d)
      JN = 0                         # total sediment diagenesis flux of PON                                          (g-N/m2/d)

      JPOP= 0                   # depositional flux of POP G1/G2/G3 from water column to bed sediments           (g-P/m2/d)
      POP2_Diagenesis= 0         # sediment disgenesis of POP G1/G2/G3 in sediment layer 2                        (g-P/m2/d)
      POP2_Burial= 0             # burial of POP G1/G2/G3 in sediment layer 2                                     (g-P/m2/d)
      JP                         # total sediment diagenesis flux of POP                                          (g-P/m2/d)

      w12= 0                        # particle mixing mass transfer coefficient                                  (m/d)
      KL12= 0                       # mass transfer velocity between the two sediment layers                     (m/d)
      KL01= 0                       # mass transfer velocity between overlying water and the aerobic layer       (m/d)
      SOD_Bed= 0                    # SedFlux sediment oxygen demand                                             (g-O2/m2/d)

      JNH4= 0                       # sediment-water flux of ammonia                           (g-N/m2/d)
      TNH41_Burial= 0               # burial of TNH41 in sediment layer 1                      (g-N/m2/d)   
      NH41_Nitrification= 0         # nitrification of TNH41 in sediment layer 1               (g-N/m2/d)
      PNH41_PNH42= 0                # mass transfer between TNH41 and TNH42 in adsorbed form   (g-N/m2/d)
      NH41_NH42= 0                  # mass transfer between TNH41 and TNH42 in dissolved form  (g-N/m2/d)
      TNH42_Burial= 0               # burial of TNH42 in sediment layer 2                      (g-N/m2/d)

      JNO3= 0                       # sediment-water flux of nitrate                           (g-N/m2/d)
      NO31_Denit= 0                 # denitrification of NO31 in sediment layer 1              (g-N/m2/d)
      NO31_NO32= 0                  # mass transfer between NO31 and NO32                      (g-N/m2/d)
      NO32_Denit= 0                 # denitrification of NO32 in sediment layer 2              (g-N/m2/d)

      JCH4= 0                       # sediment-water flux of methane                           (g-O2/m2/d)
      JCH4g= 0                      # methane loss as bubbles from sediment                    (g-O2/m2/d)
      CH4sat= 0                     # saturated concentration of methane in oxygen equivalents (mg-O2/L)
      JCc_CH4= 0                    # carbon diagenesis flux consumed for methane formation    (g-O2/m2/d)
      CH41_Oxidation= 0             # methane oxidation in sediment layer 1                    (g-O2/m2/d)

      JSO4= 0                       # sediment-water flux of sulfate                           (g-O2/m2/d)
      JCc_SO4= 0                    # carbon diagenesis flux consumed for sulfate reduction    (g-O2/m2/d)
      SO41_SO42= 0                  # mass transfer between SO41 and SO42                      (g-O2/m2/d) 

      JH2S= 0                       # sediment-water flux of sulfide                           (g-O2/m2/d)
      H2S1_Oxidation= 0             # sulfide oxidation in sediment layer 1                    (g-O2/m2/d)
      TH2S1_Burial= 0               # burial of H2S1 in sediment layer 1                       (g-O2/m2/d)
      H2S1_H2S2= 0                  # mass transfer between H2S1 and H2S2 in dissolved form    (g-O2/m2/d)
      PH2S1_PH2S2= 0                # mass transfer between H2S1 and H2S2 in adsorbed form     (g-O2/m2/d)
      TH2S2_Burial= 0               # burial of H2S2 in sediment layer 2                       (g-O2/m2/d)

      JDIC= 0                       # sediment-water flux of dissolved inorganic carbon        (g-C/m2/d)
      DIC1_CH41_Oxidation= 0        # DIC1 produced by CH41 oxidation in sediment layer 1      (g-C/m2/d)
      DIC1_NO31_Denit= 0            # DIC1 produced by NO31 denitrification in sediment layer 1(g-C/m2/d) 
      DIC1_DIC2= 0                  # mass transfer between DIC1 and DIC2 in dissolved form    (g-C/m2/d)
      DIC2_POC2_SO42= 0             # DIC2 produced by sulfate reduction in sediment layer 2   (g-C/m2/d)
      DIC2_CH42= 0                  # DIC2 produced by mathene formation in sediment layer 2   (g-C/m2/d)
      DIC2_NO32_Denit= 0            # DIC2 produced by NO32 denitrification in sediment layer 2(g-C/m2/d)

      JDIP= 0                       # sediment-water flux of phosphate                         (g-P/m2/d)
      TIP_TIP2= 0                   # settling of PIP of water column into PIP2 in layer 2     (g-P/m2/d)
      TIP1_Burial= 0                # burial of TIP1 in sediment layer 1                       (g-P/m2/d)
      PIP1_PIP2= 0                  # mass transfer between TIP1 and TPO42 in adsorbed form    (g-P/m2/d)
      DIP1_DIP2= 0                  # mass transfer between TIP1 and TPO42 in dissolved form   (g-P/m2/d)
      TIP2_Burial= 0                # burial of TIP2 in sediment layer 2                       (g-P/m2/d)


  # local variables
      TNH41, TNH42 = 0                          # total ammonia concentration in sediment layer   (mg-N/L)
      Si1, Si2 = 0                              # total silica concentration in sediment layer     (mg-Si/L)
      CSOD_CH4, CSOD_H2S, NSOD  = 0            # carbonaceous oxygen demand, nitrogenous oxygen demand (mg-O2/L)
      fd1, fd2= 0                              # fraction of inorganic matter (ammonia, phosphate) in dissolved form in sediment layer  
      fp1, fp2= 0                              # fraction of inorganic matter (ammonia, phosphate) in particulate form in sediment layer   
      FNH4= 0                                  # modification of nitrification reaction in layer 1 
      FOxna= 0                                 # nitrification attenuation due to low oxygen in layer 1 
      FOxch= 0                                 # methane oxidation attenuation due to low oxygen in layer 1 
      rondn = 32.0 / 12.0 * 10.0 / 8.0 * 12.0 / 14.0  # oxygen stoichiometric coeff for denitrification      (g-O2/g-N)
      rcdn  = 5.0 * 12.0 / (4.0 * 14.0)               # carbon stoichiometric coeff for denitrification      (g-C/g-N)
      roc   = 32.0 / 12.0                             # oxygen stoichiometric coeff for organic carbon decay (g-O2/g-C)
      ron   = 2.0 * 32.0 / 14.0                       # oxygen stoichiometric coeff for nitrification        (g-O2/g-N)
      roso4 = 2.0 * 32.0 / 98.0                       # oxygen stoichiometric coeff for sulfate SO4          (g-O2/g-SO4)
      a12_TNH4, a21_TNH4, a22_TNH4, b2_TNH4= 0           # matrix coefficient for TNH41, TNH42
      a12_NO3, a21_NO3, a22_NO3, b2_NO3= 0               # matrix coefficient for NO31, NO32
      a12_CH4, a21_CH4, a22_CH4= 0                       # matrix coefficient for CH41, CH42
      a12_SO4, a21_SO4, a22_SO4= 0                       # matrix coefficient for SO41, SO42
      a12_TH2S, a21_TH2S, a22_TH2S= 0                    # matrix coefficient for TH2S1, TH2S2
      a11, a12, b1, a21, a22, b2= 0                       
      con_nit, con_cox, con_sox= 0                        
      fds1, fps1, fds2, fps2= 0                          # dissolved and particulate fraction for H2S1 and H2S2 in layer 1 and 2
      CH42_prev, TH2S2_prev= 0       

      JCc= 0                         # carbon diagenesis flux corrected for denitrification     (g-O2/m2/d)
      CSODmax= 0                     # used for analytical soluton of methane                   (g-O2/m2/d)
      TempBen = 10.0              # critical temperature for benthic stress                  (oC)
      ISWBEN= 0                      # check if water temperature of last time step < TempBen
      FORmax= 0                     # maximum benthic stress oxygen correction coefficient
      BFORmax= 0

      CH4 = 0.0                   # methane concentration in overlying water column                (mg-O2/L)
      H2S = 0.0                   # H2S concentration in overlying water column                    (mg-O2/L)

      SO4= 0                         # sulfate concentration in overlying water column                (mg-O2/L)
      KL12SO4= 0                     # dissolved mass transfer velocity of sulfate between two layers (m/day)
      SO42_prev, HSO4_prev, TH2S1_prev= 0
      i
    '''
    #ZZ h2, KL01, KSh2s, POCr, CSS2, KsOxch should not be set to zero. 

    Dd = self.sedFlux_constants['Dd']
    Dp = self.sedFlux_constants['Dp']
    vnh41 = self.sedFlux_constants['vnh41']
    vno31= self.sedFlux_constants['vno31']
    vch41= self.sedFlux_constants['vch41']
    vh2sd= self.sedFlux_constants['vh2sd']
    vh2sp= self.sedFlux_constants['vh2sp']
    KPONG1= self.sedFlux_constants['KPONG1']
    KPONG2= self.sedFlux_constants['KPONG2']
    KPOPG1= self.sedFlux_constants['KPOPG1']
    KPOPG2= self.sedFlux_constants['KPOPG2']
    KPOCG1= self.sedFlux_constants['KPOCG1']
    KPOCG2= self.sedFlux_constants['KPOCG2']
    vno32= self.sedFlux_constants['vno32']

    KsDp= self.sedFlux_constants['KsDp'] 
    POCr= self.sedFlux_constants['POCr']
    res= self.sedFlux_constants['res']
    maxit= self.sedFlux_constants['maxit']
    kst= self.sedFlux_constants['kst']
    SO4_fresh= self.sedFlux_constants['So4_fresh']
    focm2= self.sedFlux_constants['focm2']
    d_kpo41= self.sedFlux_constants['d_kpo41']
    DOcr= self.sedFlux_constants['DOcr']
    KsOxch= self.sedFlux_constants['KsOxch']
    KSh2s= self.sedFlux_constants['KSh2s']
    KsOxna1= self.sedFlux_constants['KsOxna1']
    KsNh4= self.sedFlux_constants['KsNh4']
    Css1= self.sedFlux_constants['Css1']
    FPON1= self.sedFlux_constants['FPON1']
    FPON2= self.sedFlux_constants['FPON2']
    FPOP1= self.sedFlux_constants['FPOP1']
    FPOP2= self.sedFlux_constants['FPOP2']
    FPOC1= self.sedFlux_constants['FPOC1'] 
    FPOC2= self.sedFlux_constants['FPOC2']
    FAP1= self.sedFlux_constants['FAP1']
    FAP2= self.sedFlux_constants['FAP2']
    FAB1= self.sedFlux_constants['FAB1']
    FAB2= self.sedFlux_constants['FAB2']
    FCBOD1= self.sedFlux_constants['FCBOD1']
    FCBOD2= self.sedFlux_constants['FCBOD2']
    kdnh42= self.sedFlux_constants['kdnh42']
    kdh2s2= self.sedFlux_constants['kdh2s2']
    kdpo42= self.sedFlux_constants['kdpo42']
    KsSO4= self.sedFlux_constants['KsSO4']
    Css2= self.sedFlux_constants['Css2']
    SedFlux_solution_option= self.sedFlux_constants['SedFlux_solution_option']
    Methane_solution_option= self.sedFlux_constants['Methane_solution_option']
    POCdiagenesis_part_option= self.sedFlux_constants['POCdiagenesis_part_option']
    BFORmax= self.sedFlux_constants['BFORmax']

    vsoc=self.global_vars['vsoc']
    POC= self.global_vars['POC']
    vson=self.global_vars['vson']
    OrgN=self.global_vars['OrgN']
    vsop=self.global_vars['vsop'] 
    OrgP=self.global_vars['OrgP']
    rca=self.algae_pathways['rca'] 
    Ap=self.global_vars['Ap'] 
    vsap=self.algae_pathways['vsap']
    rna=self.algae_pathways['rna']
    rpa=self.algae_pathways['rpa']
    AbDeath=self.algae_pathways['AbDeath']
    Fw=self.Balgae_pathways['Fw']
    rcb=self.Balgae_pathways['rcb']
    rnb=self.Balgae_pathways['rnb']
    rpb=self.Balgae_pathways['rpb']
    dt=self.global_vars['dt']
    TsedC=self.global_vars['TsedC']
    vb= self.global_vars['vb']
    h2= self.global_vars['h2']
    DOX=self.global_vars['DOX']

    DOX = max(DOX, 0.01) #TODO check this, why this too? was in the subroutine that called all the subroutines

    DIC = self.global_vars['DIC']
    fdp = self.global_vars['fdp']
    TIP = self.global_vars['TIP']
    vs = self.global_vars['vs']
    Salinity = self.global_vars['Salinity']
    depth = self.global_vars['depth']
    NH4 = self.global_vars['NH4']
    NO3 = self.global_vars['NO3']
    CBOD=self.global_vars['CBOD']

    NH41 = self.global_vars_sedflux['NH41']
    NH42 = self.global_vars_sedflux['NH42']
    NO31 = self.global_vars_sedflux['NO31']
    NO32 = self.global_vars_sedflux['NO32']
    TH2S1 = self.global_vars_sedflux['TH2S1']
    TH2S2 = self.global_vars_sedflux['TH2S2']
    CH41 = self.global_vars_sedflux['CH41']
    CH42 = self.global_vars_sedflux['CH42']
    DIC1 = self.global_vars_sedflux['DIC1']
    DIC2 = self.global_vars_sedflux['DIC2']
    TIP1 = self.global_vars_sedflux['TIP1']
    TIP2 = self.global_vars_sedflux['TIP2']
    SO41 = self.global_vars_sedflux['SO41']
    SO42 = self.global_vars_sedflux['SO42']

    POC2=[0]*3        #concentration of sediment particulate organic carbon
    PON2=[0]*3        #concentration of sediment particulate organic nitrogen
    POP2=[0]*3        #concentration of sediment particulate organic phosphrous 

    POC2[1] = self.global_vars_sedflux['POC2_1']
    POC2[2] = self.global_vars_sedflux['POC2_2']
    POC2[3] = self.global_vars_sedflux['POC2_3']

    PON2[1] = self.global_vars_sedflux['PON2_1']
    PON2[2] = self.global_vars_sedflux['PON2_2']
    PON2[3] = self.global_vars_sedflux['PON2_3']

    POP2[1] = self.global_vars_sedflux['POP2_1']
    POP2[2] = self.global_vars_sedflux['POP2_2']
    POP2[3] = self.global_vars_sedflux['POP2_3']

    t=self.global_vars_sedflux['t']

    ksbod_tc = self.CBOD_pathways['ksbod_tc']

    #local variables
    roc = 32.0 / 12.0
    TempBen = 10
    roso4 = 2.0 * 32.0 / 98.0
    rondn = 32.0 / 12.0 * 10.0 / 8.0 * 12.0 / 14.0
    H2S = 0.0
    CH4 = 0.0
    ron   = 2.0 * 32.0 / 14.0
    rcdn  = 5.0 * 12.0 / (4.0 * 14.0)

    KPON_tc=[0]*3
    KPOP_tc=[0]*3
    KPOC_tc = [0]*3
  
    Dd_tc=TempCorrection(Dd, 1.08).arrhenius_correction(TsedC) 
    Dp_tc=TempCorrection(Dp, 1.117).arrhenius_correction(TsedC)

    vnh41_tc= vnh41 * 1.123 * ((TsedC-20)/2)
    vno31_tc= vno31 * 1.08 * ((TsedC-20)/2)
    vch41_tc = vch41 * 1.079 * ((TsedC-20)/2)
    vh2sd_tc = vh2sd * 1.079 * ((TsedC-20)/2)
    vh2sp_tc = vh2sp * 1.079 * ((TsedC-20)/2)

    #only uses arrhenius correction and passes two arguments (rc20 and TsedC not theta)
    KPON_tc[1] = TempCorrection(KPONG1, 1.1).arrhenius_correction(TsedC)
    KPON_tc[2] = TempCorrection(KPONG2, 1.15).arrhenius_correction(TsedC)
    KPON_tc[3] = 0
    KPOP_tc[1] = TempCorrection(KPOPG1, 1.1).arrhenius_correction(TsedC)
    KPOP_tc[2] = TempCorrection(KPOPG2, 1.15).arrhenius_correction(TsedC)
    KPOP_tc[3] = 0
    KPOC_tc[1] = TempCorrection(KPOCG1, 1.1).arrhenius_correction(TsedC)
    KPOC_tc[2] = TempCorrection(KPOCG2, 1.15).arrhenius_correction(TsedC)
    KPOC_tc[3] = 0

    vno32_tc= TempCorrection(vno32, 1.08).arrhenius_correction(TsedC)     #TODO why is this different from vno31 

    FAP3, FAB3 =0    #fraction algal/benthic algal death to G3

    #compute JPOC, JPON, JPOP total depositional flux to sediment of particulate oragnic matter
 
    JPOC=[0]*3
    JPON=[0]*3
    JPOP=[0]*3

    #Carbon settling from water (settling velocity * particulate organic carbon * fraction settling to G1/G2/G3 sediment)
    if self.global_module_choices['use_POC'] :
      JPOC[1] = vsoc * POC * FPOC1                  #(g-C/d*m^2)
      JPOC[2]=  vsoc * POC * FPOC2
      JPOC[3]=  vsoc * POC * (1 - FPOC1 - FPOC2)
    else :
      JPOC[1]= 0.0
      JPOC[2]= 0.0
      JPOC[3]= 0.0

    if self.global_module_choices['use_OrgN'] :
      JPON[1] = vson * OrgN * FPON1                  #(g-N/d*m^2)
      JPON[2] = vson * OrgN * FPON2
      JPON[3] = vson * OrgN * (1 -FPON1 - FPON2)
    else :
      JPON[1] = 0.0
      JPON[2] = 0.0
      JPON[3] = 0.0

    if self.global_module_choices['use_OrgP'] :
      JPOP[1] = vsop * OrgP * FPOP1                  #(g-P/d*m^2)
      JPOP[2] = vsop * OrgP * FPOP2
      JPOP[3] = vsop * OrgP * (1 - FPOP1 - FPOP2)
    else :
      JPOP[1] = 0.0
      JPOP[2] = 0.0
      JPOP[3] = 0.0
    
    #CBOD sedimentation 
    JPOC[1] = JPOC[1] + ksbod_tc * FCBOD1 * CBOD / roc
    JPOC[2] = JPOC[2] + ksbod_tc * FCBOD2 * CBOD / roc
    JPOC[3] = JPOC[3] + ksbod_tc * (1.0 - FCBOD1 - FCBOD2) * CBOD / roc

    #Algae settling 
    if self.global_module_choices['use_Algae'] :
      FAP3 = 1.0 - FAP1 - FAP2
      JPOC[1]  = JPOC[1] + FAP1 * rca * Ap * vsap     #g-C/d*m^2  
      JPOC[2]  = JPOC[2] + FAP2 * rca * Ap * vsap
      JPOC[3]  = JPOC[3] + FAP3 * rca * Ap * vsap
      JPON[1]  = JPON[1] + FAP1 * rna * Ap * vsap
      JPON[2]  = JPON[2] + FAP2 * rna * Ap * vsap
      JPON[3]  = JPON[3] + FAP3 * rna * Ap * vsap
      JPOP[1]  = JPOP[1] + FAP1 * rpa * Ap * vsap
      JPOP[2]  = JPOP[2] + FAP2 * rpa * Ap * vsap 
      JPOP[3]  = JPOP[3] + FAP3 * rpa * Ap * vsap

    #Benthic algae death
    if self.global_module_choices['use_BAlgae'] :
      FAB3 = 1.0 - FAB1 - FAB2
      JPOC[1] = JPOC[1] + AbDeath * (1.0 - Fw) * rcb * FAB1     #g-C/d*m^2 
      JPOC[2] = JPOC[2] + AbDeath * (1.0 - Fw) * rcb * FAB2
      JPOC[3] = JPOC[3] + AbDeath * (1.0 - Fw) * rcb * FAB3
      JPON[1] = JPON[1] + AbDeath * (1.0 - Fw) * rnb * FAB1
      JPON[2] = JPON[2] + AbDeath * (1.0 - Fw) * rnb * FAB2
      JPON[3] = JPON[3] + AbDeath * (1.0 - Fw) * rnb * FAB3
      JPOP[1] = JPOP[1] + AbDeath * (1.0 - Fw) * rpb * FAB1
      JPOP[2] = JPOP[2] + AbDeath * (1.0 - Fw) * rpb * FAB2
      JPOP[3] = JPOP[3] + AbDeath * (1.0 - Fw) * rpb * FAB3

  #compute POC2/PON2/POP2, POC2/PON2/POP2 pathways and JC/JN/JP
    POC2=[0]*3        #concentration of sediment particulate organic carbon
    PON2=[0]*3        #concentration of sediment particulate organic nitrogen
    POP2=[0]*3        #concentration of sediment particulate organic phosphrous 
    POC2_Diagenesis=[0]*3
    PON2_Diagenesis=[0]*3
    POP2_Diagenesis=[0]*3
    POC2_Burial=[0]*3
    PON2_Burial=[0]*3
    POP2_Burial=[0]*3
  
  #POC G1/G2/G3 in the second layer. Able for diagenesis and depends on depositional flux above 
    for i in (1, 3) :
      if self.sedFlux_constants['SedFlux_solution_option'] == 1 :       # steady state solution
        POC2[i] = JPOC[i] / (KPOC_tc[i] * h2 + vb)        #mg-C/L
        PON2[i] = JPON[i] / (KPON_tc[i] * h2 + vb)
        POP2[i] = JPOP[i] / (KPOP_tc[i] * h2 + vb)
      elif (self.sedFlux_constants['SedFlux_solution_option'] == 2):    #unsteady state solution
        POC2[i] = (JPOC[i] + POC2[i] * h2 / dt) / (h2 / dt + KPOC_tc[i] * h2 + vb)      #mg-C/L
        PON2[i] = (JPON[i] + PON2[i] * h2 / dt) / (h2 / dt + KPON_tc[i] * h2 + vb)
        POP2[i] = (JPOP[i] + POP2[i] * h2 / dt) / (h2 / dt + KPOP_tc[i] * h2 + vb)

      if math.isnan(POC2[i]) :
        POC2[i] = 0.0
      if math.isnan(PON2[i]):
        PON2[i] = 0.0
      if math.isnan(POP2[i]):
        POP2[i] = 0.0

      POC2[i] = max(POC2[i], 0.0) 
      PON2[i] = max(PON2[i], 0.0) 
      POP2[i] = max(POP2[i], 0.0) 
      
      POC2_Diagenesis[i] = KPOC_tc[i] * h2 * POC2[i]      # g-C/d*m^2
      PON2_Diagenesis[i] = KPON_tc[i] * h2 * PON2[i]
      POP2_Diagenesis[i] = KPOP_tc[i] * h2 * POP2[i]
      POC2_Burial[i] = vb * POC2[i]
      PON2_Burial[i] = vb * PON2[i]
      POP2_Burial[i] = vb * POP2[i]

    #Calculate sediment diagenesis flux C, N, P
    JC = POC2_Diagenesis[1] + POC2_Diagenesis[2]      # g-C/d*m^2
    JN = PON2_Diagenesis[1] + PON2_Diagenesis[2]
    JP = POP2_Diagenesis[1] + POP2_Diagenesis[2]

    # compute SOD
    SOD_old =0
    # coefficients of a quadratic equation
    ra2, ra1, ra0 =0
    #root of the quadratic equation      
    sn1, disc, r1, r2 =0     
  
    SOD_Bed = JC * roc + 1.714 * JN
    
    #Partical mixing transfer velocity: Transfer for NH4, H2S, and PIP between layer 1 and 2
    if SedFlux_solution_option == 1 :
      w12 = Dp_tc / (0.5 * h2)         
    elif SedFlux_solution_option == 2 :
      w12 = Dp_tc / (0.5 * h2) * (POC2[1] / (1000.0 * POCr * Css2)) #TODO the 1000 does not cancel with anything if that is what it is supposed to do
    if math.isnan(w12):
      w12 = 0.0
  
    #Sediment Benthic Stress: Low DO will eliminate bioturbation. Particle phase mixing coefficient is modified. 
    if SedFlux_solution_option == 1 :
      ST = 0.0                          #Benthic stress does not accumulate when temperature is lower than TempBen=10oC.
    elif SedFlux_solution_option == 2 :
    
      '''
      This was commented out in the FORTRAN
       if (ISWBEN) then
         if (TsedC >= TempBen) then
           ISWBEN   = .false
           BFORmax  = 0.0
         end if
         ST      = (ST + dt * KsDp(r) / (KsDp(r) + DOX)) / (1.0 + kst(r) * dt)
         if (isnan(ST)) ST = 0.0 
       else
         if (TsedC < TempBen) ISWBEN = .true.
         BFORmax = max(BFORmax, KsDp(r) / (KsDp(r) + DOX))
         ST      = (ST + dt * BFORmax) / (1.0 + kst(r) * dt)
         if (isnan(ST)) ST = 0.0 
       end if
      '''
  
      if TsedC < TempBen :
        ST = (ST + dt * KsDp / (KsDp + DOX)) / (1.0 + kst * dt)      
        BFORmax = 0.0
      else :  
        BFORmax = max(BFORmax, KsDp / KsDp + DOX)
        ST = (ST + dt * BFORmax) / (1.0 + kst * dt)   

      if math.isnan(ST) :
        ST = 0.0
 
    w12  = w12 * (1.0 - kst * ST)
    
    if (w12 < 0.0) :   
      w12  = 0.0

    #Dissolved and particulate phase mixing coefficient between layer 1 and layer 2
    KL12 = Dd_tc / (0.5 * h2)         
    if math.isnan(KL12):
      KL12 = 0.0

    #Dissolved constituent exchange between water and layer 1
    KL01 = SOD_Bed / DOX                
    if math.isnan(KL01) or KL01 == 0.0:
      KL01 = 1.0E-8

    #Sediment ammonium (TNH41 and TNH42) produced by decomposition of reactive G1 and G2 classes of PON in layer 2
    #Calculate dissolved and particulate fraction
    fd1  = 1.0 / (1.0 + Css1 * kdnh42)      #dissolved fraction NH4 layer 1
    fd2  = 1.0 / (1.0 + Css2 * kdnh42)      #dissolved fraction NH4 layer 2
    fp1  = 1.0 - fd1                        #particulate fraction NH4 layer 1
    fp2  = 1.0 - fd2                        #particulate fraction NH4 layer 2
    TNH41 = NH41 / fd1                      #total concentration NH4 dissolved layer 1
    TNH42 = NH42 / fd2                      #total concentration NH4 dissolved layer 2

    FOxna = DOX / (KsOxna1 * 2.0 + DOX)     #Oxygen attenuation factor for sediment nitrification 

    if math.isnan(FOxna):
      FOxna = 0.0
    
    con_nit  = vnh41_tc * vnh41_tc * FOxna * fd1    
    
    #coefficents for implicit finite difference form for TNH4 (a11, a12, a21, a22, b1, b2)

    #NH41 and NH42
    a12_TNH4 = -w12 * fp2 - KL12 * fd2
    a21_TNH4 = -w12 * fp1 - KL12 * fd1 - vb
    if SedFlux_solution_option == 1 :         #steady
      a22_TNH4 = -a12_TNH4 + vb
      b2_TNH4  = JN
    elif SedFlux_solution_option == 2:        #unsteady
      a22_TNH4 = -a12_TNH4 + vb + h2 / dt
      b2_TNH4  = JN + h2 * TNH42 / dt

    # NO31 and NO32
    a12_NO3 = -KL12
    a21_NO3 = -KL12
    if SedFlux_solution_option == 1 : 
      a22_NO3 = KL12 + vno32_tc
      b2_NO3  = 0.0
    elif SedFlux_solution_option == 2 :
      a22_NO3 = KL12 + vno32_tc + h2 / dt
      b2_NO3  = h2 * NO32 / dt

    # SO41, SO42 and TH2S1, TH2S2
    # compute water column SO4
    if (Salinity > 0.01) :
      # SO4 is estimated based on salinity for salt water
      SO4 = (20.0 + 27.0 / 190.0 * 607.445 * Salinity) * roso4    #TODO find this formula
    else:
      #SO4 is user-defined if salinity is not simulated for fresh water. 
      SO4 = SO4_fresh

    # Half-saturation method
    if POCdiagenesis_part_option == 1 :  
      SO42_prev = SO42
      TH2S2_prev = TH2S2
      
      fds1 = 1.0 / (1.0 + Css1 * kdh2s2)      #dissolved fraction layer 1
      fds2 = 1.0 / (1.0 + Css2 * kdh2s2)      #dissolved fraction layer 2
      fps1 = 1.0 - fds1                       #particulate fraction layer 1
      fps2 = 1.0 - fds2                       #particulate fraction layer 2

      con_sox = (vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc * vh2sp_tc * fps1) * DOX / 2.0 / KSh2s
      if math. isnan(con_sox):
        con_sox = (vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc * vh2sp_tc * fps1)
    
    # Sulfate reduction depth method
    elif POCdiagenesis_part_option == 2 :   
      #Set initial value for HSO4
      if (t < 1.0E-10) :
        HSO4 = math.sqrt(Dd_tc * SO4 * h2 / (max(roc * JC, 1.0E-10)))         #sulfate penetration into layer 2 from layer 1
        if (HSO4 > h2) :
          HSO4 = h2

      #SO41 and SO42
      TH2S1_prev = TH2S1
      HSO4_prev  = HSO4
      SO42_prev  = SO42
      a12_SO4    = - KL12
      a21_SO4    = KL12

      if SedFlux_solution_option == 1 :
        a22_SO4  = - KL12
      elif SedFlux_solution_option == 2 :
        a22_SO4  = - KL12 - h2 / dt

      # TH2S1 and TH2S2
      fds1 = 1.0 / (1.0 + Css1 * kdh2s2)
      fds2 = 1.0 / (1.0 + Css2 * kdh2s2)
      fps1 = 1.0 - fds1
      fps2 = 1.0 - fds2
      TH2S2_prev = TH2S2
      con_sox = (vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc * vh2sp_tc * fps1) * DOX / 2.0 / KSh2s
      if math.isnan(con_sox) :
        con_sox = (vh2sd_tc * vh2sd_tc * fds1 + vh2sp_tc * vh2sp_tc * fps1)
      
      a12_TH2S = -w12 * fps2 - KL12 * fds2
      a21_TH2S = -w12 * fps1 - KL12 * fds1 - vb
      if SedFlux_solution_option == 1 : 
        a22_TH2S = -a12_TH2S + vb
      elif SedFlux_solution_option == 2 :
        a22_TH2S = -a12_TH2S + vb + h2 / dt

    # CH41 and CH42
    CH4sat = 100.0 * (1.0 + depth / 10.0) * 1.024**(20.0 - TsedC)
    if Methane_solution_option ==2 :
      CH42_prev = CH42
      FOxch = DOX / (KsOxch * 2.0 + DOX)
      if math.isnan(FOxch) :
        FOxch = 0.0
      con_cox = vch41_tc * vch41_tc * FOxch
      a12_CH4 = -KL12
      a21_CH4 = -KL12
      if SedFlux_solution_option == 1 :
        a22_CH4 = KL12
      elif SedFlux_solution_option == 2 :
        a22_CH4 = KL12 + h2 / dt

    #compute SOD
    for i in (1, int(maxit)) :
      #TNH41 and TNH42
      if (KsNh4 > 0.0) :
        FNH4 = KsNh4 / (KsNh4 + fd1 * TNH41)
      else :
        FNH4  = 1.0

      a11 = -a21_TNH4 + con_nit * FNH4 / KL01 + KL01 * fd1
      b1  = KL01 * NH4
      TNH41, TNH42 = MatrixSolution(TNH41, TNH42, a11, a12_TNH4, b1, a21_TNH4, a22_TNH4, b2_TNH4)
      TNH41 = max(TNH41, 0.0)
      TNH42 = max(TNH42, 0.0)
      
      # NO31 and NO32
      a11 = -a21_NO3 + vno31_tc * vno31_tc / KL01 + KL01
      b1  = con_nit * FNH4 / KL01 * TNH41 + KL01 * NO3
      NO31, NO32 = MatrixSolution(NO31, NO32, a11, a12_NO3, b1, a21_NO3, a22_NO3, b2_NO3)
      NO31 = max(NO31, 0.0)
      NO32 = max(NO32, 0.0)
      
      JC_dn = rondn * (vno31_tc * vno31_tc / KL01 * NO31 + vno32_tc * NO32)
      JCc = max(roc * JC - JC_dn, 1.0E-10)
      
      # SO41, SO42 and TH2S1, TH2S2
      #coefficients for SO41, SO42, H2S1 and H2S2 equations

      bx=[0]*4
      ad=[(0,0,0,0),(0,0,0,0), (0,0,0,0), (0,0,0,0)] # 4x4
      
      #coefficients for SO41 and SO42 equations    
      g=[0]*2
      h=[(0,0),(0,0)] #2x2          
  
      # Half-saturation method
      if POCdiagenesis_part_option == 1 :
        #compute HSO4
        HSO4 = math.sqrt(2.0 * Dd_tc * SO4 * h2 / JCc) 
        if (HSO4 > h2):
          HSO4 = h2
        if HSO4 == 0.0:
          KL12SO4 = 1.0       # set a large value (1) for KL12SO4. 
        else:
          KL12SO4 = Dd_tc / (0.5 * HSO4)

        # four equations
        # 0 = bx(0) + ad(0,0) * SO41 + ad(0,1) * SO42 + ad(0,2) * TH2S1
        # 0 = bx(1) + ad(1,0) * SO41 + ad(1,1) * SO42 - JCc * SO42 / (SO42 + KsSO4)
        # 0 = bx(2) + ad(2,2) * TH2S1 + ad(2,3) * TH2S2
        # 0 = bx(3) + ad(3,2) * TH2S1 + ad(3,3) * TH2S2 + JCc * SO42 / (SO42 + KsSO4)

        bx[0] = KL01 * SO4                     
        bx[2] = KL01 * H2S                     
        ad[0][0] = - KL01 - KL12SO4
        ad[0][1] = KL12SO4
        ad[0][2] = con_sox / KL01
        ad[2][2] = - vb - fps1 * w12 - fds1 * KL01 - fds1 * KL12SO4 - con_sox / KL01
        ad[2][3] = w12 * fps2 + KL12SO4 * fds2
        
        ad[1][0] = KL12SO4
        ad[3][2] = vb + fps1 * w12 + fds1 * KL12SO4
        if SedFlux_solution_option == 1 :
          ad[1][1] = - KL12SO4
          bx[1]   = 0.0
          ad[3][3] = - vb - KL12SO4 * fds2 - w12 * fps2
          bx[3]   = 0.0
        else:
          ad[1][1] = - KL12SO4 - h2 / dt
          bx[1]   = h2 * SO42_prev / dt
          ad[3][3] = - vb - KL12SO4 * fds2 - w12 * fps2 - h2 / dt 
          bx[3]   = h2 * TH2S2_prev / dt

        #eliminate H2ST1 and H2ST2 from above equation sets, get two equations
        # 0 = g(0) + h(0,0) * SO41 + h(0,1) * SO42
        # 0 = g(1) + h(1,0) * SO41 + h(1,1) * SO42 + JCc * SO42 / (SO42 + KsSO4)
        
        g[1] = ((bx[0] * ad[2][2] - ad[0][2] * bx[2]) * ad[3][3] - bx[0] * ad[2][3] * ad[3][2]) / (ad[0][2] * ad[2][3]) + bx[3]   
        g[0] = g[1] + bx[1]  
        h[1][0] = ad[0][0] * (ad[2][2] * ad[3][3] - ad[2][3] * ad[3][2]) / (ad[0][2] * ad[2][3])
        h[1][1] = ad[0][1] * (ad[2][2] * ad[3][3] - ad[2][3] * ad[3][2]) / (ad[0][2] * ad[2][3])
        h[0][0] = h[1][0] + ad[1][0] 
        h[0][1] = h[1][1] + ad[1][1]
        
        #eliminate SO41 and get a quadratic equation of SO42
        # ra2 * SO42 * SO42 + ra1 * SO42 + ra0 = 0.0 
        ra0 = (h[0][0] * g[1] - h[1][0] * g[0]) * KsSO4
        ra1 = h[0][0] * g[1] - h[1][0] * g[0] + (h[0][0] * h[1][1] - h[0][1] * h[1][0]) * KsSO4 + h[0][0] * JCc
        ra2 = h[0][0] * h[1][1] - h[0][1] * h[1][0]
      
        #solve SO42
        sn1 = 1.0
        if (ra1 <= 0.0):
          sn1 = - 1.0
        disc = (- ra1 - sn1 * math.sqrt(ra1 * ra1 - 4.0 * ra2 * ra0)) / 2.0
        if (disc != 0.0) :
          r1 = disc / ra2
          r2 = ra0 / disc
        else:
          if (ra2 == 0.0) :
            r1 = - ra0 / ra1
            r2 = r1
          else:   #TODO what does this calculate? This is when ra0 = 0 
            r1 = -ra1 / ra2
  
        SO42 = r1
        if (SO42 < 0.0) :
          SO42 = r2
        
        #solve SO41, H2ST1, H2ST2
        SO41  = - (g[0] + h[0][1] * SO42) / h[0][0]
        TH2S1 = - (bx[0] + ad[0][0] * SO41 + ad[0][1] * SO42) / ad[0][2]
        TH2S2 = - (bx[2] + ad[2][2] * TH2S1) / ad[2][3]
        CSOD_H2S = con_sox / KL01 * TH2S1
        JCc_CH4  = JCc * KsSO4 / (SO42 + KsSO4)       #TODO find this formula. See page 170 is SO42 supposed to be on top as well?
        
      #sulfate reduction depth method
      elif POCdiagenesis_part_option == 2 :
        
        if (SO4 <= 0.1) :
          HSO4    = 0.0
          KL12SO4 = 1.0
          SO41    = SO4
          SO42    = SO4
        else:
          
          # first compute HSO4 TODO I do not know where the equation comes from: see page 159 seems to be close
          if SedFlux_solution_option == 1:
            ra2 = 2.0 * KL01 * JCc / h2
            ra1 = 2.0 * Dd_tc * JCc / h2
            ra0 = - 2.0 * Dd_tc * (KL01 * SO4 + con_sox / KL01 * TH2S1_prev)      
          elif SedFlux_solution_option == 2 :
            ra2 = 2.0 * KL01 * JCc / h2 + (KL01 * SO4 + con_sox / KL01 * TH2S1_prev) / dt
            ra1 = 2.0 * Dd_tc * JCc / h2 - 2.0 * KL01 * HSO4_prev * SO42_prev / dt
            ra0 = - 2.0 * Dd_tc * (KL01 * SO4 + con_sox / KL01 * TH2S1_prev + HSO4_prev * SO42_prev / dt)

          HSO4 = (- ra1 + math.sqrt(ra1 * ra1 - 4.0 * ra2 * ra0)) / 2.0 / ra2
          if (HSO4 > h2) :
            HSO4 = h2
          
          # TH2S1, TH2S2
          a11 = -a21_TH2S + con_sox / KL01 + KL01 * fds1
          b1  = KL01 * H2S
          if SedFlux_solution_option == 1 :  
            b2  = JCc * HSO4 / h2
          elif SedFlux_solution_option == 2:
            b2  = JCc * HSO4 / h2 + TH2S2_prev * h2 / dt

          TH2S1, TH2S2 = MatrixSolution (TH2S1, TH2S2, a11, a12_TH2S, b1, a21_TH2S, a22_TH2S, b2)
          TH2S1 = max(TH2S1, 0.0)
          TH2S2 = max(TH2S2, 0.0)

          CSOD_H2S = con_sox / KL01 * TH2S1
          JCc_CH4  = JCc * (h2 - HSO4) / h2     #TODO should it be JCc_CH4 = JCc * (1-(HSO4/H2)) page 170 in PDF
          
          # SO41, SO42
          if (HSO4 == h2) :
            KL12SO4 = KL12
            a11  = KL01 + KL12
            b1   = KL01 * SO4 + CSOD_H2S
            if SedFlux_solution_option == 1 : #TODO check thes b2
              b2 = JCc
            elif SedFlux_solution_option == 2 :
              b2 = JCc - SO42_prev * HSO4_prev / dt

            SO41, SO42 = MatrixSolution (SO41, SO42, a11, a12_SO4, b1, a21_SO4, a22_SO4, b2) 
            SO41 = max(SO41, 0.0)
            SO42 = max(SO42, 0.0)  

          else :
            KL12SO4 = Dd_tc / (0.5 * HSO4)
            SO41 = (KL01 * SO4 + CSOD_H2S) / (KL01 + KL12SO4 * 0.5)
            SO42 = SO41 / 2.0

        # TH2S1, TH2S2
        a11 = -a21_TH2S + con_sox / KL01 + KL01 * fds1
        b1  = KL01 * H2S
        if SedFlux_solution_option == 1 :
          b2  = JCc * HSO4 / h2
        elif SedFlux_solution_option == 2 :
          b2  = JCc * HSO4 / h2 + TH2S2_prev * h2 / dt

        TH2S1, TH2S2 = MatrixSolution(TH2S1, TH2S2, a11, a12_TH2S, b1, a21_TH2S, a22_TH2S, b2)
        TH2S1 = max(TH2S1, 0.0)
        TH2S2 = max(TH2S2, 0.0)
        CSOD_H2S = con_sox / KL01 * TH2S1
        JCc_CH4  = JCc * (h2 - HSO4) / h2

      # CH41 and CH42

      # analytical solutions
      if Methane_solution_option == 1 :
        CSODmax  = min(math.sqrt(2.0 * KL12 * CH4sat * JCc_CH4), JCc_CH4)
        
        if vch41_tc / KL01 < 100.0 :
          CSOD_CH4 = CSODmax * (1.0 - 2.0 / (math.exp(vch41_tc / KL01) + math.exp(-vch41_tc / KL01))) 
        else:
          CSOD_CH4 = CSODmax

      #numerical solutions
      elif Methane_solution_option == 2 : 
        a11 = KL12 + con_cox / KL01 + KL01
        b1  = KL01 * CH4
        if SedFlux_solution_option == 1 : 
          b2 = JCc_CH4
        elif SedFlux_solution_option == 2 :
          b2 = JCc_CH4 + CH42_prev * h2 / dt
    
        CH41, CH42 = MatrixSolution(CH41, CH42, a11, a12_CH4, b1, a21_CH4, a22_CH4, b2)
        
        if CH42 > CH4sat :
          CH42 = CH4sat
          CH41 = (b1 - a12_CH4 * CH42) / a11

        CH41     = max(CH41, 0.0)
        CH42     = max(CH42, 0.0)
        CSOD_CH4 = con_cox / KL01 * CH41

      #update SOD
      NSOD    = ron * con_nit * FNH4 / KL01 * TNH41     #TODO potentially missing a dn1 equation 5.37 page 162
      SOD_old = SOD_Bed
      SOD_Bed  = (CSOD_CH4 + CSOD_H2S + NSOD + SOD_Bed) / 2.0
      if (abs(SOD_Bed - SOD_old) / SOD_Bed * 100.0 < res) :
        break    
      KL01 = SOD_Bed / DOX
      if math.isnan(KL01) or KL01 == 0.0:
        KL01 = 1.0E-8

    #TODO not sure if i is defined outside of the loop
    #determine whether SOD is converged
    if i > int(maxit) :
      print('SOD iterations exceeded.')

    # inorganic species reactions and mass transfers
    kdpo41=0
    hsat=0    #depth where methane reaches saturation
    
    KL01 = SOD_Bed / DOX
      
    # pathways of TNH41/2, NO31/2, CH41/2, SO41/2, TH2S1/2, DIC1/2, TIP1/2
    #TNH41 and TNH42
    NH41= fd1 * TNH41
    NH42= fd2 * TNH42
    
    JNH4= KL01 * (NH41 - NH4)
    TNH41_Burial= vb * TNH41
    NH41_Nitrification = con_nit * FNH4 / KL01 * TNH41
    #PNH41_PNH42 = vb * (fp2 * TNH42 - fp1 * TNH41)
    NH41_NH42 = KL12 * (NH42 - NH41)
    TNH42_Burial = vb * TNH42
    
    # NO31 and NO32
    JNO3 = KL01 * (NO31 - NO3)
    NO31_Denit = vno31_tc * vno31_tc / KL01 * NO31
    NO31_NO32 = KL12 * (NO32 - NO31)
    NO32_Denit = vno32_tc * NO32
    
    # SO41 and SO42
    if POCdiagenesis_part_option == 1 :
      JCc_SO4 = JCc * SO42 / (SO42 + KsSO4)
    elif POCdiagenesis_part_option == 2 :
      JCc_SO4 = JCc * HSO4 / h2

    JSO4 = KL01 * (SO41 - SO4)
    SO41_SO42 = KL12SO4 * (SO42 - SO41)
    
    # TH2S1 and TH2S2
    H2S1 = fds1 * TH2S1
    H2S2 = fds2 * TH2S2
    
    JH2S = KL01 * (H2S1 - H2S)
    H2S1_Oxidation = con_sox / KL01 * TH2S1
    TH2S1_Burial = vb * TH2S1
    H2S1_H2S2 = KL12 * (H2S2 - H2S1)
    
    #PH2S1_PH2S2 = vb  * (TH2S2 * fps2 - TH2S1 * fps1)
    TH2S2_Burial = vb * TH2S2
    
    # CH41 and CH42
    if Methane_solution_option == 1 :
      CH41_Oxidation = CSOD_CH4
      JCH4 = CSODmax - CSOD_CH4
      JCH4g = JCc_CH4 - CSODmax
      if vch41_tc <= 0 :
        CH41 = 0.0
      else :
        CH41 = CSOD_CH4 / (vch41_tc * vch41_tc / KL01)

      #CH42 is not computed???
      CH42 = 0.0                                        
      
    elif Methane_solution_option == 2 :
      CH41_Oxidation = con_cox / KL01 * CH41
      JCH4 = KL01 * (CH41 - CH4)
      if CH42 == CH4sat :
        JCH4g = JCc_CH4 - JCH4 - CH41_Oxidation - (CH42 - CH42_prev) / dt * h2
      else:
        JCH4g = 0.0

    # DIC1 and DIC2
    a11 = KL01 + KL12
    a12 = -KL12
    b1 = KL01 * DIC * 12000.0 + CH41_Oxidation / 2.0 / roc + rcdn * NO31_Denit
    a21   = -KL12
    if SedFlux_solution_option == 1 :
      a22 = KL12
      b2  = (JCc_CH4 / 2.0 + JCc_SO4) / roc  + rcdn * NO32_Denit
    elif SedFlux_solution_option == 2 :
      a22 = KL12 + h2 / dt
      b2  = (JCc_CH4 / 2.0 + JCc_SO4) / roc + rcdn * NO32_Denit + DIC2 * h2 / dt

    DIC1, DIC2 = MatrixSolution(DIC1, DIC2, a11, a12, b1, a21, a22, b2)
    DIC1  = max(DIC1, 0.0)
    DIC2  = max(DIC2, 0.0)
    
    JDIC = KL01 * (DIC1 - DIC * 12000.0)
    DIC1_CH41_Oxidation = CH41_Oxidation / 2.0 / roc
    DIC1_NO31_Denit = rcdn * NO31_Denit
    DIC1_DIC2 = KL12 * (DIC2 - DIC1)
    DIC2_POC2_SO42 = JCc_SO4 / roc
    DIC2_CH42 = JCc_CH4 / 2.0 / roc
    DIC2_NO32_Denit = rcdn * NO32_Denit
    
    # TIP1 and TIP2
    if DOX >= DOcr : 
      kdpo41 = kdpo42 * d_kpo41
    else :
      kdpo41 = kdpo42 * d_kpo41**(DOX / DOcr)

    fd1 = 1.0 / (1.0 + Css1 * kdpo41)
    fd2 = 1.0 / (1.0 + Css2 * kdpo42)
    fp1 = 1.0 - fd1
    fp2 = 1.0 - fd2
    
    a21 = -w12 * fp1 - KL12 * fd1 - vb
    a11 = -a21 + KL01 * fd1
    a12 = -w12 * fp2 - KL12 * fd2
    b1  = KL01 * fdp * TIP

    if SedFlux_solution_option == 1 :
      a22 = -a12 + vb
      b2  = JP + TIP * (1.0 - fdp) * vs
    elif SedFlux_solution_option == 2 :
      a22 = -a12 + vb + h2 / dt
      b2  =  JP + TIP * (1.0 - fdp) * vs + h2 * TIP2 / dt

    TIP1, TIP2 = MatrixSolution(TIP1, TIP2, a11, a12, b1, a21, a22, b2)
    TIP1 = max(TIP1, 0.0)
    TIP2 = max(TIP2, 0.0)
    
    DIP1 = TIP1 * fd1
    DIP2 = TIP2 * fd2
    
    JDIP = KL01 * (DIP1 - fdp * TIP)
    TIP_TIP2 = TIP * (1.0 - fdp) * vs
    TIP1_Burial = vb * TIP1
    DIP1_DIP2 = KL12 * (DIP2 - DIP1)
    #PIP1_PIP2 = w12  * (TIP2 * fp2 - TIP1 * fp1)
    TIP2_Burial = vb * TIP2

  # solve mass balance equations by a matrix solution

  '''
# output sediment diagenesis pathways
  subroutine SedFluxPathwayOutput(na, a)
    integer  :: na
    real(R8) :: a(na)
    !
    do i = 1, 3
      if (JPOC_index(i) > 0)            a(JPOC_index(i))            = JPOC(i)
      if (JPON_index(i) > 0)            a(JPON_index(i))            = JPON(i)
      if (JPOP_index(i) > 0)            a(JPOP_index(i))            = JPOP(i)
      if (POC2_Diagenesis_index(i) > 0) a(POC2_Diagenesis_index(i)) = POC2_Diagenesis(i)
      if (PON2_Diagenesis_index(i) > 0) a(PON2_Diagenesis_index(i)) = PON2_Diagenesis(i)
      if (POP2_Diagenesis_index(i) > 0) a(POP2_Diagenesis_index(i)) = POP2_Diagenesis(i)
      if (POC2_Burial_index(i) > 0)     a(POC2_Burial_index(i))     = POC2_Burial(i)
      if (PON2_Burial_index(i) > 0)     a(PON2_Burial_index(i))     = PON2_Burial(i)
      if (POP2_Burial_index(i) > 0)     a(POP2_Burial_index(i))     = POP2_Burial(i)
    end do
    !
    if (JC_index > 0)    a(JC_index)    = JC
    if (JC_dn_index > 0) a(JC_dn_index) = JC_dn
    if (JN_index > 0)    a(JN_index)    = JN
    if (JP_index > 0)    a(JP_index)    = JP
    !
    if (w12_index > 0)      a(w12_index)      = w12
    if (KL12_index > 0)     a(KL12_index)     = KL12
    if (KL01_index > 0)     a(KL01_index)     = KL01
    if (SOD_Bed_index > 0)  a(SOD_Bed_index)  = SOD_Bed
    !
    if (JNH4_index > 0)                 a(JNH4_index)               = JNH4
    if (TNH41_Burial_index > 0)         a(TNH41_Burial_index)       = TNH41_Burial             
    if (NH41_Nitrification_index > 0)   a(NH41_Nitrification_index) = NH41_Nitrification
    !if (PNH41_PNH42_index > 0)          a(PNH41_PNH42_index)        = PNH41_PNH42
    if (NH41_NH42_index > 0)            a(NH41_NH42_index)          = NH41_NH42
    if (TNH42_Burial_index > 0)         a(TNH42_Burial_index)       = TNH42_Burial
    !
    if (JNO3_index > 0)                 a(JNO3_index)               = JNO3
    if (NO31_Denit_index > 0)           a(NO31_Denit_index)         = NO31_Denit
    if (NO31_NO32_index > 0)            a(NO31_NO32_index)          = NO31_NO32
    if (NO32_Denit_index > 0)           a(NO32_Denit_index)         = NO32_Denit
    !
    if (CH4sat_index > 0)               a(CH4sat_index)             = CH4sat
    if (JCH4_index > 0)                 a(JCH4_index)               = JCH4
    if (CH41_Oxidation_index > 0)       a(CH41_Oxidation_index)     = CH41_Oxidation
    if (JCc_CH4_index > 0)              a(JCc_CH4_index)            = JCc_CH4
    if (JCH4g_index > 0)                a(JCH4g_index)              = JCH4g
    !
    if (JSO4_index > 0)                 a(JSO4_index)               = JSO4
    if (JCc_SO4_index > 0)              a(JCc_SO4_index)            = JCc_SO4
    if (SO41_SO42_index > 0)            a(SO41_SO42_index)          = SO41_SO42
    !
    if (JH2S_index > 0)                 a(JH2S_index)               = JH2S
    if (H2S1_Oxidation_index > 0)       a(H2S1_Oxidation_index)     = H2S1_Oxidation
    if (TH2S1_Burial_index > 0)         a(TH2S1_Burial_index)       = TH2S1_Burial
    if (H2S1_H2S2_index > 0)            a(H2S1_H2S2_index)          = H2S1_H2S2
    !if (PH2S1_PH2S2_index > 0)          a(PH2S1_PH2S2_index)        = PH2S1_PH2S2
    if (TH2S2_Burial_index > 0)         a(TH2S2_Burial_index)       = TH2S2_Burial
    !
    if (JDIC_index > 0)                 a(JDIC_index)               = JDIC
    if (DIC1_CH41_Oxidation_index > 0)  a(DIC1_CH41_Oxidation_index)= DIC1_CH41_Oxidation
    if (DIC1_NO31_Denit_index > 0)      a(DIC1_NO31_Denit_index)    = DIC1_NO31_Denit
    if (DIC1_DIC2_index > 0)            a(DIC1_DIC2_index)          = DIC1_DIC2
    if (DIC2_POC2_SO42_index > 0)       a(DIC2_POC2_SO42_index)     = DIC2_POC2_SO42
    if (DIC2_CH42_index > 0)            a(DIC2_CH42_index)          = DIC2_CH42
    if (DIC2_NO32_Denit_index > 0)      a(DIC2_NO32_Denit_index)    = DIC2_NO32_Denit
    !
    if (JDIP_index > 0)                 a(JDIP_index)               = JDIP
    if (TIP_TIP2_index > 0)             a(TIP_TIP2_index)           = TIP_TIP2
    if (TIP1_Burial_index > 0)          a(TIP1_Burial_index)        = TIP1_Burial
    !if (PIP1_PIP2_index > 0)            a(PIP1_PIP2_index)          = PIP1_PIP2
    if (DIP1_DIP2_index > 0)            a(DIP1_DIP2_index)          = DIP1_DIP2
    if (TIP2_Burial_index > 0)          a(TIP2_Burial_index)        = TIP2_Burial
  ''' 

# compute sediment diagenesis derived variables

  #real(R8) :: kdpo41 TODO why is this here?
  TPOC2 = 0.0
  TPON2 = 0
  TPOP2 = 0.0
  for i in (1, 3) :
    TPOC2 = TPOC2 + POC2(i)
    TPON2 = TPON2 + PON2(i)
    TPOP2 = TPOP2 + POP2(i)

  POM2 = TPOC2 / focm2
  
  # TH2S1 and TH2S2
  H2S1 = TH2S1 / (1.0 + Css1 * kdh2s2)
  H2S2 = TH2S2 / (1.0 + Css2 * kdh2s2)
  
  # TIP1 and TIP2
  if DOX >= DOcr :
    kdpo41 = kdpo42 * d_kpo41
  else :
    kdpo41 = kdpo42 * d_kpo41**(DOX / DOcr)

  DIP1  = TIP1 / (1.0 + Css1 * kdpo41)
  DIP2  = TIP2 / (1.0 + Css2 * kdpo42)


def MatrixSolution(x1, x2, a11, a12, b1, a21, a22, b2) :
  x1=x1
  x2=x2
  a11=a11
  a12=a12
  b1=b1
  a21=a21
  a22=a22
  b2=b2

  if (a11 * a22 - a12 * a21 == 0.0) :
    print('Twod is singular: A11,A12,A21,A22')  
    print('a11, a12, a21, a22')

  x1 = (a22 * b1 - a12 * b2) / (a11 * a22 - a12 * a21)
  x2 = (a11 * b2 - a21 * b1) / (a11 * a22 - a12 * a21)
  return x1 , x2