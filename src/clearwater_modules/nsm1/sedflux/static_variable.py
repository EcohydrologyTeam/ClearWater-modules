"""
File contains static variables related to the SedFlux module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...




Variable(
    name='Awd',
    long_name='Algal Dry Weight',
    units='mg',
    description='Algal Dry Weight',
    use='static',
)

"""
       self.sedFlux_constants = OrderedDict()
        self.sedFlux_constants = {
            'Dd': 0.0025,
            'Dp': 0.00006,
            'vnh41': 0.1313,
            'vno31': 0.1,
            'vch41': 0.7,
            'vh2sd': 0.2,
            'vh2sp': 0.4,
            'KPONG1': 0.035,
            'KPONG2': 0.0018,
            'KPOPG1': 0.035,
            'KPOPG2': 0.0018,
            'KPOCG1': 0.035,
            'KPOCG2': 0.0018,
            'vno32': 0.25,

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

        roc = 32.0 / 12.0
        TempBen = 10
        roso4 = 2.0 * 32.0 / 98.0
        rondn = 32.0 / 12.0 * 10.0 / 8.0 * 12.0 / 14.0
        H2S = 0.0
        CH4 = 0.0
        ron = 2.0 * 32.0 / 14.0
        rcdn = 5.0 * 12.0 / (4.0 * 14.0)
"""