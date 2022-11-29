
'''
=======================================================================================
Nutrient Simulation Module 1 (NSM1): Initial Values
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
Last Revision Date: June 12, 2021
'''


class WQ_Variable:
    def __init__(self, name: str, value: float, index: str, units: str):
        pass
'''    
        self.name = name
        self.value = value
        self.index = index
        self.units = units


class Globals:
    def __init__(self):

        self.initial_values: dict = {
            # 'pb': 1600.0,
            # 'Cps': 1673.0,
            # 'h2': 0.1,
            # 'alphas': 0.0432,
            # 'Richardson_option': True
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

        self.globals: dict = {
            # Global variables
            'nCBOD': 1,  # Number of CBOD groups
            'nGS': 1,   # Number of suspended solids

            # Sediment Diagenesis
            'sedFlux_name': 'Sediment Diagenesis',

            # Dependent hydraulic variables, (1D and 2D) (3)
            'depth': WQ_Variable('Depth', 0.0, 0, ''),
            'velocity': WQ_Variable('Velocity', 0.0, 0, ''),
            'shear_velocity': WQ_Variable('Velocity', 0.0, 0, ''),

            # Dependent hydraulic variables, specific for 1D model (3)
            'flow': WQ_Variable('Flow', 0.0, 0, ''),
            'slope': WQ_Variable('Slope', 0.0, 0, ''),
            'topwidth': WQ_Variable('Top Width', 0.0, 0, ''),

            # Dependent meteorological variables (3)
            'q_solar': WQ_Variable('Solar Radiation', 0.0, 0, ''),
            'wind_speed': WQ_Variable('Wind Speed', 0.0, 0, ''),
            'pressure_atm': WQ_Variable('Atmospheric Pressure', 0.0, 1, ''),

            # Dependent water quality variables (4)
            'TwaterC': WQ_Variable('Water Temperature', 0.0, 0, ''),
            'TsedC': WQ_Variable('Sediment Temperature', 0.0, 0, ''),
            'salinity': WQ_Variable('Salinity', 0.0, 0, ''),
            'Suspended_Solids_list': [],  # Define below

            # State variables and kinetic rates (16)
            'Ap': WQ_Variable('Ap', 0.0, 0, 'µg-Chla/L'),
            'dApdt': WQ_Variable('dApdt', 0.0, 0, '(µg-Chla/L)/day'),
            'Ab': WQ_Variable('Ap', 0.0, 0, 'g-Ab/m2'),
            'dAbdt': WQ_Variable('dAbdt', 0.0, 0, '(g-Ab/m2)/day'),
            'NH4': WQ_Variable('NH4', 0.0, 0, 'mg-N/L'),
            'dNH4dt': WQ_Variable('dNH4dt', 0.0, 0, '(mg-N/L)/day'),
            'NO3': WQ_Variable('NO3', 0.0, 0, 'mg-N/L'),
            'dNO3dt': WQ_Variable('dNO3dt', 0.0, 0, '(mg-N/L)/day'),
            'OrgN': WQ_Variable('OrgN', 0.0, 0, 'mg-N/L'),
            'dOrgNdt': WQ_Variable('dOrgNdt', 0.0, 0, '(mg-N/L)/day'),
            'TIP': WQ_Variable('TIP', 0.0, 0, 'mg-P/L'),
            'dTIPdt': WQ_Variable('dTIPdt', 0.0, 0, '(mg-P/L)/day'),
            'OrgP': WQ_Variable('OrgP', 0.0, 0, 'mg-P/L'),
            'dOrgPdt': WQ_Variable('dOrgPdt', 0.0, 0, '(mg-P/L)/day'),
            'POC': WQ_Variable('POC', 0.0, 0, 'mg-C/L'),
            'dPOCdt': WQ_Variable('dPOCdt', 0.0, 0, '(mg-C/L)/day'),
            'DOC': WQ_Variable('DOC', 0.0, 0, 'mg-C/L'),
            'dDOCdt': WQ_Variable('dDOCdt', 0.0, 0, '(mg-C/L)/day'),
            'DIC': WQ_Variable('DIC', 0.0, 0, 'mol-C/L'),
            'dDICdt': WQ_Variable('dDICdt', 0.0, 0, '(mol-C/L)/day'),
            'CBOD_list': [],  # Define below
            'dCBODdt_list': [],  # Define below
            'DOX': WQ_Variable('DOX', 0.0, 0, 'mg-O/L'),
            'dDOXdt': WQ_Variable('dDOXdt', 0.0, 0, '(mg-O/L)/day'),
            'N2': WQ_Variable('N2', 0.0, 0, 'mg-N/L'),
            'dN2dt': WQ_Variable('dN2dt', 0.0, 0, '(mg-N/L)/day'),
            'Alk': WQ_Variable('Alk', 0.0, 0, 'mg-CaCO3/L'),
            'dAlkdt': WQ_Variable('dAlkdt', 0.0, 0, '(mg-CaCO3/L)/day'),
            'PX': WQ_Variable('PX', 0.0, 0, 'cfu/100mL'),
            'dPXdt': WQ_Variable('dPXdt', 0.0, 0, '(cfu/100mL)/day'),
            'POM2': WQ_Variable('POM2', 0.0, 0, 'mg/L'),
            'dPOM2dt': WQ_Variable('dPOM2dt', 0.0, 0, '(mg/L)/day'),

            # Sediment diagenesis state variables (19)

            # Sediment layer 1
            'NH41': WQ_Variable('NH41', 0.0, 0, 'mg-N/L'),
            'dNH41dt': WQ_Variable('dNH41dt', 0.0, 0, '(mg-N/L)/day'),
            'NO31': WQ_Variable('NO31', 0.0, 0, 'mg-N/L'),
            'dNO31dt': WQ_Variable('dNO31dt', 0.0, 0, '(mg-N/L)/day'),
            'SO41': WQ_Variable('SO41', 0.0, 0, 'mg-O2/L'),
            'dSO41dt': WQ_Variable('dSO41dt', 0.0, 0, '(mg-O2/L)/day'),
            'TH2S1': WQ_Variable('TH2S1', 0.0, 0, 'mg-O2/L'),
            'dTH2S1dt': WQ_Variable('dTH2S1dt', 0.0, 0, '(mg-O2/L)/day'),
            'CH41': WQ_Variable('CH41', 0.0, 0, 'mg-O2/L'),
            'dCH41dt': WQ_Variable('dCH41dt', 0.0, 0, '(mg-O2/L)/day'),
            'DIC1': WQ_Variable('DIC1', 0.0, 0, 'mg-C/L'),
            'dDIC1dt': WQ_Variable('dDIC1dt', 0.0, 0, '(mg-C/L)/day'),
            'TIP1': WQ_Variable('TIP1', 0.0, 0, 'mg-P/L'),
            'dTIP1dt': WQ_Variable('dTIP1dt', 0.0, 0, '(mg-P/L)/day'),

            # Sediment layer 2
            'NH42': WQ_Variable('NH42', 0.0, 0, 'mg-N/L'),
            'dNH42dt': WQ_Variable('dNH42dt', 0.0, 0, '(mg-N/L)/day'),
            'NO32': WQ_Variable('NO32', 0.0, 0, 'mg-N/L'),
            'dNO32dt': WQ_Variable('dNO32dt', 0.0, 0, '(mg-N/L)/day'),
            'SO42': WQ_Variable('SO42', 0.0, 0, 'mg-O2/L'),
            'dSO42dt': WQ_Variable('dSO42dt', 0.0, 0, '(mg-O2/L)/day'),
            'TH2S2': WQ_Variable('TH2S2', 0.0, 0, 'mg-O2/L'),
            'dTH2S2': WQ_Variable('dTH2S2', 0.0, 0, '(mg-O2/L)/day'),
            'CH42': WQ_Variable('CH42', 0.0, 0, 'mg-O2/L'),
            'dCH42dt': WQ_Variable('dCH42dt', 0.0, 0, '(mg-O2/L)/day'),
            'DIC2': WQ_Variable('DIC2', 0.0, 0, 'mg-C/L'),
            'dDIC2dt': WQ_Variable('dDIC2dt', 0.0, 0, '(mg-C/L)/day'),
            'TIP2': WQ_Variable('TIP2', 0.0, 0, 'mg-P/L'),
            'dTIP2dt': WQ_Variable('dTIP2dt', 0.0, 0, '(mg-P/L)/day'),
            'ST': WQ_Variable('ST', 0.0, 0, 'day'),
            'dSTdt': WQ_Variable('dSTdt', 0.0, 0, ''),
            'HSO4': WQ_Variable('HSO4', 0.0, 0, 'm'),
            'dSTdt': WQ_Variable('dSTdt', 0.0, 0, 'm/day'),
            'POC2_list': [],  # Define below
            'dPOC2dt_list': [],  # Define below
            'PON2_list': [],  # Define below
            'dPON2dt_list': [],  # Define below
            'POP2_list': [],  # Define below
            'dPOP2dt_list': [],  # Define below
            # Derived variables (17 + 7)
            'Apd': WQ_Variable('Apd', 0.0, 0, ''),  # Algal biomass
            # Benthic Chlorophyll-a concentration
            'Chlb': WQ_Variable('Chlb', 0.0, 0, ''),
            # Dissolved Inorganic Nitrogen concentration
            'DIN': WQ_Variable('DIN', 0.0, 0, ''),
            # Total Organice Nitrogen concentration
            'TON': WQ_Variable('TON', 0.0, 0, ''),
            # Total Kjeldahl Nitrogen concentration
            'TKN': WQ_Variable('TKN', 0.0, 0, ''),
            'TN': WQ_Variable('TN', 0.0, 0, ''),  # Total Nitrogen concentration
            # Dissolved Inorganic Phosphorous concentration
            'DIP': WQ_Variable('DIP', 0.0, 0, ''),
            # Dissolved Organic Phosphorous concentration
            'TOP': WQ_Variable('TOP', 0.0, 0, ''),
            'TP': WQ_Variable('TP', 0.0, 0, ''),  # Total Phosphorous concentration
            # Total Organic Carbon concentration
            'TOC': WQ_Variable('TOC', 0.0, 0, ''),
            # Particulate Organic Matter concentration
            'POM': WQ_Variable('POM', 0.0, 0, 'mg/L'),
            # Total Suspended Solids concentration
            'TSS': WQ_Variable('TSS', 0.0, 0, ''),
            # 5-day CBOD concentration
            'CBOD5day': WQ_Variable('CBOD5day', 0.0, 0, ''),
            'TDG': WQ_Variable('TDG', 0.0, 0, ''),  # Total Dissolved Gases
            # Light attenuation coefficient
            'lambda': WQ_Variable('lambda', 0.0, 0, ''),
            'ka': WQ_Variable('ka', 0.0, 0, ''),  # O2 reaeration rate
            'pH': WQ_Variable('pH', 0.0, 0, ''),

            # Sediment layers
            'H2S1': WQ_Variable('H2S1', 0.0, 0, ''),
            'DIP1': WQ_Variable('DIP1', 0.0, 0, ''),
            'H2S2': WQ_Variable('H2S2', 0.0, 0, ''),
            'DIP2': WQ_Variable('DIP2', 0.0, 0, ''),
            # Note: the key and variable name differ
            'TPOC2': WQ_Variable('POC2', 0.0, 0, ''),
            # Note: the key and variable name differ
            'TPON2': WQ_Variable('PON2', 0.0, 0, ''),
            # Note: the key and variable name differ
            'TPOP2': WQ_Variable('POP2', 0.0, 0, ''),
        }

        for i in range(global_vals['nGs']):
            global_vals['Suspended_Solids_list'].append(
                WQ_Variable('Suspended Solid', 0.0, 0, ''))

        for i in range(global_vals['nCBOD']):
            global_vals['CBOD_list'].append(WQ_Variable('CBOD', 0.0, 0, 'mg-O/L'))
            global_vals['dCBODdt_list'].append(
                WQ_Variable('dCBODdt', 0.0, 0, '(mg-O/L)/day'))

        for i in range(3):
            # G1/G2/G3 -- what does this mean??
            global_vals['POC2_list'].append(
                WQ_Variable('Sediment_POC', 0.0, 0, 'mg-C/L'))
            global_vals['dPOC2dt_list'].append(
                WQ_Variable('dSediment_POCdt', 0.0, 0, '(mg-C/L)/day'))
            global_vals['PON2_list'].append(
                WQ_Variable('Sediment_PON', 0.0, 0, 'mg-N/L'))
            global_vals['dPON2dt_list'].append(
                WQ_Variable('dSediment_PONdt', 0.0, 0, '(mg-N/L)/day'))
            global_vals['POP2_list'].append(
                WQ_Variable('Sediment_POP', 0.0, 0, 'mg-P/L'))
            global_vals['dPOP2dt_list'].append(
                WQ_Variable('dSediment_POPdt', 0.0, 0, '(mg-P/L)/day'))

        # WQ groups (parameters and pathways)
        self.wq_groups: dict = {
            'groupGlobal': 'NSMI Global',
            'groupAlgae': 'Algae',
            'groupBenthicAlgae': 'Benthic Algae',
            'groupNitrogen': 'Nitrogen Cycle',
            'groupPhosphorus': 'Phosphorus Cycle',
            'groupCarbon': 'Carbon Cycle',
            'groupCBOD': 'Carbonaceous BOD',
            'groupDO': 'Dissolved Oxygen',
            'groupN2': 'Nitrogen Gas',
            'groupAlkalinity': 'Alkalinity',
            'groupPathogen': 'Pathogen',
            'groupPOM': 'Particulate Organic Matter',
            'groupSedFlux': 'Sediment Diagenesis'
        }

        # Parameter questions
        self.parameter_questions: dict = {
            'paramO2HydroReaeration': 'Flow O2 reaeration',
            'paramO2WindReaeration': 'Wind O2 reaeration',
            'paramPHSolution': 'pH solution',
            'paramApGrowth': 'Ap growth',
            'paramApLight': 'Ap light',
            'paramAbGrowth': 'Ab growth',
            'paramAbLight': 'Ab light',
            'paramSedFlux': 'SedFlux solution',
            'paramCH4Solution': 'CH4 solution',
            'paramPOCDiagPart': 'POC diagenesis partition'
        }

        # Parameter names
        self.parameter_names: list = [

            # Light extinction
            'lambda0', 'lambdas', 'lambdam', 'lambda1', 'lambda2',

            # Solid partition and settling velocity
            'kdpo4', 'vs', 'vson', 'vsop', 'vsoc',

            # POM stoichiomety
            'focm',

            # SOD and ka
            'SOD', 'KsSod', 'kah', 'kaw',

            # Algae
            'AWd', 'AWc', 'AWn', 'AWp', 'AWa', 'mu_max', 'krp',
            'kdp', 'vsap', 'KL', 'KsN', 'KsP', 'PN', 'Fpocp',

            # Benthic algae
            'BWd', 'BWc', 'BWn', 'BWp', 'BWa', 'mub_max', 'krb',
            'kdb', 'KLb', 'KsNb', 'KsPb', 'KSb', 'PNb', 'Fpocb',
            'Fb_name', 'Fw_name',

            # Nitrogen
            'kon', 'knit', 'rnh4', 'kdnit', 'KsOxdn', 'vno3',

            # Phosphorus
            'kop', 'rpo4',

            # Carbon
            'kpoc', 'kdoc', 'KsOxmc', 'Fco2', 'pco2',

            # CBOD (3 * nCBOD)
            'kbod', 'KsOxbod', 'ksbod',

            # DO - no parameter, merged to global group

            # Pathogen
            'kdx', 'apx', 'vx',

            # POM group
            'kpom2', 'h2', 'vb',

            # SedFlux group

            # General parameters
            'Dd', 'Dp', 'KsDp', 'POCr', 'SO4_fresh',
            'kst', 'focm2', 'res', 'maxit',
            # 'h2', 'vb',

            # Sediment layer 1
            'vnh41', 'vno31', 'vch41', 'vh2sd', 'vh2sp', 'KsOxna1',
            'KsNh4', 'KsOxch', 'KSh2s', 'DOcr', 'd_kpo41', 'Css1',

            # Sediment layer 2
            'vno32', 'kdnh42', 'kdh2s2', 'kdpo42', 'KsSO4', 'Css2',
            'FPON1', 'FPON2', 'FPOP1', 'FPOP2', 'FPOC1', 'FPOC2',
            'FAP1', 'FAP2', 'FAB1', 'FAB2', 'FCBOD1', 'FCBOD2',
            'KPONG1', 'KPONG2', 'KPOPG1', 'KPOPG2', 'KPOCG1', 'KPOCG2',
        ]

        # Pathways
        self.pathways: dict = {

            # Algae
            'Ap_growth': 'ApGrowth',
            'Ap_respiration': 'ApRespiration',
            'Ap_death': 'ApDeath',
            'Ap_settling': 'ApSettling',
            'FL_name': 'FL',
            'FN_name': 'FN',
            'FP_name': 'FP',

            # Benthic algae
            'Ab_growth': 'AbGrowth',
            'Ab_respiration': 'AbRespiration',
            'Ab_death': 'AbDeath',
            'FLb_name': 'FLb',
            'FNb_name': 'FNb',
            'FPb_name': 'FPb',
            'FSb_name': 'FSb',

            # Nitrogen
            'Ap_OrgN': 'ApDeath_OrgN',
            'Ap_NH4': 'NH4_ApRespiration',
            'NH4_Ap': 'NH4_ApGrowth',
            'NO3_Ap': 'NO3_ApGrowth',

            'OrgN_Bed': 'OrgN_Settling',
            'OrgN_NH4': 'OrgN_NH4_Decay',
            'Bed_NH4': 'NH4fromBed',
            'NH4_NO3': 'NH4_Nitrification',
            'NO3_Denit_name': 'NO3_Denit',
            'NO3_Bed': 'NO3_BedDenit',

            'Ab_OrgN': 'AbDeath_OrgN',
            'Ab_NH4': 'NH4_AbRespiration',
            'NH4_Ab': 'NH4_AbGrowth',
            'NO3_Ab': 'NO3_AbGrowth',

            # Phosphorous
            'Ap_OrgP': 'ApDeath_OrgP',
            'Ap_DIP': 'DIP_ApRespiration',
            'DIP_Ap': 'DIP_ApGrowth',

            'OrgP_Bed': 'OrgP_Settling',
            'OrgP_DIP': 'OrgP_DIP_Decay',
            'TIP_Bed': 'TIP_Settling',
            'Bed_DIP': 'DIPfromBed',

            'Ab_OrgP': 'AbDeath_OrgP',
            'Ab_DIP': 'DIP_AbRespiration',
            'DIP_Ab': 'DIP_AbGrowth',

            # Carbon
            'Ap_POC': 'ApDeath_POC',
            'Ap_DOC': 'ApDeath_DOC',
            'Ap_DIC': 'ApRespiration_DIC',
            'DIC_Ap': 'DIC_ApGrowth',

            'POC_Bed': 'POC_Settling',
            'POC_DOC': 'POC_DOC_Hydrolysis',
            'DOC_DIC': 'DOC_DIC_Oxidation',
            'CBOD_DIC': 'CBOD_DIC_Oxidation',
            'Atm_DIC': 'DIC_Reaeration',
            'Bed_DIC': 'DICfromBed',

            'Ab_POC': 'AbDeath_POC',
            'Ab_DOC': 'AbDeath_DOC',
            'Ab_DIC': 'AbRespiration_DIC',
            'DIC_Ab': 'DIC_AbGrowth',

            # CBOD (2 * nCBOD)
            'CBOD_Oxidation_name': 'CBOD_Oxidation',
            'CBOD_Sediment_name': 'CBOD_Sediment',

            # Oxygen
            'Atm_O2': 'O2_Reaeration',
            'O2sat_name': 'O2sat',
            'Ap_O2': 'O2_ApGrowth',
            'O2_Ap': 'O2_ApRespiration',
            'Ab_O2': 'O2_AbGrowth',
            'O2_Ab': 'O2_AbRespiration',
            'O2_Nitrification_name': 'O2_Nitrification',
            'O2_DOC': 'O2_DOC_Oxidation',
            'O2_CBOD': 'O2_CBOD_Oxidation',
            'O2_SOD_name': 'O2_SOD',

            # N2
            'Atm_N2': 'N2_Reaeration',
            'N2_sat': 'N2sat',

            # Alkalinity
            'Alk_Ap': 'Alk_ApGrowth',
            'Ap_Alk': 'Alk_ApRespiration',
            'Alk_Ab': 'Alk_AbGrowth',
            'Ab_Alk': 'Alk_AbRespiration',
            'Alk_Nitrification_name': 'Alk_Nitrification',
            'Alk_Denit_name': 'Alk_Denit',

            # Pathogen
            'PX_death': 'PathogenDeath',
            'PX_decay': 'PathogenDecay',
            'PX_settling': 'PathogenSettling',

            # POM
            'Ap_POM2': 'ApSettling_POM2',
            'Ab_POM2': 'AbDeath_POM2',
            'POC_POM2': 'POCSettling_POM2',
            'POM2_Dissolution_name': 'POM2_Dissolution',
            'POM2_Burial_name': 'POM2_Burial',

            # Sediment diagenesis
            'JPOC_name': 'JPOC',
            'JPON_name': 'JPON',
            'JPOP_name': 'JPOP',
            'POC2_Diagenesis_name': 'POC2_Diagenesis',
            'PON2_Diagenesis_name': 'PON2_Diagenesis',
            'POP2_Diagenesis_name': 'POP2_Diagenesis',
            'POC2_Burial_name': 'POC2_Burial',
            'PON2_Burial_name': 'PON2_Burial',
            'POP2_Burial_name': 'POP2_Burial',
            'JC_name': 'JC',
            'JC_dn_name': 'JC_dn',
            'JN_name': 'JN',
            'JP_name': 'JP',

            'w12_name': 'w12',
            'KL12_name': 'KL12',
            'KL01_name': 'KL01',
            'SOD_Bed_name': 'SOD_Bed',

            'NH41_NH4': 'JNH4',
            'TNH41_Burial_name': 'TNH41_Burial',
            'NH41_NO31': 'NH41_Nitrification',
            'NH41_NH42_name': 'NH41_NH42',
            # 'PNH41_PNH42_name'  : 'PNH41_PNH42',
            'TNH42_Burial_name': 'TNH42_Burial',

            'NO31_NO3': 'JNO3',
            'NO31_Denit_name': 'NO31_Denit',
            'NO31_NO32_name': 'NO31_NO32',
            'NO32_Denit_name': 'NO32_Denit',

            'POC2_CH42': 'JCc_CH4',
            'CH4sat_name': 'CH4sat',
            'CH41_CH4': 'JCH4',
            'CH41_CSOD': 'CH41_Oxidation',
            'CH42_Gas': 'JCH4g',

            'SO41_SO4': 'JSO4',
            'JCc_SO4_name': 'JCc_SO4',
            'SO41_SO42_name': 'SO41_SO42',

            'H2S1_H2S': 'JH2S',
            'H2S1_CSOD': 'H2S1_Oxidation',
            'TH2S1_Burial_name': 'TH2S1_Burial',
            'H2S1_H2S2_name': 'H2S1_H2S2',
            # 'PH2S1_PH2S2_name'  : 'PH2S1_PH2S2',
            'TH2S3_Burial_name': 'TH2S2_Burial',

            'DIC1_DIC': 'JDIC',
            'CH41_DIC1': 'DIC1_CH41_Oxidation',
            'NO31_DIC1': 'DIC1_NO31_Denit',
            'DIC1_DIC2_name': 'DIC1_DIC2',
            'POC2_DIC2': 'DIC2_POC2_SO42',
            'CH42_DIC2_name': 'DIC2_CH42',
            'NO32_DIC2': 'DIC2_NO32_Denit',

            'DIP1_DIP': 'JDIP',
            'TIP_TIP2_name': 'TIP_TIP2',
            'TIP1_Burial_name': 'TIP1_Burial',
            'DIP1_DIP2_name': 'DIP1_DIP2',
            # 'PIP1_PIP2_name' : 'PIP1_PIP2',
            'TIP2_Burial_name': 'TIP2_Burial',
        }
'''