"""Constants for the NSM1. The TypedDicts allow for updating upon module init"""
from typing import (
    TypedDict,
)


class AlgaeStaticVariables(TypedDict):
    AWd: float                          
    AWc: float                               
    AWn: float                          
    AWp: float                          
    AWa: float                          
    KL: float                           
    KsN: float                          
    KsP: float                          
    mu_max_20: float                       
    kdp_20: float                          
    krp_20: float                          
    vsap: float                         
    growth_rate_option: int             
    light_limitation_option: int 

DEFAULT_ALGAE = AlgaeStaticVariables(
    AWd = 100,
    AWc= 40,
    AWn= 7.2,
    AWp= 1,
    AWa= 1000,
    KL= 10,
    KsN= 0.04,
    KsP= 0.0012,
    mu_max_20= 1,
    kdp_20= 0.15,
    krp_20= 0.2,
    vsap= 0.15,
    growth_rate_option = 1,
    light_limitation_option = 1,
)

class AlkalinityStaticVariables(TypedDict):
    r_alkaa: float
    r_alkan: float
    r_alkn: float
    r_alkden: float
    r_alkba: float
    r_alkbn: float
      
DEFAULT_ALKALINITY = AlkalinityStaticVariables(
    r_alkaa = 14.0 / 106.0 / 12.0 / 1000.0,
    r_alkan= 18.0 / 106.0 / 12.0 / 1000.0,
    r_alkn = 2.0 / 14.0 / 1000.0,
    r_alkden = 4.0 / 14.0 / 1000.0,
    r_alkba = 14.0 / 106.0 / 12.0 / 1000.0,
    r_alkbn =18.0 / 106.0 / 12.0 / 1000.0
   
)

class BalgaeStaticVariables(TypedDict):
    BWd: float
    BWc: float
    BWn: float
    BWp: float
    BWa: float

    KLb: float
    KsNb: float
    KsPb: float
    Ksb: float
    mub_max_20: float
    krb_20: float
    kdb_20: float                         
    b_growth_rate_option: float
    b_light_limitation_option: float
    Fw: float
    Fb: float

DEFAULT_BALGAE = BalgaeStaticVariables(
    BWd = 100,
    BWc= 40,
    BWn=7.2,
    BWp= 1,
    BWa= 3500,

    KLb= 10,
    KsNb= 0.25,
    KsPb=0.125,
    Ksb=10,
    mub_max_20=0.4,
    krb_20=0.2,
    kdb_20=0.3,
    b_growth_rate_option=1,
    b_light_limitation_option=1,
    Fw=0.9,
    Fb=0.9
)


class NitrogenStaticVariables(TypedDict):
    KNR: float
    knit_20: float
    kon_20: float
    kdnit_20: float
    rnh4_20: float
    vno3_20: float
    KsOxdn: float
    PN: float
    PNb: float

DEFAULT_NITROGEN = NitrogenStaticVariables(
    KNR= 0.6 ,
    knit_20= 0.1,
    kon_20=0.1,
    kdnit_20=0.002,
    rnh4_20=0,
    vno3_20=0,
    KsOxdn=0.1,
    PN=0.5,
    PNb=0.5
)


class CarbonStaticVariables(TypedDict):
    f_pocp: float
    kdoc_20: float
    f_pocb: float
    kpoc_20: float
    KsOxmc: float
    pCO2: float
    FCO2: float
    roc: float

DEFAULT_CARBON = CarbonStaticVariables(
    f_pocp = 0.9,
    kdoc_20= 0.01,
    f_pocb=0.9,
    kpoc_20= 0.005,
    KsOxmc=1.0,
    pCO2 = 383.0,
    FCO2 = 0.2,
    roc = 32.0/12.0
)

class CBODStaticVariables(TypedDict):
    KsOxbod: float
    kbod_20: float
    ksbod_20: float

DEFAULT_CBOD = CBODStaticVariables(
    KsOxbod = 0.5,
    kbod_20 =  0.12,
    ksbod_20 = 0.0
)

class DOXStaticVariables(TypedDict):
    ron : float
    KsSOD : float

DEFAULT_DOX = DOXStaticVariables(
    ron = 2.0 * 32.0 / 14.0,
    KsSOD =1,
)

class N2StaticVariables(TypedDict):
    pass

DEFAULT_N2 = N2StaticVariables(

)


class POMStaticVariables(TypedDict):
    kpom_20: float

DEFAULT_POM = POMStaticVariables(
    kpom_20 = 0.1
)


class PathogenStaticVariables(TypedDict):
    kdx_20: float
    apx: float
    vx: float                  

DEFAULT_PATHOGEN = PathogenStaticVariables(
    kdx_20=0.8,
    apx=1,
    vx=1
)

class PhosphorusStaticVariables(TypedDict):
    kop_20: float
    rpo4_20: float
    kdpo4: float

DEFAULT_PHOSPHORUS = PhosphorusStaticVariables(
    kop_20 = 0.1,
    rpo4_20 =0,
    kdpo4 = 0.0,
)


class GlobalParameters(TypedDict):
    use_NH4 : bool
    use_NO3 : bool  
    use_OrgN: bool
    use_OrgP: bool                 
    use_TIP : bool   
    use_SedFlux: bool
    use_POC: bool
    use_DOC: bool
    use_DOX: bool
    use_DIC: bool
    use_Algae: bool
    use_Balgae: bool
    use_N2: bool
    use_Pathogen: bool
    use_Alk: bool
    use_POM: bool

 
DEFAULT_GLOBALPARAMETERS = GlobalParameters(
    use_NH4= True,
    use_NO3= True, 
    use_OrgN= True,
    use_OrgP = True,
    use_TIP= True,  
    use_SedFlux= False,
    use_POC = True,
    use_DOC = True,
    use_DOX= True,
    use_DIC= True,
    use_Algae= True,
    use_Balgae= True,
    use_N2 = True,
    use_Pathogen = True,
    use_Alk = True,
    use_POM = True         
)

class GlobalVars(TypedDict):
    vson: float
    vsoc: float
    vsop: float
    vs: float
    SOD_20: float
    SOD_theta: float
    vb: float
    fcom: float
    kaw_20_user: float
    kah_20_user: float
    hydraulic_reaeration_option: int
    wind_reaeration_option: int
    timestep: float
    depth: float
    TwaterC: float
    theta: float
    velocity: float
    flow: float
    topwidth: float
    slope: float
    shear_velocity: float
    pressure_atm: float
    wind_speed: float
    q_solar: float
    Solid: int
    lambda0: float
    lambda1: float
    lambda2: float
    lambdas: float
    lambdam: float
    Fr_PAR: float        

              

DEFAULT_GLOBALVARS = GlobalVars(
    vson = 0.01,
    vsoc = 0.01,
    vsop = 999,
    vs = 999,
    SOD_20 = 999,
    SOD_theta = 999,
    vb = 0.01,
    fcom = 0.4,
    kaw_20_user = 999,
    kah_20_user = 999,
    hydraulic_reaeration_option = 2,
    wind_reaeration_option = 2,  
    timestep = 86400,    #TODO Dynamic or static?
    depth = 1.5,         #TODO Dynamic or static?
    TwaterC = 20,
    theta = 1.047,
    velocity = 1,
    flow = 2,
    topwidth = 1,
    slope = 2,
    shear_velocity = 4,
    pressure_atm = 2,
    wind_speed = 4,
    q_solar = 500,
    Solid = 1,
    lambda0 = .02,
    lambda1 = .0088,
    lambda2 = .054,
    lambdas = .052,
    lambdam = .0174, 
    Fr_PAR = .47   
)
