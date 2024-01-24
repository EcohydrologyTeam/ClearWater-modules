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
    light_limitation_option = 1
)


class BalgaeStaticVariables(TypedDict):
    BWd: float
    BWc: float
    Bwn: float
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
    Bwn=7.2,
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


class N2StaticVariables(TypedDict):
    pass                  

DEFAULT_N2 = N2StaticVariables(
)


class PhosphorusStaticVariables(TypedDict):
    kop_20: float
    rpo4_20: float

DEFAULT_PHOSPHORUS = PhosphorusStaticVariables(
    kop_20 = 0.1,
    rpo4_20 = 0
)


class POMStaticVariables(TypedDict):
    kpom_20: float

DEFAULT_POM = POMStaticVariables(
    kpom_20 = 0.1
)


class CBODStaticVariables(TypedDict):
    kbod_20: float
    ksbod_20: float
    ksOxbod: float

DEFAULT_CBOD = CBODStaticVariables(
    kbod_20 = 0.12,
    ksbod_20 = 0,
    ksOxbod = 0.5
)


class CarbonStaticVariables(TypedDict):
    F_pocp: float
    kdoc_20: float
    F_pocb: float
    kpoc_20: float
    K_sOxmc: float
    pCO2: float
    FCO2: float
    
DEFAULT_CARBON = CarbonStaticVariables(
    F_pocp = 0.9,
    kdoc_20 = 0.01,
    F_pocb = 0.9,
    kpoc_20 = 0.005,
    K_sOxmc = 1,
    pCO2 = 383,
    FCO2 = 0.2
)


class DOXStaticVariables(TypedDict):
    ...

DEFAULT_DOX = DOXStaticVariables(
    
)


class PathogenStaticVariables(TypedDict):
    kdx: float
    apx: float
    vx: float                  

DEFAULT_PATHOGEN = PathogenStaticVariables(
    kdx=0.8,
    apx=1,
    vx=1
)

class AlkalinityStaticVariables(TypedDict):
    r_alkaa: float
    r_alkan: float
    r_alkn: float
    r_alkden: float
    r_alkba: float
    r_alkbn: float
    
DEFAULT_ALKALINITY = AlkalinityStaticVariables(
    r_alkaa = 1,
    r_alkan = 1,
    r_alkn = 1,
    r_alkden = 1,
    r_alkba = 1,
    r_alkbn = 1 
)

class GlobalParameters(TypedDict):
    use_NH4 : bool
    use_NO3 : bool  
    use_OrgN: bool                 
    use_TIP : bool   
    use_SedFlux: bool
    use_DOX: bool
    use_Algae: bool
    use_Balgae: bool
    use_TIP: bool 
 
DEFAULT_GLOBALPARAMETERS = GlobalParameters(
    use_NH4= True,
    use_NO3= True, 
    use_OrgN= True,
    use_TIP= True,  
    use_SedFlux= True,
    use_DOX= True,
    use_Algae= True,
    use_Balgae= True,
    use_TIP= True
)

class GlobalVars(TypedDict):
    L : float #lambda
    fdp: float
    PAR: float
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

DEFAULT_GLOBALVARS = GlobalVars(
    L = 1, #lambda
    fdp= 0.5,
    PAR= 1, 
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
    hydraulic_reaeration_option = 999,
    wind_reaeration_option = 999
)

