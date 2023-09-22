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
    mu_max: float                       
    kdp: float                          
    krp: float                          
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
    mu_max= 1,
    kdp= 0.15,
    krp= 0.2,
    vsap= 0.15,
    growth_rate_option = 1,
    light_limitation_option = 1
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
    mub_max: float
    krb: float
    kdb: float                         
    b_growth_rate_option: float
    b_light_limitation_option: float
    Fw: float
    Fb: float

DEFAULT_BALGAE = BalgaeStaticVariables(
    BWd = 100
    BWc= 40
    Bwn=7.2
    BWp= 1
    BWa= 3500

    KLb= 10
    KsNb= 0.25
    KsPb=0.125
    Ksb=10
    mub_max=0.4
    krb=0.2
    kdb=0.3
    b_growth_rate_option=1
    b_light_limitation_option=1
    Fw=0.9
    Fb=0.9 
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
    use_TIP= True, 
)

class GlobalVars(TypedDict):
    L : float #lambda
    fdp: float
    PAR: float
    vson: float                     

DEFAULT_GLOBALVARS = GlobalVars(
    L = 1, #lambda
    fdp= 0.5,
    PAR= 1, 
    vson = 0.01,
)
