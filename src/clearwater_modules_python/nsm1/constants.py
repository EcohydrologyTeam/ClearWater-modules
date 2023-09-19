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

class GlobalParameters(TypedDict):
    use_NH4 : bool
    use_NO3 : bool                   
    use_TIP : bool    
 
DEFAULT_GLOBALPARAMETERS = GlobalParameters(
    use_NH4= True,
    use_NO3= True,                   
    use_TIP= True  
)

class GlobalVars(TypedDict):
    depth: float
    L : float #lambda
    fdp: float
    PAR: float                     

DEFAULT_GLOBALVARS = GlobalVars(
    depth= 1,
    L = 1, #lambda
    fdp= 0.5,
    PAR= 1, 
)
