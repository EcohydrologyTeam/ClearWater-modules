"""Constants for the NSM1. The TypedDicts allow for updating upon module init"""
from typing import (
    TypedDict,
)


class Algae(TypedDict):
    AWd: float                          #Algal Dry [mg]
    AWc: float                          #Carbon [mg]      
    AWn: float                          #Nitrogen [mg]
    AWp: float                          #Phosphorus [mg]
    AWa: float                          #Algal Chlorophyl [mg]

    KL: float                           #Light liminting constant for algal growth [W/m^2]
    KsN: float                          #Half saturation N limiting constant for algal growth [mg-N/L]
    KsP: float                          #Half saturation P limiting constant for algal growth [mg-P/L]
    mu_max: float                       #Max algal growh [1/d]
    kdp: float                          #Algal mortality rate [1/d]
    krp: float                          #Algal respiration rate [1/d]
    vsap: float                         #Algal settling velocity [m/d]
    growth_rate_option: int             #Algal growth rate option 1) multiplicative, 2) Limiting Nutrient, 3) Harmonic Mean Option
    light_limitation_option: int        #Algal light limitation 1) half saturation, 2) Smith model, 3) Steele model

DEFAULT_Algae = Algae(
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

class global_parameters(TypedDict):
    use_NH4 : bool
    use_NO3 : bool                   
    use_TIP : bool    
 
DEFAULT_global_parameters = global_parameters(
    use_NH4= True,
    use_NO3= True,                   
    use_TIP= True  
)

class global_vars(TypedDict):
    depth: float
    L : float #lambda
    fdp: float
    PAR: float                     

DEFAULT_global_vars = global_vars(
    depth= 1,
    L = 1, #lambda
    fdp= 0.5,
    PAR= 1, 
)
