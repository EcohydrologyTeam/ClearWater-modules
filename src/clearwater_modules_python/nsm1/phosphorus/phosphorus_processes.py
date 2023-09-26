"""
File contains process to calculate new phosphorus concentration and associated dependent variables
"""
"""
(Global) module_choices T/F
'use_Algae'
'use_BAlgae'
'use_OrgP'
'use_TIP'
'use_SedFlux'

(Global) global_vars
'Ap'       Algae concentration                              [ug/Chla/L]
'TwaterC'  Water temperature                                [C]
'depth'    Depth from water surface                         [m]
'OrgP'     Organic phosphorus                               [mg-P/L]
'TIP'      Total inorganic phosphorus                       [mg-P/L]
'vs'       Sediment settling velocity                       [m/d]
'fdp'      fraction P dissolved                             [unitless]

phosphorus_constant_changes
'kop'       Decay rate or orgnaic P to DIP                  [1/d]
'rpo4'      Benthic sediment release rate of DIP            [g-P/m2*d]

from Algae
'rPa'           AlgalP : Chla ratio                              [mg-P/ugChla]
'ApGrowth'      Algal growth rate                                [ug-chla/L/d]
'ApDeath'       Algal death rate                                 [ug-chla/L/d]
'ApRespiration' AlgalRespiration rate                            [ug-chla/L/d]

from Benthic Algae
'rpb'           Benthic Algal P: Benthic Algal Dry Weight        [mg-P/mg-D]
'AbGrowth'      Benthic Algal growth rate                        [g/m^2*d]
'AbDeath'       Benthic Algal death rate                         [g/m^2*d]
'AbRespiration' Benthic Algal respiration rate                   [g/m^2*d]

from SedFlux
'JDIP'     Sediment water flux of phosphate                   [g-P/m^2*d]      

        Organic Phosphorus  (mgP/day)
    
        dOrgP/dt =    Algae_OrgP              (Algae -> OrgP)
                    - OrgP Decay              (OrgP -> DIP)	
                    - OrgP Settling           (OrgP -> bed) 
                    + BenthicAlgalDeath       (Benthic Algae -> OrgP)    	
        
"""
import math
from clearwater_modules_python.shared.processes import arrhenius_correction
import numba

@numba.njit
def kop_tc(
    TwaterC : float,
    kop_20: float
) -> float :

    """Calculate kop_tc: Decay rate of organic P to DIP temperature correction (1/d).

    Args:
        TwaterC: Water temperature (C)
        kop_20: Decay rate of organic P to DIP at 20C (1/d)
    """

    return arrhenius_correction(TwaterC, kop_20, 1.047)

@numba.njit
def rop4_tc(
    TwaterC : float,
    rop4_20: float
) -> float :

    """Calculate rop4_tc: Benthic sediment release rate of DIP temperature correction(g-P/m2/d).

    Args:
        TwaterC: Water temperature (C)
        kop_20: Benthic sediment release rate of DIP at 20C (1/d)
    """

    return arrhenius_correction(TwaterC, rop4_20, 1.074)

@numba.njit
def OrgP_DIP_decay(
    kop_tc : float,
    OrgP: float,
) -> float :
    
    """Calculate OrgP_DIP: organic phosphorus decay to dissolve inorganic phosphorus (mg-P/L/d).

    Args:
        kop_tc: Decay rate of organic P to DIP temperature correction (1/d)
        OrgP: Organic phosphorus concentration (mg-P/L)
    """        
    return kop_tc * OrgP

@numba.njit
def OrgP_Settling(
    vsop : float,
    depth: float,
    OrgP: float,
) -> float :
    
    """Calculate OrgP_Settling: organic phosphorus settling to sediment (mg-P/L/d).

    Args:
        vsop: Organic phosphorus settling velocity (m/d)
        depth: water depth (m)
        OrgP: Organic phosphorus concentration (mg-P/L)
    """        
    return (vsop / depth) * OrgP

@numba.njit
def ApDeath_OrgP(
    rpa : float,
    ApDeath: float,
) -> float :
    
    """Calculate ApDeath_OrgP: Algal death turning into organic phosphorus  (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApDeath: Algal death rate (ug-Chla/L/d)

    """        

    return rpa * ApDeath

@numba.njit
def AbDeath_OrgP(
    rpb : float,
    AbDeath: float,
    Fw: float,
    Fb: float,
    depth: float
) -> float :
    
    """Calculate AbDeath_OrgP: Benthic algal death turning into organic phosphorus (mg-P/L/d).

    Args:
        rpb : Benthic algal P: Benthic algal dry (mg-P/mg-D)
        AbDeath: Benthic algal death rate (g/m^2/d)
        Fw: Fraction benthic algal death to water column (unitless)
        Fb: Fraction bottom area avalible for benthic algae (unitless)
        depth: water depth (m)

    """        

    return (rpb * Fw *Fb * AbDeath) / depth      

@numba.njit
def dOrgPdt(
    ApDeath_OrgP : float,
    AbDeath_OrgP: float,
    OrgP_DIP_decay: float,
    OrgP_Settling: float,
    use_OrgP: bool,
    use_Algae: bool,
    use_Balgae: bool,
) -> float :
    """Calculate dOrgPdt: change in organic phosphorus concentration (mg-P/L/d).

    Args:
        ApDeath_OrgP: Algal death turns into organic phosphrous 
        AbDeath_OrgP: Benthic algal death turns into organic phosphrous
        OrgP_DIP_decay: Organic phosphrous decaying into dissolve inorganic phosphrous
        OrgP_Settling: Organic phosphrous settling into sediment
        use_OrgP: true/false to use organic phosphorus module (true/false)
        use_Algae: true/false to use algae module (true/false)
        use_Balgae: true/false to use benthic algae module (true/false)
    """     
    dOrgPdt=0
    if use_OrgP:
        dOrgPdt=dOrgPdt-OrgP_DIP_decay-OrgP_Settling

        if use_Algae:
            dOrgPdt=dOrgPdt+ApDeath_OrgP
        if use_Balgae:
            dOrgPdt=dOrgPdt+AbDeath_OrgP

    return dOrgPdt
    

@numba.njit
def DIPfromBed_SedFlux(
    JDIP: float,
    depth:float
) -> float :
    """Calculate DIPfromBed_SedFlux: Dissolved Organic Phosphorus coming from Bed calculated using SedFlux modules (mg-P/L/d).

    Args:
        JDIP: Sediment-water flux of phosphate (g-P/m^2/d)
        depth: water depth (m)
    """    
    return JDIP / depth

@numba.njit
def DIPfromBed_NoSedFlux(
    rpo4_tc: float,
    depth:float
) -> float :
    """Calculate DIPfromBed_NoSedFlux: Dissolved Organic Phosphorus coming from Bed calculated without SedFlux modules (mg-P/L/d).

    Args:
        rpo4_tc: Benthic sediment release rate of DIP with temperature correction (g-P/m^2/d)
        depth: water depth (m)
    """             
    return rpo4_tc / depth

@numba.njit
def TIP_Settling(
    vs: float,
    depth: float,
    fdp: float,
    TIP: float
) -> float :

    """Calculate TIP_Settling: Total inorganic phosphorus settling from water to bed (mg-P/L/d).

    Args:
        vs: Sediment settling velocity (m/d)
        depth: water depth (m)
        fdp: Fraction phosphorus dissolved (unitless)
        TIP: Total inorganic phosphorus water concentration (mg-P/L)
    """             
    return vs / depth * (1.0 - fdp) * TIP

@numba.njit
def OrgP_DIP_decay(
    kop_tc: float,
    OrgP: float,

) -> float :
    """Calculate OrgP_DIP_decay: Total organic phosphorus decaying to dissolved inorganic phosphrous (mg-P/L/d).

    Args:
        kop_tc: Decay rate of organic phosphorus to dissolved inorganic phosphorus with temperature correction (1/d)
        OrgP: Total organic phosphorus water concentration (mg-P/L)
    """  
    return kop_tc * OrgP

@numba.njit
def DIP_ApRespiration(
    rpa: float,
    ApRespiration: float,

) -> float :
    """Calculate DIP_ApRespiration: Dissolved inorganic phosphorus released from algal respiration (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApRespiration: Algal respiration rate (ug-Chla/L/d)
    """ 
    return rpa * ApRespiration

@numba.njit
def DIP_ApGrowth(
    rpa: float,
    ApGrowth: float,

) -> float :
    """Calculate DIP_ApGrowth: Dissolved inorganic phosphorus consumed for algal growth (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApGrowth: Algal growth rate (ug-Chla/L/d)
    """ 
    return rpa * ApGrowth

@numba.njit
def DIP_AbRespiration(
    rpb: float,
    AbRespiration: float,

) -> float :
    """Calculate DIP_AbRespiration: Dissolved inorganic phosphorus released for benthic algal respiration (mg-P/L/d).

    Args:
        rpb: Benthic algal P : Benthic algal dry ratio (mg-P/mg-D)
        AbRespiration: Benthic algal respiration rate (g/m^2/d)
    """     
    return rpb * AbRespiration

@numba.njit
def DIP_AbGrowth(
    rpb: float,
    AbGrowth: float,
    Fb: float,
    depth: float

) -> float :
    """Calculate DIP_AbGrowth: Dissolved inorganic phosphorus consumed for benthic algal growth (mg-P/L/d).

    Args:
        rpb: Benthic algal P : Benthic algal dry ratio (mg-P/mg-D)
        AbGrowth: Benthic algal growth rate (g/m^2/d)
        Fb: Fraction of bottom area available for benthic algal (unitless)
        depth: water depth (m)
    """     
    return rpb * Fb * AbGrowth / depth

@numba.njit
def dTIPdt(
    OrgP_DIP_decay: float,
    TIP_Settling: float,
    DIPfromBed_NoSedFlux: float,
    DIPfromBed_SedFlux: float,
    DIP_ApRespiration: float,
    DIP_ApGrowth: float,
    DIP_AbRespiration: float,
    DIP_AbGrowth: float,
    use_TIP: bool, 
    use_SedFlux: bool,
    use_OrgP: bool,
    use_Algae: bool,
    use_Balgae: bool,
) -> float :    
    
    """Calculate dTIPdt: Change in dissolved inorganic phosphorus water concentration (mg-P/L/d).

    Args:
        OrgP_DIP_decay: Total organic phosphorus decaying to dissolved inorganic phosphrous (mg-P/L/d),
        TIP_Settling: Total inorganic phosphorus settling from water to bed (mg-P/L/d),
        DIPfromBed_NoSedFlux: Dissolved Organic Phosphorus coming from Bed calculated without SedFlux modules (mg-P/L/d),
        DIPfromBed_SedFlux: Dissolved Organic Phosphorus coming from Bed calculated using SedFlux modules (mg-P/L/d),
        DIP_ApRespiration: Dissolved inorganic phosphorus released from algal respiration (mg-P/L/d),
        DIP_ApGrowth: Dissolved inorganic phosphorus consumed for algal growth (mg-P/L/d),
        DIP_AbRespiration: Dissolved inorganic phosphorus released for benthic algal respiration (mg-P/L/d),
        DIP_AbGrowth: Dissolved inorganic phosphorus consumed for benthic algal growth (mg-P/L/d),
        use_TIP: true/false to use total inorganic phosphorus module (true/false), 
        use_SedFlux: true/false to use sediment flux module (true/false)
        use_OrgP: true/false to use organic phosphorus module (true/false),
        use_Algae: true/false to use algae module (true/false),
        use_Balgae: true/false to use benthic algae module (true/false),

    dTIP/dt =     OrgP Decay                (OrgP -> DIP)
                - DIP AlgalUptake           (DIP -> Floating Algae)
                - DIP BenthicAlgae Uptake   (DIP -> Floating Algae)	
                - TIP Settling              (TIP -> bed)
                + DIP From Benthos          (Benthos -> DIP) 
    """
    dTIPdt = 0
    if use_TIP:
        dTIPdt = dTIPdt - TIP_Settling
        if use_SedFlux:
            dTIPdt = dTIPdt + DIPfromBed_SedFlux
        else: 
            dTIPdt = dTIPdt + DIPfromBed_NoSedFlux
        if use_OrgP:
            dTIPdt = dTIPdt + OrgP_DIP_decay 
        if use_Algae:
            dTIPdt = dTIPdt + DIP_ApRespiration - DIP_ApGrowth
        if use_Balgae:
            dTIPdt = dTIPdt + DIP_AbRespiration - DIP_AbGrowth 

    return dTIPdt

@numba.njit
def TIP_new(
    TIP: float,
    dTIPdt: float,

) -> float :
    """Calculate TIP_new: New total inorganic phosphorus (mg-P/L).

    Args:
        dTIPdt: Change in total inorganic phosphorus (mg-P/L/d)
        TIP: Total inorganic phosphorus water concentration (mg-P/L),

    """     
    return TIP+dTIPdt
@numba.njit
def OrgP_new(
    OrgP: float,
    dOrgPdt: float,

) -> float :
    """Calculate OrgP_new: New total organic phosphorus (mg-P/L).

    Args:
        dOrgPdt: Change in total organic phosphorus (mg-P/L/d)
        OrgP: Total organic phosphorus water concentration (mg-P/L),

    """     
    return OrgP+dOrgPdt

@numba.njit
def TOP(
    use_OrgP: bool,
    OrgP_new: float,
    use_Algae: bool,
    rpa: float,
    Ap: float

) -> float :
    """Calculate TOP: Total organic phosphorus (mg-P/L).

    Args:
        use_OrgP: true/false to use organic phosphorus module (true/false),
        OrgP_new: New organic phosphorus water concentration (mg-P/L),
        use_Algae: true/false to use algae module (true/false),
        rpa: Algal P: Chla ratio (mg-P/ug-Chla),
        Ap: Algal water concentration (ug-Chla/L)
    """     
    TOP = 0.0
    if use_OrgP:
        TOP = TOP + OrgP_new
    if use_Algae:
        TOP = TOP + rpa * Ap
    return TOP


@numba.njit
def TP(
    use_TIP: bool,
    TOP: float,
    TIP_new: float

) -> float :
    """Calculate TP: Total phosphorus (mg-P/L).

    Args:
        use_TIP: true/false to use total inorganic phosphorus module (true/false),
        TIP_new: New total inorganic phosphorus water concentration (mg-P/L),
        TOP: Total organic phosphorus water concentration (mg-P/L)
    """  
    TP = TOP
    if use_TIP:
        TP = TP + TIP_new

@numba.njit
def DIP(
    fdp: float,
    TIP_new: float

) -> float :
    """Calculate DIP: Dissolve inorganich phosphorus (mg-P/L).

    Args:
        fdp: fraction P dissolved
        TIP_new: New total inorganic phosphorus water concentration (mg-P/L),
    """
    return TIP_new * fdp

