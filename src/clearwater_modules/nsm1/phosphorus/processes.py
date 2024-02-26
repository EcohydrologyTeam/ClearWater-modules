"""
File contains process to calculate new phosphorus concentration and associated dependent variables
"""
import numba
import xarray as xr
from clearwater_modules.shared.processes import arrhenius_correction
import math


@numba.njit
def fdp(
    use_TIP: bool,
    Solid : xr.DataArray,
    kdop4: xr.DataArray
) -> xr.DataArray :

    """Calculate kop_tc: Decay rate of organic P to DIP temperature correction (1/d).

    Args:
        use_TIP: true/false use total inorganic phosphrous,
        Solid : #TODO define this
        kdop4: solid partitioning coeff. of PO4 (L/kg)
    """
  
    return xr.where(use_TIP, 1/(1+kdop4 * Solid/0.000001), 0)

@numba.njit
def kop_tc(
    TwaterC : xr.DataArray,
    kop_20: xr.DataArray
) -> xr.DataArray :

    """Calculate kop_tc: Decay rate of organic P to DIP temperature correction (1/d).

    Args:
        TwaterC: Water temperature (C)
        kop_20: Decay rate of organic P to DIP at 20C (1/d)
    """

    return arrhenius_correction(TwaterC, kop_20, 1.047)

@numba.njit
def rpo4_tc(
    TwaterC : xr.DataArray,
    rpo4_20: xr.DataArray
) -> xr.DataArray :

    """Calculate rpo4_tc: Benthic sediment release rate of DIP temperature correction(g-P/m2/d).

    Args:
        TwaterC: Water temperature (C)
        kop_20: Benthic sediment release rate of DIP at 20C (1/d)
    """

    return arrhenius_correction(TwaterC, rpo4_20, 1.074)

@numba.njit
def OrgP_DIP_decay(
    kop_tc : xr.DataArray,
    OrgP: xr.DataArray,
    use_OrgP: bool,
) -> xr.DataArray :
    
    """Calculate OrgP_DIP: organic phosphorus decay to dissolve inorganic phosphorus (mg-P/L/d).

    Args:
        kop_tc: Decay rate of organic P to DIP temperature correction (1/d)
        OrgP: Organic phosphorus concentration (mg-P/L)
        use_OrgP: true/false use organic phosphorus (t/f)
    """        
    return xr.where(use_OrgP,kop_tc * OrgP,0)

@numba.njit
def OrgP_Settling(
    vsop : xr.DataArray,
    depth: xr.DataArray,
    OrgP: xr.DataArray,
) -> xr.DataArray :
    
    """Calculate OrgP_Settling: organic phosphorus settling to sediment (mg-P/L/d).

    Args:
        vsop: Organic phosphorus settling velocity (m/d)
        depth: water depth (m)
        OrgP: Organic phosphorus concentration (mg-P/L)
    """        
    return (vsop / depth) * OrgP

@numba.njit
def ApDeath_OrgP(
    rpa : xr.DataArray,
    ApDeath: xr.DataArray,
    use_Algae: bool,
) -> xr.DataArray :
    
    """Calculate ApDeath_OrgP: Algal death turning into organic phosphorus  (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApDeath: Algal death rate (ug-Chla/L/d)
        use_Algae: true/false to use algae module (T/F)

    """        

    return xr.where(use_Algae, rpa * ApDeath,0)

@numba.njit
def AbDeath_OrgP(
    rpb : xr.DataArray,
    AbDeath: xr.DataArray,
    Fw: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: bool
) -> xr.DataArray :
    
    """Calculate AbDeath_OrgP: Benthic algal death turning into organic phosphorus (mg-P/L/d).

    Args:
        rpb : Benthic algal P: Benthic algal dry (mg-P/mg-D)
        AbDeath: Benthic algal death rate (g/m^2/d)
        Fw: Fraction benthic algal death to water column (unitless)
        Fb: Fraction bottom area avalible for benthic algae (unitless)
        depth: water depth (m)
        use_Balgae: true/false use benthic algae module (t/f)

    """        

    return xr.where(use_Balgae, (rpb * Fw *Fb * AbDeath) / depth,0)      

@numba.njit
def dOrgPdt(
    ApDeath_OrgP : xr.DataArray,
    AbDeath_OrgP: xr.DataArray,
    OrgP_DIP_decay: xr.DataArray,
    OrgP_Settling: xr.DataArray,
    use_OrgP: bool,
) -> xr.DataArray :
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

    return xr.where(use_OrgP, -OrgP_DIP_decay-OrgP_Settling + ApDeath_OrgP + AbDeath_OrgP, 0)
    
#TODO will this be a problem if use_SedFlux is False
@numba.njit
def DIPfromBed_SedFlux(
    use_SedFlux: bool,
    JDIP: xr.DataArray,
    depth:xr.DataArray,
    rpo4_tc: xr.DataArray,
) -> xr.DataArray :
    """Calculate DIPfromBed_SedFlux: Dissolved Organic Phosphorus coming from Bed calculated using SedFlux modules (mg-P/L/d).

    Args:
        use_SedFlux: true/false to use the sediment flux module (unitless)
        JDIP: Sediment-water flux of phosphate (g-P/m^2/d)
        depth: water depth (m)
        rpo4_tc: Benthic sediment release rate of DIP temperature correction(g-P/m2/d)
    """    
    return xr.where(use_SedFlux, JDIP / depth, rpo4_tc/depth)

#TODO calcuate fdp?
@numba.njit
def TIP_Settling(
    vs: xr.DataArray,
    depth: xr.DataArray,
    fdp: xr.DataArray,
    TIP: xr.DataArray
) -> xr.DataArray :

    """Calculate TIP_Settling: Total inorganic phosphorus settling from water to bed (mg-P/L/d).

    Args:
        vs: Sediment settling velocity (m/d)
        depth: water depth (m)
        fdp: Fraction phosphorus dissolved (unitless)
        TIP: Total inorganic phosphorus water concentration (mg-P/L)
    """             
    return vs / depth * (1.0 - fdp) * TIP

@numba.njit
def DIP_ApRespiration(
    rpa: xr.DataArray,
    ApRespiration: xr.DataArray,
    use_Algae: bool

) -> xr.DataArray :
    """Calculate DIP_ApRespiration: Dissolved inorganic phosphorus released from algal respiration (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApRespiration: Algal respiration rate (ug-Chla/L/d)
        use_Algae: true/false to use algae module (t/f)
    """ 
    return xr.where(use_Algae, rpa * ApRespiration,0)

@numba.njit
def DIP_ApGrowth(
    rpa: xr.DataArray,
    ApGrowth: xr.DataArray,
    use_Algae: bool

) -> xr.DataArray :
    """Calculate DIP_ApGrowth: Dissolved inorganic phosphorus consumed for algal growth (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApGrowth: Algal growth rate (ug-Chla/L/d)
        use_Algae: true/false to use algae module (t/f)
    """ 
    return xr.where(use_Algae, rpa * ApGrowth,0)

@numba.njit
def DIP_AbRespiration(
    rpb: xr.DataArray,
    AbRespiration: xr.DataArray,
    use_Balgae: bool

) -> xr.DataArray :
    """Calculate DIP_AbRespiration: Dissolved inorganic phosphorus released for benthic algal respiration (mg-P/L/d).

    Args:
        rpb: Benthic algal P : Benthic algal dry ratio (mg-P/mg-D)
        AbRespiration: Benthic algal respiration rate (g/m^2/d)
        use_Blgae: true/false to use benthic algae module (t/f)        
    """     
    return xr.where(use_Balgae, rpb * AbRespiration,0)

@numba.njit
def DIP_AbGrowth(
    rpb: xr.DataArray,
    AbGrowth: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: bool

) -> xr.DataArray :
    """Calculate DIP_AbGrowth: Dissolved inorganic phosphorus consumed for benthic algal growth (mg-P/L/d).

    Args:
        rpb: Benthic algal P : Benthic algal dry ratio (mg-P/mg-D)
        AbGrowth: Benthic algal growth rate (g/m^2/d)
        Fb: Fraction of bottom area available for benthic algal (unitless)
        depth: water depth (m)
        use_Balgae: true/false to use benthic algae module (t/f) 
    """     
    return xr.where(use_Balgae, rpb * Fb * AbGrowth / depth,0)

@numba.njit
def dTIPdt(
    OrgP_DIP_decay: xr.DataArray,
    TIP_Settling: xr.DataArray,
    DIPfromBed_SedFlux: xr.DataArray,
    DIP_ApRespiration: xr.DataArray,
    DIP_ApGrowth: xr.DataArray,
    DIP_AbRespiration: xr.DataArray,
    DIP_AbGrowth: xr.DataArray,
    use_TIP: bool, 

) -> xr.DataArray :    
    
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


    dTIP/dt =     OrgP Decay                (OrgP -> DIP)
                - DIP AlgalUptake           (DIP -> xr.DataArraying Algae)
                - DIP BenthicAlgae Uptake   (DIP -> xr.DataArraying Algae)	
                - TIP Settling              (TIP -> bed)
                + DIP From Benthos          (Benthos -> DIP) 
    """

    return xr.where(use_TIP, - TIP_Settling + DIPfromBed_SedFlux + OrgP_DIP_decay + DIP_ApRespiration - DIP_ApGrowth + DIP_AbRespiration - DIP_AbGrowth, 0)


@numba.njit
def TIP(
    TIP: xr.DataArray,
    dTIPdt: xr.DataArray,
    timestep: xr.DataArray

) -> xr.DataArray :
    """Calculate TIP: New total inorganic phosphorus (mg-P/L).

    Args:
        dTIPdt: Change in total inorganic phosphorus (mg-P/L/d)
        TIP: Total inorganic phosphorus water concentration (mg-P/L),
        timestep: current iteration timestep (d)
    """     
    return TIP+dTIPdt*timestep

@numba.njit
def OrgP(
    OrgP: xr.DataArray,
    dOrgPdt: xr.DataArray,
    timestep: xr.DataArray

) -> xr.DataArray :
    """Calculate OrgP: New total organic phosphorus (mg-P/L).

    Args:
        dOrgPdt: Change in total organic phosphorus (mg-P/L/d)
        OrgP: Total organic phosphorus water concentration (mg-P/L),
        timestep: current iteration timestep (d)
    """     
    return OrgP+dOrgPdt*timestep

@numba.njit
def TOP(
    use_OrgP: bool,
    OrgP: xr.DataArray,
    use_Algae: bool,
    rpa: xr.DataArray,
    Ap: xr.DataArray

) -> xr.DataArray :
    """Calculate TOP: Total organic phosphorus (mg-P/L).

    Args:
        use_OrgP: true/false to use organic phosphorus module (true/false),
        OrgP: New organic phosphorus water concentration (mg-P/L),
        use_Algae: true/false to use algae module (true/false),
        rpa: Algal P: Chla ratio (mg-P/ug-Chla),
        Ap: Algal water concentration (ug-Chla/L)
    """     
    TOP = 0.0
    TOP = xr.where(use_OrgP, TOP + OrgP,TOP)
    TOP = xr.where(use_Algae, TOP + rpa*Ap, TOP)

    return TOP


@numba.njit
def TP(
    use_TIP: bool,
    TOP: xr.DataArray,
    TIP: xr.DataArray

) -> xr.DataArray :
    """Calculate TP: Total phosphorus (mg-P/L).

    Args:
        use_TIP: true/false to use total inorganic phosphorus module (true/false),
        TIP: New total inorganic phosphorus water concentration (mg-P/L),
        TOP: Total organic phosphorus water concentration (mg-P/L)
    """  
    TP = TOP
    TP = xr.where(use_TIP,TP + TIP,TP)

@numba.njit
def DIP(
    fdp: xr.DataArray,
    TIP: xr.DataArray

) -> xr.DataArray :
    """Calculate DIP: Dissolve inorganich phosphorus (mg-P/L).

    Args:
        fdp: fraction P dissolved
        TIP: New total inorganic phosphorus water concentration (mg-P/L),
    """
    return TIP * fdp

