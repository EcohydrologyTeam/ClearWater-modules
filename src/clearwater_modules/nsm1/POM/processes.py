import numba
import math
from clearwater_modules.shared.processes import (
    arrhenius_correction
)
import xarray as xr
from clearwater_modules.nsm1.POM import dynamic_variables
from clearwater_modules.nsm1.POM import static_variables
from clearwater_modules.nsm1 import static_variables_global
from clearwater_modules.nsm1 import dynamic_variables_global
from clearwater_modules.nsm1 import state_variables

@numba.njit
def kpom_tc(
    TwaterC: float,
    kpom_20: float,
) -> float:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        TwaterC: Water temperature in Celsius
        kpom_20: POM dissolution rate at 20 degrees Celsius (1/d)
    """
    return arrhenius_correction(TwaterC, kpom_20, 1.047)


def POM_algal_settling(
    Ap: xr.DataArray,
    vsap: xr.DataArray,
    rda: xr.DataArray,
    depth: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculates the particulate organic matter concentration change due to algal mortality
    
    Args:
        Ap: Algae concentration (mg/L)
        vsap: Algal settling velocity (m/d)
        rda: Ratio of algal biomass to chlorophyll-a 
        depth: Depth of water in computation cell (m)
        use_Algae: Option to consider algal kinetics  
    """
    da: xr.DataArray = xr.where(use_Algae == True, vsap * Ap * rda / depth, 0)

    return da


@numba.njit
def POM_dissolution(
    POM: xr.DataArray,
    kpom_tc: xr.DataArray
) -> xr.DataArray:
    """Calculates the particulate organic matter concentration change due to POM dissolution

    Args:
        POM: Concentration of particulate organic matter (mg/L)
        kpom_tc: POM dissolution rate corrected for temperature (1/d)
    """

    return POM * kpom_tc


def POM_POC_settling(
    POC: xr.DataArray,
    vsoc: xr.DataArray,
    depth: xr.DataArray,
    fcom: xr.DataArray,
    use_POC: xr.DataArray
) -> xr.DataArray:
    """Calculates particulate organic matter concentration change due to POM settling
    
    Args:
        POC: Concentration of particulate organic carbon (mg/L)
        vsoc: POC settling velocity (m/d)
        depth: Depth of water (m)
        fcom: Fraction of carbon in organic matter (mg-C/mg-D) 
        use_POC: Option to consider particulate organic carbon
    """
    da: xr.DataArray = xr.where(use_POC == True, vsoc * POC / depth / fcom, 0)
    
    return da


def POM_benthic_algae_mortality(
    Ab: xr.DataArray,
    kdb_tc: xr.DataArray,
    Fb: xr.DataArray,
    Fw: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculates particulate organic matter concentration change due to benthic algae mortality
    
    Args:
        Ab: Benthic algae concentration (mg/L)
        kdb_tc: Benthic algae death rate (1/d)
        Fb: Fraction of bottom area available for benthic algae growth
        Fw: Fraction of benthic algae mortality into water column
        depth: Depth of water in computation cell (m)
        use_Balgae: Option for considering benthic algae in DOC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, Ab * kdb_tc * Fb * (1 - Fw) / depth, 0)

    return da


@numba.njit
def POM_burial(
    vb: xr.DataArray,
    POM: xr.DataArray,
    depth: xr.DataArray
) -> xr.DataArray:
    """Calculates particulate organic matter concentration change due to POM burial in the sediments
    
    Args:
        vb: Velocity of burial (m/d)
        POM: POM concentration (mg/L)
        depth: Depth of water in computation cell (m)
    """
    return vb * POM / depth


@numba.njit
def dPOMdt(
    POM_algal_settling: xr.DataArray,
    POM_dissolution: xr.DataArray,
    POM_POC_settling: xr.DataArray,
    POM_benthic_algae_mortality: xr.DataArray,
    POM_burial: xr.DataArray,
) -> xr.DataArray:
    """Calculates the concentration change of POM for one timestep

    Args:
        POM_algal_settling: POM concentration change due to algal settling (mg/L/d)
        POM_dissolution: POM concentration change due to dissolution (mg/L/d)
        POM_POC_settling: POM concentration change due to POC settling (mg/L/d)
        POM_benthic_algae_mortality: POM concentration change due to benthic algae mortality (mg/L/d)
        POM_burial: POM concentration change due to burial (mg/L/d)
    """
    return POM_algal_settling - POM_dissolution + POM_POC_settling + POM_benthic_algae_mortality - POM_burial


@numba.njit
def POM(
    dPOMdt: xr.DataArray,
    POM: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Computes updated particulate organic matter concentration (mg/L)
    
    Args:
        dPOMdt: Change in POM concentration over timestep (mg/L/d)
        POM: POM concentration from previous timestep (mg/L)
        timestep: Current iteration timestep (d)
    """
    return POM + dPOMdt * timestep


    
    

