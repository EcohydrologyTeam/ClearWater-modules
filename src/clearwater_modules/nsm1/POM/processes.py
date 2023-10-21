import numpy as np
import numba
import math
from clearwater_modules.shared.processes import (
    arrhenius_correction
)
import xarray as xr
from clearwater_modules.nsm1.alkalinity import dynamic_variables
from clearwater_modules.nsm1.alkalinity import static_variables
from clearwater_modules.nsm1 import static_variables_global
from clearwater_modules.nsm1 import dynamic_variables_global
from clearwater_modules.nsm1 import state_variables


@numba.njit
def kdp_T(
    water_temp_c: float,
    kdp_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kdp_20: Algae mortality rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient for kdp
    """
    return arrhenius_correction(water_temp_c, kdp_20, theta)


@numba.njit
def kdb_T(
    water_temp_c: float,
    kdb_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kdb_20: Benthic algae mortality rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient for kdb
    """
    return arrhenius_correction(water_temp_c, kdb_20, theta)


@numba.njit
def kpom_T(
    water_temp_c: float,
    kpom_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kpom_20: POM dissolution rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient for kpom
    """
    return arrhenius_correction(water_temp_c, kpom_20, theta)

def POM_algal_mortality(
    Ap: xr.DataArray,
    kdp_T: xr.DataArray,
    rda: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculates the particulate organic matter concentration change due to algal mortality
    
    Args:
        Ap: 
        kdp_T:
        rda:
        use_Algae:  
    """
    da: xr.DataArray = xr.where(use_Algae == True, kdp_T * Ap * rda, 0)

    return da


@numba.njit
def POM_dissolution(
    POM: xr.DataArray,
    kpom_T: xr.DataArray
) -> xr.DataArray:
    """Calculates the particulate organic matter concentration change due to POM dissolution

    Args:
        POM: Concentration of particulate organic matter (mg/L)
        kpom_T: POM dissolution rate corrected for temperature (1/d)
    """

    return POM * kpom_T


@numba.njit
def POM_settling(
    POM: xr.DataArray,
    vsom: xr.DataArray,
    depth: xr.DataArray
) -> xr.DataArray:
    """Calculates particulate organic matter concentration change due to POM settling
    
    Args:
        POM: Concentration of particulate organic matter (mg/L)
        vsom: POM settling velocity (m/d)
        depth: Depth of water (m)
    """

    return vsom / depth * POM


def POM_benthic_algae_mortality(
    Ab: xr.DataArray,
    kdb_T: xr.DataArray,
    Fb: xr.DataArray,
    Fw: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculates particulate organic matter concentration change due to benthic algae mortality
    
    Args:
        Ab: Benthic algae concentration (mg/L)
        kdb_T: Benthic algae death rate (1/d)
        Fb: Fraction of bottom area available for benthic algae growth
        Fw: Fraction of benthic algae mortality into water column
        depth: Depth of water in computation cell (m)
        use_Balgae: Option for considering benthic algae in DOC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, Ab * kdb_T * Fb * (1 - Fw) / depth, 0)

    return da


    
    

