"""
File contains process to calculate new CBOD concentration and associated dependent variables
"""

import numba
import xarray as xr
from clearwater_modules.shared.processes import arrhenius_correction
import math

@numba.njit
def kbod_tc(
    TwaterC: xr.DataArray,
    kbod_20: xr.DataArray,
) -> xr.DataArray:
    """Calculate the temperature adjusted CBOD oxidation rate (1/d)

    Args:
        TwaterC: water temperature in Celsius
        kbod_20: CBOD oxidation rate at 20 degrees Celsius (1/d)
    """

    kbod_tc = arrhenius_correction(TwaterC, kbod_20, 1.047)
    return kbod_tc


@numba.njit
def ksbod_tc(
    TwaterC: xr.DataArray,
    ksbod_20: xr.DataArray,
) -> xr.DataArray:
    """Calculate the temperature adjusted CBOD sedimentation rate (m/d)

    Args:
        TwaterC: water temperature in Celsius
        ksbod_20: CBOD sedimentation rate at 20 degrees Celsius (m/d)
    """

    ksbod_tc = arrhenius_correction(TwaterC, ksbod_20, 1.024)
    return ksbod_tc



def CBOD_oxidation(
    DOX: xr.DataArray,
    CBOD: xr.DataArray,
    kbod_tc: xr.DataArray,
    KsOxbod: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Calculates CBOD oxidation

    Args:
        DOX: Dissolved oxygen concentration (mg-O2/L)
        CBOD: Carbonaceous biochemical oxygen demand (mg-O2/L)
        kbod_tc: Temperature adjusted CBOD oxidation rate (1/d)
        KsOxbod: Half-saturation oxygen attenuation for CBOD oxidation (mg-O2/L)
        use_DOX: Option to consider DOX concentration in calculation of CBOD oxidation
    """
    da: xr.DataArray = xr.where(use_DOX == True, (DOX / (KsOxbod + DOX)) * kbod_tc * CBOD, kbod_tc * CBOD)
    
    return da


@numba.njit
def CBOD_sedimentation(
    CBOD: xr.DataArray,
    ksbod_tc: xr.DataArray
) -> xr.DataArray:
    """Calculates CBOD sedimentation for each group

    Args:
        CBOD: CBOD concentration (mg-O2/L)
        ksbod_tc: Temperature adjusted sedimentation rate (m/d)
    """
    
    CBOD_sedimentation = CBOD * ksbod_tc
    return CBOD_sedimentation


@numba.njit
def dCBODdt(
    CBOD_oxidation: xr.DataArray,
    CBOD_sedimentation: xr.DataArray
) -> xr.DataArray:
    """Computes change in each CBOD group for a given timestep

    Args:
        CBOD_oxidation: CBOD concentration change due to oxidation (mg/L/d)
        CBOD_sedimentation: CBOD concentration change due to sedimentation (mg/L/d)
    """
    return - CBOD_oxidation - CBOD_sedimentation


@numba.njit
def CBOD(
    CBOD: xr.DataArray,
    dCBODdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculates new CBOD concentration for next timestep

    Args:
        CBOD: CBOD concentration from previous timestep (mg/L)
        dCBODdt: CBOD concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return CBOD + dCBODdt * timestep
