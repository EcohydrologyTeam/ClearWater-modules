import numpy as np
import numba
import xarray as xr
from clearwater_modules.shared.processes import (
    arrhenius_correction,
)


@numba.njit
def kbod_T(
    water_temp_c: xr.DataArray,
    kbod_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted CBOD oxidation rate (1/d)

    Args:
        water_temp_c: water temperature in Celsius
        kbod_20: CBOD oxidation rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient
    """

    kbod_T = arrhenius_correction(water_temp_c, kbod_20, theta)
    return kbod_T


@numba.njit
def ksbod_T(
    water_temp_c: xr.DataArray,
    ksbod_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted CBOD sedimentation rate (m/d)

    Args:
        water_temp_c: water temperature in Celsius
        ksbod_20: CBOD sedimentation rate at 20 degrees Celsius (m/d)
        theta: Arrhenius coefficient
    """

    ksbod_T = arrhenius_correction(water_temp_c, ksbod_20, theta)
    return ksbod_T



def CBOD_oxidation(
    DOX: xr.DataArray,
    CBOD: xr.DataArray,
    kbod_T: xr.DataArray,
    KsOxbod: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Calculates CBOD oxidation

    Args:
        DOX: Dissolved oxygen concentration (mg-O2/L)
        CBOD: Carbonaceous biochemical oxygen demand (mg-O2/L)
        kbod_T: Temperature adjusted CBOD oxidation rate (1/d)
        KsOxbod: Half-saturation oxygen attenuation for CBOD oxidation (mg-O2/L)
        use_DOX: Option to consider DOX concentration in calculation of CBOD oxidation
    """
    da: xr.DataArray = xr.where(use_DOX == True, (DOX / (KsOxbod + DOX)) * kbod_T * CBOD, kbod_T * CBOD)
    
    return da


@numba.njit
def CBOD_sedimentation(
    CBOD: xr.DataArray,
    ksbod_T: xr.DataArray
) -> xr.DataArray:
    """Calculates CBOD sedimentation for each group

    Args:
        CBOD: CBOD concentration (mg-O2/L)
        ksbod_T: Temperature adjusted sedimentation rate (m/d)
    """
    
    CBOD_sedimentation = CBOD * ksbod_T
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
def CBOD_new(
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
