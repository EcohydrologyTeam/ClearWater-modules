import numpy as np
import numba
import xarray as xr
from clearwater_modules.shared.processes import (
    arrhenius_correction,
)

# TODO: Figure out how to handle multiple CBOD groups... right now using np.append...xr.concat? Or build the multiple groups outside of modules

@numba.njit
def kbod_i_T(
    water_temp_c: xr.DataArray,
    kbod_i_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted CBOD oxidation rates for each group (1/d)

    Args:
        water_temp_c: water temperature in Celsius
        kbod_i_20: CBOD oxidation rate at 20 degrees Celsius for each CBOD group (1/d)
        theta: Arrhenius coefficient
    """
    kbod_i_T = xr.DataArray([])
    for i in kbod_i_20:
        np.append(kbod_i_T, arrhenius_correction(water_temp_c, i, theta))
    return kbod_i_T


@numba.njit
def ksbod_i_T(
    water_temp_c: xr.DataArray,
    ksbod_i_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted CBOD sedimentation rates for each group (m/d)

    Args:
        water_temp_c: water temperature in Celsius
        ksbod_i_20: CBOD sedimentation rate at 20 degrees Celsius for each CBOD group (m/d)
        theta: Arrhenius coefficient
    """
    ksbod_i_T = xr.DataArray([])
    for i in ksbod_i_20:
        np.append(ksbod_i_T, arrhenius_correction(water_temp_c, i, theta))
    return ksbod_i_T


@numba.njit
def CBOD_oxidation(
    DOX: xr.DataArray,
    CBOD: xr.DataArray,
    kbod_i_T: xr.DataArray,
    KsOxbod_i: xr.DataArray,
    use_DOX:xr.DataArray
) -> xr.DataArray:
    """Calculates CBOD oxidation for each group

    Args:
        DOX: Dissolved oxygen concentration (mg-O2/L)
        CBOD: Carbonaceous biochemical oxygen demand for each CBOD group (mg-O2/L, array)
        kbod_i_T: Temperature adjusted CBOD oxidation rate for each CBOD group (1/d, array)
        KsOxbod_i: Half-saturation oxygen attenuation for CBOD oxidation for each CBOD group (mg-O2/L, array)
        use_DOX: Option to consider DOX concentration in calculation of CBOD oxidation
    """
    nCBOD = len(CBOD)
    CBOD_ox = xr.DataArray([])

    if use_DOX:
        for i in nCBOD:
            np.append(
                CBOD_ox, (DOX / (KsOxbod_i[i] + DOX)) * kbod_i_T[i] * CBOD[i])
    else:
        for i in nCBOD:
            np.append(kbod_i_T[i] * CBOD[i])

    return CBOD_ox


@numba.njit
def CBOD_sedimentation(
    CBOD: xr.DataArray,
    ksbod_i_T: xr.DataArray
) -> xr.DataArray:
    """Calculates CBOD sedimentation for each group

    Args:
        CBOD: CBOD concentration for each CBOD group (mg-O2/L)
        ksbod_i_T: Temperature adjusted sedimentation rate for each CBOD group (m/d, array)
    """
    nCBOD = len(CBOD)
    CBOD_sedimentation = xr.DataArray([])
    for i in nCBOD:
        np.append(CBOD_sedimentation, CBOD[i] * ksbod_i_T[i])
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
