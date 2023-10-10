import numpy as np
import numba
from clearwater_modules.shared.processes import (
    arrhenius_correction,
)


@numba.njit
def kbod_i_T(
    water_temp_c: float,
    kbod_i_20: np.array,
    theta: float
) -> np.array:
    """Calculate the temperature adjusted CBOD oxidation rates for each group (1/d)

    Args:
        water_temp_c: water temperature in Celsius
        kbod_i_20: CBOD oxidation rate at 20 degrees Celsius for each CBOD group (1/d)
        theta: Arrhenius coefficient
    """
    kbod_i_T = np.array([])
    for i in kbod_i_20:
        np.append(kbod_i_T, arrhenius_correction(water_temp_c, i, theta))
    return kbod_i_T


@numba.njit
def ksbod_i_T(
    water_temp_c: float,
    ksbod_i_20: np.array,
    theta: float
) -> np.array:
    """Calculate the temperature adjusted CBOD sedimentation rates for each group (m/d)

    Args:
        water_temp_c: water temperature in Celsius
        ksbod_i_20: CBOD sedimentation rate at 20 degrees Celsius for each CBOD group (m/d)
        theta: Arrhenius coefficient
    """
    ksbod_i_T = np.array([])
    for i in ksbod_i_20:
        np.append(ksbod_i_T, arrhenius_correction(water_temp_c, i, theta))
    return ksbod_i_T


@numba.njit
def CBOD_oxidation(
    DOX: float,
    CBOD: np.array,
    kbod_i_T: np.array,
    KsOxbod_i: np.array,
    use_DOX: bool
) -> np.array:
    """Calculates CBOD oxidation for each group

    Args:
        DOX: Dissolved oxygen concentration (mg-O2/L)
        CBOD: Carbonaceous biochemical oxygen demand for each CBOD group (mg-O2/L, array)
        kbod_i_T: Temperature adjusted CBOD oxidation rate for each CBOD group (1/d, array)
        KsOxbod_i: Half-saturation oxygen attenuation for CBOD oxidation for each CBOD group (mg-O2/L, array)
        use_DOX: Option to consider DOX concentration in calculation of CBOD oxidation
    """
    nCBOD = len(CBOD)
    CBOD_ox = np.array([])

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
    CBOD: np.array,
    ksbod_i_T: np.array
) -> np.array:
    """Calculates CBOD sedimentation for each group

    Args:
        CBOD: CBOD concentration for each CBOD group (mg-O2/L)
        ksbod_i_T: Temperature adjusted sedimentation rate for each CBOD group (m/d, array)
    """
    nCBOD = len(CBOD)
    CBOD_sedimentation = np.array([])
    for i in nCBOD:
        np.append(CBOD_sedimentation, CBOD[i] * ksbod_i_T[i])
    return CBOD_sedimentation


@numba.njit
def dCBODdt(
    CBOD_oxidation: np.array,
    CBOD_sedimentation: np.array
) -> np.array:
    """Computes change in each CBOD group for a given timestep

    Args:
        CBOD_oxidation: CBOD concentration change due to oxidation (mg/L/d)
        CBOD_sedimentation: CBOD concentration change due to sedimentation (mg/L/d)
    """
    return - CBOD_oxidation - CBOD_sedimentation


@numba.njit
def CBOD_new(
    CBOD: np.array,
    dCBODdt: np.array,
    timestep: float
) -> np.array:
    """Calculates new CBOD concentration for next timestep

    Args:
        CBOD: CBOD concentration from previous timestep (mg/L)
        dCBODdt: CBOD concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return CBOD + dCBODdt * timestep
