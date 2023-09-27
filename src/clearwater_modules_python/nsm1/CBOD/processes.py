import numpy as np
from clearwater_modules_python.shared.processes import (
    arrhenius_correction,
)

def kbod_i_T(
    water_temp_c: float,
    kbod_i_20: np.array,
    theta: float
) -> np.array:
    """Calculate the temperature adjusted CBOD oxidation rates for each group (/d)

    Args:
        water_temp_c: water temperature in Celsius
        kbod_i_20: CBOD oxidation rate at 20 degrees Celsius for each CBOD group
        theta: Arrhenius coefficient
    """
    kbod_i_T = np.array([])
    for i in kbod_i_20:
        np.append(kbod_i_T,arrhenius_correction(water_temp_c, i, theta))        
    return kbod_i_T

def ksbod_i_T(
    water_temp_c: float,
    ksbod_i_20: np.array,
    theta: float
) -> np.array:
    """Calculate the temperature adjusted CBOD sedimentation rates for each group (m/d)

    Args:
        water_temp_c: water temperature in Celsius
        ksbod_i_20: CBOD sedimentation rate at 20 degrees Celsius for each CBOD group
        theta: Arrhenius coefficient
    """
    ksbod_i_T = np.array([])
    for i in ksbod_i_20:
        np.append(ksbod_i_T,arrhenius_correction(water_temp_c, i, theta))        
    return ksbod_i_T

def CBOD_oxidation(
    DOX: float,
    CBOD_i: np.array,
    kbod_i_T: np.array,
    KsOxbod_i: np.array,
    use_DOX: bool
) -> np.array:
    """Calculates CBOD oxidation for each group
    
    Args:
        DOX:
        CBOD_i:
        kbod_i_T:
        KsOxbod_i:
        use_DOX:
    """
    nCBOD = len(CBOD_i)
    CBOD_ox = np.array([])

    if use_DOX:
        for i in nCBOD:
            np.append(CBOD_ox, (DOX / (KsOxbod_i[i] + DOX)) * kbod_i_T[i] * CBOD_i[i])
    else:
        for i in nCBOD:
            np.append(kbod_i_T[i] * CBOD_i[i])
    
    return CBOD_ox

def CBOD_sedimentation(
    CBOD_i: np.array,
    ksbod_i_T: np.array
) -> np.array:
    """Calculates CBOD sedimentation for each group
    
    Args:
        CBOD_i:
        ksbod_i_T:
    """
    nCBOD = len(CBOD_i)
    CBOD_sedimentation = np.array([])
    for i in nCBOD:
        np.append(CBOD_sedimentation, CBOD_i[i] * ksbod_i_T[i])
    return CBOD_sedimentation

def CBOD_change(
    CBOD_oxidation: np.array,
    CBOD_sedimentation: np.array
) -> np.array:
    """Computes change in each CBOD group for a given timestep
    
    Args:
        CBOD_oxidation:
        CBOD_sedimentation:
    """
    return - CBOD_oxidation - CBOD_sedimentation

def update_CBOD(
    CBOD_i: np.array,
    dCBOD_idt: np.array
) -> np.array:
    """Calculates new CBOD concentration for next timestep
    
    Args:
        CBOD_i:
        dCBOD_idt:
    """
    return CBOD_i + dCBOD_idt
