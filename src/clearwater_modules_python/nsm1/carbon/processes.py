import numpy as np
from clearwater_modules_python.shared.processes import (
    arrhenius_correction,
)        


def compute_h(
    surface_area: float,
    volume: float
) -> float:
    """Computes average water depth in a cell

    Args:
        surface_area:
        volume:
    """
    return volume / surface_area

def k_poc_T(
    water_temp_c: float,
    k_poc_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        water_temp_c: water temperature in Celsius
        k_poc_20: POC hydrolysis rate at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, k_poc_20, theta)

def POC_hydrolysis(
    k_poc_T: float,
    POC: float,
) -> float:
    """Calculate the POC lost to hydrolysis for a given timestep

    Args:
        k_poc_T: POC hydrolysis rate at given water temperature
        POC: POC concentration
    """
    return k_poc_T * POC

def POC_settling(
    vsoc: float,
    h: float,
    POC: float
) -> float:
    """Calculate the POC lost due to settling for a given timestep

    Args:
        vsoc: POC settling velocity
        h: water depth of cell
        POC: POC concentration
    """
    return vsoc / h * POC

def POC_algal_mortality(
    f_pocp: float,
    kdp_T: float,
    rca: float,
    Ap: float,
) -> float:
    """Calculate the POC created due to algal mortality

    Args:
        f_pocp: fraction of algal mortality into POC
        kdp_T: algal death rate at water temperature
        rca: algal C to chlorophyll-a ratio
        Ap: algae concentration
    """
    return f_pocp * kdp_T * rca * Ap

def POC_benthic_algae_mortality(
    h: float,
    F_pocb: float,
    kdb_T: float,
    rcb: float,
    Ab: float,
    Fb: float,
    Fw: float
) -> float:
    """Calculate the POC created due to benthic algae mortality

    Args: 
        h: water depth in cell
        F_pocb: fraction of benthic algal mortality into POC
        kdb_T: benthic algae death rate
        rcb: benthic algae C to biomass weight ratio
        Ab: benthic algae concentration
        Fb:
        Fw:
    """
    return (1/h) * F_pocb * kdb_T * rcb * Ab * Fb * Fw

def POC_change(
    POC_settling: float,
    POC_hydrolysis: float,
    POC_algal_mortality: float,
    POC_benthic_algae_mortality: float 
) -> float:
    """Calculate the change in POC

    Args:
        POC_settling:
        POC_hydrolysis:
        POC_algal_mortality:
        POC_benthic_algae_mortality:
    """
    return POC_algal_mortality + POC_benthic_algae_mortality - POC_settling - POC_hydrolysis

def update_POC(
    POC: float,
    dPOCdt: float,
) -> float:
    """Calculate the POC concentration at the next time step

    Args:
        POC:
        dPOCdt:
    """
    return POC + dPOCdt

        #DOC###################################################################

def DOC_algal_mortality(
        f_pocp: float,
        kdp_T: float,
        rca: float,
        Ap: float,
) -> float:
    """Calculate the DOC created due to algal mortality

    Args:
        f_pocp: fraction of algal mortality into POC
        kdp_T: algal death rate at water temperature
        rca: algal C to chlorophyll-a ratio
        Ap: algae concentration
    """
    return (1 - f_pocp) * kdp_T * rca * Ap

def DOC_benthic_algae_mortality(
    h: float,
    F_pocb: float,
    kdb_T: float,
    rcb: float,
    Ab: float,
    Fb: float,
    Fw: float
) -> float:
    """Calculate the DOC created due to benthic algae mortality

    Args: 
        h: water depth in cell
        F_pocb: fraction of benthic algal mortality into POC
        kdb_T: benthic algae death rate
        rcb: benthic algae C to biomass weight ratio
        Ab: benthic algae concentration
        Fb:
        Fw:
    """
    return (1/h) * (1 - F_pocb) * kdb_T * rcb * Ab * Fb * Fw

def DOC_oxidation(
    DOX: float,
    KsOxmc: float,
    kdoc_T: float,
    DOC: float 
) -> float:
    """Calculates the DOC lost due to oxidation

    Args:
        DOX: concentration of dissolved oxygen
        KsOxmc: half saturation oxygen attenuation constant for DOC oxidation rate 
        kdoc_T: DOC oxidation rate
        DOC: concentration of dissolved organic carbon
    """
    return DOX / (KsOxmc + DOX) * kdoc_T * DOC

def DOC_change(
    DOC_oxidation: float,
    POC_hydrolysis: float,
    DOC_algal_mortality: float,
    DOC_benthic_algae_mortality: float
) -> float:
    """Calculates the change in DOC

    Args:
        DOC_oxidation:
        POC_hydrolysis:
        DOC_algal_mortality: 
        DOC_benthic_algae_mortality: 
    """
    return POC_hydrolysis + DOC_algal_mortality + DOC_benthic_algae_mortality - DOC_oxidation

def update_DOC(
    DOC: float,
    dDOCdt: float,
) -> float:
    """Calculate the DOC concentration at the next time step

    Args:
        DOC:
        dDOCdt:
    """
    return DOC + dDOCdt


        #DIC####################################################################

def Henrys_k(
    water_temp_k: float
)-> float:
    """Calculates the temperature dependent Henry's coefficient

    Args:
        water_temp_k: water temperature in kelvin
    """
    return 10**(2385.73 / water_temp_k + .0152642 * water_temp_k - 14.0184)

def kac_T(
    water_temp_c: float,
    kac_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        water_temp_c: water temperature in Celsius
        k_poc_20: POC hydrolysis rate at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kac_20, theta)


def Atmospheric_CO2_reaeration(
    kac_T: float,
    K_H: float,
    pCO2: float,
    FCO2: float,
    DIC: float
) -> float:
    """Calculates the atmospheric input of CO2 into the waterbody

    Args:
        kac_T:
        K_H:
        pCO2:
        FCO2:
        DIC:
    """
    return 12 * kac_T * (10**-3 * K_H * pCO2 - 10**3 * FCO2 * DIC)

def DIC_algal_respiration(
    ApRespiration: float,
    rca: float
) -> float:
    """Calculates DIC flux due to algal respiration

    Args:
        ApRespiration: 
        rca: 
    """
    return ApRespiration * rca 
    
def DIC_algal_photosynthesis(
    ApGrowth: float,
    rca: float
) -> float:
    """Calculates DIC flux due to algal growth

    Args:
        ApGrowth:
        rca:
    """
    return ApGrowth * rca

def DIC_benthic_algae_respiration(
    AbRespiration: float,
    rcb: float,
    Fb: float,
    h: float
) -> float:
    """Calculates DIC flux due to benthic algae respiration

    Args:
        AbRespiration: 
        rcb:
        Fb:
        h:
    """
    return AbRespiration * rcb * Fb * (1/h) 
    
def DIC_benthic_algae_photosynthesis(
    AbGrowth: float,
    rcb: float,
    Fb: float,
    h: float
) -> float:
    """Calculates DIC flux due to benthic algae growth

    Args:
        AbGrowth:
        rcb:
        Fb:
        h:
    """
    return AbGrowth * rcb * Fb * (1/h)

def DIC_CBOD_oxidation(
    DOX: float,
    CBOD_i: np.array,
    roc: float,
    kbod_i_T: np.array,
    KsOxbod_i: np.array
) -> float:
    """Calculates CBOD oxidation
    
    Args:
        DOX:
        CBOD_i
        roc:
        kbod_i_T:
        KsOxbod_i:
    """
    nCBOD = len(CBOD_i)
    CBOD_ox = 0
    for i in nCBOD:
        CBOD_ox = CBOD_ox + (DOX / (KsOxbod_i[i] + DOX)) * kbod_i_T[i] * CBOD_i[i]
    return CBOD_ox / roc

def DIC_sed_release(
    SOD_T: float,
    roc: float,
    h: float
) -> float:
    """Computes the sediment release of DIC

    Args:
        SOD_T:
        roc:
        h:
    """

def DIC_change(
    Atm_CO2_reaeration: float,
    DIC_algal_respiration: float,
    DIC_algal_photosynthesis: float,
    DIC_benthic_algae_respiration: float,
    DIC_benthic_algae_photosynthesis: float,
    DIC_CBOD_oxidation: float,
    DIC_sed_release: float
) -> float:
    """Calculates the change in DIC

    Args:
        Atm_CO2_reaeration:
        DIC_algal_respiration:
        DIC_algal_photosynthesis:
        DIC_benthic_algae_respiration:
        DIC_benthic_algae_photosynthesis:
        DIC_CBOD_oxidation:
        DIC_sed_release: 
    """
    return Atm_CO2_reaeration + DIC_algal_respiration + DIC_benthic_algae_respiration + DIC_CBOD_oxidation + DIC_sed_release - DIC_algal_photosynthesis - DIC_benthic_algae_photosynthesis

def update_DIC(
    DIC: float,
    dDICdt: float,
) -> float:
    """Calculate the DIC concentration at the next time step

    Args:
        DIC:
        dDICdt:
    """
    return DIC + dDICdt

