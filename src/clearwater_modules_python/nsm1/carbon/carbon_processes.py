import numpy as np
import numba
from clearwater_modules_python.shared.processes import (
    arrhenius_correction
)
from clearwater_modules_python.nsm1.carbon import dynamic_variables_carbon
from clearwater_modules_python.nsm1.carbon import static_variables_carbon
from clearwater_modules_python.nsm1 import static_variables_global
#from clearwater_modules_python.nsm1 import dynamic_variables_global
from clearwater_modules_python.nsm1 import state_variables


@numba.njit
def kpoc_T(
    water_temp_c: float,
    kpoc_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kpoc_20: POC hydrolysis rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient for kpoc
    """
    return arrhenius_correction(water_temp_c, kpoc_20, theta)

@numba.njit
def POC_hydrolysis(
    kpoc_T: float,
    POC: float,
) -> float:
    """Calculate the POC concentration change due to hydrolysis for a given timestep

    Args:
        kpoc_T: POC hydrolysis rate at given water temperature (1/d)
        POC: POC concentration (mg/L)
    """
    return kpoc_T * POC

@numba.njit
def POC_settling(
    vsoc: float,
    depth: float,
    POC: float
) -> float:
    """Calculate the POC concentration change due to settling for a given timestep

    Args:
        vsoc: POC settling velocity (m/d)
        depth: Water depth of cell (m)
        POC: POC concentration (mg/L)
    """
    return vsoc / depth * POC

@numba.njit
def POC_algal_mortality(
    f_pocp: float,
    kdp_T: float,
    rca: float,
    Ap: float,
    use_Algae: bool
) -> float:
    """Calculate the POC concentration change due to algal mortality

    Args:
        f_pocp: Fraction of algal mortality into POC
        kdp_T: Algal death rate at water temperature (1/d)
        rca: Algal C to chlorophyll-a ratio (mg-C/ugChla)
        Ap: Algae concentration (mg/L)
        use_Algae: Option for considering algae in POC budget (boolean)
    """
    if use_Algae:
        return f_pocp * kdp_T * rca * Ap
    else:
        return 0

@numba.njit
def POC_benthic_algae_mortality(
    depth: float,
    F_pocb: float,
    kdb_T: float,
    rcb: float,
    Ab: float,
    Fb: float,
    Fw: float,
    use_Balgae: bool
) -> float:
    """Calculate the POC concentration change due to benthic algae mortality

    Args: 
        depth: Water depth in cell (m)
        F_pocb: Fraction of benthic algal mortality into POC
        kdb_T: Benthic algae death rate (1/d)
        rcb: Benthic algae C to biomass weight ratio (mg-C/mg-D)
        Ab: Benthic algae concentration (mg/L)
        Fb: Fraction of bottom area available for benthic algae growth
        Fw: Fraction of benthic algae mortality into water column
        use_Balgae: Option for considering benthic algae in POC budget (boolean)
    """
    if use_Balgae:
        return (1 / depth) * F_pocb * kdb_T * rcb * Ab * Fb * Fw
    else:
        return 0

@numba.njit
def POC_change(
    POC_settling: float,
    POC_hydrolysis: float,
    POC_algal_mortality: float,
    POC_benthic_algae_mortality: float 
) -> float:
    """Calculate the change in POC concentration

    Args:
        POC_settling: Concentration change of POC due to settling (mg/L/d)
        POC_hydrolysis: Concentration change of POC due to hydrolysis (mg/L/d)
        POC_algal_mortality: Concentration change of POC due to algal mortality (mg/L/d)
        POC_benthic_algae_mortality: Concentration change of POC due to benthic algae mortality (mg/L/d)
    """
    return POC_algal_mortality + POC_benthic_algae_mortality - POC_settling - POC_hydrolysis

@numba.njit
def update_POC(
    POC: float,
    dPOCdt: float,
    timestep: float
) -> float:
    """Calculate the POC concentration at the next time step

    Args:
        POC: Concentration of POC from previous timestep (mg/L)
        dPOCdt: POC concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return POC + dPOCdt * timestep

        #DOC###################################################################

@numba.njit
def DOC_algal_mortality(
        f_pocp: float,
        kdp_T: float,
        rca: float,
        Ap: float,
        use_Algae: bool
) -> float:
    """Calculate the DOC concentration change due to algal mortality

    Args:
        f_pocp: Fraction of algal mortality into POC 
        kdp_T: Algal death rate at water temperature (1/d) 
        rca: Algal C to chlorophyll-a ratio (mg-C/ug-Chla)
        Ap: Algae concentration (mg/L)
        use_Algae: Option for considering algae in DOC budget (boolean)
    """
    return (1 - f_pocp) * kdp_T * rca * Ap

@numba.njit
def DOC_benthic_algae_mortality(
    depth: float,
    F_pocb: float,
    kdb_T: float,
    rcb: float,
    Ab: float,
    Fb: float,
    Fw: float,
    use_Balgae: bool
) -> float:
    """Calculate the DOC concentration change due to benthic algae mortality

    Args: 
        depth: Water depth in cell (m)
        F_pocb: Fraction of benthic algal mortality into POC
        kdb_T: Benthic algae death rate (1/d)
        rcb: Benthic algae C to biomass weight ratio (mg-C/mg-D)
        Ab: Benthic algae concentration (mg/L)
        Fb: Fraction of bottom area available for benthic algae growth
        Fw: Fraction of benthic algae mortality into water column
        use_Balgae: Option for considering benthic algae in DOC budget (boolean)
    """
    return (1 / depth) * (1 - F_pocb) * kdb_T * rcb * Ab * Fb * Fw

@numba.njit
def kdoc_T(
    water_temp_c: float,
    kdoc_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted DOC oxidation rate (1/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kdoc_20: DOC oxidation rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kdoc_20, theta)

@numba.njit
def kpom_T(
    water_temp_c: float,
    kpom_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted POM dissolution rate (1/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kpom_20: POM dissolution rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kpom_20, theta)

@numba.njit
def DOC_POM_dissolution(
    kpom_T: float,
    fcom: float,
    POM: float,
)-> float:
    """Calculates the DOC concentration change due to dissolution of particulate organic matter

    Args:
        kpom_T: POM dissolution rate (1/d)
        fcom: Fraction of carbon in organic matter (mg-C/mg-D)
        POM: Particulate organic matter concentration (mg/L)
    """
    return kpom_T * fcom * POM

@numba.njit
def DOC_denitrification(
    DOX: float,
    NO3: float,
    kdnit_T: float,
    KsOxdn: float,
    use_DOX: bool
) -> float:
    """Calculates the DOC concentration change due to denitrification
    
    Args:
        DOX: Dissolved oxygen concentration (mg/L)
        NO3: Nitrate concentration (mg/L)
        kdnit_T: Denitrification rate (1/d)
        KsOxdn: Half-saturation oxygen inhibition constant for denitrification (mg-O2/L) 
        use_DOX: Option for considering dissolved oxygen concentration in DOC denitrification calculation (boolean) 
    """
    if use_DOX:
        return (5 / 4) * (12 / 14) * (1 - DOX / KsOxdn + DOX) * kdnit_T * NO3
    else:
        return  (5 / 4) * (12 / 14) * kdnit_T * NO3 

@numba.njit
def DOC_oxidation(
    DOX: float,
    KsOxmc: float,
    kdoc_T: float,
    DOC: float,
    use_DOX: bool 
) -> float:
    """Calculates the DOC concentration change due to oxidation

    Args:
        DOX: Concentration of dissolved oxygen (mg/L)
        KsOxmc: Half saturation oxygen attenuation constant for DOC oxidation rate (mg-O2/L)
        kdoc_T: DOC oxidation rate (1/d)
        DOC: Concentration of dissolved organic carbon (mg/L)
        use_DOX: Option for considering dissolved oxygen concentration in DOC oxidation calculation (boolean)
    """
    if use_DOX:
        return DOX / (KsOxmc + DOX) * kdoc_T * DOC
    else:
        return kdoc_T * DOC

@numba.njit
def DOC_change(
    DOC_oxidation: float,
    POC_hydrolysis: float,
    DOC_algal_mortality: float,
    DOC_benthic_algae_mortality: float,
    DOC_POM_dissolution: float,
    DOC_denitrification: float
) -> float:
    """Calculates the change in DOC concentration

    Args:
        POC_hydrolysis: DOC concentration change due to POC hydrolysis (mg/L/d)
        DOC_POM_dissolution: DOC concentration change due to POM dissolution (mg/L/d)
        DOC_denitrification: DOC concentration change due to DOC denitrification (mg/L/d)
        DOC_algal_mortality: DOC concentration change due to algal mortality (mg/L/d)
        DOC_benthic_algae_mortality: DOC concentration change due to benthic algae mortality (mg/L/d)
        DOC_oxidation: DOC concentration change due to DOC oxidation (mg/L/d)
    """
    return POC_hydrolysis + DOC_POM_dissolution - DOC_denitrification + DOC_algal_mortality + DOC_benthic_algae_mortality - DOC_oxidation

@numba.njit
def update_DOC(
    DOC: float,
    dDOCdt: float,
    timestep: float
) -> float:
    """Calculate the DOC concentration at the next time step

    Args:
        DOC: Dissolved organic carbon concentration from previous timestep (mg/L)
        dDOCdt: Dissolved organic carbon concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return DOC + dDOCdt * timestep


        #DIC####################################################################

@numba.njit
def Henrys_k(
    water_temp_c: float
)-> float:
    """Calculates the temperature dependent Henry's coefficient (mol/L/atm)

    Args:
        water_temp_c: Water temperature in celsius
    """
    return 10**(2385.73 / (water_temp_c + 273.15) + .0152642 * (water_temp_c + 273.15) - 14.0184)

@numba.njit
def kac_T(
    water_temp_c: float,
    kac_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted CO2 reaeration rate (1/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kac_20: CO2 reaeration rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kac_20, theta)

@numba.njit
def Atmospheric_CO2_reaeration(
    kac_T: float,
    K_H: float,
    pCO2: float,
    FCO2: float,
    DIC: float
) -> float:
    """Calculates the atmospheric input of CO2 into the waterbody

    Args:
        kac_T: CO2 reaeration rate adjusted for temperature (1/d)
        K_H: Henry's Law constant (mol/L/atm)
        pCO2: Partial pressure of CO2 in the atmosphere (ppm)
        FCO2: Fraction of CO2 in total inorganic carbon
        DIC: Dissolved inorganic carbon concentration (mg/L)
    """
    return 12 * kac_T * (10**-3 * K_H * pCO2 - 10**3 * FCO2 * DIC)

@numba.njit
def DIC_algal_respiration(
    ApRespiration: float,
    rca: float,
    use_Algae: bool
) -> float:
    """Calculates DIC concentration change due to algal respiration

    Args:
        ApRespiration: Algae respiration calculated in algae module (ug-Chla/L/d)
        rca: Ratio of carbon to chlorophyll-a (mg-C/ug-Chla)
        use_Algae: Option to consider algae in the DIC budget (boolean)
    """
    if use_Algae:
        return ApRespiration * rca
    else:
        return 0 

@numba.njit
def DIC_algal_photosynthesis(
    ApGrowth: float,
    rca: float,
    use_Algae: bool
) -> float:
    """Calculates DIC concentration change due to algal photosynthesis

    Args:
        ApGrowth: Algal photosynthesis calculated in algae module (ug-Chla/L/d)
        rca: Ratio of carbon to chlorophyll-a (mg-C/ug-Chla)
        use_Algae: Option to consider algae in the DIC budget (boolean)
    """
    if use_Algae:
        return ApGrowth * rca
    else:
        return 0

@numba.njit
def DIC_benthic_algae_respiration(
    AbRespiration: float,
    rcb: float,
    Fb: float,
    depth: float,
    use_Balgae: bool
) -> float:
    """Calculates DIC flux due to benthic algae respiration

    Args:
        AbRespiration: Benthic algae respiration calculated in benthic algae module (g/m2/d)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Depth of water (m)
        use_Balgae: Option to consider benthic algae in the DIC budget (boolean)
    """
    if use_Balgae:
        return AbRespiration * rcb * Fb * (1 / depth)
    else:
        return 0

@numba.njit
def DIC_benthic_algae_photosynthesis(
    AbGrowth: float,
    rcb: float,
    Fb: float,
    depth: float,
    use_Balgae: bool
) -> float:
    """Calculates DIC flux due to benthic algae growth

    Args:
        AbGrowth: Benthic algae photosynthesis calculated in the benthic algae module (g/m2/d)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Depth of water (m)
        use_Balgae: Option to consider benthic algae in the DIC budget (boolean)
    """
    if use_Balgae:
        return AbGrowth * rcb * Fb * (1 / depth)
    else:
        return 0

@numba.njit
def DIC_CBOD_oxidation(
    DOX: float,
    CBOD: np.array,
    roc: float,
    kbod_i_T: np.array,
    KsOxbod_i: np.array,
    use_DOX: bool
) -> float:
    """Calculates DIC concentration change due to CBOD oxidation
    
    Args:
        DOX: Dissolved oxygen concentration (mg/L)
        CBOD: Carbonaceous biochemical oxygen demand concentration for multiple groups (mg/L, array) 
        roc: Ratio of O2 to carbon for carbon oxidation (mg-O2/mg-C)
        kbod_i_T: CBOD oxidation rate (1/d, array)
        KsOxbod_i: Half saturation oxygen attenuation constant for CBOD oxidation (mg-O2/L, array)
        use_DOX: Option to consider dissolved oxygen in CBOD oxidation calculation (boolean)
    """
    nCBOD = len(CBOD)
    CBOD_ox = 0

    if use_DOX:
        for i in nCBOD:
            CBOD_ox = CBOD_ox + (DOX / (KsOxbod_i[i] + DOX)) * kbod_i_T[i] * CBOD[i]
        return CBOD_ox / roc
    else:
        for i in nCBOD:
            CBOD_ox = CBOD_ox + CBOD[i] * kbod_i_T[i]

@numba.njit
def DIC_sed_release(
    SOD_tc: float,
    roc: float,
    depth: float,
    JDIC: float,
    use_SedFlux: bool
) -> float:
    """Computes the sediment release of DIC

    Args:
        SOD_T: Sediment oxygen demand adjusted for water temperature (mg-O2/L/d)
        roc: Ratio of O2 to carbon for carbon oxidation (mg-O2/mg-C)
        depth: Water depth (m)
        JDIC: Sediment-water flux of dissolved inorganic carbon (g-C/m2/d)
        use_SedFlux: Option to consider full sediment flux budget in DIC sediment contribution (bool)
    """
    if use_SedFlux:
        return JDIC / depth 
    else:
        return SOD_tc / roc / depth

@numba.njit
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
        Atm_CO2_reaeration: DIC concentration change due to atmospheric CO2 reaeration (mg/L/d)
        DIC_algal_respiration: DIC concentration change due to algal respiration (mg/L/d)
        DIC_algal_photosynthesis: DIC concentration change due to algal photosynthesis (mg/L/d)
        DIC_benthic_algae_respiration: DIC concentration change due to benthic algae respiration (mg/L/d)
        DIC_benthic_algae_photosynthesis: DIC concentration change due to benthic algae photosynthesis (mg/L/d)
        DIC_CBOD_oxidation: DIC concentration change due to CBOD oxidation (mg/L/d)
        DIC_sed_release: DIC concentration change due to sediment release (mg/L/d)
    """
    return Atm_CO2_reaeration + DIC_algal_respiration + DIC_benthic_algae_respiration + DIC_CBOD_oxidation + DIC_sed_release - DIC_algal_photosynthesis - DIC_benthic_algae_photosynthesis

@numba.njit
def update_DIC(
    DIC: float,
    dDICdt: float,
    timestep: float
) -> float:
    """Calculate the DIC concentration at the next time step

    Args:
        DIC: Concentration of DIC from previous timestep (mg/L)
        dDICdt: Change in concentration of DIC for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return DIC + dDICdt * timestep

