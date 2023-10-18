import numpy as np
import numba
import xarray as xr
from clearwater_modules.shared.processes import (
    arrhenius_correction
)
from clearwater_modules.nsm1.carbon import dynamic_variables
from clearwater_modules.nsm1.carbon import static_variables
from clearwater_modules.nsm1 import static_variables_global
# from clearwater_modules.nsm1 import dynamic_variables_global
from clearwater_modules.nsm1 import state_variables


@numba.njit
def kpoc_T(
    water_temp_c: xr.DataArray,
    kpoc_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kpoc_20: POC hydrolysis rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient for kpoc
    """
    return arrhenius_correction(water_temp_c, kpoc_20, theta)


@numba.njit
def POC_hydrolysis(
    kpoc_T: xr.DataArray,
    POC: xr.DataArray,
) -> xr.DataArray:
    """Calculate the POC concentration change due to hydrolysis for a given timestep

    Args:
        kpoc_T: POC hydrolysis rate at given water temperature (1/d)
        POC: POC concentration (mg/L)
    """
    return kpoc_T * POC


@numba.njit
def POC_settling(
    vsoc: xr.DataArray,
    depth: xr.DataArray,
    POC: xr.DataArray
) -> xr.DataArray:
    """Calculate the POC concentration change due to settling for a given timestep

    Args:
        vsoc: POC settling velocity (m/d)
        depth: Water depth of cell (m)
        POC: POC concentration (mg/L)
    """
    return vsoc / depth * POC


def POC_algal_mortality(
    f_pocp: xr.DataArray,
    kdp_T: xr.DataArray,
    rca: xr.DataArray,
    Ap: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculate the POC concentration change due to algal mortality

    Args:
        f_pocp: Fraction of algal mortality into POC
        kdp_T: Algal death rate at water temperature (1/d)
        rca: Algal C to chlorophyll-a ratio (mg-C/ugChla)
        Ap: Algae concentration (mg/L)
        use_Algae: Option for considering algae in POC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Algae == True, f_pocp * kdp_T * rca * Ap, 0)

    return da


def POC_benthic_algae_mortality(
    depth: xr.DataArray,
    F_pocb: xr.DataArray,
    kdb_T: xr.DataArray,
    rcb: xr.DataArray,
    Ab: xr.DataArray,
    Fb: xr.DataArray,
    Fw: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
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
    da: xr.DataArray = xr.where(use_Balgae == True, (1 / depth) * F_pocb * kdb_T * rcb * Ab * Fb * Fw, 0)

    return da

@numba.njit
def dPOCdt(
    POC_settling: xr.DataArray,
    POC_hydrolysis: xr.DataArray,
    POC_algal_mortality: xr.DataArray,
    POC_benthic_algae_mortality: xr.DataArray
) -> xr.DataArray:
    """Calculate the change in POC concentration

    Args:
        POC_settling: Concentration change of POC due to settling (mg/L/d)
        POC_hydrolysis: Concentration change of POC due to hydrolysis (mg/L/d)
        POC_algal_mortality: Concentration change of POC due to algal mortality (mg/L/d)
        POC_benthic_algae_mortality: Concentration change of POC due to benthic algae mortality (mg/L/d)
    """
    return POC_algal_mortality + POC_benthic_algae_mortality - POC_settling - POC_hydrolysis


@numba.njit
def POC_new(
    POC: xr.DataArray,
    dPOCdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculate the POC concentration at the next time step

    Args:
        POC: Concentration of POC from previous timestep (mg/L)
        dPOCdt: POC concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return POC + dPOCdt * timestep


def DOC_algal_mortality(
        f_pocp: xr.DataArray,
        kdp_T: xr.DataArray,
        rca: xr.DataArray,
        Ap: xr.DataArray,
        use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculate the DOC concentration change due to algal mortality

    Args:
        f_pocp: Fraction of algal mortality into POC 
        kdp_T: Algal death rate at water temperature (1/d) 
        rca: Algal C to chlorophyll-a ratio (mg-C/ug-Chla)
        Ap: Algae concentration (mg/L)
        use_Algae: Option for considering algae in DOC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Algae == True, (1 - f_pocp) * kdp_T * rca * Ap, 0)

    return da


def DOC_benthic_algae_mortality(
    depth: xr.DataArray,
    F_pocb: xr.DataArray,
    kdb_T: xr.DataArray,
    rcb: xr.DataArray,
    Ab: xr.DataArray,
    Fb: xr.DataArray,
    Fw: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
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
    da: xr.DataArray = xr.where(use_Balgae == True, (1 / depth) * (1 - F_pocb) * kdb_T * rcb * Ab * Fb * Fw, 0)

    return da


@numba.njit
def kdoc_T(
    water_temp_c: xr.DataArray,
    kdoc_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted DOC oxidation rate (1/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kdoc_20: DOC oxidation rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kdoc_20, theta)


def DOC_oxidation(
    DOX: xr.DataArray,
    KsOxmc: xr.DataArray,
    kdoc_T: xr.DataArray,
    DOC: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Calculates the DOC concentration change due to oxidation

    Args:
        DOX: Concentration of dissolved oxygen (mg/L)
        KsOxmc: Half saturation oxygen attenuation constant for DOC oxidation rate (mg-O2/L)
        kdoc_T: DOC oxidation rate (1/d)
        DOC: Concentration of dissolved organic carbon (mg/L)
        use_DOX: Option for considering dissolved oxygen concentration in DOC oxidation calculation (boolean)
    """
    da: xr.DataArray = xr.where(use_DOX == True, DOX / (KsOxmc + DOX) * kdoc_T * DOC, kdoc_T * DOC)

    return da


@numba.njit
def dDOCdt(
    DOC_oxidation: xr.DataArray,
    POC_hydrolysis: xr.DataArray,
    DOC_algal_mortality: xr.DataArray,
    DOC_benthic_algae_mortality: xr.DataArray
) -> xr.DataArray:
    """Calculates the change in DOC concentration

    Args:
        POC_hydrolysis: DOC concentration change due to POC hydrolysis (mg/L/d)
        DOC_POM_dissolution: DOC concentration change due to POM dissolution (mg/L/d)
        DOC_denitrification: DOC concentration change due to DOC denitrification (mg/L/d)
        DOC_algal_mortality: DOC concentration change due to algal mortality (mg/L/d)
        DOC_benthic_algae_mortality: DOC concentration change due to benthic algae mortality (mg/L/d)
        DOC_oxidation: DOC concentration change due to DOC oxidation (mg/L/d)
    """
    return POC_hydrolysis + DOC_algal_mortality + DOC_benthic_algae_mortality - DOC_oxidation


@numba.njit
def DOC_new(
    DOC: xr.DataArray,
    dDOCdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculate the DOC concentration at the next time step

    Args:
        DOC: Dissolved organic carbon concentration from previous timestep (mg/L)
        dDOCdt: Dissolved organic carbon concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return DOC + dDOCdt * timestep


@numba.njit
def Henrys_k(
    water_temp_c: xr.DataArray
) -> xr.DataArray:
    """Calculates the temperature dependent Henry's coefficient (mol/L/atm)

    Args:
        water_temp_c: Water temperature in celsius
    """
    return 10**(2385.73 / (water_temp_c + 273.15) + .0152642 * (water_temp_c + 273.15) - 14.0184)


@numba.njit
def kac_T(
    water_temp_c: xr.DataArray,
    kac_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted CO2 reaeration rate (1/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kac_20: CO2 reaeration rate at 20 degrees Celsius (1/d)
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kac_20, theta)


@numba.njit
def Atmospheric_CO2_reaeration(
    kac_T: xr.DataArray,
    K_H: xr.DataArray,
    pCO2: xr.DataArray,
    FCO2: xr.DataArray,
    DIC: xr.DataArray
) -> xr.DataArray:
    """Calculates the atmospheric input of CO2 into the waterbody

    Args:
        kac_T: CO2 reaeration rate adjusted for temperature (1/d)
        K_H: Henry's Law constant (mol/L/atm)
        pCO2: Partial pressure of CO2 in the atmosphere (ppm)
        FCO2: Fraction of CO2 in total inorganic carbon
        DIC: Dissolved inorganic carbon concentration (mg/L)
    """
    return 12 * kac_T * (10**-3 * K_H * pCO2 - 10**3 * FCO2 * DIC)


def DIC_algal_respiration(
    ApRespiration: xr.DataArray,
    rca: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC concentration change due to algal respiration

    Args:
        ApRespiration: Algae respiration calculated in algae module (ug-Chla/L/d)
        rca: Ratio of carbon to chlorophyll-a (mg-C/ug-Chla)
        use_Algae: Option to consider algae in the DIC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Algae == True, ApRespiration * rca, 0)

    return da


def DIC_algal_photosynthesis(
    ApGrowth: xr.DataArray,
    rca: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC concentration change due to algal photosynthesis

    Args:
        ApGrowth: Algal photosynthesis calculated in algae module (ug-Chla/L/d)
        rca: Ratio of carbon to chlorophyll-a (mg-C/ug-Chla)
        use_Algae: Option to consider algae in the DIC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Algae == True, ApGrowth * rca, 0)

    return da


def DIC_benthic_algae_respiration(
    AbRespiration: xr.DataArray,
    rcb: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC flux due to benthic algae respiration

    Args:
        AbRespiration: Benthic algae respiration calculated in benthic algae module (g/m2/d)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Depth of water (m)
        use_Balgae: Option to consider benthic algae in the DIC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, AbRespiration * rcb * Fb * (1 / depth), 0)

    return da


def DIC_benthic_algae_photosynthesis(
    AbGrowth: xr.DataArray,
    rcb: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC flux due to benthic algae growth

    Args:
        AbGrowth: Benthic algae photosynthesis calculated in the benthic algae module (g/m2/d)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Depth of water (m)
        use_Balgae: Option to consider benthic algae in the DIC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, AbGrowth * rcb * Fb * (1 / depth), 0)

    return da


def DIC_CBOD_oxidation(
    DOX: xr.DataArray,
    CBOD: xr.DataArray,
    roc: xr.DataArray,
    kbod_T: xr.DataArray, #imported from CBOD module
    KsOxbod: xr.DataArray, #imported from CBOD module
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC concentration change due to CBOD oxidation

    Args:
        DOX: Dissolved oxygen concentration (mg/L)
        CBOD: Carbonaceous biochemical oxygen demand concentration (mg/L) 
        roc: Ratio of O2 to carbon for carbon oxidation (mg-O2/mg-C)
        kbod_T: CBOD oxidation rate (1/d)
        KsOxbod: Half saturation oxygen attenuation constant for CBOD oxidation (mg-O2/L)
        use_DOX: Option to consider dissolved oxygen in CBOD oxidation calculation (boolean)
    """
    
    da: xr.DataArray = xr.where(use_DOX == True, (1 / roc) * (DOX / (KsOxbod + DOX)) * kbod_T * CBOD, CBOD * kbod_T)

    return da



def DIC_sed_release(
    SOD_tc: xr.DataArray,
    roc: xr.DataArray,
    depth: xr.DataArray,
    JDIC: xr.DataArray,
    use_SedFlux: xr.DataArray
) -> xr.DataArray:
    """Computes the sediment release of DIC

    Args:
        SOD_T: Sediment oxygen demand adjusted for water temperature (mg-O2/L/d)
        roc: Ratio of O2 to carbon for carbon oxidation (mg-O2/mg-C)
        depth: Water depth (m)
        JDIC: Sediment-water flux of dissolved inorganic carbon (g-C/m2/d)
        use_SedFlux: Option to consider full sediment flux budget in DIC sediment contribution (bool)
    """
    da: xr.DataArray = xr.where(use_SedFlux == True, JDIC / depth, SOD_tc / roc / depth)

    return da

@numba.njit
def dDICdt(
    Atm_CO2_reaeration: xr.DataArray,
    DIC_algal_respiration: xr.DataArray,
    DIC_algal_photosynthesis: xr.DataArray,
    DIC_benthic_algae_respiration: xr.DataArray,
    DIC_benthic_algae_photosynthesis: xr.DataArray,
    DIC_DOC_oxidation: xr.DataArray,
    DIC_CBOD_oxidation: xr.DataArray,
    DIC_sed_release: xr.DataArray
) -> xr.DataArray:
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
    return Atm_CO2_reaeration + DIC_algal_respiration + DIC_benthic_algae_respiration + DIC_DOC_oxidation + DIC_CBOD_oxidation + DIC_sed_release - DIC_algal_photosynthesis - DIC_benthic_algae_photosynthesis


@numba.njit
def DIC_new(
    DIC: xr.DataArray,
    dDICdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculate the DIC concentration at the next time step

    Args:
        DIC: Concentration of DIC from previous timestep (mg/L)
        dDICdt: Change in concentration of DIC for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return DIC + dDICdt * timestep
