import numpy as np
import numba
import xarray as xr
from clearwater_modules.shared.processes import (
    arrhenius_correction,
)


#TODO: make sure np.exp will work here...
@numba.njit
def pwv(
    t_water_k: xr.DataArray
) -> xr.DataArray:
    """Calculate partial pressure of water vapor

    Args:
        t_water_k: Water temperature kelvin
    """
    return np.exp(11.8571 - 3840.70 / t_water_k - 216961 / t_water_k ** 2)


@numba.njit
def DOs_atm_alpha(
    t_water_k: xr.DataArray
) -> xr.DataArray:
    """Calculate DO saturation atmospheric correction coefficient

    Args:
        t_water_k: Water temperature kelvin
    """
    return .000975 - 1.426 * 10 ** -5 * t_water_k + 6.436 * 10 ** -8 * t_water_k ** 2


@numba.njit
def DOX_sat(
    t_water_k: xr.DataArray,
    patm: xr.DataArray,
    pwv: xr.DataArray,
    DOs_atm_alpha: xr.DataArray
) -> xr.DataArray:
    """Calculate DO saturation value

    Args:
        t_water_k: Water temperature kelvin
        patm: Atmospheric pressure (atm)
        pwv: Patrial pressure of water vapor (atm)
        DOs_atm_alpha: DO saturation atmospheric correction coefficient
    """
    DOX_sat_uncorrected = np.exp(-139.34410 + 1.575701 * 10 ** 5 / t_water_k - 6.642308 * 10 ** 7 / t_water_k ** 2
                                 + 1.243800 * 10 ** 10 / t_water_k - 8.621949 * 10 ** 11 / t_water_k)

    DOX_sat_corrected = DOX_sat_uncorrected * patm * \
        (1 - pwv / patm) * (1 - DOs_atm_alpha * patm) / \
        ((1 - pwv) * (1 - DOs_atm_alpha))
    return DOX_sat_corrected


@numba.njit
def Atm_O2_reaeration(
    ka_T: xr.DataArray,
    DOX_sat: xr.DataArray,
    DOX: xr.DataArray
) -> xr.DataArray:
    """Compute the atmospheric O2 reaeration flux

    Args: 
        ka_T: Oxygen reaeration rate adjusted for temperature (1/d)
        DOX_sat: Dissolved oxygen saturation concentration (mg/L)
        DOX: Dissolved oxygen concentration (mg/L)
    """
    return ka_T * (DOX_sat - DOX)


def DOX_ApGrowth(
    ApGrowth: xr.DataArray,
    rca: xr.DataArray,
    roc: xr.DataArray,
    ApUptakeFr_NH4: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Compute DOX flux due to algal photosynthesis

    Args:
        ApGrowth: Algae photosynthesis, calculated in the algae module (ug-Chla/L/d)
        rca: Ratio of algal carbon to chlorophyll-a (mg-C/ug-Chla)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
        ApUptakeFr_NH4: Fraction of actual algal uptake that is from the ammonia pool, calculated in nitrogen module 
    """
    da: xr.DataArray = xr.where(use_Algae == True, ApGrowth * rca * roc * (138 / 106 - 32 * ApUptakeFr_NH4 / 106), 0)

    return da


def DOX_ApRespiration(
    ApRespiration: xr.DataArray,
    rca: xr.DataArray,
    roc: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Compute DOX flux due to algal photosynthesis

    Args:
        ApRespiration: algae respiration, calculated in the algae module
        rca: Ratio of algal carbon to chlorophyll-a (mg-C/ug-Chla)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C) 
    """
    da: xr.DataArray = xr.where(use_Algae == True, ApRespiration * rca * roc, 0)

    return da


def DOX_Nitrification(
    KNR: xr.DataArray,
    DOX: xr.DataArray,
    ron: xr.DataArray,
    knit_tc: xr.DataArray,
    NH4: xr.DataArray,
    use_NH4: xr.DataArray
) -> xr.DataArray:
    """Compute DOX flux due to nitrification of ammonia

    Args:
        KNR: Oxygen inhibition factor for nitrification (mg-O2/L)
        DOX: Dissolved oxygen concentration (mg/L)
        ron: Ratio of oxygen to nitrogen for nitrificiation (mg-O2/mg-N)
        knit_tc: Nitrification rate of NH4 to NO3 (1/d)
        NH4: Ammonia/ammonium concentration
    """
    da: xr.DataArray = xr.where(use_NH4 == True, (1.0 - np.exp(-KNR * DOX)) * ron * knit_tc * NH4, 0)

    return da


def DOX_DOC_Oxidation(
    DOC_Oxidation: xr.DataArray,
    roc: xr.DataArray,
    use_DOC: xr.DataArray
) -> xr.DataArray:
    """Computes dissolved oxygen flux due to oxidation of dissolved organic carbon

    Args:
        DOC_Oxidation: Dissolved organic carbon oxidation, calculated in carbon module (mg/L/d)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
    """
    da: xr.DataArray = xr.where(use_DOC == True, roc * DOC_Oxidation, 0)

    return da


@numba.njit
def DOX_CBOD_Oxidation(
    DIC_CBOD_Oxidation: xr.DataArray,
    roc: xr.DataArray
) -> xr.DataArray:
    """Compute dissolved oxygen flux due to CBOD oxidation

    Args:
        DIC_CBOD_Oxidation: Carbonaceous biochemical oxygen demand oxidation, calculated in CBOD module (mg/L/d)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
    """
    return DIC_CBOD_Oxidation * roc


def DOX_AbGrowth(
    AbUptakeFr_NH4: xr.DataArray,
    roc: xr.DataArray,
    rcb: xr.DataArray,
    AbGrowth: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_BAlgae: xr.DataArray
) -> xr.DataArray:
    """Compute dissolved oxygen flux due to benthic algae growth

    Args:
        AbUptakeFr_NH4: Fraction of actual benthic algal uptake that is from the ammonia pool, calculated in nitrogen module
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        AbGrowth: Benthic algae photosynthesis, calculated in benthic algae module (mg/L/d)
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Water depth (m)
        use_BAlgae: Option to consider benthic algae in the DOX budget
    """
    da: xr.DataArray = xr.where(use_BAlgae == True, (138 / 106 - 32 / 106 * AbUptakeFr_NH4) * roc * rcb * AbGrowth * Fb / depth, 0)

    return da


def DOX_AbRespiration(
    roc: xr.DataArray,
    rcb: xr.DataArray,
    AbRespiration: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_BAlgae: xr.DataArray
) -> xr.DataArray:
    """Compute dissolved oxygen flux due to benthic algae respiration

    Args:
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        AbRespiration: Benthic algae respiration, calculated in the benthic algae module
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Water depth (m)
        use_BAlgae: Option to consider benthic algae in the DOX budget
    """

    da: xr.DataArray = xr.where(use_BAlgae == True, roc * rcb * AbRespiration * Fb / depth, 0)

    return da


def DOX_SOD(
    SOD_Bed: xr.DataArray,
    depth: xr.DataArray,
    SOD_tc: xr.DataArray,
    use_SedFlux: xr.DataArray
) -> xr.DataArray:
    """Compute dissolved oxygen flux due to sediment oxygen demand

    Args:
        SOD_Bed: Sediment oxygen demand if calculated using the SedFlux module (mg-O2/m2)
        depth: Water depth (m)
        SOD_tc: Sediment oxygen demand not considering the SedFlux budget (mg-O2/m2)
        use_SedFlux: Option to consider sediment flux in DOX budget (boolean)
    """

    da: xr.DataArray = xr.where(use_SedFlux == 1, SOD_Bed / depth, SOD_tc / depth)

    return da

@numba.njit
def dDOXdt(
    Atm_O2_reaeration: xr.DataArray,
    DOX_ApGrowth: xr.DataArray,
    DOX_ApRespiration: xr.DataArray,
    DOX_Nitrification: xr.DataArray,
    DOX_DOC_Oxidation: xr.DataArray,
    DOX_CBOD_Oxidation: xr.DataArray,
    DOX_AbGrowth: xr.DataArray,
    DOX_AbRespiration: xr.DataArray,
    DOX_SOD: xr.DataArray
) -> xr.DataArray:
    """Compute change in dissolved oxygen concentration for one timestep

    Args:
        Atm_O2_reaeration: DOX concentration change due to atmospheric O2 reaeration (mg/L/d)
        DOX_ApGrowth: DOX concentration change due to algal photosynthesis (mg/L/d)
        DOX_ApRespiration: DOX concentration change due to algal respiration (mg/L/d)
        DOX_Nitrification: DOX concentration change due to nitrification (mg/L/d)
        DOX_DOC_Oxidation: DOX concentration change due to DOC oxidation (mg/L/d)
        DOX_CBOD_Oxidation: DOX concentration change due to CBOD oxidation (mg/L/d)
        DOX_AbGrowth: DOX concentration change due to benthic algae photosynthesis (mg/L/d)
        DOX_AbRespiration: DOX concentration change due to benthic algae respiration (mg/L/d)
        DOX_SOD: DOX concentration change due to sediment oxygen demand (mg/L/d)
    """
    return Atm_O2_reaeration + DOX_ApGrowth - DOX_ApRespiration - DOX_Nitrification - DOX_DOC_Oxidation - DOX_CBOD_Oxidation + DOX_AbGrowth - DOX_AbRespiration - DOX_SOD


@numba.njit
def DOX_new(
    DOX: xr.DataArray,
    dDOXdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Computes updated dissolved oxygen concentration

    Args:
        DOX: Dissolved oxygen concentration from previous timestep
        dDOXdt: Change in dissolved oxygen concentration over timestep
        timestep: Current iteration timestep (d)
    """
    return DOX + dDOXdt * timestep
