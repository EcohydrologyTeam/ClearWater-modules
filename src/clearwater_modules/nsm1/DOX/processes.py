"""
File contains dynamic variables related to the DOX module
"""

import numba
import xarray as xr
from clearwater_modules.shared.processes import arrhenius_correction
import math


#TODO: make sure np.exp will work here...
@numba.njit
def pwv(
    TwaterK: xr.DataArray
) -> xr.DataArray:
    """Calculate partial pressure of water vapor

    Args:
        TwaterK: Water temperature kelvin
    """
    return np.exp(11.8571 - 3840.70 / TwaterK - 216961 / TwaterK ** 2)


@numba.njit
def DOs_atm_alpha(
    TwaterK: xr.DataArray
) -> xr.DataArray:
    """Calculate DO saturation atmospheric correction coefficient

    Args:
        TwaterK: Water temperature kelvin
    """
    return .000975 - 1.426 * 10 ** -5 * TwaterK + 6.436 * 10 ** -8 * TwaterK ** 2


@numba.njit
def DOX_sat(
    TwaterK: xr.DataArray,
    pressure_atm: xr.DataArray,
    pwv: xr.DataArray,
    DOs_atm_alpha: xr.DataArray
) -> xr.DataArray:
    """Calculate DO saturation value

    Args:
        TwaterK: Water temperature kelvin
        pressure_atm: Atmospheric pressure (atm)
        pwv: Patrial pressure of water vapor (atm)
        DOs_atm_alpha: DO saturation atmospheric correction coefficient
    """
    DOX_sat_uncorrected = np.exp(-139.34410 + 1.575701 * 10 ** 5 / TwaterK - 6.642308 * 10 ** 7 / TwaterK ** 2
                                 + 1.243800 * 10 ** 10 / TwaterK - 8.621949 * 10 ** 11 / TwaterK)

    DOX_sat_corrected = DOX_sat_uncorrected * pressure_atm * \
        (1 - pwv / pressure_atm) * (1 - DOs_atm_alpha * pressure_atm) / \
        ((1 - pwv) * (1 - DOs_atm_alpha))
    return DOX_sat_corrected

#TODO potentially move kah_20 to a global file; found in global parameters
def kah_20(
    kah_20_user: xr.DataArray,
    hydraulic_reaeration_option: xr.DataArray,
    velocity: xr.DataArray,
    depth: xr.DataArray,
    flow: xr.DataArray,
    topwidth: xr.DataArray,
    slope: xr.DataArray,
    shear_velocity: xr.DataArray
) -> xr.DataArray:
    """Calculate hydraulic oxygen reaeration rate based on flow parameters in different cells

    Args:
        kah_20_user: User defined O2 reaeration rate at 20 degrees (1/d)
        hydraulic_reaeration_option: Integer value which selects method for computing O2 reaeration rate 
        velocity: Average water velocity in cell (m/s)
        depth: Average water depth in cell (m)
        flow: Average flow rate in cell (m3/s)
        topwidth: Average topwidth of cell (m)
        slope: Average slope of bottom surface 
        shear_velocity: Average shear velocity on bottom surface (m/s)
    """

    da: xr.DataArray = xr.where(hydraulic_reaeration_option == 1, kah_20_user,
                        xr.where(hydraulic_reaeration_option == 2, (3.93 * velocity**0.5) / (depth**1.5),
                        xr.where(hydraulic_reaeration_option == 3, (5.32 * velocity**0.67) / (depth**1.85),
                        xr.where(hydraulic_reaeration_option == 4, (5.026 * velocity) / (depth**1.67),
                        xr.where(hydraulic_reaeration_option == 5, xr.where(depth < 0.61, (5.32 * velocity**0.67) / (depth**1.85), xr.where(depth > 0.61, (3.93 * velocity**0.5) / (depth**1.5), (5.026 * velocity) / (depth**1.67))),
                        xr.where(hydraulic_reaeration_option == 6, xr.where(flow < 0.556, 517 * (velocity * slope)**0.524 * flow**-0.242, 596 * (velocity * slope)**0.528 * flow**-0.136),
                        xr.where(hydraulic_reaeration_option == 7, xr.where(flow < 0.556, 88 * (velocity * slope)**0.313 * depth**-0.353, 142 * (velocity * slope)**0.333 * depth**-0.66 * topwidth**-0.243),
                        xr.where(hydraulic_reaeration_option == 8, xr.where(flow < 0.425, 31183 * velocity * slope, 15308 * velocity * slope),
                        xr.where(hydraulic_reaeration_option == 9, 2.16 * (1 + 9 * (velocity / (9.81 * depth)**0.5)**0.25) * shear_velocity / depth, -9999
                                 )))))))))
    return da


@numba.njit
def kah_tc(
    TwaterC: xr.DataArray,
    kah_20: xr.DataArray,
) -> xr.DataArray:
    """Calculate the temperature adjusted hydraulic oxygen reaeration rate (/d)

    Args:
        TwaterC: Water temperature in Celsius
        kah_20: Hydraulic oxygen reaeration rate at 20 degrees Celsius
    """
    return arrhenius_correction(TwaterC, kah_20, 1.024)


def kaw_20(
    kaw_20_user: xr.DataArray,
    wind_speed: xr.DataArray,
    wind_reaeration_option: xr.DataArray
) -> xr.DataArray:
    """Calculate the wind oxygen reaeration velocity (m/d) based on wind speed, r stands for regional

    Args:
        kaw_20_user: User defined wind oxygen reaeration velocity at 20 degrees C (m/d)
        wind_speed: Wind speed at 10 meters above the water surface (m/s)
        wind_reaeration_option: Integer value which selects method for computing wind oxygen reaeration velocity
    """
    Uw10 = wind_speed * (10 / 2)**0.143

    da: xr.DataArray = xr.where(wind_reaeration_option == 1, kaw_20_user,
                        xr.where(wind_reaeration_option == 2, 0.864 * Uw10,
                        xr.where(wind_reaeration_option == 3, xr.where(Uw10 <= 3.5, 0.2 * Uw10, 0.057 * Uw10**2),
                        xr.where(wind_reaeration_option == 4, 0.728 * Uw10**0.5 - 0.317 * Uw10 + 0.0372 * Uw10**2,
                        xr.where(wind_reaeration_option == 5, 0.0986 * Uw10**1.64,
                        xr.where(wind_reaeration_option == 6, 0.5 + 0.05 * Uw10**2,
                        xr.where(wind_reaeration_option == 7, xr.where(Uw10 <= 5.5, 0.362 * Uw10**0.5, 0.0277 * Uw10**2),
                        xr.where(wind_reaeration_option == 8, 0.64 + 0.128 * Uw10**2,
                        xr.where(wind_reaeration_option == 9, xr.where(Uw10 <= 4.1, 0.156 * Uw10**0.63, 0.0269 * Uw10**1.9),
                        xr.where(wind_reaeration_option == 10, 0.0276 * Uw10**2,
                        xr.where(wind_reaeration_option == 11, 0.0432 * Uw10**2,
                        xr.where(wind_reaeration_option == 12, 0.319 * Uw10,
                        xr.where(wind_reaeration_option == 13, xr.where(Uw10 < 1.6, 0.398, 0.155 * Uw10**2), -9999
                                 )))))))))))))
    
    return da


@numba.njit
def kaw_tc(
    TwaterC: xr.DataArray,
    kaw_20: xr.DataArray,
) -> xr.DataArray:
    """Calculate the temperature adjusted wind oxygen reaeration velocity (m/d)

    Args:
        TwaterC: Water temperature in Celsius
        kaw_20: Wind oxygen reaeration velocity at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(TwaterC, kaw_20, 1.024)


@numba.njit
def ka_tc(
    kah_tc: xr.DataArray,
    kaw_tc: xr.DataArray,
    depth: xr.DataArray
) -> xr.DataArray:
    """Compute the oxygen reaeration rate, adjusted for temperature (1/d)

    Args:
        kah_tc: Oxygen reaeration rate adjusted for temperature (1/d)
        kaw_tc: Wind oxygen reaeration velocity adjusted for temperature (m/d)
        depth: Average water depth in cell (m)
    """
    return kaw_tc / depth + kah_tc


@numba.njit
def Atm_O2_reaeration(
    ka_tc: xr.DataArray,
    DOX_sat: xr.DataArray,
    DOX: xr.DataArray
) -> xr.DataArray:
    """Compute the atmospheric O2 reaeration flux

    Args: 
        ka_tc: Oxygen reaeration rate adjusted for temperature (1/d)
        DOX_sat: Dissolved oxygen saturation concentration (mg/L)
        DOX: Dissolved oxygen concentration (mg/L)
    """
    return ka_tc * (DOX_sat - DOX)


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
        ApUptakeFr_NH4: Algae preference for ammonia 
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
        AbUptakeFr_NH4: Fraction of actual benthic algal uptake that is form the ammonia pool, calculated in nitrogen module
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

#TODO potentially move to global parameter 
def SOD_tc(
    SOD_20: xr.DataArray,
    TwaterC: xr.DataArray,
    DOX: xr.DataArray,
    KsSOD: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Compute the sediment oxygen demand corrected by temperature and dissolved oxygen concentration

    Args:
        SOD_20: Sediment oxygen demand at 20 degrees celsius (mg-O2/m2)
        TwaterC: Water temperature in degrees C
        KsSod: Half saturation oxygen attenuation constant for SOD (mg-O/L)
        use_DOX: Option to consider DOX concentration in water in calculation of sediment oxygen demand
    """
    SOD_tc = arrhenius_correction(TwaterC, SOD_20, 1.060)

    da: xr.DataArray = xr.where(use_DOX == True, SOD_tc * DOX / (DOX + KsSOD), SOD_tc)

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
def DOX(
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
