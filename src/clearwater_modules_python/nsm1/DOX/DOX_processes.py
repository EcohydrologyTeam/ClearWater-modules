import numpy as np
import numba
from clearwater_modules.shared.processes import (
    arrhenius_correction,
)


@numba.njit
def pwv(
    t_water_k: float
) -> float:
    """Calculate partial pressure of water vapor

    Args:
        t_water_k: Water temperature kelvin
    """
    return np.exp(11.8571 - 3840.70 / t_water_k - 216961 / t_water_k ** 2)


@numba.njit
def DOs_atm_alpha(
    t_water_k: float
) -> float:
    """Calculate DO saturation atmospheric correction coefficient

    Args:
        t_water_k: Water temperature kelvin
    """
    return .000975 - 1.426 * 10 ** -5 * t_water_k + 6.436 * 10 ** -8 * t_water_k ** 2


@numba.njit
def DOX_sat(
    t_water_k: float,
    patm: float,
    pwv: float,
    DOs_atm_alpha: float
) -> float:
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
def kah_20(
    kah_20_user: float,
    hydraulic_reaeration_option: int,
    velocity: float,
    depth: float,
    flow: float,
    topwidth: float,
    slope: float,
    shear_velocity: float
) -> float:
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
    if hydraulic_reaeration_option == 1:
        return kah_20_user

    elif hydraulic_reaeration_option == 2:
        return (3.93 * velocity**0.5) / (depth**1.5)

    elif hydraulic_reaeration_option == 3:
        return (5.32 * velocity**0.67) / (depth**1.85)

    elif hydraulic_reaeration_option == 4:
        return (5.026 * velocity) / (depth**1.67)

    elif hydraulic_reaeration_option == 5:
        if depth < 0.61:
            return (5.32 * velocity**0.67) / (depth**1.85)
        elif depth > 0.61:
            return (3.93 * velocity**0.5) / (depth**1.5)
        else:
            return (5.026 * velocity) / (depth**1.67)

    elif hydraulic_reaeration_option == 6:
        if flow < 0.556:
            return 517 * (velocity * slope)**0.524 * flow**-0.242
        elif flow >= 0.556:
            return 596 * (velocity * slope)**0.528 * flow**-0.136

    elif hydraulic_reaeration_option == 7:
        if flow < 0.556:
            return 88 * (velocity * slope)**0.313 * depth**-0.353
        elif flow >= 0.556:
            return 142 * (velocity * slope)**0.333 * depth**-0.66 * topwidth**-0.243

    elif hydraulic_reaeration_option == 8:
        if flow < 0.425:
            return 31183 * velocity * slope
        elif flow >= 0.425:
            return 15308 * velocity * slope

    elif hydraulic_reaeration_option == 9:
        froude = velocity / (9.81 * depth)**0.5
        return 2.16 * (1 + 9 * froude**0.25) * shear_velocity / depth

    else:
        raise ValueError('Hydraulic Reaeration Option not properly selected')


@numba.njit
def kah_T(
    water_temp_c: float,
    kah_20: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted hydraulic oxygen reaeration rate (/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kah_20: Hydraulic oxygen reaeration rate at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kah_20, theta)


@numba.njit
def kaw_20(
    kaw_20_user: float,
    wind_speed: float,
    wind_reaeration_option: int
) -> float:
    """Calculate the wind oxygen reaeration velocity (m/d) based on wind speed, r stands for regional

    Args:
        kaw_20_user: User defined wind oxygen reaeration velocity at 20 degrees C (m/d)
        wind_speed: Wind speed at 10 meters above the water surface (m/s)
        wind_reaeration_option: Integer value which selects method for computing wind oxygen reaeration velocity
    """
    Uw10 = wind_speed * (10 / 2)**0.143

    if wind_reaeration_option == 1:
        return kaw_20_user

    elif wind_reaeration_option == 2:
        return 0.864 * Uw10

    elif wind_reaeration_option == 3:
        if Uw10 <= 3.5:
            return 0.2 * Uw10
        else:
            return 0.057 * Uw10**2

    elif wind_reaeration_option == 4:
        return 0.728 * Uw10**0.5 - 0.317 * Uw10 + 0.0372 * Uw10**2

    elif wind_reaeration_option == 5:
        return 0.0986 * Uw10**1.64

    elif wind_reaeration_option == 6:
        return 0.5 + 0.05 * Uw10**2

    elif wind_reaeration_option == 7:
        if Uw10 <= 5.5:
            return 0.362 * Uw10**0.5
        else:
            return 0.0277 * Uw10**2

    elif wind_reaeration_option == 8:
        return 0.64 + 0.128 * Uw10**2

    elif wind_reaeration_option == 9:
        if Uw10 <= 4.1:
            return 0.156 * Uw10**0.63
        else:
            return 0.0269 * Uw10**1.9

    elif wind_reaeration_option == 10:
        return 0.0276 * Uw10**2

    elif wind_reaeration_option == 11:
        return 0.0432 * Uw10**2

    elif wind_reaeration_option == 12:
        return 0.319 * Uw10

    elif wind_reaeration_option == 13:
        if Uw10 < 1.6:
            return 0.398
        else:
            return 0.155 * Uw10**2

    else:
        raise ValueError('Wind Reaeration Option not properly selected')


@numba.njit
def kaw_T(
    water_temp_c: float,
    kaw_20_r: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted wind oxygen reaeration velocity (m/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kaw_20: Wind oxygen reaeration velocity at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kaw_20, theta)


@numba.njit
def ka_T(
    kah_T: float,
    kaw_T: float,
    depth: float
) -> float:
    """Compute the oxygen reaeration rate, adjusted for temperature (1/d)

    Args:
        kah_T: Oxygen reaeration rate adjusted for temperature (1/d)
        kaw_T: Wind oxygen reaeration velocity adjusted for temperature (m/d)
        depth: Average water depth in cell (m)
    """
    return kaw_T / depth + kah_T


@numba.njit
def Atm_O2_reaeration(
    ka_T: float,
    DOX_sat: float,
    DOX: float
) -> float:
    """Compute the atmospheric O2 reaeration flux

    Args: 
        ka_T: Oxygen reaeration rate adjusted for temperature (1/d)
        DOX_sat: Dissolved oxygen saturation concentration (mg/L)
        DOX: Dissolved oxygen concentration (mg/L)
    """
    return ka_T * (DOX_sat - DOX)


@numba.njit
def DOX_ApGrowth(
    ApGrowth: float,
    rca: float,
    roc: float,
    F1: float,
    use_Algae: bool
) -> float:
    """Compute DOX flux due to algal photosynthesis

    Args:
        ApGrowth: Algae photosynthesis, calculated in the algae module (ug-Chla/L/d)
        rca: Ratio of algal carbon to chlorophyll-a (mg-C/ug-Chla)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
        F1: Algae preference for ammonia 
    """
    if use_Algae:
        return ApGrowth * rca * roc * (138 / 106 - 32 * F1 / 106)
    else:
        return 0


@numba.njit
def DOX_ApRespiration(
    ApRespiration: float,
    rca: float,
    roc: float,
    use_Algae: bool
) -> float:
    """Compute DOX flux due to algal photosynthesis

    Args:
        ApRespiration: algae respiration, calculated in the algae module
        rca: Ratio of algal carbon to chlorophyll-a (mg-C/ug-Chla)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C) 
    """
    if use_Algae:
        return ApRespiration * rca * roc
    else:
        return 0


@numba.njit
def DOX_Nitrification(
    KNR: float,
    DOX: float,
    ron: float,
    knit_tc: float,
    NH4: float,
    use_NH4: bool
) -> float:
    """Compute DOX flux due to nitrification of ammonia

    Args:
        KNR: Oxygen inhibition factor for nitrification (mg-O2/L)
        DOX: Dissolved oxygen concentration (mg/L)
        ron: Ratio of oxygen to nitrogen for nitrificiation (mg-O2/mg-N)
        knit_tc: Nitrification rate of NH4 to NO3 (1/d)
        NH4: Ammonia/ammonium concentration
    """
    if use_NH4:
        return (1.0 - np.exp(-KNR * DOX)) * ron * knit_tc * NH4
    else:
        return 0


@numba.njit
def DOX_DOC_Oxidation(
    DOC_Oxidation: float,
    roc: float,
    use_DOC: bool
) -> float:
    """Computes dissolved oxygen flux due to oxidation of dissolved organic carbon

    Args:
        DOC_Oxidation: Dissolved organic carbon oxidation, calculated in carbon module (mg/L/d)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
    """
    if use_DOC:
        return roc * DOC_Oxidation
    else:
        return 0


@numba.njit
def DOX_CBOD_Oxidation(
    DIC_CBOD_Oxidation: float,
    roc: float
) -> float:
    """Compute dissolved oxygen flux due to CBOD oxidation

    Args:
        DIC_CBOD_Oxidation: Carbonaceous biochemical oxygen demand oxidation, calculated in CBOD module (mg/L/d)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
    """
    return DIC_CBOD_Oxidation * roc


@numba.njit
def DOX_AbGrowth(
    AbUptakeFr_NH4: float,
    roc: float,
    rcb: float,
    AbGrowth: float,
    Fb: float,
    depth: float,
    use_BAlgae: bool
) -> float:
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
    if use_BAlgae:
        return (138 / 106 - 32 / 106 * AbUptakeFr_NH4) * roc * rcb * AbGrowth * Fb / depth
    else:
        return 0


@numba.njit
def DOX_AbRespiration(
    roc: float,
    rcb: float,
    AbRespiration: float,
    Fb: float,
    depth: float,
    use_BAlgae: bool
) -> float:
    """Compute dissolved oxygen flux due to benthic algae respiration

    Args:
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        AbRespiration: Benthic algae respiration, calculated in the benthic algae module
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Water depth (m)
        use_BAlgae: Option to consider benthic algae in the DOX budget
    """
    if use_BAlgae:
        return roc * rcb * AbRespiration * Fb / depth
    else:
        return 0


@numba.njit
def SOD_tc(
    SOD_20: float,
    t_water_C: float,
    theta: float,
    DOX: float,
    KsSOD: float,
    use_DOX: bool
) -> float:
    """Compute the sediment oxygen demand corrected by temperature and dissolved oxygen concentration

    Args:
        SOD_20: Sediment oxygen demand at 20 degrees celsius (mg-O2/m2)
        t_water_C: Water temperature in degrees C
        theta: Arrhenius coefficient
        use_DOX: Option to consider DOX concentration in water in calculation of sediment oxygen demand
    """
    SOD_tc = arrhenius_correction(t_water_C, SOD_20, theta)
    if use_DOX:
        return SOD_tc * DOX / (DOX + KsSOD)
    else:
        return SOD_tc


@numba.njit
def DOX_SOD(
    SOD_Bed: float,
    depth: float,
    SOD_tc: float,
    use_SedFlux: bool
) -> float:
    """Compute dissolved oxygen flux due to sediment oxygen demand

    Args:
        SOD_Bed: Sediment oxygen demand if calculated using the SedFlux module (mg-O2/m2)
        depth: Water depth (m)
        SOD_tc: Sediment oxygen demand not considering the SedFlux budget (mg-O2/m2)
        use_SedFlux: Option to consider sediment flux in DOX budget (boolean)
    """
    if use_SedFlux:
        return SOD_Bed / depth
    else:
        return SOD_tc / depth


@numba.njit
def dDOXdt(
    Atm_O2_reaeration: float,
    DOX_ApGrowth: float,
    DOX_ApRespiration: float,
    DOX_Nitrification: float,
    DOX_DOC_Oxidation: float,
    DOX_CBOD_Oxidation: float,
    DOX_AbGrowth: float,
    DOX_AbRespiration: float,
    DOX_SOD: float
) -> float:
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
    DOX: float,
    dDOXdt: float,
    timestep: float
) -> float:
    """Computes updated dissolved oxygen concentration

    Args:
        DOX: Dissolved oxygen concentration from previous timestep
        dDOXdt: Change in dissolved oxygen concentration over timestep
        timestep: Current iteration timestep (d)
    """
    return DOX + dDOXdt * timestep
