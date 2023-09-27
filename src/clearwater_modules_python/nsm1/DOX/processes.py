import numpy as np
from clearwater_modules_python.shared.processes import (
    arrhenius_correction,
)     

def pwv(
    t_water_k: float
) -> float:
    """Calculate partial pressure of water vapor

    Args:
        t_water_k: water temperature kelvin
    """
    return np.exp(11.8571 - 3840.70 / t_water_k - 216961 / t_water_k **2)

def DOs_atm_alpha(
    t_water_k: float
) -> float:
    """Calculate DO saturation atmospheric correction coefficient
    
    Args:
        t_water_k: water temperature kelvin
    """
    return .000975 - 1.426 * 10 ** -5 * t_water_k + 6.436 * 10 ** -8 * t_water_k ** 2

def DOX_sat(
    t_water_k: float,
    patm: float,
    pwv: float,
    DOs_atm_alpha: float
) -> float:
    """Calculate DO saturation value
    
    Args:
        t_water_k: water temperature kelvin
        patm:
        pwv:
        DOs_atm_alpha:
    """
    DOX_sat_uncorrected = np.exp(-139.34410 + 1.575701 * 10 ** 5 / t_water_k - 6.642308 * 10 ** 7 / t_water_k ** 2 
                  + 1.243800 * 10 ** 10 / t_water_k - 8.621949 * 10 ** 11 / t_water_k)

    DOX_sat_corrected = DOX_sat_uncorrected * patm * (1 - pwv / patm) * (1 - DOs_atm_alpha * patm) / ((1 - pwv) * (1 - DOs_atm_alpha))
    return DOX_sat_corrected

def kah_20_r(
    kah_20: float,
    hydraulic_reaeration_option: int,
    velocity: float,
    depth: float,
    flow: float,
    topwidth: float,
    slope: float,
    shear_velocity: float
) -> float:
    """Calculate hydraulic oxygen reaeration rate based on flow parameters in different cells, r stands for regional
    
    Args:
        kah_20:
        hydraulic_reaeration_option:
        velocity:
        depth:
        flow:
        topwidth:
        slope:
        shear_velocity:
    """
    if hydraulic_reaeration_option == 1:
        return kah_20
    
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

def kah_T_r(
    water_temp_c: float,
    kah_20_r: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted hydraulic oxygen reaeration rate (/d)

    Args:
        water_temp_c: water temperature in Celsius
        kah_20: hydraulic oxygen reaeration rate at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kah_20_r, theta)

def kaw_20_r(
    kaw_20: float,
    wind_speed: float,
    wind_reaeration_option: int
) -> float:
    """Calculate the wind oxygen reaeration velocity (m/d) based on wind speed, r stands for regional
    
    Args:
        kaw_20:
        wind_speed:
        wind_reaeration_option:
    """
    Uw10 = wind_speed * (10 / 2)**0.143
    
    if wind_reaeration_option == 1:
        return kaw_20
    
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
        if Uw10 <=  4.1:
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

def kaw_T_r(
    water_temp_c: float,
    kaw_20_r: float,
    theta: float
) -> float:
    """Calculate the temperature adjusted wind oxygen reaeration velocity (m/d)

    Args:
        water_temp_c: water temperature in Celsius
        kaw_20: wind oxygen reaeration velocity at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kaw_20_r, theta)

def ka_T(
    kah_T_r: float,
    kaw_T_r: float,
    h: float
) -> float:
    """Compute the oxygen reaeration rate, adjusted for temperature (/d)

    Args:
        ka_T_r:
        kaw_T_r:
        h:
    """
    return kaw_T_r / h + kah_T_r

def Atm_O2_reaeration(
    ka_T: float,
    DOX_sat: float,
    DOX: float
) -> float:
    """Compute the atmospheric O2 reaeration flux
    
    Args: 
        ka_T:
        DOX_sat:
        DOX:
    """
    return ka_T * (DOX_sat - DOX)

def DOX_ApGrowth(
    ApGrowth: float,
    rca: float,
    roc: float,
    F1: float,
    use_Algae: bool
) -> float:
    """Compute DOX flux due to algal photosynthesis

    Args:
        ApGrowth: algae growth computed in the algae module
        rca: ratio of algal carbon to chlorophyll-a
        roc: stoichiometric ratio of oxygen to carbon
        F1: algae preference for ammonia 
    """
    if use_Algae:
        return ApGrowth * rca * roc * (138 / 106 - 32 * F1 / 106)
    else:
        return 0

def DOX_ApRespiration(
    ApRespiration: float,
    rca: float,
    roc: float,
    use_Algae: bool
) -> float:
    """Compute DOX flux due to algal photosynthesis

    Args:
        ApRespiration: algae respiration computed in the algae module
        rca: ratio of algal carbon to chlorophyll-a
        roc: stoichiometric ratio of oxygen to carbon 
    """
    if use_Algae:
        return ApRespiration * rca * roc
    else:
        return 0

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
        KNR:
        DOX:
        ron:
        knit_tc:
        NH4:
    """
    if use_NH4:
        return  1.0 - np.exp(-KNR * DOX) * ron * knit_tc * NH4
    else:
        return 0

def DOX_DOC_Oxidation(
    DOC_Oxidation: float,
    roc: float,
    use_DOC: bool
) -> float:
    """Computes dissolved oxygen flux due to oxidation of dissolved organic carbon
    
    Args:
        DOC_Oxidation: 
        roc: 
    """
    if use_DOC:
        return roc * DOC_Oxidation
    else:
        return 0

def DOX_CBOD_Oxidation(
    DIC_CBOD_Oxidation: float,
    roc: float
) -> float:
    """Compute dissolved oxygen flux due to CBOD oxidation
    
    Args:
        DIC_CBOD_Oxidation:
        roc:
    """
    return DIC_CBOD_Oxidation * roc

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
        AbUptakeFr_NH4:
        roc:
        rcb:
        AbGrowth:
        Fb:
        depth:
        use_BAlgae:
    """
    if use_BAlgae:
        return (138 / 106 - 32 / 106 * AbUptakeFr_NH4) * roc * rcb * AbGrowth * Fb / depth 
    else:
        return 0
    
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
        roc:
        rcb:
        AbRespiration:
        Fb:
        depth:
        use_BAlgae:
    """
    if use_BAlgae:
        return roc * rcb * AbRespiration * Fb / depth
    else:
        return 0
    
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
        SOD_20:
        t_water_C:
        theta:
        use_DOX:
    """
    SOD_tc = arrhenius_correction(t_water_C, SOD_20, theta)
    if use_DOX:
        return SOD_tc * DOX / (DOX + KsSOD)
    else:
        return SOD_tc

def DOX_SOD(
    SOD_Bed: float,
    depth: float,
    SOD_tc: float,
    use_SedFlux: bool
) -> float:
    """Compute dissolved oxygen flux due to sediment oxygen demand
    
    Args:
        SOD_Bed:
        depth:
        SOD_tc:
        use_SedFlux:
    """
    if use_SedFlux:
        return SOD_Bed / depth
    else:
        return SOD_tc / depth
    
def DOX_change(
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
        Atm_O2_reaeration:
        DOX_ApGrowth:
        DOX_ApRespiration:
        DOX_Nitrification:
        DOX_DOC_Oxidation:
        DOX_CBOD_Oxidation:
        DOX_AbGrowth:
        DOX_AbRespiration:
        DOX_SOD: 
    """
    return Atm_O2_reaeration + DOX_ApGrowth - DOX_ApRespiration - DOX_Nitrification - DOX_DOC_Oxidation - DOX_CBOD_Oxidation + DOX_AbGrowth - DOX_AbRespiration - DOX_SOD

def update_DOX(
    DOX: float,
    dDOXdt: float
) -> float:
    """Computes updated dissolved oxygen concentration

    Args:
        DOX: Dissolved oxygen concentration from previous timestep
        dDOXdt: Change in dissolved oxygen concentration over timestep
    """
    return DOX + dDOXdt