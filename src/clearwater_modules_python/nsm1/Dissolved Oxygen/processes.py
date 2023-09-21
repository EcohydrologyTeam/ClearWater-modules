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