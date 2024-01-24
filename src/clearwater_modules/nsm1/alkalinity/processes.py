import numba
import math
from clearwater_modules.shared.processes import (
    arrhenius_correction
)
import xarray as xr
from clearwater_modules.nsm1.alkalinity import dynamic_variables
from clearwater_modules.nsm1.alkalinity import static_variables
from clearwater_modules.nsm1 import static_variables_global
from clearwater_modules.nsm1 import dynamic_variables_global
from clearwater_modules.nsm1 import state_variables



def Alk_denitrification(
    DOX: xr.DataArray,
    NO3: xr.DataArray,
    kdnit_tc: xr.DataArray,
    KsOxdn: xr.DataArray,
    r_alkden: xr.DataArray,
    use_NO3: xr.DataArray,
    use_DOX: xr.DataArray        
) -> xr.DataArray:
    """Calculate the alkalinity concentration change due to denitrification of nitrate
    
    Args:
        DOX: Concentration of dissolved oxygen (mg/L)
        NO3: Concentration of nitrate (mg/L)
        kdnit_T: Denitrification rate corrected for temperature (1/d)
        KsOxdn: Half-saturation oxygen inhibition constant for denitrification (mg-O2/L)
        ralkden: Ratio translating NO3 denitrification into Alk (eq/mg-N)
        use_NO3: Option to use nitrate
        use_DOX: Option to use dissolved oxygen 
    """
    da: xr.DataArray = xr.where(use_NO3 == True, xr.where(use_DOX == True, r_alkden * (1.0 - (DOX / (DOX + KsOxdn))) * kdnit_tc * NO3, r_alkden * kdnit_tc * NO3), 0)

    return da


def Alk_nitrification(
    DOX: xr.DataArray,
    NH4: xr.DataArray,
    knit_tc: xr.DataArray,
    KNR: xr.DataArray,
    r_alkn: xr.DataArray,
    use_NH4: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Calculate the alkalinity concentration change due to nitrification of ammonium

    Args:
        DOX: Concentration of dissolved oxygen (mg/L)
        NH4: Concentration of ammonia/ammonium (mg/L)
        knit_T: Nitrification rate corrected for temperature (1/d)
        KNR: Oxygen inhibition factor for nitrification (mg-O2/L)
        r_alkn: Ratio translating NH4 nitrification into Alk (eq/mg-N)
        use_NH4: Option to use ammonium
        use_DOX: Option to use dissolved oxygen
    """
    da: xr.DataArray = xr.where(use_NH4 == True, xr.where(use_DOX == True, r_alkn * (1 - math.exp(-KNR * DOX)) * knit_tc * NH4, knit_tc * NH4), 0)

    return da


def Alk_algal_growth(
    ApGrowth: xr.DataArray,
    r_alkaa: xr.DataArray,
    r_alkan: xr.DataArray,
    ApUptakeFr_NH4: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculate the alkalinity concentration change due to algal growth

    Args:
        ApGrowth: Algal photosynthesis calculated in algae module (ug-Chla/L/d)
        r_alkaa: Ratio translating algal growth into Alk if NH4 is the N source (eq/ug-Chla)
        r_alkan: Ratio translating algal growth into Alk if NO3 is the N source (eq/ug-Chla)
        ApUptakeFr_NH4: Preference fraction of algal N uptake from NH4
        use_Algae: Option to use algae
    """
    da: xr.DataArray = xr.where(use_Algae == True, (r_alkaa * ApUptakeFr_NH4 - r_alkan * (1 - ApUptakeFr_NH4)) * ApGrowth, 0)

    return da


def Alk_algal_respiration(
    ApRespiration: xr.DataArray,
    r_alkaa: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculate the alkalinity concentration change due to algal respiration

    Args:
        ApRespiration: Algae respiration calculated in algae module (ug-Chla/L/d)
        r_alkaa: Ratio translating algal growth into Alk if NH4 is the N source (eq/ug-Chla)
        use_Algae: Option to use algae
    """
    da: xr.DataArray = xr.where(use_Algae == True, ApRespiration * r_alkaa, 0)

    return da


def Alk_benthic_algae_growth(
    AbGrowth: xr.DataArray,
    depth: xr.DataArray,
    r_alkba: xr.DataArray,
    r_alkbn: xr.DataArray,
    AbUptakeFr_NH4: xr.DataArray,
    Fb: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculate the alkalinity concentration change due to algal growth

    Args:
        ApGrowth: Algal photosynthesis calculated in algae module (ug-Chla/L/d)
        depth: Depth of water (m)
        r_alkaa: Ratio translating algal growth into Alk if NH4 is the N source (eq/ug-Chla)
        r_alkan: Ratio translating algal growth into Alk if NO3 is the N source (eq/ug-Chla)
        AbUptakeFr_NH4: Preference fraction of benthic algae N uptake from NH4
        Fb: Fraction of bottom area available for benthic algae growth
        use_Balgae: Option to use benthic algae
    """
    da: xr.DataArray = xr.where(use_Balgae == True, (1 / depth) * (r_alkba * AbUptakeFr_NH4 - r_alkbn * (1 - AbUptakeFr_NH4)) * AbGrowth * Fb, 0)

    return da


def Alk_benthic_algae_respiration(
    AbRespiration: xr.DataArray,
    depth: xr.DataArray,
    r_alkba: xr.DataArray,
    Fb: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculate the alkalinity concentration change due to algal respiration

    Args:
        ApRespiration: Algae respiration calculated in algae module (ug-Chla/L/d)
        r_alkaa: Ratio translating algal growth into Alk if NH4 is the N source (eq/ug-Chla)
        Fb: Fraction of bottom area available for benthic algae growth
        use_Balgae: Option to use betnhic algae
    """
    da: xr.DataArray = xr.where(use_Balgae == True, (1 / depth) * r_alkba * AbRespiration * Fb, 0)

    return da


@numba.njit
def dAlkdt(
    Alk_denitrification: xr.DataArray,
    Alk_nitrification: xr.DataArray,
    Alk_algal_growth: xr.DataArray,
    Alk_algal_respiration: xr.DataArray,
    Alk_benthic_algae_growth: xr.DataArray,
    Alk_benthic_algae_respiration: xr.DataArray
) -> xr.DataArray:
    """Computes the change in alkalinity for timestep

    Args:
        Alk_denitrification: xr.DataArray,
        Alk_nitrification: xr.DataArray,
        Alk_algal_growth: xr.DataArray,
        Alk_algal_respiration: xr.DataArray,
        Alk_benthic_algae_growth: xr.DataArray,
        Alk_benthic_algae_respiration: xr.DataArray
    """
    return Alk_denitrification - Alk_nitrification - Alk_algal_growth + Alk_algal_respiration - Alk_benthic_algae_growth + Alk_benthic_algae_respiration


@numba.njit
def Alk_new(
    Alk: xr.DataArray,
    dAlkdt: xr.DataArray,
    timestep: xr.DataArray,
) -> xr.DataArray:
    """Computes the alkalinity concentration at the next timestep

    Args:
        Alk: Concentration of alkalinity from previous timestep (mg/L)
        dAlkdt: Change in concentration of alkalinity for current timestep (mg/L/d)
        timestep: Current iteration timestep (d)
    """
    return Alk + dAlkdt * timestep