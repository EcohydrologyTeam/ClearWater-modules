"""
File contains process to calculate new benthic algae biomass concentration and associated dependent variables
"""

import math
from clearwater_modules.shared.processes import arrhenius_correction
import numba
import xarray as xr


@numba.njit
def mub_max_tc(
    mub_max_20: xr.DataArray,
    TwaterC: xr.DataArray
) -> xr.DataArray:
    """Calculate mub_max_tc: Maximum benthic algal growth rate with temperature correction (1/d).

    Args:
        mu_max_20: Maximum benthic algal growth rate at 20C (1/d)
        TwaterC: Water temperature (C)
    """
    return arrhenius_correction(TwaterC, mub_max_20, 1.047)


@numba.njit
def krb_tc(
    krb_20: xr.DataArray,
    TwaterC: xr.DataArray
) -> xr.DataArray:
    """Calculate krb_tc: Benthic algae respiration rate with temperature correction (1/d).

    Args:
        krb_20: Benthic algae respiration rate at 20C (1/d)
        TwaterC: Water temperature (C)
    """
    return arrhenius_correction(TwaterC, krb_20, 1.06)


@numba.njit
def kdb_tc(
    kdb_20: xr.DataArray,
    TwaterC: xr.DataArray
) -> xr.DataArray:
    """Calculate kdb_tc: Benthic algae mortality rate with temperature correction (1/d).

    Args:
        kdb_20: Benthic algae mortality rate at 20C (1/d)
        TwaterC: Water temperature (C)
    """
    return arrhenius_correction(TwaterC, kdb_20, 1.047)


@numba.njit
def rnb(
    BWn: xr.DataArray,
    BWd: xr.DataArray
) -> xr.DataArray:
    """Calculate rnb (mg-N/mg-D).

    Args:
        BWn: Benthic algae nitrogen (unitless)
        BWd: Benthic algae dry weight (unitless)
    """
    return BWn/BWd


@numba.njit
def rpb(
    BWp: xr.DataArray,
    BWd: xr.DataArray
) -> xr.DataArray:
    """Calculate rpd: Benthic algae phosphorus to dry weight ratio (mg-P/mg-D).

    Args:
        BWp: Benthic algae phosphorus (mg-P)
        BWd: Benthic algae dry weight (mg-D)
    """
    return BWp/BWd


@numba.njit
def rcb(
    BWc: xr.DataArray,
    BWd: xr.DataArray
) -> xr.DataArray:
    """Calculate rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D).

    Args:
        BWc: Benthic algae carbon (mg-C)
        BWd: Benthic algae dry weight (mg-D)
    """
    return BWc/BWd


@numba.njit
def rab(
    BWa: xr.DataArray,
    BWd: xr.DataArray
) -> xr.DataArray:
    """Calculate rab: Benthic algae chlorophyll-a to dry weight ratio (ug-Chla-a/mg-D).

    Args:
        BWa: Benthic algae chlorophyll-a (ug-Chla-a)
        BWd: Benthic algae dry weight (mg-D)
    """
    return BWa/BWd


@numba.njit
def FLb(
    L: xr.DataArray,
    depth: xr.DataArray,
    Ab: xr.DataArray,
    PAR: xr.DataArray,
    b_light_limitation_option: xr.DataArray,
    KLb: xr.DataArray,

) -> xr.DataArray:
    """Calculate Benthic algal light limitation: FLb (unitless).

    Args:
        L: Lambda light attenuation coefficient (unitless)
        depth: Water depth (m)
        Ab: Benthic algae Concentration (g/m^2)
        PAR: Surface light intensity (W/m^2)
        b_light_limitation_option: Benthic light limitation option 1) Half-saturation, 2) Smith model, 3) Steele model (unitless)
        KLb: Light limitation constant for benthic algal growth (W/m^2) 
    """

    # Note that KENT is defined differently here than it was for the algal equations.
    # The equations are different, this expression is more convenient here.
    KEXT = math.exp(-L*depth)

    FLb = xr.where(Ab <= 0.0 or KEXT <= 0.0 or PAR <= 0.0, 0.0,
          xr.where(b_light_limitation_option == 1, PAR * KEXT / (KLb + PAR * KEXT),
          xr.where(b_light_limitation_option == 2, PAR * KEXT / ((KLb**2.0 + (PAR * KEXT)**2.0)**0.5),
          xr.where(b_light_limitation_option == 3, 
          xr.where(abs(KLb) < 1.0E-10, 0.0, PAR * KEXT / KLb * math.exp(1.0 - PAR * KEXT / KLb)), "NaN"
          ))))
    FLb = xr.where(FLb > 1.0, 1.0,
          xr.where(FLb < 0.0, 0.0, FLb))

    return FLb


@numba.njit
def FNb(
    use_NH4: bool,
    use_NO3: bool,
    NH4: xr.DataArray,
    NO3: xr.DataArray,
    KsNb: xr.DataArray,

) -> xr.DataArray:
    """Calculate Benthic algae nitrogen limitation: FNb (unitless).

    Args:
        use_NH4: Use NH4 module true or false (true/false)
        use_NO3: Use NO3 module true or false (true/false)
        NH4: Ammonium concentration (mg-N/L)
        NO3: Nitrate concentration (mg-N/L)
        KsNb: Michaelis-Menton half-saturation constant relating inorganic N to benthic algal growth (mg-N/L)
    """
    FNb = xr.where(use_NH4 or use_NO3, (NH4 + NO3) / (KsNb + NH4 + NO3),1)
    FNb = xr.where(math.isnan(FNb),0.0,
          xr.where(FNb < 1.0, 1, FNb))

    return FNb


@numba.njit
def FPb(
    fdp: xr.DataArray,
    TIP: xr.DataArray,
    use_TIP: bool,
    KsPb: xr.DataArray
) -> xr.DataArray:
    """Calculate benthic algae phosphorous limitation: FPb (unitless).

    Args:
        use_TIP: Use Total Inorganic Phosphorus module true or false (true/false)
        TIP: Total Inorganic Phosphorus concentration (mg-P/L)
        KsPb: Michaelis-Menton half-saturation constant relating inorganic P to benthic algal growth (mg-P/L)
        fdp: Fraction P dissolved (unitless)
    """

    FPb = xr.where(use_TIP, fdp * TIP / (KsPb + fdp * TIP),1.0)
    FPb = xr.where(math.isnan(FPb),0.0,
          xr.where(FPb > 1.0, 1.0, FPb))

    return FPb


@numba.njit
def FSb(
    Ab: xr.DataArray,
    Ksb: xr.DataArray,

) -> xr.DataArray:
    """Calculate benthic density attenuation (unitless)

    Args:
        Ab: Benthic algae concentration (g/m^2)
        Ksb: Half-saturation density constant for benthic algae growth (g-D/m^2)

    """

    FSb = 1.0 - (Ab / (Ab + Ksb))
    FSb = xr.where(math.isnan(FSb), 1.0, 
          xr.where(FSb > 1.0, 1.0, FSb))
    
    return FSb


@numba.njit
def mub(
    mub_max_tc: xr.DataArray,
    b_growth_rate_option: int,
    FLb: xr.DataArray,
    FPb: xr.DataArray,
    FNb: xr.DataArray,
    FSb: xr.DataArray,

) -> xr.DataArray:
    """Calculate benthic algae specific growth rate (1/d)

    Args:
        mub_max_tc: Maximum benthic algal growth rate with temperature correction (1/d)
        b_growth_rate_option: Benthic Algal growth rate with three options 1) Multiplicative, 2) Limiting Nutrient
        FLb: Benethic algal light limitation (unitless)
        FPb: Benthic algae phosphorous limitation (unitless)
        FNb: Benthic algae nitrogen limitation (unitless)
        FSb: Benthic density attenuation (unitless),
    """

    # Benthic Local Specific Growth Rate
    return xr.where(b_growth_rate_option == 1, mub_max_tc * FLb * FPb * FNb * FSb, mub_max_tc * FLb * FSb * min(FPb, FNb))



@numba.njit
def AbGrowth(
    mub: xr.DataArray,
    Ab: xr.DataArray
) -> xr.DataArray:
    """Calculate Benthic algal growth (g/m^2/d)

    Args:
        mub: Benthic algae growth rate (1/d)
        Ab: Benthic algae concentration (g/m^2)
    """

    return mub * Ab


@numba.njit
def AbRespiration(
    krb_tc: xr.DataArray,
    Ab: xr.DataArray
) -> xr.DataArray:
    """Calculate benthic algal Respiration (g/m^2/d)

    Args:
        krb_tc: Benthic algal respiration rate temperature corrected (1/d)
        Ab: Algae concentration (g/m^2)
    """
    return krb_tc * Ab


@numba.njit
def AbDeath(
    kdb_tc: xr.DataArray,
    Ab: xr.DataArray
) -> xr.DataArray:
    """Calculate benthic algae death (g/m^2/d)

    Args:
        kdb_tc: Benthic algae death rate temperature corrected (1/d)
        Ab: Benthic algae concentration (g/m^2)
    """

    return kdb_tc * Ab


@numba.njit
def dAbdt(
    AbGrowth: xr.DataArray,
    AbRespiration: xr.DataArray,
    AbDeath: xr.DataArray

) -> xr.DataArray:
    """Calculate change in benthic algae concentration (g/m^2/d)

    Args:
        AbGrowth: Benthic algae growth rate (g/m^2/d)
        AbRespiration: Benthic algae respiration rate (g/m^2/d)
        AbDeath: Benthic algae death rate (g/m^2/d)
    """
    return AbGrowth - AbRespiration - AbDeath


@numba.njit
def Chlb(
    rab: xr.DataArray,
    Ab: xr.DataArray,

) -> xr.DataArray:
    """Calculate chlorophyll-a concentration (mg-Chla/m^2)

    Args:
        rab: Balgae Chla to Dry ratio (mg-D/ug-Chla)
        Ab: Benthic algae concentration (g/m^2)
    """

    return rab * Ab
