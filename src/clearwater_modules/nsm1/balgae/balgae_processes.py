"""
File contains process to calculate new benthic algae biomass concentration and associated dependent variables
"""

import math
from clearwater_modules.shared.processes import arrhenius_correction
import numba


@numba.njit
def mub_max_tc(
    mub_max_20: float,
    TwaterC: float
) -> float:
    """Calculate mub_max_tc: Maximum benthic algal growth rate with temperature correction (1/d).

    Args:
        mu_max_20: Maximum benthic algal growth rate at 20C (1/d)
        TwaterC: Water temperature (C)
    """
    return arrhenius_correction(TwaterC, mub_max_20, 1.047)


@numba.njit
def krb_tc(
    krb_20: float,
    TwaterC: float
) -> float:
    """Calculate krb_tc: Benthic algae respiration rate with temperature correction (1/d).

    Args:
        krb_20: Benthic algae respiration rate at 20C (1/d)
        TwaterC: Water temperature (C)
    """
    return arrhenius_correction(TwaterC, krb_20, 1.06)


@numba.njit
def kdb_tc(
    kdb_20: float,
    TwaterC: float
) -> float:
    """Calculate kdb_tc: Benthic algae mortality rate with temperature correction (1/d).

    Args:
        kdb_20: Benthic algae mortality rate at 20C (1/d)
        TwaterC: Water temperature (C)
    """
    return arrhenius_correction(TwaterC, kdb_20, 1.047)


@numba.njit
def rnb(
    BWn: float,
    BWd: float
) -> float:
    """Calculate rnb (mg-N/mg-D).

    Args:
        BWn: Benthic algae nitrogen (unitless)
        BWd: Benthic algae dry weight (unitless)
    """
    return BWn/BWd


@numba.njit
def rpb(
    BWp: float,
    BWd: float
) -> float:
    """Calculate rpd: Benthic algae phosphorus to dry weight ratio (mg-P/mg-D).

    Args:
        BWp: Benthic algae phosphorus (mg-P)
        BWd: Benthic algae dry weight (mg-D)
    """
    return BWp/BWd


@numba.njit
def rcb(
    BWc: float,
    BWd: float
) -> float:
    """Calculate rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D).

    Args:
        BWc: Benthic algae carbon (mg-C)
        BWd: Benthic algae dry weight (mg-D)
    """
    return BWc/BWd


@numba.njit
def rab(
    BWa: float,
    BWd: float
) -> float:
    """Calculate rab: Benthic algae chlorophyll-a to dry weight ratio (ug-Chla-a/mg-D).

    Args:
        BWa: Benthic algae chlorophyll-a (ug-Chla-a)
        BWd: Benthic algae dry weight (mg-D)
    """
    return BWa/BWd


@numba.njit
def FLb(
    L: float,
    depth: float,
    Ab: float,
    PAR: float,
    b_light_limitation_option: float,
    KLb: float,

) -> float:
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

    if Ab <= 0.0 or KEXT <= 0.0 or PAR <= 0.0:
        # After sunset, no growth
        FLb = 0.0
    elif b_light_limitation_option == 1:
        # Use half-saturation formulation
        FLb = PAR * KEXT / (KLb + PAR * KEXT)
    elif b_light_limitation_option == 2:
        # Use Smith's equation
        FLb = PAR * KEXT / ((KLb**2.0 + (PAR * KEXT)**2.0)**0.5)
    elif b_light_limitation_option == 3:
        # Use Steele's equation
        if abs(KLb) < 1.0E-10:
            FLb = 0.0
    else:
        FLb = PAR * KEXT / KLb * math.exp(1.0 - PAR * KEXT / KLb)

    # Limit the benthic light limitation factor to between 0.0 and 1.0
    if FLb > 1.0:
        FLb = 1.0
    if FLb < 0.0:
        FLb = 0.0

    return FLb


@numba.njit
def FNb(
    use_NH4: bool,
    use_NO3: bool,
    NH4: float,
    NO3: float,
    KsNb: float,

) -> float:
    """Calculate Benthic algae nitrogen limitation: FNb (unitless).

    Args:
        use_NH4: Use NH4 module true or false (true/false)
        use_NO3: Use NO3 module true or false (true/false)
        NH4: Ammonium concentration (mg-N/L)
        NO3: Nitrate concentration (mg-N/L)
        KsNb: Michaelis-Menton half-saturation constant relating inorganic N to benthic algal growth (mg-N/L)
    """

    if use_NH4 or use_NO3:
        FNb = (NH4 + NO3) / (KsNb + NH4 + NO3)
        if math.isnan(FNb):
            FNb = 0.0
        if FNb > 1.0:
            FNb = 1.0
    else:
        FNb = 1.0

    return FNb


@numba.njit
def FPb(
    fdp: float,
    TIP: float,
    use_TIP: bool,
    KsPb: float
) -> float:
    """Calculate benthic algae phosphorous limitation: FPb (unitless).

    Args:
        use_TIP: Use Total Inorganic Phosphorus module true or false (true/false)
        TIP: Total Inorganic Phosphorus concentration (mg-P/L)
        KsPb: Michaelis-Menton half-saturation constant relating inorganic P to benthic algal growth (mg-P/L)
        fdp: Fraction P dissolved (unitless)
    """

    if use_TIP:
        FPb = fdp * TIP / (KsPb + fdp * TIP)
        if math.isnan(FPb):
            FPb = 0.0
        if FPb > 1.0:
            FPb = 1.0
    else:
        FPb = 1.0

    return FPb


@numba.njit
def FSb(
    Ab: float,
    Ksb: float,

) -> float:
    """Calculate benthic density attenuation (unitless)

    Args:
        Ab: Benthic algae concentration (g/m^2)
        Ksb: Half-saturation density constant for benthic algae growth (g-D/m^2)

    """

    FSb = 1.0 - (Ab / (Ab + Ksb))
    if math.isnan(FSb):
        FSb = 1.0
    if FSb > 1.0:
        FSb = 1.0

    return FSb


@numba.njit
def mub(
    Ab: float,
    mub_max_tc: float,
    b_growth_rate_option: int,
    FLb: float,
    FPb: float,
    FNb: float,
    FSb: float,

) -> float:
    """Calculate benthic algae specific growth rate (1/d)

    Args:
        Ab: Benthic algae concentration (g/m^2)
        mub_max_tc: Maximum benthic algal growth rate with temperature correction (1/d)
        b_growth_rate_option: Benthic Algal growth rate with three options 1) Multiplicative, 2) Limiting Nutrient
        FLb: Benethic algal light limitation (unitless)
        FPb: Benthic algae phosphorous limitation (unitless)
        FNb: Benthic algae nitrogen limitation (unitless)
        FSb: Benthic density attenuation (unitless),
    """

    # Benthic Local Specific Growth Rate
    if b_growth_rate_option == 1:
        # (a) Multiplicative (day-1)
        mub = mub_max_tc * FLb * FPb * FNb * FSb
    elif b_growth_rate_option == 2:
        # (b) Limiting nutrient (day-1)
        mub = mub_max_tc * FLb * FSb * min(FPb, FNb)

    return mub


@numba.njit
def AbGrowth(
    mub: float,
    Ab: float
) -> float:
    """Calculate Benthic algal growth (g/m^2/d)

    Args:
        mub: Benthic algae growth rate (1/d)
        Ab: Benthic algae concentration (g/m^2)
    """

    return mub * Ab


@numba.njit
def AbRespiration(
    krb_tc: float,
    Ab: float
) -> float:
    """Calculate benthic algal Respiration (g/m^2/d)

    Args:
        krb_tc: Benthic algal respiration rate temperature corrected (1/d)
        Ab: Algae concentration (g/m^2)
    """
    return krb_tc * Ab


@numba.njit
def AbDeath(
    kdb_tc: float,
    Ab: float
) -> float:
    """Calculate benthic algae death (g/m^2/d)

    Args:
        kdb_tc: Benthic algae death rate temperature corrected (1/d)
        Ab: Benthic algae concentration (g/m^2)
    """

    return kdb_tc * Ab


@numba.njit
def dAbdt(
    AbGrowth: float,
    AbRespiration: float,
    AbDeath: float

) -> float:
    """Calculate change in benthic algae concentration (g/m^2/d)

    Args:
        AbGrowth: Benthic algae growth rate (g/m^2/d)
        AbRespiration: Benthic algae respiration rate (g/m^2/d)
        AbDeath: Benthic algae death rate (g/m^2/d)
    """
    return AbGrowth - AbRespiration - AbDeath


@numba.njit
def Chlb(
    rab: float,
    Ab: float,

) -> float:
    """Calculate chlorophyll-a concentration (mg-Chla/m^2)

    Args:
        rab: Balgae Chla to Dry ratio (mg-D/ug-Chla)
        Ab: Benthic algae concentration (g/m^2)
    """

    return rab * Ab
