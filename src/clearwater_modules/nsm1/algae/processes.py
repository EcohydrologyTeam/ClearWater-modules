"""
File contains process to calculate new algae biomass concentration and associated dependent variables
"""
import numba
import xarray as xr
import numpy as np
from clearwater_modules.shared.processes import arrhenius_correction
import math


@numba.njit
def rna(
    AWn: xr.DataArray,
    AWa: xr.DataArray
) -> xr.DataArray:
    """Calculate rna (mg-N/ug-Chla).

    Args:
        AWn: Nitrogen Weight (mg)
        AWa: Algal Chlorophyll (ug-Chla)
    """
    return AWn/AWa


@numba.njit
def rpa(
    AWp: xr.DataArray,
    AWa: xr.DataArray
) -> xr.DataArray:
    """Calculate rpa (mg-P/ug-Chla).

    Args:
        AWp: Phosphorus Weight (mg)
        AWa: Algal Chlorophyll (ug-Chla)
    """

    return AWp/AWa


@numba.njit
def rca(
    AWc: xr.DataArray,
    AWa: xr.DataArray
) -> xr.DataArray:
    """Calculate rca (mg-C/ug-Chla).

    Args:
        AWc: Carbon Weight (mg)
        AWa: Algal Chlorophyll (ug-Chla)
    """
    return AWc/AWa


@numba.njit
def rda(
    AWd: xr.DataArray,
    AWa: xr.DataArray
) -> xr.DataArray:
    """Calculate rda (mg-D/ug-Chla).

    Args:
        AWd: Dry Algal Weight (mg)
        AWa: Algal Chlorophyll (ug-Chla)
    """
    return AWd/AWa


@numba.njit
def mu_max_tc(
    TwaterC: xr.DataArray,
    mu_max_20: xr.DataArray
) -> xr.DataArray:
    """Calculate mu_max_tc (1/d).

    Args:
        TwaterC: Water temperature (C)
        mu_max: Max Algae growth (1/d)
    """

    return arrhenius_correction(TwaterC, mu_max_20, 1.047)


@numba.njit
def krp_tc(
    TwaterC: xr.DataArray,
    krp_20: xr.DataArray
) -> xr.DataArray:
    """Calculate krp_tc (1/d).

    Args:
        TwaterC: Water temperature (C)
        krp: Algal respiration rate at 20 degree (1/d)
    """

    return arrhenius_correction(TwaterC, krp_20, 1.047)


@numba.njit
def kdp_tc(
    TwaterC: xr.DataArray,
    kdp_20: xr.DataArray
) -> xr.DataArray:
    """Calculate kdp_tc (1/d).

    Args:
        TwaterC: Water temperature (C)
        kdp: Algal death rate at 20 degree (1/d)
    """

    return arrhenius_correction(TwaterC, kdp_20, 1.047)


@numba.njit
def FL(
    L: xr.DataArray,
    depth: xr.DataArray,
    Ap: xr.DataArray,
    PAR: xr.DataArray,
    light_limitation_option: int,
    KL: int,
) -> xr.DataArray:
    """Calculate Algal light limitation: FL (unitless).

    Args:
        L: Lambda light attenuation  coefficient (unitless)
        depth: Water depth (m)
        Ap: Algae Concentration (mg-Chla/L)
        PAR: Surface light intensity (W/m^2)
        light_limitation_option: Algal light limitation  option 1) Half-saturation, 2) Smith model, 3) Steele model (unitless)
        KL: Light limitation  constant for algal growth (W/m^2)
    """

    KEXT = L * depth


    FL = xr.where(Ap <= 0.0 or KEXT <= 0.0 or PAR <= 0.0, 0, 
         xr.where(light_limitation_option==1, (1.0 / KEXT) * math.log((KL + PAR) /(KL + PAR * math.exp(-KEXT))),
         xr.where(light_limitation_option==2,
         xr.where(abs(KL)< 0.0000000001, 1, (1.0 / KEXT) * math.log( (PAR / KL + ((1.0 + (PAR / KL)**2.0)**0.5)) / (PAR * math.exp(-KEXT) / KL + ((1.0 + (PAR * math.exp(-KEXT) / KL)**2.0)**0.5)))), 
         xr.where(light_limitation_option==3,
         xr.where(abs(KL)< 0.0000000001,0,(2.718/KEXT) * (math.exp(-PAR/KL * math.exp(-KEXT)) - math.exp(-PAR/KL))), "NaN"))))
          
    
    FL= xr.where(FL > 1.0, 1.0,
        xr.where(FL<0.0, 0.0, FL))

    return FL


@numba.njit
def FN(
    use_NH4: bool,
    use_NO3: bool,
    NH4: xr.DataArray,
    NO3: xr.DataArray,
    KsN: xr.DataArray,

) -> xr.DataArray:
    """Calculate Algal nitrogen limitation: FN (unitless).

    Args:
        use_NH4: Use NH4 module true or false (true/false)
        use_NO3: Use NO3 module true or false (true/false)
        NH4: Ammonium concentration (mg-N/L)
        NO3: Nitrate concentration (mg-N/L)
        KsN: Michaelis-Menton half-saturation constant relating inorganic N to algal growth (mg-N/L)
    """

    FN = xr.where(use_NH4 or use_NO3, (NH4 + NO3) / (KsN + NH4 + NO3), 1)
    FN = xr.where(math.isnan(FN), 0,
         xr.where(FN>1.0,1.0,FN))

    return FN


@numba.njit
def FP(
    fdp: xr.DataArray,
    TIP: xr.DataArray,
    use_TIP: bool,
    KsP: xr.DataArray
) -> xr.DataArray:
    """Calculate Algal phosphorous limitation: FP (unitless).

    Args:
        use_TIP: Use Total Inorganic Phosphorus module true or false (true/false)
        TIP: Total Inorganic Phosphorus concentration (mg-P/L)
        KsP: Michaelis-Menton half-saturation constant relating inorganic P to algal growth (mg-P/L)
        fdp: Fraction P dissolved (unitless)
    """

    FP = xr.where(use_TIP, fdp * TIP / (KsP + fdp * TIP), 1.0)
    FP = xr.where(math.isnan(FP), 0, 
         xr.where(FP>1.0, 1, FP))

    return FP


@numba.njit
def mu(
    mu_max_tc: xr.DataArray,
    growth_rate_option: int,
    FL: xr.DataArray,
    FP: xr.DataArray,
    FN: xr.DataArray

) -> xr.DataArray:
    """Calculate Algal growth rate with three options 1) Multiplicative, 2) Limiting nutrient, 3) Harmonic Mean (1/d)

    Args:
        mu_max_tc: Max algae growth temperature corrected (1/d)
        growth_rate_option: Algal growth rate with options 1) Multiplicative, 2) Limiting nutrient, 3) Harmonic Mean (unitless)
        FL: Algae light limitation factor (unitless)
        FP: Algae phosphorus limitation factor (unitless)
        FN: Algae nitrogen limitation factor (unitless)
    """

    return xr.where(growth_rate_option == 1, mu_max_tc * FL * FP * FN,
           xr.where(growth_rate_option == 2, mu_max_tc * FL * min(FP, FN),
           xr.where(growth_rate_option == 3,
           xr.where(FN==0.0 or FP==0.0, 0.0, mu_max_tc * FL * 2.0 / (1.0 / FN + 1.0 / FP)), "NaN")))


@numba.njit
def ApGrowth(
    mu: xr.DataArray,
    Ap: xr.DataArray
) -> xr.DataArray:
    """Calculate Algal growth (ug-Chla/L/d)

    Args:
        mu: Algal growth rate (1/d)
        Ap: Algae concentration (ug-Chla/L)
    """

    return mu * Ap


@numba.njit
def ApRespiration(
    krp_tc: xr.DataArray,
    Ap: xr.DataArray
) -> xr.DataArray:
    """Calculate Algal Respiration (ug-Chla/L/d)

    Args:
        krp_tc: Algal respiration rate temperature corrected (1/d)
        Ap: Algae concentration (ug-Chla/L)
    """

    return krp_tc * Ap


@numba.njit
def ApDeath(
    kdp_tc: xr.DataArray,
    Ap: xr.DataArray
) -> xr.DataArray:
    """Calculate Algal death (ug-Chla/L/d)

    Args:
        kdp_tc: Algal death rate temperature corrected (1/d)
        Ap: Algae concentration (ug-Chla/L)
    """
    return kdp_tc * Ap


@numba.njit
def ApSettling(
    vsap: xr.DataArray,
    Ap: xr.DataArray,
    depth: xr.DataArray
) -> xr.DataArray:
    """Calculate Algal setting rate (ug-Chla/L/d)

    Args:
        vsap: Algal settling velocity (m/d)
        Ap: Algae concentration (ug-Chla/L)
        depth: Depth from Water Surface (m)s
    """
    return vsap / depth * Ap


@numba.njit
def dApdt(
    ApGrowth: xr.DataArray,
    ApRespiration: xr.DataArray,
    ApDeath: xr.DataArray,
    ApSettling: xr.DataArray
) -> xr.DataArray:
    """Calculate change in algae biomass concentration (ug-Chla/L/d)

    Args:
        ApGrowth: Algal growth (ug-Chla/L/d)
        ApRespiration: Algal respiration (ug-Chla/L/d)
        ApDeath: Algal death (ug-Chla/L/d)
        ApSettling: Algal settling (ug-Chla/L/d)
    """

    return ApGrowth - ApRespiration - ApDeath - ApSettling


@numba.njit
def Ap(
    Ap: xr.DataArray,
    dApdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculate new algae concentration (ug-Chla/L)

    Args:
        Ap: Initial algae biomass concentration (ug-Chla/L)
        dApdt: Change in algae biomass concentration (ug-Chla/L/d)
        timestep: current iteration timestep (d)
    """
    return Ap + dApdt*timestep
