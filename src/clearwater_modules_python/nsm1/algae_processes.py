# Changes: lambda is not calcualted
# Some options such as light calculation method was not changeable in FORTRAN but is now changeable

import math
from clearwater_modules_python.shared.processes import arrhenius_correction
import numba

@numba.njit
def rna(
    AWn : float,
    AWa: float
) -> float :
    return AWn/AWa

@numba.njit
def rpa(
    AWp : float,
    AWa: float
) -> float :
    return AWp/AWa

@numba.njit
def rca(
    AWc : float,
    AWa: float
) -> float :
    return AWc/AWa

@numba.njit
def rda(
    AWd : float,
    AWa: float
) -> float :
    return AWd/AWa

@numba.njit
def mu_max_tc(
    TwaterC : float,
    mu_max: float
) -> float :
    return arrhenius_correction(TwaterC, mu_max, 1.047)

@numba.njit
def krp_tc(
    TwaterC : float,
    krp: float
) -> float :
    return arrhenius_correction(TwaterC, krp, 1.047)

@numba.njit
def kdp_tc(
    TwaterC : float,
    kdp: float
) -> float :
    return arrhenius_correction(TwaterC, kdp, 1.047)

@numba.njit
def FL(
    L : float,
    depth: float,
    Ap: float,
    PAR: float,
    light_limitation_option: int,
    KL : float,
) -> float :

    # Depth averaged light function
    # lambda is light attenuation coefficient (1/m). depth is depth from water surface (m) [unitless] TODO: other depth not initalized method
    KEXT = L * depth
    sqrt1 = 0.0
    sqrt2 = 0.0
    # (1) Algal light limitation (FL)

    if (Ap <= 0.0 or KEXT <= 0.0 or PAR <= 0.0):
        # After sunset or if there is no algae present
        # light limiting factor for algal growth [unitless]
        FL = 0.0
    elif light_limitation_option == 1:
        # Half-saturation formulation
        FL = (1.0 / KEXT) * math.log((KL + PAR) /
                                        (KL + PAR * math.exp(-KEXT)))
    elif light_limitation_option == 2:
        # Smith's model
        if abs(KL) < 1.0E-10:
            FL = 1.0
        else:
            sqrt1 = (1.0 + (PAR / KL)**2.0)**0.5
            sqrt2 = (1.0 + (PAR * math.exp(-KEXT) / KL)**2.0)**0.5
            FL = (1.0 / KEXT) * math.log((PAR / KL + sqrt1) /
                                            (PAR * math.exp(-KEXT) / KL + sqrt2))
    elif light_limitation_option == 3:
        # Steele's model
        if abs(KL) < 1.0E-10:
            FL = 0.0
        else:
            FL = (2.718/KEXT) * (math.exp(-PAR/KL *
                                            math.exp(-KEXT)) - math.exp(-PAR/KL))

    # Limit factor to between 0.0 and 1.0.
    # This should never happen, but it would be a mess if it did.
    if FL > 1.0:
        FL = 1.0
    if FL < 0.0:
        FL = 0.0

    return FL

@numba.njit
def FN(
    use_NH4 : bool,
    use_NO3 : bool,
    NH4 : float,
    NO3 : float,
    KsN : float,

) -> float :

# (2) Algal nitrogen limitation (FN)
# KsN = Michaelis-Menton half-saturation constant (mg N/L) relating inorganic N to algal growth
    if use_NH4 or use_NO3:
        FN = (NH4 + NO3) / (KsN + NH4 + NO3)
        if math.isnan(FN):
            FN = 0.0
        if FN > 1.0:
            FN = 1.0
    else:
        FN = 1.0
    
    return FN

@numba.njit
def FP(
    fdp : float,
    TIP: float,
    use_TIP : bool,
    KsP : float
) -> float :

# (3) Algal phosphorous limitation (FP)
# KsP = Michaelis-Menton half-saturation constant (mg-P/L) relating inorganic P to algal growth
    if use_TIP:

        FP = fdp * TIP / (KsP + fdp * TIP)
        if math.isnan(FP):
            FP = 0.0
        if FP > 1.0:
            FP = 1.0
    else:
        FP = 1.0
    
    return FP

@numba.njit
def mu(
    mu_max_tc : float,
    growth_rate_option: int,
    FL: float,
    FP: float,
    FN: float

) -> float :

# Algal growth rate with three options
# (a) Multiplicative (b) Limiting nutrient (c) Harmonic Mean
    if growth_rate_option == 1:
        # (a) Multiplicative (day-1)
        mu = mu_max_tc * FL * FP * FN                   # [1/d]
    elif growth_rate_option == 2:
        # (b) Limiting nutrient (day-1)
        mu = mu_max_tc * FL * min(FP, FN)
    elif growth_rate_option == 3:
        # (c) Harmonic Mean Option (day-1)
        if FN == 0.0 or FP == 0.0:
            mu = 0.0
        else:
            mu = mu_max_tc * FL * 2.0 / (1.0 / FN + 1.0 / FP)
    
    return mu

@numba.njit
def ApGrowth(
    mu : float,
    Ap: float
) -> float :

    return  mu * Ap                   

@numba.njit
def ApRespiration(
    krp_tc : float,
    Ap: float
) -> float :

    return krp_tc * Ap         

@numba.njit
def ApDeath(
    kdp_tc : float,
    Ap: float
) -> float :

    return kdp_tc * Ap                

@numba.njit
def ApSettling(
    vsap : float,
    Ap: float,
    depth: float
) -> float :

    return vsap / depth * Ap

@numba.njit
def dApdt(
    ApGrowth: float,
    ApRespiration: float,
    ApDeath : float,
    ApSettling : float
) -> float :

    return ApGrowth - ApRespiration - ApDeath - ApSettling
