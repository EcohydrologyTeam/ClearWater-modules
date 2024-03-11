import warnings
import numba
import numpy as np
import xarray as xr
import math

############################################ From shared processes
@numba.njit
def celsius_to_kelvin(tempc: xr.DataArray) -> xr.DataArray:
    return tempc + 273.16


@numba.njit
def kelvin_to_celsius(tempk: xr.DataArray) -> xr.DataArray:
    return tempk - 273.16

@numba.njit
def arrhenius_correction(
    TwaterC: xr.DataArray,
    rc20: xr.DataArray,
    theta: xr.DataArray,
) -> xr.DataArray:
    """
    Computes an adjusted kinetics reaction rate coefficient for the specified water 
    temperature using the van't Hoff form of the Arrhenius equation

    Parameters
    ----------
    TwaterC : xr.DataArray
        Water temperature in degrees Celsius
    rc20 : xr.DataArray
        Kinetics reaction (decay) coefficient at 20 degrees Celsius
    theta : xr.DataArray
        Temperature correction factor

    Returns
    ----------
    float
        Adjusted kinetics rate for the specified water temperature
    """
    return rc20 * theta**(TwaterC - 20.0)


@numba.njit
def TwaterK(
    TwaterC : xr.DataArray,
) -> xr.DataArray :
    """Calculate temperature in kelvin (K)
    Args:
        TwaterC: water temperature celcius (C)
    """
    return celsius_to_kelvin(TwaterC)


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
    condlist = [hydraulic_reaeration_option == 1, 
                hydraulic_reaeration_option == 2, 
                hydraulic_reaeration_option == 3, 
                hydraulic_reaeration_option == 4, 
                hydraulic_reaeration_option == 5 & depth < 0.61,
                hydraulic_reaeration_option == 5 & depth > 0.61,
                hydraulic_reaeration_option == 5 & depth == 0.61, 
                hydraulic_reaeration_option == 6 & flow < 0.556,
                hydraulic_reaeration_option == 6 & flow >= 0.556,
                hydraulic_reaeration_option == 7 & flow < 0.556, 
                hydraulic_reaeration_option == 7 & flow >= 0.556, 
                hydraulic_reaeration_option == 8 & flow < 0.425, 
                hydraulic_reaeration_option == 8 & flow >= 0.425,
                hydraulic_reaeration_option == 9]
    
    choicelist = [kah_20_user, 
                  (3.93 * velocity**0.5) / (depth**1.5), 
                  (5.32 * velocity**0.67) / (depth**1.85), 
                  (5.026 * velocity) / (depth**1.67), 
                  (5.32 * velocity**0.67) / (depth**1.85),
                  (3.93 * velocity**0.5) / (depth**1.5),
                  (5.026 * velocity) / (depth**1.67),
                  517 * (velocity * slope)**0.524 * flow**-0.242,
                  596 * (velocity * slope)**0.528 * flow**-0.136,
                  88 * (velocity * slope)**0.313 * depth**-0.353,
                  142 * (velocity * slope)**0.333 * depth**-0.66 * topwidth**-0.243, 
                  31183 * velocity * slope,
                  15308 * velocity * slope,
                  2.16 * (1 + 9 * (velocity / (9.81 * depth)**0.5)**0.25) * shear_velocity / depth
                  ]

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
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted hydraulic oxygen reaeration rate (/d)

    Args:
        TwaterC: Water temperature in Celsius
        kah_20: Hydraulic oxygen reaeration rate at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(TwaterC, kah_20, theta)


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
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted wind oxygen reaeration velocity (m/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kaw_20: Wind oxygen reaeration velocity at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(TwaterC, kaw_20, theta)


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

def SOD_tc(
    SOD_20: xr.DataArray,
    TwaterC: xr.DataArray,
    SOD_theta: xr.DataArray,
    DOX: xr.DataArray,
    KsSOD: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Compute the sediment oxygen demand corrected by temperature and dissolved oxygen concentration

    Args:
        SOD_20: Sediment oxygen demand at 20 degrees celsius (mg-O2/m2)
        TwaterC: Water temperature in degrees C
        theta: Arrhenius coefficient
        use_DOX: Option to consider DOX concentration in water in calculation of sediment oxygen demand
    """
    SOD_tc = arrhenius_correction(TwaterC, SOD_20, SOD_theta)

    da: xr.DataArray = xr.where(use_DOX == True, SOD_tc * DOX / (DOX + KsSOD), SOD_tc)

    return da

def L(
    lambda0: xr.DataArray,
    lambda1: xr.DataArray,
    lambda2: xr.DataArray,
    lambdas: xr.DataArray,
    lambdam: xr.DataArray,
    Solid: xr.DataArray,
    POC: xr.DataArray,
    fcom: xr.DataArray,
    use_Algae: xr.DataArray,
    use_POC: xr.DataArray,
    Ap: xr.DataArray,

) -> xr.DataArray:
    """Compute L: lambda: light extinction coefficient (unitless)

    Args:
        lambda0: background portion (1/m)
        lambda1: linear self shading (1/m/(ug Chla/L))
        lambda2: non-linear (unitless),
        lambdas: ISS portion (L/mg/m),
        lambdam: POM portion (L/mg/m)
        Solid: #TODO define this
        POC: particulate organic carbon (mg-C/L)
        fcom: ratio of carbon to organic matter (mg-C/mg-D)
        use_Algae: true/false use algae module (t/f)
        use_POC: true/falseo use particulate organic carbon module (t/f)
        Ap: algae concentration (ug-Chla/L)
    """
    L = lambda0 + lambdas * Solid

    L: xr.DataArray = xr.where (use_POC, L+lambdam*POC/fcom, L)
    L: xr.DataArray = xr.where (use_Algae, L+lambda1*Ap + lambda2*Ap**0.66667, L)

    return L

def PAR(
    use_Algae : bool,
    use_Balgae: bool,
    q_solar: xr.DataArray,
    Fr_PAR: xr.DataArray,
) -> xr.DataArray :
    """Calculate temperature in kelvin (K)
    Args:
        use_Algae : true/false use algae module (t/f)
        use_Balgae: true/falsoe use balgae module (t/f)
        q_solar: solar radiation (1/d),
        Fr_PAR: fraction of solar radiation within the PAR of the spectrum
    """
    return xr.where (use_Algae or use_Balgae, q_solar * Fr_PAR, 0)


def fdp(
    use_TIP: bool,
    Solid : xr.DataArray,
    kdpo4: xr.DataArray
) -> xr.DataArray :

    """Calculate kop_tc: Decay rate of organic P to DIP temperature correction (1/d).

    Args:
        use_TIP: true/false use total inorganic phosphrous,
        Solid : #TODO define this
        kdpo4: solid partitioning coeff. of PO4 (L/kg)
    """
  
    return xr.where(use_TIP, 1/(1+kdpo4 * Solid/0.000001), 0)

############################################ From algae
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


def FL(
    L: xr.DataArray,
    depth: xr.DataArray,
    Ap: xr.DataArray,
    PAR: xr.DataArray,
    light_limitation_option: int,
    KL: xr.DataArray,
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
         xr.where(abs(KL)<0.0000000001, 1, (1.0 / KEXT) * math.log( (PAR / KL + ((1.0 + (PAR / KL)**2.0)**0.5)) / (PAR * math.exp(-KEXT) / KL + ((1.0 + (PAR * math.exp(-KEXT) / KL)**2.0)**0.5)))), 
         xr.where(light_limitation_option==3,
         xr.where(abs(KL)<0.0000000001,0,(2.718/KEXT) * (math.exp(-PAR/KL * math.exp(-KEXT)) - math.exp(-PAR/KL))), "NaN"))))
          
    
    FL= xr.where(FL > 1.0, 1.0,
        xr.where(FL<0.0, 0.0, FL))

    return FL


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
    return Ap + dApdt * timestep

############################################ From benthic algae
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
def Ab(
    Ab: xr.DataArray,
    dAbdt: xr.DataArray,
    timestep: xr.DataArray,

) -> xr.DataArray:
    """Calculate Ab: New concentration benthic algae (mg-N/L)

    Args:
        Ab: Concentration of benthic algae (mg-N/L)
        dAbdt: Change in Ab (mg-N/L/d)
        timestep: current iteration timestep (d)

    """

    return Ab + dAbdt * timestep

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


############################################ From nitrogen

@numba.njit
def knit_tc(
    TwaterC: xr.DataArray,
    knit_20: xr.DataArray
) -> xr.DataArray:
    """Calculate knit_tc: Nitrification rate ammonia decay NH4 to NO3 temperature correction (1/d). #TODO only if use_NH4 = true

    Args:
        TwaterC: Water temperature (C)
        knit_20: Nitrification rate ammonia decay (1/d)
    """

    return arrhenius_correction(TwaterC, knit_20, 1.083)


@numba.njit
def rnh4_tc(
    TwaterC: xr.DataArray,
    rnh4_20: xr.DataArray
) -> xr.DataArray:
    """Calculate rnh4_tc: Sediment release rate of NH4 temperature correction(1/d). #TODO only if use_sedflux = true

    Args:
        TwaterC: Water temperature (C)
        rnh4_20: Sediment release rate of NH4 (1/d)
    """

    return arrhenius_correction(TwaterC, rnh4_20, 1.074)


@numba.njit
def vno3_tc(
    TwaterC: xr.DataArray,
    vno3_20: xr.DataArray
) -> xr.DataArray:
    """Calculate vno3_tc: Sediment denitrification velocity temperature correction (m/d). #TODO only if use_sedflux = true

    Args:
        TwaterC: Water temperature (C)
        vno3_20: Sedimet release rate of NO3 (1/d)
    """

    return arrhenius_correction(TwaterC, vno3_20, 1.08)


@numba.njit
def kon_tc(
    TwaterC: xr.DataArray,
    kon_20: xr.DataArray
) -> xr.DataArray:
    """Calculate kon_tc: Decay rate of OrgN to NH4 temperature correction(1/d). #TODO only if use_OrgN = true

    Args:
        TwaterC: Water temperature (C)
        kon_20: Decay rate of OrgN to NH4 (1/d)
    """

    return arrhenius_correction(TwaterC, kon_20, 1.074)


@numba.njit
def kdnit_tc(
    TwaterC: xr.DataArray,
    kdnit_20: xr.DataArray
) -> xr.DataArray:
    """Calculate kdnit_tc: Denitrification rate temperature correction (1/d). #TODO only if use_NO3 = true

    Args:
        TwaterC: Water temperature (C)
        kdnit_20: Denitrification rate (1/d)
    """

    return arrhenius_correction(TwaterC, kdnit_20, 1.045)


def ApUptakeFr_NH4(
    use_NH4: bool,
    use_NO3: bool,
    use_Algae: bool,
    PN: xr.DataArray,
    NH4: xr.DataArray,
    NO3: xr.DataArray,

) -> xr.DataArray:
    """Calculate ApUptakeFr_NH4: 

    Args:
        use_NH4: use ammonium module (unitless)
        use_NO3: use nitrate module (unitless)
        use_Algae: use algae module (unitless)
        PN: NH4 preference factor algae (unitless)
        NH4: Ammonium water concentration (mg-N/L)
        NO3: Nitrate water concentration (mg-N/L)
    """
    ApUptakeFr_NH4 = 0

    # set value of UptakeFr_NH4/NO3 for special conditions
    ApUptakeFr_NH4 = xr.where(use_NH4 and not use_NO3, 1.0,
                     xr.where(not use_NH4 and use_NO3, 0.0,
                     xr.where(not use_NH4 and not use_NO3, 0.5,
                     xr.where(use_Algae and use_NH4 and use_NO3, PN * NH4 / (PN * NH4 + (1.0 - PN) * NO3), "NaN"))))
    
    # Check for case when NH4 and NO3 are very small.  If so, force uptake_fractions appropriately.
    ApUptakeFr_NH4 = xr.where(math.isnan(ApUptakeFr_NH4),PN,ApUptakeFr_NH4)

    return ApUptakeFr_NH4


@numba.njit
def ApUptakeFr_NO3(
    ApUptakeFr_NH4: xr.DataArray
) -> xr.DataArray:
    """Calculate ApUptakeFr_NO3: Fraction of actual xr.DataArraying algal uptake from nitrate pool (unitless)

    Args:
        ApUptakeFr_NH4: Fraction of actual xr.DataArraying algal uptake from ammonia pool
    """

    return 1 - ApUptakeFr_NH4


def AbUptakeFr_NH4(
    use_NH4: bool,
    use_NO3: bool,
    use_Balgae: bool,
    PNb: xr.DataArray,
    NH4: xr.DataArray,
    NO3: xr.DataArray,

) -> xr.DataArray:
    """Calculate AbUptakeFr_NH4: Fraction of actual benthic algal uptake from ammonia pool

    Args:
        use_NH4: use ammonium module (unitless)
        use_NO3: use nitrate module (unitless)
        use_Balgae: use benthic algae module (unitless)
        PNb: NH4 preference factor benthic algae (unitless)
        NH4: Ammonium water concentration (mg-N/L)
        NO3: Nitrate water concentration (mg-N/L)
    """
    AbUptakeFr_NH4 = 0
    AbUptakeFr_NH4 = xr.where(use_NH4 and not use_NO3, 1.0,
                     xr.where(not use_NH4 and use_NO3, 0.0,
                     xr.where(not use_NH4 and not use_NO3, 0.5,
                     xr.where(use_Balgae and use_NH4 and use_NO3, (PNb * NH4) / (PNb * NH4 + (1.0 - PNb) * NO3), "NaN"))))
    
    AbUptakeFr_NH4 = xr.where(math.isnan(AbUptakeFr_NH4),PNb,AbUptakeFr_NH4)

    return AbUptakeFr_NH4


@numba.njit
def AbUptakeFr_NO3(
    AbUptakeFr_NH4: xr.DataArray
) -> xr.DataArray:
    """Calculate AbUptakeFr_NO3: Fraction of actual benthic algal uptake from nitrate pool (unitless)

    Args:
        AbUptakeFr_NH4: Fraction of actual benthic algal uptake from ammonia pool
    """

    return 1 - AbUptakeFr_NH4

def OrgN_NH4_Decay(
    kon_tc: xr.DataArray,
    OrgN: xr.DataArray,
    use_OrgN: bool
) -> xr.DataArray:
    """Calculate OrgN_NH4: OrgN -> NH4 (mg-N/L/d)

    Args:
        kon_tc: Decay rate of organic nitrogen to nitrate with temperature correction (1/d),
        OrgN: Concentration of organic nitrogen (mg-N/L)
        use_OrgN: true/false use organic nitrogen (t/f)
    """

    return xr.where(use_OrgN, kon_tc * OrgN,0)

@numba.njit
def OrgN_Settling(
    vson: xr.DataArray,
    depth: xr.DataArray,
    OrgN: xr.DataArray,
) -> xr.DataArray:
    """Calculate OrgN_Settling: OrgN -> bed (mg-N/L/d)

    Args:
        vson: Organic nitrogen settling velocity (m/d)
        depth: water depth (m)
    """

    return vson / depth * OrgN

def ApDeath_OrgN(
    use_Algae: bool,
    rna: xr.DataArray,
    ApDeath: xr.DataArray,
) -> xr.DataArray:
    """Calculate ApDeath_OrgN: Algae -> OrgN (mg-N/L/d)

    Args:
        use_Algae: true/false to use algae module (unitless)
        rna: Algal N: Chla ratio (mg-N/ug-Chla)
        ApDeath: Algal death rate (ug-Chla/L/d)
    """

    return xr.where(use_Algae, rna * ApDeath, 0.0)

def AbDeath_OrgN(
    use_Balgae: bool,
    rnb: xr.DataArray,
    Fw: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    AbDeath: xr.DataArray,
) -> xr.DataArray:
    """Calculate ApDeath_OrgN: Algae -> OrgN (mg-N/L/d)

    Args:
        use_Balgae: true/false to use benthic algae module (unitless),
        rnb: Benthic algal N: Benthic Algal Dry Weight (mg-N/mg-D)
        Fw: Fraction benthic algae mortality into water column (unitless)
        Fb: Fraction of bottom area for benthic algae (unitless)
        depth: water depth (m)
        AbDeath: Benthic algal death rate (g/m^2/d)
    """

    return xr.where(use_Balgae, rnb * Fw * Fb * AbDeath / depth, 0.0)

def dOrgNdt(
    use_OrgN: bool,
    ApDeath_OrgN: xr.DataArray,
    AbDeath_OrgN: xr.DataArray,
    OrgN_NH4_Decay: xr.DataArray,
    OrgN_Settling: xr.DataArray,

) -> xr.DataArray:
    """Calculate dOrgNdt: Change in Organic Nitrogen (mg-N/L/d)

    Args:
        use_OrgN: true/false to use organic nitrogen module (unitless)
        ApDeath_OrgN: Algae -> OrgN (mg-N/L/d)
        AbDeath_OrgN: Benthic Algae -> OrgN (mg-N/L/d)
        OrgN_NH4_Decay: OrgN -> NH4 (mg-N/L/d)
        OrgN_Settling: OrgN -> bed (mg-N/L/d)

    """

    return xr.where(use_OrgN, ApDeath_OrgN + AbDeath_OrgN - OrgN_NH4_Decay - OrgN_Settling,0)

@numba.njit
def OrgN(
    OrgN: xr.DataArray,
    dOrgNdt: xr.DataArray,
    timestep: xr.DataArray,

) -> xr.DataArray:
    """Calculate OrgN: New concentration OrgN (mg-N/L)

    Args:
        OrgN: Concentration of organic nitrogen (mg-N/L)
        dOrgNdt: Change in Organic Nitrogen (mg-N/L/d)
        timestep: current iteration timestep (d)

    """

    return OrgN + dOrgNdt * timestep

def NitrificationInhibition(
    use_DOX: bool,
    KNR: xr.DataArray,
    DOX: xr.DataArray,

) -> xr.DataArray:
    """Calculate NitrificationInhibition: Nitrification Inhibitation (limits nitrification under low DO conditions)

    Args:
        KNR: Oxygen inhabitation factor for nitrification (mg-O2/L),
        DOX: Dissolved oxygen concentration (mg-O2/L),
        use_DOX: true/false to use dissolve oxygen module (unitless),

    """

    return xr.where (use_DOX, 1.0 - math.exp(-KNR * DOX), 1.0)

def NH4_Nitrification(
    NitrificationInhibition: xr.DataArray,
    NH4: xr.DataArray,
    knit_tc: xr.DataArray,
    use_NH4: xr.DataArray

) -> xr.DataArray:
    """Calculate NH4_Nitrification: NH4 -> NO3  Nitrification  (mg-N/L/day)

    Args:
        NitrificationInhibition: Nitrification Inhibitation (limits nitrification under low DO conditions)
        knit_tc: Nitrification rate ammonia decay NH4 to NO3 temperature correction (1/d).
        NH4: Ammonium concentration (mg-N/L),
    """

    return xr.where(use_NH4,NitrificationInhibition * knit_tc * NH4,0)

@numba.njit
def NH4fromBed(
    depth: xr.DataArray,
    rnh4_tc: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH4fromBed: bed ->  NH4 (diffusion)    (mg-N/L/day)

    Args:
        depth: water depth (m),
        rnh4_tc: Sediment release rate of NH4 temperature correction(1/d).

    """

    return rnh4_tc / depth

def NH4_ApRespiration(
    use_Algae: bool,
    ApRespiration: xr.DataArray,
    rna: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH4_ApRespiration: Floating algae -> NH4      (mg-N/L/day)   

    Args:
        use_Algae: true/false to use algae module (unitless),
        rna: Algal N: Chla ratio (mg-N/ug-Chla),
        ApRespiration: Algal respiration rate (ug-Chla/L/d),
    """

    return xr.where (use_Algae, rna * ApRespiration, 0.0)

def NH4_ApGrowth(
    use_Algae: bool,
    ApGrowth: xr.DataArray,
    rna: xr.DataArray,
    ApUptakeFr_NH4: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH4_ApGrowth: NH4 -> Floating algae      (mg-N/L/day)

    Args:
        use_Algae: true/false to use algae module (unitless),
        rna: Algal N: Chla ratio (mg-N/ug-Chla),
        ApGrowth: Algal growth rate (ug-Chla/L/d),
        ApUptakeFr_NH4: Fraction of actual xr.DataArraying algal uptake from ammonia pool
    """

    return xr.where(use_Algae, ApUptakeFr_NH4 * rna * ApGrowth, 0.0)

def NH4_AbRespiration(
    use_Balgae: bool,
    rnb: xr.DataArray,
    AbRespiration: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH4_AbRespiration: Benthic algae -> NH4       (mg-N/L/day) 

    Args:
        use_Balgae: true/false to use benthic algae module (unitless),
        rnb: xr.DataArray,
        AbRespiration: Benthic algal respiration rate (g/m^2/d),
    """
    # TODO changed the calculation for respiration from the inital FORTRAN due to conflict with the reference guide

    return xr.where(use_Balgae, rnb * AbRespiration, 0.0 )

def NH4_AbGrowth(
    use_Balgae: bool,
    rnb: xr.DataArray,
    AbGrowth: xr.DataArray,
    AbUptakeFr_NH4: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH4_AbGrowth: NH4 -> Benthic Algae       (g-N/L/day)

    Args:
        use_Balgae: true/false to use benthic algae module (unitless),
        rnb: xr.DataArray,
        AbGrowth: Benthic alga growth rate (g/m^2/d),  
        depth: water depth (m),
        Fb: Fraction of bottom area for benthic algae (unitless),
        AbUptakeFr_NH4: Fraction of actual benthic algal uptake from ammonia pool
    """

    return xr.where(use_Balgae,(AbUptakeFr_NH4 * rnb * Fb * AbGrowth) / depth, 0.0 )

def dNH4dt(
    use_NH4: bool,
    NH4_Nitrification: xr.DataArray,
    NH4fromBed: xr.DataArray, 
    NH4_ApRespiration: xr.DataArray, 
    NH4_ApGrowth: xr.DataArray,
    NH4_AbRespiration: xr.DataArray,
    NH4_AbGrowth: xr.DataArray,


) -> xr.DataArray:
    """Calculate dNH4dt: Change in Ammonium (mg-N/L)

    Args:
        use_OrgN: true/false to use organic nitrogen module (unitless),
        use_NH4: true/false to use ammonium module (unitless),
        NH4_Nitrification: NH4 -> NO3  Nitrification  (mg-N/L/day)
        NH4fromBed: bed ->  NH4 (diffusion)    (mg-N/L/day)
        NH4_ApRespiration: Floating algae -> NH4      (mg-N/L/day)   
        NH4_ApGrowth: NH4 -> Floating algae      (mg-N/L/day)
        NH4_AbRespiration: Benthic algae -> NH4       (mg-N/L/day) 
        NH4_AbGrowth: NH4 -> Benthic Algae       (g-N/L/day)

    Ammonia Nitrogen (NH4)                 (mg-N/day*L)
    dNH4/dt   =    OrgN_NH4_Decay   
                  (OrgN -> NH4)  
                    - NH4 Oxidation         (NH4 -> NO3)
                    - NH4AlgalUptake        (NH4 -> xr.DataArraying Algae)
                    + Benthos NH4           (Benthos -> NH4)
                    - Benthic Algae Uptake  (NH4 -> Benthic Algae)		

    """

    return xr.where(use_NH4, OrgN_NH4_Decay - NH4_Nitrification + NH4fromBed + NH4_ApRespiration - NH4_ApGrowth + NH4_AbRespiration - NH4_AbGrowth, 0.0)

@numba.njit
def NH4(
    NH4: xr.DataArray,
    dNH4dt: xr.DataArray,
    timestep: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH4: New concentration NH4 (mg-N/L)

    Args:
        NH4: Concentration of NH4 (mg-N/L)
        dNH4dt: Change in NH4 (mg-N/L/d)
        timestep: current iteration timestep (d)

    """

    return NH4 + dNH4dt * timestep

def NO3_Denit(
    use_DOX: bool,
    DOX: xr.DataArray,
    KsOxdn: xr.DataArray,
    kdnit_tc: xr.DataArray,
    NO3: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO3_Denit:  NO3 -> Loss                (mg-N/L/day)

    Args:
        use_DOX: true/false to use dissolve oxygen module (unitless),
        KsOxdn: Half-saturation oxygen inhibition constant for denitrification (mg-O2/L)
        DOX: Dissolved oxygen concentration (mg-O2/L),
        NO3: Nitrate concentration (mg-N/L),
        kdnit_tc: Denitrification rate temperature correction (1/d)

    """
    return xr.where(use_DOX,xr.where(math.isnan((1.0 - (DOX / (DOX + KsOxdn))) * kdnit_tc * NO3),kdnit_tc * NO3,(1.0 - (DOX / (DOX + KsOxdn))) * kdnit_tc * NO3),0.0)

@numba.njit
def NO3_BedDenit(
    depth: xr.DataArray,
    vno3_tc: xr.DataArray,
    NO3: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO3_BedDenit: Sediment denitrification   (mg-N/L/day) 

    Args:
        depth: water depth (m),
        NO3: Nitrate concentration (mg-N/L)
        vno3_tc: Sediment denitrification velocity temperature correction (m/d)

    """
    
    return vno3_tc * NO3 / depth

def NO3_ApGrowth(
    use_Algae: bool,
    ApUptakeFr_NO3: xr.DataArray,
    rna: xr.DataArray,
    ApGrowth: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO3_ApGrowth: NO3 -> Floating algae      (mg-N/L/day)

    Args:
        use_Algae: true/false to use algae module (unitless),
        rna: Algal N: Chla ratio (mg-N/ug-Chla),
        ApGrowth: Algal growth rate (ug-Chla/L/d),
        ApUptakeFr_NO3: Fraction of actual algal uptake from nitrate pool (unitless)


    """

    return xr.where(use_Algae, ApUptakeFr_NO3 * rna * ApGrowth, 0.0)

def NO3_AbGrowth(
    use_Balgae: bool,
    AbUptakeFr_NO3: xr.DataArray,
    rnb: xr.DataArray,
    Fb: xr.DataArray,
    AbGrowth: xr.DataArray,
    depth: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO3_AbGrowth: NO3 -> Benthic Algae       (g-N/L/day)

    Args:
        use_Balgae: true/false to use benthic algae module (unitless),
        depth: water depth (m),
        rnb: Benthic algal N: Benthic Algal Dry Weight (mg-N/mg-D),
        Fb: Fraction of bottom area for benthic algae (unitless),
        AbGrowth: Benthic alga growth rate (g/m^2/d),
        AbUptakeFr_NO3: Fraction of actual benthic algal uptake from nitrate pool (unitless)
    """

    return xr.where(use_Balgae, (AbUptakeFr_NO3 * rnb * Fb * AbGrowth) / depth, 0.0)


def dNO3dt(
    use_NO3: bool,
    NH4_Nitrification: xr.DataArray,
    NO3_Denit: xr.DataArray,
    NO3_BedDenit: xr.DataArray,
    NO3_ApGrowth: xr.DataArray,
    NO3_AbGrowth: xr.DataArray,

) -> xr.DataArray:
    """Calculate dNO3dt: Change in nitrate (mg-N/L)

    Args:
        use_NH4: true/false to use ammonium module (unitless),
        use_NO3: true/false to use nitrate module (unitless),
        NH4_Nitrification: NH4 -> NO3  Nitrification  (mg-N/L/day)
        NO3_Denit: NO3 -> Loss                (mg-N/L/day),
        NO3_BedDenit: Sediment denitrification   (mg-N/L/day) 
        NO3_ApGrowth: NO3 -> Floating algae      (mg-N/L/day)
        NO3_AbGrowth: NO3 -> Benthic Algae       (g-N/L/day)

    Nitrite Nitrogen  (NO3)                       (mg-N/day*L)
    dNO3/dt  =      NH4 Oxidation                 (NH4 -> NO3) 
                    - NO3 Sediment Denitrification
                    - NO3 Algal Uptake            (NO3-> xr.DataArraying Algae) 
                    - NO3 Benthic Algal Uptake    (NO3-> Benthic  Algae) 	

    """


    return xr.where(use_NO3, NH4_Nitrification - NO3_Denit - NO3_BedDenit - NO3_ApGrowth - NO3_AbGrowth ,0)

@numba.njit
def NO3(
    NO3: xr.DataArray,
    dNO3dt: xr.DataArray,
    timestep: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO3: New concentration NO# (mg-N/L)

    Args:
        NO3: Concentration of NO3 (mg-N/L)
        dNO3dt: Change in NO3(mg-N/L/d)
        timestep: current iteration timestep (d)

    """

    return NO3 + dNO3dt * timestep


def DIN(
    use_NH4: bool,
    use_NO3: bool,
    NH4: xr.DataArray,
    NO3: xr.DataArray,

) -> xr.DataArray:
    """Calculate DIN: Dissolve inorganic nitrogen (mg-N/L)

    Args:
        use_NH4: true/false to use ammonium module (unitless),
        use_NO3: true/false to use nitrate module (unitless),
        NH4: Ammonium concentration (mg-N/L), 
        NO3: Nitrate concentration (mg-N/L), 
    """
    DIN = 0.0
    DIN = xr.where(use_NH4, DIN + NH4,DIN)
    DIN = xr.where(use_NO3, DIN + NO3, DIN)

    return DIN


def TON(
    use_OrgN: bool,
    use_Algae: bool,
    OrgN: xr.DataArray,
    rna: xr.DataArray,
    Ap: xr.DataArray

) -> xr.DataArray:
    """Calculate TON: Total organic nitrogen (mg-N/L)

    Args:
        use_OrgN: true/false to use organic nitrogen module (unitless),
        use_Algae: true/false to use algae module (unitless),
        OrgN: Organic nitrogen concentration (mg-N/L), 
        rna: Algal N: Chla ratio (mg-N/ug-Chla),
        Ap: Algae water concentration (ug-Chla/L)

    """
    TON = 0.0
    TON = xr.where(use_OrgN, TON + OrgN, TON)
    TON = xr.where(use_Algae, TON + rna * Ap, TON)
    
    return TON


def TKN(
    use_NH4: bool,
    NH4: xr.DataArray,
    TON: xr.DataArray

) -> xr.DataArray:
    """Calculate TKN: Total kjeldhl (mg-N/L)

    Args:
        use_NH4: true/false to use organic nitrogen module (unitless),
        NH4: Ammonium concentration (mg-N/L)
        TON: Total organic nitrogen (mg-N/L)
    """
    TKN = 0.0
    TKN = xr.where(use_NH4, TKN + NH4, TKN)

    return TKN + TON


@numba.njit
def TN(
    DIN: xr.DataArray,
    TON: xr.DataArray,

) -> xr.DataArray:
    """Calculate TN: Total nitrogen (mg-N/L)

    Args:
        DIN: Dissolve inorganic nitrogen (mg-N/L)
        TON: Total organic nitrogen (mg-N/L)
    """

    return DIN + TON

################################### From phosphorus
@numba.njit
def kop_tc(
    TwaterC : xr.DataArray,
    kop_20: xr.DataArray
) -> xr.DataArray :

    """Calculate kop_tc: Decay rate of organic P to DIP temperature correction (1/d).

    Args:
        TwaterC: Water temperature (C)
        kop_20: Decay rate of organic P to DIP at 20C (1/d)
    """

    return arrhenius_correction(TwaterC, kop_20, 1.047)

@numba.njit
def rpo4_tc(
    TwaterC : xr.DataArray,
    rpo4_20: xr.DataArray
) -> xr.DataArray :

    """Calculate rpo4_tc: Benthic sediment release rate of DIP temperature correction(g-P/m2/d).

    Args:
        TwaterC: Water temperature (C)
        kop_20: Benthic sediment release rate of DIP at 20C (1/d)
    """

    return arrhenius_correction(TwaterC, rpo4_20, 1.074)

def OrgP_DIP_decay(
    kop_tc : xr.DataArray,
    OrgP: xr.DataArray,
    use_OrgP: bool,
) -> xr.DataArray :
    
    """Calculate OrgP_DIP: organic phosphorus decay to dissolve inorganic phosphorus (mg-P/L/d).

    Args:
        kop_tc: Decay rate of organic P to DIP temperature correction (1/d)
        OrgP: Organic phosphorus concentration (mg-P/L)
        use_OrgP: true/false use organic phosphorus (t/f)
    """        
    return xr.where(use_OrgP,kop_tc * OrgP,0)

@numba.njit
def OrgP_Settling(
    vsop : xr.DataArray,
    depth: xr.DataArray,
    OrgP: xr.DataArray,
) -> xr.DataArray :
    
    """Calculate OrgP_Settling: organic phosphorus settling to sediment (mg-P/L/d).

    Args:
        vsop: Organic phosphorus settling velocity (m/d)
        depth: water depth (m)
        OrgP: Organic phosphorus concentration (mg-P/L)
    """        
    return (vsop / depth) * OrgP

def ApDeath_OrgP(
    rpa : xr.DataArray,
    ApDeath: xr.DataArray,
    use_Algae: bool,
) -> xr.DataArray :
    
    """Calculate ApDeath_OrgP: Algal death turning into organic phosphorus  (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApDeath: Algal death rate (ug-Chla/L/d)
        use_Algae: true/false to use algae module (T/F)

    """        

    return xr.where(use_Algae, rpa * ApDeath,0)

def AbDeath_OrgP(
    rpb : xr.DataArray,
    AbDeath: xr.DataArray,
    Fw: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: bool
) -> xr.DataArray :
    
    """Calculate AbDeath_OrgP: Benthic algal death turning into organic phosphorus (mg-P/L/d).

    Args:
        rpb : Benthic algal P: Benthic algal dry (mg-P/mg-D)
        AbDeath: Benthic algal death rate (g/m^2/d)
        Fw: Fraction benthic algal death to water column (unitless)
        Fb: Fraction bottom area avalible for benthic algae (unitless)
        depth: water depth (m)
        use_Balgae: true/false use benthic algae module (t/f)

    """        

    return xr.where(use_Balgae, (rpb * Fw *Fb * AbDeath) / depth,0)      

def dOrgPdt(
    ApDeath_OrgP : xr.DataArray,
    AbDeath_OrgP: xr.DataArray,
    OrgP_DIP_decay: xr.DataArray,
    OrgP_Settling: xr.DataArray,
    use_OrgP: bool,
) -> xr.DataArray :
    """Calculate dOrgPdt: change in organic phosphorus concentration (mg-P/L/d).

    Args:
        ApDeath_OrgP: Algal death turns into organic phosphrous 
        AbDeath_OrgP: Benthic algal death turns into organic phosphrous
        OrgP_DIP_decay: Organic phosphrous decaying into dissolve inorganic phosphrous
        OrgP_Settling: Organic phosphrous settling into sediment
        use_OrgP: true/false to use organic phosphorus module (true/false)
        use_Algae: true/false to use algae module (true/false)
        use_Balgae: true/false to use benthic algae module (true/false)
    """     

    return xr.where(use_OrgP, -OrgP_DIP_decay-OrgP_Settling + ApDeath_OrgP + AbDeath_OrgP, 0)

@numba.njit
def DIPfromBed(
    depth:xr.DataArray,
    rpo4_tc: xr.DataArray,
) -> xr.DataArray :
    """Calculate DIPfromBed: Dissolved Organic Phosphorus coming from Bed calculated without a SedFlux module (mg-P/L/d).

    Args:
        depth: water depth (m)
        rpo4_tc: Benthic sediment release rate of DIP temperature correction(g-P/m2/d)
    """    
    return rpo4_tc / depth

#TODO calcuate fdp?
@numba.njit
def TIP_Settling(
    vs: xr.DataArray,
    depth: xr.DataArray,
    fdp: xr.DataArray,
    TIP: xr.DataArray
) -> xr.DataArray :

    """Calculate TIP_Settling: Total inorganic phosphorus settling from water to bed (mg-P/L/d).

    Args:
        vs: Sediment settling velocity (m/d)
        depth: water depth (m)
        fdp: Fraction phosphorus dissolved (unitless)
        TIP: Total inorganic phosphorus water concentration (mg-P/L)
    """             
    return vs / depth * (1.0 - fdp) * TIP

def DIP_ApRespiration(
    rpa: xr.DataArray,
    ApRespiration: xr.DataArray,
    use_Algae: bool

) -> xr.DataArray :
    """Calculate DIP_ApRespiration: Dissolved inorganic phosphorus released from algal respiration (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApRespiration: Algal respiration rate (ug-Chla/L/d)
        use_Algae: true/false to use algae module (t/f)
    """ 
    return xr.where(use_Algae, rpa * ApRespiration,0)

def DIP_ApGrowth(
    rpa: xr.DataArray,
    ApGrowth: xr.DataArray,
    use_Algae: bool

) -> xr.DataArray :
    """Calculate DIP_ApGrowth: Dissolved inorganic phosphorus consumed for algal growth (mg-P/L/d).

    Args:
        rpa: Algal P : Chla ratio (mg-P/ug-Chla)
        ApGrowth: Algal growth rate (ug-Chla/L/d)
        use_Algae: true/false to use algae module (t/f)
    """ 
    return xr.where(use_Algae, rpa * ApGrowth,0)

def DIP_AbRespiration(
    rpb: xr.DataArray,
    AbRespiration: xr.DataArray,
    use_Balgae: bool

) -> xr.DataArray :
    """Calculate DIP_AbRespiration: Dissolved inorganic phosphorus released for benthic algal respiration (mg-P/L/d).

    Args:
        rpb: Benthic algal P : Benthic algal dry ratio (mg-P/mg-D)
        AbRespiration: Benthic algal respiration rate (g/m^2/d)
        use_Blgae: true/false to use benthic algae module (t/f)        
    """     
    return xr.where(use_Balgae, rpb * AbRespiration,0)

def DIP_AbGrowth(
    rpb: xr.DataArray,
    AbGrowth: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: bool

) -> xr.DataArray :
    """Calculate DIP_AbGrowth: Dissolved inorganic phosphorus consumed for benthic algal growth (mg-P/L/d).

    Args:
        rpb: Benthic algal P : Benthic algal dry ratio (mg-P/mg-D)
        AbGrowth: Benthic algal growth rate (g/m^2/d)
        Fb: Fraction of bottom area available for benthic algal (unitless)
        depth: water depth (m)
        use_Balgae: true/false to use benthic algae module (t/f) 
    """     
    return xr.where(use_Balgae, rpb * Fb * AbGrowth / depth,0)

def dTIPdt(
    OrgP_DIP_decay: xr.DataArray,
    TIP_Settling: xr.DataArray,
    DIPfromBed: xr.DataArray,
    DIP_ApRespiration: xr.DataArray,
    DIP_ApGrowth: xr.DataArray,
    DIP_AbRespiration: xr.DataArray,
    DIP_AbGrowth: xr.DataArray,
    use_TIP: bool, 

) -> xr.DataArray :    
    
    """Calculate dTIPdt: Change in dissolved inorganic phosphorus water concentration (mg-P/L/d).

    Args:
        OrgP_DIP_decay: Total organic phosphorus decaying to dissolved inorganic phosphrous (mg-P/L/d),
        TIP_Settling: Total inorganic phosphorus settling from water to bed (mg-P/L/d),
        DIPfromBed_NoSedFlux: Dissolved Organic Phosphorus coming from Bed calculated without SedFlux modules (mg-P/L/d),
        DIPfromBed_SedFlux: Dissolved Organic Phosphorus coming from Bed calculated using SedFlux modules (mg-P/L/d),
        DIP_ApRespiration: Dissolved inorganic phosphorus released from algal respiration (mg-P/L/d),
        DIP_ApGrowth: Dissolved inorganic phosphorus consumed for algal growth (mg-P/L/d),
        DIP_AbRespiration: Dissolved inorganic phosphorus released for benthic algal respiration (mg-P/L/d),
        DIP_AbGrowth: Dissolved inorganic phosphorus consumed for benthic algal growth (mg-P/L/d),
        use_TIP: true/false to use total inorganic phosphorus module (true/false), 


    dTIP/dt =     OrgP Decay                (OrgP -> DIP)
                - DIP AlgalUptake           (DIP -> xr.DataArraying Algae)
                - DIP BenthicAlgae Uptake   (DIP -> xr.DataArraying Algae)	
                - TIP Settling              (TIP -> bed)
                + DIP From Benthos          (Benthos -> DIP) 
    """

    return xr.where(use_TIP, - TIP_Settling + DIPfromBed + OrgP_DIP_decay + DIP_ApRespiration - DIP_ApGrowth + DIP_AbRespiration - DIP_AbGrowth, 0)


@numba.njit
def TIP(
    TIP: xr.DataArray,
    dTIPdt: xr.DataArray,
    timestep: xr.DataArray

) -> xr.DataArray :
    """Calculate TIP: New total inorganic phosphorus (mg-P/L).

    Args:
        dTIPdt: Change in total inorganic phosphorus (mg-P/L/d)
        TIP: Total inorganic phosphorus water concentration (mg-P/L),
        timestep: current iteration timestep (d)
    """     
    return TIP + dTIPdt * timestep

@numba.njit
def OrgP(
    OrgP: xr.DataArray,
    dOrgPdt: xr.DataArray,
    timestep: xr.DataArray

) -> xr.DataArray :
    """Calculate OrgP: New total organic phosphorus (mg-P/L).

    Args:
        dOrgPdt: Change in total organic phosphorus (mg-P/L/d)
        OrgP: Total organic phosphorus water concentration (mg-P/L),
        timestep: current iteration timestep (d)
    """     
    return OrgP + dOrgPdt * timestep

def TOP(
    use_OrgP: bool,
    OrgP: xr.DataArray,
    use_Algae: bool,
    rpa: xr.DataArray,
    Ap: xr.DataArray

) -> xr.DataArray :
    """Calculate TOP: Total organic phosphorus (mg-P/L).

    Args:
        use_OrgP: true/false to use organic phosphorus module (true/false),
        OrgP: New organic phosphorus water concentration (mg-P/L),
        use_Algae: true/false to use algae module (true/false),
        rpa: Algal P: Chla ratio (mg-P/ug-Chla),
        Ap: Algal water concentration (ug-Chla/L)
    """     
    TOP = 0.0
    TOP = xr.where(use_OrgP, TOP + OrgP,TOP)
    TOP = xr.where(use_Algae, TOP + rpa*Ap, TOP)

    return TOP


def TP(
    use_TIP: bool,
    TOP: xr.DataArray,
    TIP: xr.DataArray

) -> xr.DataArray :
    """Calculate TP: Total phosphorus (mg-P/L).

    Args:
        use_TIP: true/false to use total inorganic phosphorus module (true/false),
        TIP: New total inorganic phosphorus water concentration (mg-P/L),
        TOP: Total organic phosphorus water concentration (mg-P/L)
    """  
    TP = TOP
    TP = xr.where(use_TIP,TP + TIP,TP)

@numba.njit
def DIP(
    fdp: xr.DataArray,
    TIP: xr.DataArray

) -> xr.DataArray :
    """Calculate DIP: Dissolve inorganich phosphorus (mg-P/L).

    Args:
        fdp: fraction P dissolved
        TIP: New total inorganic phosphorus water concentration (mg-P/L),
    """
    return TIP * fdp


################################### From POM
@numba.njit
def kpom_tc(
    TwaterC: float,
    kpom_20: float,
) -> float:
    """Calculate the temperature adjusted POM dissolution rate (1/d)

    Args:
        TwaterC: Water temperature in Celsius
        kpom_20: POM dissolution rate at 20 degrees Celsius (1/d)
    """
    return arrhenius_correction(TwaterC, kpom_20, 1.047)


def POM_algal_settling(
    Ap: xr.DataArray,
    vsap: xr.DataArray,
    rda: xr.DataArray,
    depth: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculates the particulate organic matter concentration change due to algal mortality
    
    Args:
        Ap: Algae concentration (mg/L)
        vsap: Algal settling velocity (m/d)
        rda: Ratio of algal biomass to chlorophyll-a 
        depth: Depth of water in computation cell (m)
        use_Algae: Option to consider algal kinetics  
    """
    da: xr.DataArray = xr.where(use_Algae == True, vsap * Ap * rda / depth, 0)

    return da


@numba.njit
def POM_dissolution(
    POM: xr.DataArray,
    kpom_tc: xr.DataArray
) -> xr.DataArray:
    """Calculates the particulate organic matter concentration change due to POM dissolution

    Args:
        POM: Concentration of particulate organic matter (mg/L)
        kpom_tc: POM dissolution rate corrected for temperature (1/d)
    """

    return POM * kpom_tc


def POM_POC_settling(
    POC: xr.DataArray,
    vsoc: xr.DataArray,
    depth: xr.DataArray,
    fcom: xr.DataArray,
    use_POC: xr.DataArray
) -> xr.DataArray:
    """Calculates particulate organic matter concentration change due to POM settling
    
    Args:
        POC: Concentration of particulate organic carbon (mg/L)
        vsoc: POC settling velocity (m/d)
        depth: Depth of water (m)
        fcom: Fraction of carbon in organic matter (mg-C/mg-D) 
        use_POC: Option to consider particulate organic carbon
    """
    da: xr.DataArray = xr.where(use_POC == True, vsoc * POC / depth / fcom, 0)
    
    return da


def POM_benthic_algae_mortality(
    Ab: xr.DataArray,
    kdb_tc: xr.DataArray,
    Fb: xr.DataArray,
    Fw: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculates particulate organic matter concentration change due to benthic algae mortality
    
    Args:
        Ab: Benthic algae concentration (mg/L)
        kdb_tc: Benthic algae death rate (1/d)
        Fb: Fraction of bottom area available for benthic algae growth
        Fw: Fraction of benthic algae mortality into water column
        depth: Depth of water in computation cell (m)
        use_Balgae: Option for considering benthic algae in DOC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, Ab * kdb_tc * Fb * (1 - Fw) / depth, 0)

    return da


@numba.njit
def POM_burial(
    vb: xr.DataArray,
    POM: xr.DataArray,
    depth: xr.DataArray
) -> xr.DataArray:
    """Calculates particulate organic matter concentration change due to POM burial in the sediments
    
    Args:
        vb: Velocity of burial (m/d)
        POM: POM concentration (mg/L)
        depth: Depth of water in computation cell (m)
    """
    return vb * POM / depth


@numba.njit
def dPOMdt(
    POM_algal_settling: xr.DataArray,
    POM_dissolution: xr.DataArray,
    POM_POC_settling: xr.DataArray,
    POM_benthic_algae_mortality: xr.DataArray,
    POM_burial: xr.DataArray,
) -> xr.DataArray:
    """Calculates the concentration change of POM for one timestep

    Args:
        POM_algal_settling: POM concentration change due to algal settling (mg/L/d)
        POM_dissolution: POM concentration change due to dissolution (mg/L/d)
        POM_POC_settling: POM concentration change due to POC settling (mg/L/d)
        POM_benthic_algae_mortality: POM concentration change due to benthic algae mortality (mg/L/d)
        POM_burial: POM concentration change due to burial (mg/L/d)
    """
    return POM_algal_settling - POM_dissolution + POM_POC_settling + POM_benthic_algae_mortality - POM_burial


@numba.njit
def POM(
    dPOMdt: xr.DataArray,
    POM: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Computes updated particulate organic matter concentration (mg/L)
    
    Args:
        dPOMdt: Change in POM concentration over timestep (mg/L/d)
        POM: POM concentration from previous timestep (mg/L)
        timestep: Current iteration timestep (d)
    """
    return POM + dPOMdt * timestep


################################## From CBOD
@numba.njit
def kbod_tc(
    TwaterC: xr.DataArray,
    kbod_20: xr.DataArray,
) -> xr.DataArray:
    """Calculate the temperature adjusted CBOD oxidation rate (1/d)

    Args:
        TwaterC: water temperature in Celsius
        kbod_20: CBOD oxidation rate at 20 degrees Celsius (1/d)
    """

    kbod_tc = arrhenius_correction(TwaterC, kbod_20, 1.047)
    return kbod_tc


@numba.njit
def ksbod_tc(
    TwaterC: xr.DataArray,
    ksbod_20: xr.DataArray,
) -> xr.DataArray:
    """Calculate the temperature adjusted CBOD sedimentation rate (m/d)

    Args:
        TwaterC: water temperature in Celsius
        ksbod_20: CBOD sedimentation rate at 20 degrees Celsius (m/d)
    """

    ksbod_tc = arrhenius_correction(TwaterC, ksbod_20, 1.024)
    return ksbod_tc



def CBOD_oxidation(
    DOX: xr.DataArray,
    CBOD: xr.DataArray,
    kbod_tc: xr.DataArray,
    KsOxbod: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Calculates CBOD oxidation

    Args:
        DOX: Dissolved oxygen concentration (mg-O2/L)
        CBOD: Carbonaceous biochemical oxygen demand (mg-O2/L)
        kbod_tc: Temperature adjusted CBOD oxidation rate (1/d)
        KsOxbod: Half-saturation oxygen attenuation for CBOD oxidation (mg-O2/L)
        use_DOX: Option to consider DOX concentration in calculation of CBOD oxidation
    """
    da: xr.DataArray = xr.where(use_DOX == True, (DOX / (KsOxbod + DOX)) * kbod_tc * CBOD, kbod_tc * CBOD)
    
    return da


@numba.njit
def CBOD_sedimentation(
    CBOD: xr.DataArray,
    ksbod_tc: xr.DataArray
) -> xr.DataArray:
    """Calculates CBOD sedimentation for each group

    Args:
        CBOD: CBOD concentration (mg-O2/L)
        ksbod_tc: Temperature adjusted sedimentation rate (m/d)
    """
    
    CBOD_sedimentation = CBOD * ksbod_tc
    return CBOD_sedimentation


@numba.njit
def dCBODdt(
    CBOD_oxidation: xr.DataArray,
    CBOD_sedimentation: xr.DataArray
) -> xr.DataArray:
    """Computes change in each CBOD group for a given timestep

    Args:
        CBOD_oxidation: CBOD concentration change due to oxidation (mg/L/d)
        CBOD_sedimentation: CBOD concentration change due to sedimentation (mg/L/d)
    """
    return - CBOD_oxidation - CBOD_sedimentation


@numba.njit
def CBOD(
    CBOD: xr.DataArray,
    dCBODdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculates new CBOD concentration for next timestep

    Args:
        CBOD: CBOD concentration from previous timestep (mg/L)
        dCBODdt: CBOD concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return CBOD + dCBODdt * timestep

############################### From Carbon

@numba.njit
def kpoc_tc(
    TwaterC: xr.DataArray,
    kpoc_20: xr.DataArray,
) -> xr.DataArray:
    """Calculate the temperature adjusted POC hydrolysis rate (/d)

    Args:
        TwaterC: Water temperature in Celsius
        kpoc_20: POC hydrolysis rate at 20 degrees Celsius (1/d)
    """
    return arrhenius_correction(TwaterC, kpoc_20, 1.047)


@numba.njit
def POC_hydrolysis(
    kpoc_tc: xr.DataArray,
    POC: xr.DataArray,
) -> xr.DataArray:
    """Calculate the POC concentration change due to hydrolysis for a given timestep

    Args:
        kpoc_tc: POC hydrolysis rate at given water temperature (1/d)
        POC: POC concentration (mg/L)
    """
    return kpoc_tc * POC


@numba.njit
def POC_settling(
    vsoc: xr.DataArray,
    depth: xr.DataArray,
    POC: xr.DataArray
) -> xr.DataArray:
    """Calculate the POC concentration change due to settling for a given timestep

    Args:
        vsoc: POC settling velocity (m/d)
        depth: Water depth of cell (m)
        POC: POC concentration (mg/L)
    """
    return vsoc / depth * POC


def POC_algal_mortality(
    f_pocp: xr.DataArray,
    kdp_tc: xr.DataArray,
    rca: xr.DataArray,
    Ap: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculate the POC concentration change due to algal mortality

    Args:
        f_pocp: Fraction of algal mortality into POC
        kdp_tc: Algal death rate at water temperature (1/d)
        rca: Algal C to chlorophyll-a ratio (mg-C/ugChla)
        Ap: Algae concentration (mg/L)
        use_Algae: Option for considering algae in POC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Algae == True, f_pocp * kdp_tc * rca * Ap, 0)

    return da


def POC_benthic_algae_mortality(
    depth: xr.DataArray,
    f_pocb: xr.DataArray,
    kdb_tc: xr.DataArray,
    rcb: xr.DataArray,
    Ab: xr.DataArray,
    Fb: xr.DataArray,
    Fw: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculate the POC concentration change due to benthic algae mortality

    Args: 
        depth: Water depth in cell (m)
        f_pocb: Fraction of benthic algal mortality into POC
        kdb_tc: Benthic algae death rate (1/d)
        rcb: Benthic algae C to biomass weight ratio (mg-C/mg-D)
        Ab: Benthic algae concentration (mg/L)
        Fb: Fraction of bottom area available for benthic algae growth
        Fw: Fraction of benthic algae mortality into water column
        use_Balgae: Option for considering benthic algae in POC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, (1 / depth) * f_pocb * kdb_tc * rcb * Ab * Fb * Fw, 0)

    return da

@numba.njit
def dPOCdt(
    POC_settling: xr.DataArray,
    POC_hydrolysis: xr.DataArray,
    POC_algal_mortality: xr.DataArray,
    POC_benthic_algae_mortality: xr.DataArray
) -> xr.DataArray:
    """Calculate the change in POC concentration

    Args:
        POC_settling: Concentration change of POC due to settling (mg/L/d)
        POC_hydrolysis: Concentration change of POC due to hydrolysis (mg/L/d)
        POC_algal_mortality: Concentration change of POC due to algal mortality (mg/L/d)
        POC_benthic_algae_mortality: Concentration change of POC due to benthic algae mortality (mg/L/d)
    """
    return POC_algal_mortality + POC_benthic_algae_mortality - POC_settling - POC_hydrolysis


@numba.njit
def POC(
    POC: xr.DataArray,
    dPOCdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculate the POC concentration at the next time step

    Args:
        POC: Concentration of POC from previous timestep (mg/L)
        dPOCdt: POC concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return POC + dPOCdt * timestep


def DOC_algal_mortality(
    f_pocp: xr.DataArray,
    kdp_tc: xr.DataArray,
    rca: xr.DataArray,
    Ap: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculate the DOC concentration change due to algal mortality

    Args:
        f_pocp: Fraction of algal mortality into POC 
        kdp_tc: Algal death rate at water temperature (1/d) 
        rca: Algal C to chlorophyll-a ratio (mg-C/ug-Chla)
        Ap: Algae concentration (mg/L)
        use_Algae: Option for considering algae in DOC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Algae == True, (1 - f_pocp) * kdp_tc * rca * Ap, 0)

    return da


def DOC_benthic_algae_mortality(
    depth: xr.DataArray,
    f_pocb: xr.DataArray,
    kdb_tc: xr.DataArray,
    rcb: xr.DataArray,
    Ab: xr.DataArray,
    Fb: xr.DataArray,
    Fw: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculate the DOC concentration change due to benthic algae mortality

    Args: 
        depth: Water depth in cell (m)
        F_pocb: Fraction of benthic algal mortality into POC
        kdb_tc: Benthic algae death rate (1/d)
        rcb: Benthic algae C to biomass weight ratio (mg-C/mg-D)
        Ab: Benthic algae concentration (mg/L)
        Fb: Fraction of bottom area available for benthic algae growth
        Fw: Fraction of benthic algae mortality into water column
        use_Balgae: Option for considering benthic algae in DOC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, (1 / depth) * (1 - f_pocb) * kdb_tc * rcb * Ab * Fb * Fw, 0)

    return da


@numba.njit
def kdoc_tc(
    TwaterC: xr.DataArray,
    kdoc_20: xr.DataArray,
) -> xr.DataArray:
    """Calculate the temperature adjusted DOC oxidation rate (1/d)

    Args:
        TwaterC: Water temperature in Celsius
        kdoc_20: DOC oxidation rate at 20 degrees Celsius (1/d)
    """
    return arrhenius_correction(TwaterC, kdoc_20, 1.047)


def DOC_DIC_oxidation(
    DOX: xr.DataArray,
    KsOxmc: xr.DataArray,
    kdoc_tc: xr.DataArray,
    DOC: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Calculates the DOC concentration change due to oxidation

    Args:
        DOX: Concentration of dissolved oxygen (mg/L)
        KsOxmc: Half saturation oxygen attenuation constant for DOC oxidation rate (mg-O2/L)
        kdoc_tc: DOC oxidation rate (1/d)
        DOC: Concentration of dissolved organic carbon (mg/L)
        use_DOX: Option for considering dissolved oxygen concentration in DOC oxidation calculation (boolean)
    """
    da: xr.DataArray = xr.where(use_DOX == True, DOX / (KsOxmc + DOX) * kdoc_tc * DOC, kdoc_tc * DOC)

    return da


@numba.njit
def dDOCdt(
    DOC_DIC_oxidation: xr.DataArray,
    POC_hydrolysis: xr.DataArray,
    DOC_algal_mortality: xr.DataArray,
    DOC_benthic_algae_mortality: xr.DataArray
) -> xr.DataArray:
    """Calculates the change in DOC concentration

    Args:
        POC_hydrolysis: DOC concentration change due to POC hydrolysis (mg/L/d)
        DOC_POM_dissolution: DOC concentration change due to POM dissolution (mg/L/d)
        DOC_denitrification: DOC concentration change due to DOC denitrification (mg/L/d)
        DOC_algal_mortality: DOC concentration change due to algal mortality (mg/L/d)
        DOC_benthic_algae_mortality: DOC concentration change due to benthic algae mortality (mg/L/d)
        DOC_oxidation: DOC concentration change due to DOC oxidation (mg/L/d)
    """
    return POC_hydrolysis + DOC_algal_mortality + DOC_benthic_algae_mortality - DOC_DIC_oxidation


@numba.njit
def DOC(
    DOC: xr.DataArray,
    dDOCdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculate the DOC concentration at the next time step

    Args:
        DOC: Dissolved organic carbon concentration from previous timestep (mg/L)
        dDOCdt: Dissolved organic carbon concentration change for current timestep (mg/L/d)
        timestep: current iteration timestep (d)
    """
    return DOC + dDOCdt * timestep


@numba.njit
def Henrys_k(
    TwaterC: xr.DataArray
) -> xr.DataArray:
    """Calculates the temperature dependent Henry's coefficient (mol/L/atm)

    Args:
        TwaterC: Water temperature in celsius
    """
    return 10**(2385.73 / (TwaterC + 273.15) + .0152642 * (TwaterC + 273.15) - 14.0184)

@numba.njit
def Atmospheric_CO2_reaeration(
    ka_tc: xr.DataArray,
    K_H: xr.DataArray,
    pCO2: xr.DataArray,
    FCO2: xr.DataArray,
    DIC: xr.DataArray
) -> xr.DataArray:
    """Calculates the atmospheric input of CO2 into the waterbody

    Args:
        ka_tc: CO2 reaeration rate adjusted for temperature, same as O2 reaeration rate (1/d)
        K_H: Henry's Law constant (mol/L/atm)
        pCO2: Partial pressure of CO2 in the atmosphere (ppm)
        FCO2: Fraction of CO2 in total inorganic carbon
        DIC: Dissolved inorganic carbon concentration (mg/L)
    """
    return 12 * ka_tc * (10**-3 * K_H * pCO2 - 10**3 * FCO2 * DIC)


def DIC_algal_respiration(
    ApRespiration: xr.DataArray,
    rca: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC concentration change due to algal respiration

    Args:
        ApRespiration: Algae respiration calculated in algae module (ug-Chla/L/d)
        rca: Ratio of carbon to chlorophyll-a (mg-C/ug-Chla)
        use_Algae: Option to consider algae in the DIC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Algae == True, ApRespiration * rca, 0)

    return da


def DIC_algal_photosynthesis(
    ApGrowth: xr.DataArray,
    rca: xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC concentration change due to algal photosynthesis

    Args:
        ApGrowth: Algal photosynthesis calculated in algae module (ug-Chla/L/d)
        rca: Ratio of carbon to chlorophyll-a (mg-C/ug-Chla)
        use_Algae: Option to consider algae in the DIC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Algae == True, ApGrowth * rca, 0)

    return da


def DIC_benthic_algae_respiration(
    AbRespiration: xr.DataArray,
    rcb: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC flux due to benthic algae respiration

    Args:
        AbRespiration: Benthic algae respiration calculated in benthic algae module (g/m2/d)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Depth of water (m)
        use_Balgae: Option to consider benthic algae in the DIC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, AbRespiration * rcb * Fb * (1 / depth), 0)

    return da


def DIC_benthic_algae_photosynthesis(
    AbGrowth: xr.DataArray,
    rcb: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC flux due to benthic algae growth

    Args:
        AbGrowth: Benthic algae photosynthesis calculated in the benthic algae module (g/m2/d)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Depth of water (m)
        use_Balgae: Option to consider benthic algae in the DIC budget (boolean)
    """
    da: xr.DataArray = xr.where(use_Balgae == True, AbGrowth * rcb * Fb * (1 / depth), 0)

    return da


def DIC_CBOD_oxidation(
    DOX: xr.DataArray,
    CBOD: xr.DataArray,
    roc: xr.DataArray,
    kbod_tc: xr.DataArray, #imported from CBOD module
    KsOxbod: xr.DataArray, #imported from CBOD module
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Calculates DIC concentration change due to CBOD oxidation

    Args:
        DOX: Dissolved oxygen concentration (mg/L)
        CBOD: Carbonaceous biochemical oxygen demand concentration (mg/L) 
        roc: Ratio of O2 to carbon for carbon oxidation (mg-O2/mg-C)
        kbod_tc: CBOD oxidation rate (1/d)
        KsOxbod: Half saturation oxygen attenuation constant for CBOD oxidation (mg-O2/L)
        use_DOX: Option to consider dissolved oxygen in CBOD oxidation calculation (boolean)
    """
    
    da: xr.DataArray = xr.where(use_DOX == True, (1 / roc) * (DOX / (KsOxbod + DOX)) * kbod_tc * CBOD, CBOD * kbod_tc)

    return da


def DIC_sed_release(
    SOD_tc: xr.DataArray,
    roc: xr.DataArray,
    depth: xr.DataArray,

) -> xr.DataArray:
    """Computes the sediment release of DIC

    Args:
        SOD_tc: Sediment oxygen demand adjusted for water temperature (mg-O2/L/d)
        roc: Ratio of O2 to carbon for carbon oxidation (mg-O2/mg-C)
        depth: Water depth (m)
    """
    return SOD_tc / roc / depth


@numba.njit
def dDICdt(
    Atm_CO2_reaeration: xr.DataArray,
    DIC_algal_respiration: xr.DataArray,
    DIC_algal_photosynthesis: xr.DataArray,
    DIC_benthic_algae_respiration: xr.DataArray,
    DIC_benthic_algae_photosynthesis: xr.DataArray,
    DIC_CBOD_oxidation: xr.DataArray,
    DIC_sed_release: xr.DataArray
) -> xr.DataArray:
    """Calculates the change in DIC

    Args:
        Atm_CO2_reaeration: DIC concentration change due to atmospheric CO2 reaeration (mg/L/d)
        DIC_algal_respiration: DIC concentration change due to algal respiration (mg/L/d)
        DIC_algal_photosynthesis: DIC concentration change due to algal photosynthesis (mg/L/d)
        DIC_benthic_algae_respiration: DIC concentration change due to benthic algae respiration (mg/L/d)
        DIC_benthic_algae_photosynthesis: DIC concentration change due to benthic algae photosynthesis (mg/L/d)
        DIC_CBOD_oxidation: DIC concentration change due to CBOD oxidation (mg/L/d)
        DIC_sed_release: DIC concentration change due to sediment release (mg/L/d)
    """
    return Atm_CO2_reaeration + DIC_algal_respiration + DIC_benthic_algae_respiration + DIC_CBOD_oxidation + DIC_sed_release - DIC_algal_photosynthesis - DIC_benthic_algae_photosynthesis


@numba.njit
def DIC(
    DIC: xr.DataArray,
    dDICdt: xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray:
    """Calculate the DIC concentration at the next time step

    Args:
        DIC: Concentration of DIC from previous timestep (mg/L)
        dDICdt: Change in concentration of DIC for current timestep (mg/L/d)
        timestep: Current iteration timestep (d)
    """
    return DIC + dDICdt * timestep


######################################## From DOX


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


def DOX_DOC_oxidation(
    DOC_DIC_oxidation: xr.DataArray,
    roc: xr.DataArray,
    use_DOC: xr.DataArray
) -> xr.DataArray:
    """Computes dissolved oxygen flux due to oxidation of dissolved organic carbon

    Args:
        DOC_DIC_Oxidation: Dissolved organic carbon oxidation, calculated in carbon module (mg/L/d)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
    """
    da: xr.DataArray = xr.where(use_DOC == True, roc * DOC_DIC_oxidation, 0)

    return da


@numba.njit
def DOX_CBOD_oxidation(
    DIC_CBOD_oxidation: xr.DataArray,
    roc: xr.DataArray
) -> xr.DataArray:
    """Compute dissolved oxygen flux due to CBOD oxidation

    Args:
        DIC_CBOD_Oxidation: Carbonaceous biochemical oxygen demand oxidation, calculated in CBOD module (mg/L/d)
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
    """
    return DIC_CBOD_oxidation * roc


def DOX_AbGrowth(
    AbUptakeFr_NH4: xr.DataArray,
    roc: xr.DataArray,
    rcb: xr.DataArray,
    AbGrowth: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Compute dissolved oxygen flux due to benthic algae growth

    Args:
        AbUptakeFr_NH4: Fraction of actual benthic algal uptake that is from the ammonia pool, calculated in nitrogen module
        roc: Ratio of oxygen to carbon for carbon oxidation (mg-O2/mg-C)
        rcb: Benthic algae carbon to dry weight ratio (mg-C/mg-D)
        AbGrowth: Benthic algae photosynthesis, calculated in benthic algae module (mg/L/d)
        Fb: Fraction of bottom area available for benthic algae growth
        depth: Water depth (m)
        use_Balgae: Option to consider benthic algae in the DOX budget
    """
    da: xr.DataArray = xr.where(use_Balgae == True, (138 / 106 - 32 / 106 * AbUptakeFr_NH4) * roc * rcb * AbGrowth * Fb / depth, 0)

    return da


def DOX_AbRespiration(
    roc: xr.DataArray,
    rcb: xr.DataArray,
    AbRespiration: xr.DataArray,
    Fb: xr.DataArray,
    depth: xr.DataArray,
    use_Balgae: xr.DataArray
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

    da: xr.DataArray = xr.where(use_Balgae == True, roc * rcb * AbRespiration * Fb / depth, 0)

    return da


def DOX_SOD(
    depth: xr.DataArray,
    SOD_tc: xr.DataArray
) -> xr.DataArray:
    """Compute dissolved oxygen flux due to sediment oxygen demand

    Args:
        depth: Water depth (m)
        SOD_tc: Sediment oxygen demand not considering the SedFlux budget (mg-O2/m2)
    """

    return SOD_tc / depth

@numba.njit
def dDOXdt(
    Atm_O2_reaeration: xr.DataArray,
    DOX_ApGrowth: xr.DataArray,
    DOX_ApRespiration: xr.DataArray,
    DOX_Nitrification: xr.DataArray,
    DOX_DOC_oxidation: xr.DataArray,
    DOX_CBOD_oxidation: xr.DataArray,
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
        DOX_DOC_oxidation: DOX concentration change due to DOC oxidation (mg/L/d)
        DOX_CBOD_oxidation: DOX concentration change due to CBOD oxidation (mg/L/d)
        DOX_AbGrowth: DOX concentration change due to benthic algae photosynthesis (mg/L/d)
        DOX_AbRespiration: DOX concentration change due to benthic algae respiration (mg/L/d)
        DOX_SOD: DOX concentration change due to sediment oxygen demand (mg/L/d)
    """
    return Atm_O2_reaeration + DOX_ApGrowth - DOX_ApRespiration - DOX_Nitrification - DOX_DOC_oxidation - DOX_CBOD_oxidation + DOX_AbGrowth - DOX_AbRespiration - DOX_SOD


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

######################################### From pathogen


@numba.njit
def kdx_tc(
    TwaterC : xr.DataArray,
    kdx_20: xr.DataArray
) -> xr.DataArray :

    """Calculate kdx_tc: pathogen death rate (1/d).

    Args:
        TwaterC: Water temperature (C)
        kdx_20: Pathogen death rate at 20 degree (1/d)
    """

    return arrhenius_correction(TwaterC, kdx_20, 1.07)

@numba.njit
def PathogenDeath(
    kdx_tc : xr.DataArray,
    PX: xr.DataArray
) -> xr.DataArray :

    """Calculate PathogenDeath: pathogen natural death (cfu/100mL/d).

    Args:
      kdx_tc: pathogen death rate with temperature correction (1/d),
      PX: pathogen concentration (cfu/100mL)

    """
    return kdx_tc * PX

@numba.njit
def PathogenDecay(
    apx: xr.DataArray,
    q_solar: xr.DataArray,
    L: xr.DataArray,
    depth: xr.DataArray,
    PX: xr.DataArray
) -> xr.DataArray :

    """Calculate PathogenDecay: pathogen death due to light (cfu/100mL/d).

    Args:
      apx: light efficiency factor for pathogen decay,
      q_solar: Incident short-wave solar radiation (W/m2),
      L: lambda (1/m),
      depth: water depth (m),
      PX: Pathogen concentration (cfu/100mL)

    """
    return apx * q_solar / (L * depth) * (1 - math.exp(-L * depth)) * PX

@numba.njit
def PathogenSettling(
    vx: xr.DataArray,
    depth: xr.DataArray,
    PX: xr.DataArray
) -> xr.DataArray :

    """Calculate PathogenSettling: pathogen settling (cfu/100mL/d).

    Args:
      vx: pathogen net settling velocity (m)
      depth: water depth (m),
      PX: Pathogen concentration (cfu/100mL)
    """
    return vx/depth*PX

@numba.njit
def dPXdt(
    PathogenDeath: xr.DataArray,
    PathogenDecay: xr.DataArray,
    PathogenSettling: xr.DataArray

) -> xr.DataArray :

    """Calculate dPXdt: change in pathogen concentration (cfu/100mL/d).

    Args:
      PathogenSettling: pathogen settling (cfu/100mL/d)
      PathogenDecay: pathogen death due to light (cfu/100mL/d)
      PathogenDeath: pathogen natural death (cfu/100mL/d)

    """
    return -PathogenDeath - PathogenDecay - PathogenSettling

@numba.njit
def PX(
    PX:xr.DataArray,
    dPXdt: xr.DataArray,
    timestep: xr.DataArray

) -> xr.DataArray :

    """Calculate PX: New pathogen concentration (cfu/100mL).

    Args:
      dPXdt: change in pathogen concentration (cfu/100mL/d)
      PX: Pathogen concentration (cfu/100mL)
      timestep: Current iteration timestep (d)
    """
    return PX + timestep * dPXdt


##################################### From alkalinity

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
        kdnit_tc: Denitrification rate corrected for temperature (1/d)
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
        knit_tc: Nitrification rate corrected for temperature (1/d)
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
    ApUptakeFr_NH4 : xr.DataArray,
    use_Algae: xr.DataArray
) -> xr.DataArray:
    """Calculate the alkalinity concentration change due to algal growth

    Args:
        ApGrowth: Algal photosynthesis calculated in algae module (ug-Chla/L/d)
        r_alkaa: Ratio translating algal growth into Alk if NH4 is the N source (eq/ug-Chla)
        r_alkan: Ratio translating algal growth into Alk if NO3 is the N source (eq/ug-Chla)
        ApUptakeFr_NH4 : Preference fraction of algal N uptake from NH4
        use_Algae: Option to use algae
    """
    da: xr.DataArray = xr.where(use_Algae == True, (r_alkaa * ApUptakeFr_NH4  - r_alkan * (1 - ApUptakeFr_NH4 )) * ApGrowth, 0)

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
    AbUptakeFr_NH4 : xr.DataArray,
    Fb: xr.DataArray,
    use_Balgae: xr.DataArray
) -> xr.DataArray:
    """Calculate the alkalinity concentration change due to algal growth

    Args:
        ApGrowth: Algal photosynthesis calculated in algae module (ug-Chla/L/d)
        depth: Depth of water (m)
        r_alkaa: Ratio translating algal growth into Alk if NH4 is the N source (eq/ug-Chla)
        r_alkan: Ratio translating algal growth into Alk if NO3 is the N source (eq/ug-Chla)
        AbUptakeFr_NH4 : Preference fraction of benthic algae N uptake from NH4
        Fb: Fraction of bottom area available for benthic algae growth
        use_Balgae: Option to use benthic algae
    """
    da: xr.DataArray = xr.where(use_Balgae == True, (1 / depth) *(r_alkba * AbUptakeFr_NH4  - r_alkbn * (1 - AbUptakeFr_NH4 )) * AbGrowth * Fb, 0)

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
def Alk(
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

##################################### From N2


@numba.njit
def KHN2_tc(
    TwaterK : xr.DataArray,
) -> xr.DataArray :
    
    """Calculate Henry's law constant (mol/L/atm)
    
    Constant values found on NIST

    Args:
        TwaterK: water temperature kelvin (K)
        Henry's law constant for solubility in water at 298.15K: 0.00065 (mol/(kg*bar))
        Temperature dependence constant: 1300 (K) 
        Reference temperature: 298.15 (K) 
    """

    return 0.00065 * math.exp(1300.0 * (1.0 / TwaterK - 1 / 298.15))   
        
@numba.njit
def P_wv(
    TwaterK : xr.DataArray,
) -> xr.DataArray :
        
    """Calculate partial pressure water vapor (atm)

    Constant values found in documentation

    Args:
        TwaterK: water temperature kelvin (K)

    """
    return math.exp(11.8571  - (3840.70 / TwaterK) - (216961.0 / (TwaterK**2)))
  

def N2sat(
    KHN2_tc : xr.DataArray,
    pressure_atm: xr.DataArray,
    P_wv: xr.DataArray
) -> xr.DataArray:
    
    """Calculate N2 at saturation f(Twater and atm pressure) (mg-N/L)

    Args:
        KHN2_tc: Henry's law constant (mol/L/atm)
        pressure_atm: atmosphric pressure in atm (atm)
        P_wv: Partial pressure of water vapor (atm)
    """
        
    N2sat = 2.8E+4 * KHN2_tc * 0.79 * (pressure_atm - P_wv)  
    N2sat = xr.where(N2sat < 0.0,0.0,N2sat) #Trap saturation concentration to ensure never negative

    return N2sat

@numba.njit    
def dN2dt(
    ka_tc : xr.DataArray, 
    N2sat : xr.DataArray,
    N2: xr.DataArray,
) -> xr.DataArray: 
    
    """Calculate change in N2 air concentration (mg-N/L/d)

    Args:
        ka_tc: Oxygen re-aeration rate (1/d)
        N2sat: N2 at saturation f(Twater and atm pressure) (mg-N/L)
        N2: Nitrogen concentration air (mg-N/L)
    """
        
    return 1.034 * ka_tc * (N2sat - N2)

@numba.njit    
def N2(
    N2: xr.DataArray,
    dN2dt : xr.DataArray,
    timestep: xr.DataArray
) -> xr.DataArray: 
    
    """Calculate change in N2 air concentration (mg-N/L/d)

    Args:
        N2: Nitrogen concentration air (mg-N/L)
        dN2dt: Change in nitrogen concentration air
        timestep: Current iteration timestep (d)
    """
        
    return N2 + dN2dt * timestep

def TDG(
    N2: xr.DataArray,
    N2sat : xr.DataArray,
    DOX: xr.DataArray,
    DOX_sat: xr.DataArray,
    use_DOX: bool,
) -> xr.DataArray: 
    
    """Calculate total dissolved gas (%)

    Args:
        N2: Nitrogen concentration air (mg-N/L)
        N2sat: N2 at saturation f(Twater and atm pressure) (mg-N/L)
        DOX: Dissolved oxygen concentration (mg-O2/L)
        DOX_sat: O2 at saturation f(Twater and atm pressure) (mg-O2/L)
        use_DOX: true/false use dissolved oxygen module (true/false)
    """

    return xr.where(use_DOX,(79.0 * N2 / N2sat) + (21.0 * DOX / DOX_sat), N2/N2sat) 