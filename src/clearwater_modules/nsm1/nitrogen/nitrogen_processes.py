"""
File contains process to calculate nitrogen species concentration and associated dependent variables
"""
import math
from clearwater_modules.shared.processes import arrhenius_correction
import numba
import xarray as xr

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


@numba.njit
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
    if use_NH4 and not use_NO3:
        ApUptakeFr_NH4 = 1.0
    if not use_NH4 and use_NO3:
        ApUptakeFr_NH4 = 0.0
    if not use_NH4 and not use_NO3:
        ApUptakeFr_NH4 = 0.5

    #  Calculating Nitrogen Kinetics
    if use_Algae and use_NH4 and use_NO3:
        ApUptakeFr_NH4 = PN * NH4 / (PN * NH4 + (1.0 - PN) * NO3)
    # Check for case when NH4 and NO3 are very small.  If so, force uptake_fractions appropriately.
        if math.isnan(ApUptakeFr_NH4):
            ApUptakeFr_NH4 = PN

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


@numba.njit
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

    # set value of UptakeFr_NH4/NO3 for special conditions
    if use_NH4 and not use_NO3:
        AbUptakeFr_NH4 = 1.0
    if not use_NH4 and use_NO3:
        AbUptakeFr_NH4 = 0.0
    if not use_NH4 and not use_NO3:
        AbUptakeFr_NH4 = 0.5

    # Check for benthic and recompute if necessary
    if use_Balgae and use_NH4 and use_NO3:
        AbUptakeFr_NH4 = (PNb * NH4) / (PNb * NH4 + (1.0 - PNb) * NO3)
        AbUptakeFr_NO3 = 1 - AbUptakeFr_NH4

    # Check if NH4 and NO3 are very small.  If so, force uptake_fractions appropriately.
        if math.isnan(AbUptakeFr_NH4):
            AbUptakeFr_NH4 = PNb

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


@numba.njit
def dOrgNdt(
    use_OrgN: bool,
    use_Algae: bool,
    use_Balgae: bool,
    kon_tc: xr.DataArray,
    OrgN: xr.DataArray,
    vson: xr.DataArray,
    depth: xr.DataArray,
    rna: xr.DataArray,
    rnb: xr.DataArray,
    Fw: xr.DataArray,
    Fb: xr.DataArray,
    ApDeath: xr.DataArray,
    AbDeath: xr.DataArray,

) -> xr.DataArray:
    """Calculate dOrgNdt: Change in Organic Nitrogen (mg-N/L)

    Args:
        use_OrgN: true/false to use organic nitrogen module (unitless)
        use_Algae: true/false to use algae module (unitless)
        use_Balgae: true/false to use benthic algae module (unitless),
        kon_tc: Decay rate of organic nitrogen to nitrate with temperature correction (1/d),
        OrgN: Concentration of organic nitrogen (mg-N/L)
        vson: Organic nitrogen settling velocity (m/d)
        depth: water depth (m)
        rna: Algal N: Chla ratio (mg-N/ug-Chla)
        rnb: Benthic algal N: Benthic Algal Dry Weight (mg-N/mg-D)
        Fw: Fraction benthic algae mortality into water column (unitless)
        Fb: Fraction of bottom area for benthic algae (unitless)
        ApDeath: Algal death rate (ug-Chla/L/d)
        AbDeath: Benthic algal death rate (g/m^2/d)

    Organic Nitrogen                     (mg-N/d*L)
    dOrgN/dt =   Algae_OrgN              (xr.DataArraying Algae -> OrgN)
                - OrgN_NH4_Decay         (OrgN -> NH4)
                - OrgN_Settling          (OrgN -> bed)
                + Benthic Death          (Benthic Algae -> OrgN)	

    """
    if use_OrgN:
        OrgN_NH4_Decay = kon_tc * OrgN
        OrgN_Settling = vson / depth * OrgN

        if use_Algae:
            ApDeath_OrgN = rna * ApDeath
        else:
            ApDeath_OrgN = 0.0

        if use_Balgae:
            AbDeath_OrgN = rnb * Fw * Fb * AbDeath / depth
        else:
            AbDeath_OrgN = 0.0

        dOrgNdt = ApDeath_OrgN + AbDeath_OrgN - OrgN_NH4_Decay - OrgN_Settling
    else:
        dOrgNdt = 0

    return dOrgNdt


@numba.njit
def dNH4dt(
    use_OrgN: bool,
    use_Algae: bool,
    use_Balgae: bool,
    use_DOX: bool,
    use_SedFlux: bool,
    use_NH4: bool,
    depth: xr.DataArray,
    rna: xr.DataArray,
    rnb: xr.DataArray,

    Fb: xr.DataArray,
    KNR: xr.DataArray,
    DOX: xr.DataArray,
    NH4: xr.DataArray,
    JNH4: xr.DataArray,

    ApRespiration: xr.DataArray,
    ApGrowth: xr.DataArray,
    AbRespiration: xr.DataArray,
    AbGrowth: xr.DataArray,

) -> xr.DataArray:
    """Calculate dNH4dt: Change in Ammonium (mg-N/L)

    Args:
        use_OrgN: true/false to use organic nitrogen module (unitless),
        use_Algae: true/false to use algae module (unitless),
        use_Balgae: true/false to use benthic algae module (unitless),
        use_DOX: true/false to use dissolve oxygen module (unitless),
        use_SedFlux: true/false to use sediment flux module (unitless),
        use_NH4: true/false to use ammonium module (unitless),

        depth: water depth (m),
        rna: Algal N: Chla ratio (mg-N/ug-Chla),
        rnb: Benthic algal N: Benthic Algal Dry Weight (mg-N/mg-D),

        Fb: Fraction of bottom area for benthic algae (unitless),
        KNR: Oxygen inhabitation factor for nitrification (mg-O2/L),
        DOX: Dissolved oxygen concentration (mg-O2/L),
        NH4: Ammonium concentration (mg-N/L),
        JNH4: Sediment water flux of ammonium (g-N/m^2/d),

        ApRespiration: Algal respiration rate (ug-Chla/L/d),
        ApGrowth: Algal growth rate (ug-Chla/L/d),
        AbRespiration: Benthic algal respiration rate (g/m^2/d),
        AbGrowth: Benthic alga growth rate (g/m^2/d),    

    Ammonia Nitrogen (NH4)                 (mg-N/day*L)
    dNH4/dt   =    OrgN_NH4_Decay          (OrgN -> NH4)  
                    - NH4 Oxidation         (NH4 -> NO3)
                    - NH4AlgalUptake        (NH4 -> xr.DataArraying Algae)
                    + Benthos NH4           (Benthos -> NH4)
                    - Benthic Algae Uptake  (NH4 -> Benthic Algae)		

    """
    NitrificationInhibition = 0

    if use_NH4:
        if use_DOX:
            NitrificationInhibition = 1.0 - math.exp(-KNR * DOX)
        else:
            NitrificationInhibition = 1.0

        NH4_Nitrification = NitrificationInhibition * knit_tc * NH4

        if use_SedFlux:
            NH4fromBed = JNH4 / depth
        else:
            NH4fromBed = rnh4_tc / depth

        if use_Algae:
            NH4_ApRespiration = rna * ApRespiration
            NH4_ApGrowth = ApUptakeFr_NH4 * rna * ApGrowth
        else:
            NH4_ApRespiration = 0.0
            NH4_ApGrowth = 0.0

        if use_Balgae:
            # TODO changed the calculation for respiration from the inital FORTRAN due to conflict with the reference guide
            NH4_AbRespiration = rnb * AbRespiration
            NH4_AbGrowth = (AbUptakeFr_NH4 * rnb * Fb * AbGrowth) / depth
        else:
            NH4_AbRespiration = 0.0
            NH4_AbGrowth = 0.0

        if not use_OrgN:
            OrgN_NH4_Decay = 0.0

        dNH4dt = OrgN_NH4_Decay - NH4_Nitrification + NH4fromBed + \
            NH4_ApRespiration - NH4_ApGrowth + NH4_AbRespiration - NH4_AbGrowth
    else:
        dNH4dt = 0

    return dNH4dt


@numba.njit
def dNO3dt(
    use_Algae: bool,
    use_Balgae: bool,
    use_DOX: bool,
    use_SedFlux: bool,
    use_NH4: bool,
    use_NO3: bool,

    depth: xr.DataArray,
    rna: xr.DataArray,
    rnb: xr.DataArray,
    KsOxdn: xr.DataArray,

    Fb: xr.DataArray,
    DOX: xr.DataArray,
    NO3: xr.DataArray,
    JNO3: xr.DataArray,

    ApGrowth: xr.DataArray,
    AbGrowth: xr.DataArray,

) -> xr.DataArray:
    """Calculate dNO3dt: Change in nitrate (mg-N/L)

    Args:
        use_Algae: true/false to use algae module (unitless),
        use_Balgae: true/false to use benthic algae module (unitless),
        use_DOX: true/false to use dissolve oxygen module (unitless),
        use_SedFlux: true/false to use sediment flux module (unitless),
        use_NH4: true/false to use ammonium module (unitless),
        use_NO3: true/false to use nitrate module (unitless),

        depth: water depth (m),
        rna: Algal N: Chla ratio (mg-N/ug-Chla),
        rnb: Benthic algal N: Benthic Algal Dry Weight (mg-N/mg-D),
        KsOxdn: Half-saturation oxygen inhibition constant for denitrification (mg-O2/L)

        Fb: Fraction of bottom area for benthic algae (unitless),
        DOX: Dissolved oxygen concentration (mg-O2/L),
        NO3: Nitrate concentration (mg-N/L),
        JNO3: Sediment water flux of nitrate (g-N/m^2/d),

        ApGrowth: Algal growth rate (ug-Chla/L/d),
        AbGrowth: Benthic alga growth rate (g/m^2/d),    

    Nitrite Nitrogen  (NO3)                       (mg-N/day*L)
    dNO3/dt  =      NH4 Oxidation                 (NH4 -> NO3) 
                    - NO3 Sediment Denitrification
                    - NO3 Algal Uptake            (NO3-> xr.DataArraying Algae) 
                    - NO3 Benthic Algal Uptake    (NO3-> Benthic  Algae) 	

    """

    if use_NO3:
        if use_DOX:
            NO3_Denit = (1.0 - (DOX / (DOX + KsOxdn))) * kdnit_tc * NO3
            if math.isnan(NO3_Denit):
                NO3_Denit = kdnit_tc * NO3
        else:
            NO3_Denit = 0.0

        if use_SedFlux:
            NO3_BedDenit = JNO3 / depth
        else:
            NO3_BedDenit = vno3_tc * NO3 / depth

        if use_Algae:
            NO3_ApGrowth = ApUptakeFr_NO3 * rna * ApGrowth
        else:
            NO3_ApGrowth = 0.0

        if use_Balgae:
            NO3_AbGrowth = (AbUptakeFr_NO3 * rnb * Fb * AbGrowth) / depth
        else:
            NO3_AbGrowth = 0.0

        if not use_NH4:
            NH4_Nitrification = 0.0

        dNO3dt = NH4_Nitrification - NO3_Denit - \
            NO3_BedDenit - NO3_ApGrowth - NO3_AbGrowth
    else:
        dNO3dt = 0

    return dNO3dt


@numba.njit
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
    if use_NH4:
        DIN = DIN + NH4

    if use_NO3:
        DIN = DIN + NO3

    return DIN


@numba.njit
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
    if use_OrgN:
        TON = TON + OrgN

    if use_Algae:
        TON = TON + rna * Ap

    return TON


@numba.njit
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
    if use_NH4:
        TKN = TKN + NH4

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
