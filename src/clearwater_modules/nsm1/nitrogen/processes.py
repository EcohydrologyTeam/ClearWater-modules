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

@numba.njit
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

@numba.njit
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

@numba.njit
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

@numba.njit
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

    return OrgN + dOrgNdt*timestep

@numba.njit
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

@numba.njit
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
    use_SedFlux: bool,
    JNH4: xr.DataArray,
    depth: xr.DataArray,
    rnh4_tc: xr.DataArray,

) -> xr.DataArray:
    """Calculate NH4fromBed: bed ->  NH4 (diffusion)    (mg-N/L/day)

    Args:
        use_SedFlux: true/false to use sediment flux module (unitless),
        depth: water depth (m),
        JNH4: Sediment water flux of ammonium (g-N/m^2/d),
        rnh4_tc: Sediment release rate of NH4 temperature correction(1/d).

    """

    return xr.where(use_SedFlux, JNH4 / depth, rnh4_tc / depth)

@numba.njit
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

@numba.njit
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

@numba.njit
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

@numba.njit
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

@numba.njit
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

    return NH4 + dNH4dt*timestep

@numba.njit
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
    use_SedFlux: bool,
    JNO3: xr.DataArray,
    depth: xr.DataArray,
    vno3_tc: xr.DataArray,
    NO3: xr.DataArray,

) -> xr.DataArray:
    """Calculate NO3_BedDenit: Sediment denitrification   (mg-N/L/day) 

    Args:
        use_SedFlux: true/false to use sediment flux module (unitless),
        depth: water depth (m),
        NO3: Nitrate concentration (mg-N/L),
        JNO3: Sediment water flux of nitrate (g-N/m^2/d),
        vno3_tc: Sediment denitrification velocity temperature correction (m/d)

    """
    
    return xr.where(use_SedFlux, JNO3 / depth,vno3_tc * NO3 / depth)

@numba.njit
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

@numba.njit
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


@numba.njit
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

    return NO3 + dNO3dt*timestep


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
    DIN = xr.where(use_NH4, DIN + NH4,DIN)
    DIN = xr.where(use_NO3, DIN + NO3, DIN)

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
    TON = xr.where(use_OrgN, TON + OrgN, TON)
    TON = xr.where(use_Algae, TON + rna * Ap, TON)
    
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
