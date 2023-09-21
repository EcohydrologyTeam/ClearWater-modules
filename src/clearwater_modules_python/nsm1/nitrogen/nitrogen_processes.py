"""
File containes process to calculate nitrogen species concentration and associated dependent variables
"""
import math
from clearwater_modules_python.shared.processes import arrhenius_correction
import numba


# local variables
# Nitrification Inhibitation (limits nitrification under low DO conditions)
NitrificationInhibition = 0
# fraction of actual floating algal uptake that is from ammonia pool
ApUptakeFr_NH4 = 0
# fraction of actual floating algal uptake that is from nitrate pool
ApUptakeFr_NO3 = 0
# fraction of actual benthic algal uptake that is from nitrate pool
AbUptakeFr_NO3 = 0
# fraction of actual benthic algal uptake that is from ammonia pool
AbUptakeFr_NH4 = 0

@numba.njit
def knit_tc(
    TwaterC : float,
    knit_20: float
) -> float :

    """Calculate knit_tc (1/d). #TODO only if use_NH4 = true

    Args:
        TwaterC: Water temperature (C)
        knit_20: Nitrification rate ammonia decay (1/d)
    """

    return arrhenius_correction(TwaterC, knit_20, 1.083)

@numba.njit
def rnh4_tc(
    TwaterC : float,
    rnh4_20: float
) -> float :

    """Calculate rnh4_tc (1/d). #TODO only if use_sedflux = true

    Args:
        TwaterC: Water temperature (C)
        rnh4_20: Sedimet release rate of NH4 (1/d)
    """

    return arrhenius_correction(TwaterC, rnh4_20, 1.074)

@numba.njit
def vno3_tc(
    TwaterC : float,
    vno3_20: float
) -> float :

    """Calculate vno3_tc (1/d). #TODO only if use_sedflux = true

    Args:
        TwaterC: Water temperature (C)
        vno3_20: Sedimet release rate of NO3 (1/d)
    """

    return arrhenius_correction(TwaterC, vno3_20, 1.08)

@numba.njit
def kon_tc(
    TwaterC : float,
    kon_20: float
) -> float :

    """Calculate kon_tc (1/d). #TODO only if use_OrgN = true

    Args:
        TwaterC: Water temperature (C)
        kon_20: Decay rate of OrgN to NH4 (1/d)
    """

    return arrhenius_correction(TwaterC, kon_20, 1.074)

@numba.njit
def kdnit_tc(
    TwaterC : float,
    kdnit_20: float
) -> float :

    """Calculate kdnit_tc (1/d). #TODO only if use_NO3 = true

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
    PN: float,
    NH4: float,
    NO3: float,

) -> float :

    """Calculate ApUptakeFr_NH4: Fraction of actual floating algal uptake from ammonia pool

    Args:
        use_NH4: use ammonium module (unitless)
        use_NO3: use nitrate module (unitless)
        use_Algae: use algae module (unitless)
        PN: NH4 preference factor algae (unitless)
        NH4: Ammonium water concentration (mg-N/L)
        NO3: Nitrate water concentration (mg-N/L)
    """

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
    ApUptakeFr_NH4: float
) -> float :

    """Calculate ApUptakeFr_NO3: Fraction of actual floating algal uptake from nitrate pool (unitless)

    Args:
        ApUptakeFr_NH4: Fraction of actual floating algal uptake from ammonia pool
    """

    return 1- ApUptakeFr_NH4

@numba.njit
def AbUptakeFr_NH4(
    use_NH4: bool,
    use_NO3: bool,
    use_Balgae: bool,
    PNb: float,
    NH4: float,
    NO3: float,

) -> float :

    """Calculate AbUptakeFr_NH4: Fraction of actual benthic algal uptake from ammonia pool

    Args:
        use_NH4: use ammonium module (unitless)
        use_NO3: use nitrate module (unitless)
        use_Balgae: use benthic algae module (unitless)
        PNb: NH4 preference factor benthic algae (unitless)
        NH4: Ammonium water concentration (mg-N/L)
        NO3: Nitrate water concentration (mg-N/L)
    """

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
    AbUptakeFr_NH4: float
) -> float :

    """Calculate AbUptakeFr_NO3: Fraction of actual benthic algal uptake from nitrate pool (unitless)

    Args:
        AbUptakeFr_NH4: Fraction of actual benthic algal uptake from ammonia pool
    """

    return 1- AbUptakeFr_NH4

"""
Organic Nitrogen                     (mg-N/d*L)
dOrgN/dt =   Algae_OrgN              (Floating Algae -> OrgN)
            - OrgN_NH4_Decay        (OrgN -> NH4)
            - OrgN_Settling         (OrgN -> bed)
            + Benthic Death         (Benthic Algae -> OrgN)														 
"""
@numba.njit
def dOrgNdt(
    use_OrgN: bool,
    use_Algae: bool,
    use_Balgae: bool,
    kon_tc: float,
    OrgN: float,
    vson: float,
    depth: float,
    rna: float,
    rnb: float,
    Fw: float,
    Fb: float,
    ApDeath: float,
    AbDeath: float,

) -> float :

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

"""
Ammonia Nitrogen (NH4)                 (mg-N/day*L)
dNH4/dt   =    OrgN_NH4_Decay          (OrgN -> NH4)  
                - NH4 Oxidation         (NH4 -> NO3)
                - NH4AlgalUptake        (NH4 -> Floating Algae)
                + Benthos NH4           (Benthos -> NH4)
                - Benthic Algae Uptake  (NH4 -> Benthic Algae)														 
"""
# Compute nitrification inhibition coefficient used to retard oxidation
# rate in case of low dissolved oxygen.  The coefficient should range between zero and one.

# Modify Nitrogren Inhibition Factor. If the function is disabled, make the function linear at DO concentrations
# greater than zero, and shut off nitrogen oxidation completely once DO is depleted
if self.global_module_choices['use_NH4']:  # TODO this looks different
    if self.global_module_choices['use_DOX']:
        NitrificationInhibition = 1.0 - \
            math.exp(-self.nitrogen_constant['KNR']
                        * self.global_vars['DOX'])
    else:
        NitrificationInhibition = 1.0

    NH4_Nitrification = NitrificationInhibition * \
        knit_tc * self.global_vars['NH4']

    if self.global_module_choices['use_SedFlux']:
        NH4fromBed = self.sedFlux_pathways['JNH4'] / \
            self.global_vars['depth']
    else:
        NH4fromBed = rnh4_tc / self.global_vars['depth']

    if self.global_module_choices['use_Algae']:
        NH4_ApRespiration = self.algae_pathways['rna'] * \
            self.algae_pathways['ApRespiration']
        NH4_ApGrowth = ApUptakeFr_NH4 * \
            self.algae_pathways['rna'] * \
            self.algae_pathways['ApGrowth']
    else:
        NH4_ApRespiration = 0.0
        NH4_ApGrowth = 0.0

    if self.global_module_choices['use_BAlgae']:
        # TODO changed the calculation for respiration from the inital FORTRAN due to conflict with the reference guide
        NH4_AbRespiration = self.Balgae_pathways['rnb'] * \
            self.Balgae_pathways['AbRespiration']
        NH4_AbGrowth = (AbUptakeFr_NH4 * self.Balgae_pathways['rnb'] * self.Balgae_pathways['Fb']
                        * self.Balgae_pathways['AbGrowth']) / self.global_vars['depth']
    else:
        NH4_AbRespiration = 0.0
        NH4_AbGrowth = 0.0

    if not self.global_module_choices['use_OrgN']:
        OrgN_NH4_Decay = 0.0

    dNH4dt = OrgN_NH4_Decay - NH4_Nitrification + NH4fromBed + \
        NH4_ApRespiration - NH4_ApGrowth + NH4_AbRespiration - NH4_AbGrowth
else:
    dNH4dt = 0

"""
Nitrite Nitrogen  (NO3)                       (mg-N/day*L)
dNO3/dt  =      NH4 Oxidation                 (NH4 -> NO3) 
                - NO3 Sediment Denitrification
                - NO3 Algal Uptake            (NO3-> Floating Algae) 
                - NO3 Benthic Algal Uptake    (NO3-> Benthic  Algae) 
"""
if self.global_module_choices['use_NO3']:
    if self.global_module_choices['use_DOX']:
        NO3_Denit = (1.0 - (self.global_vars['DOX'] / (self.global_vars['DOX'] +
                        self.nitrogen_constant['KsOxdn']))) * kdnit_tc * self.global_vars['NO3']
        if math.isnan(NO3_Denit):
            NO3_Denit = kdnit_tc * self.global_vars['NO3']
    else:
        NO3_Denit = 0.0

    if self.global_module_choices['use_SedFlux']:
        NO3_BedDenit = self.sedFlux_pathways['JNO3'] / \
            self.global_vars['depth']
    else:
        NO3_BedDenit = vno3_tc * \
            self.global_vars['NO3'] / self.global_vars['depth']

    if self.global_module_choices['use_Algae']:
        NO3_ApGrowth = ApUptakeFr_NO3 * \
            self.algae_pathways['rna'] * \
            self.algae_pathways['ApGrowth']
    else:
        NO3_ApGrowth = 0.0

    if self.global_module_choices['use_BAlgae']:
        NO3_AbGrowth = (AbUptakeFr_NO3 * self.Balgae_pathways['rnb'] * self.Balgae_pathways['Fb']
                        * self.Balgae_pathways['AbGrowth']) / self.global_vars['depth']
    else:
        NO3_AbGrowth = 0.0

    if not self.global_module_choices['use_NH4']:
        NH4_Nitrification = 0.0

    dNO3dt = NH4_Nitrification - NO3_Denit - \
        NO3_BedDenit - NO3_ApGrowth - NO3_AbGrowth
else:
    dNO3dt = 0

# Derived variables calculations
DIN = 0.0
TON = 0.0
TKN = 0.0
if self.global_module_choices['use_NH4']:
    DIN = DIN + self.global_vars['NH4']
    TKN = TKN + self.global_vars['NH4']

if self.global_module_choices['use_NO3']:
    DIN = DIN + self.global_vars['NO3']

if self.global_module_choices['use_OrgN']:
    TON = TON + self.global_vars['OrgN']

if self.global_module_choices['use_Algae']:
    TON = TON + self.algae_pathways['rna'] * self.global_vars['Ap']

TKN = TKN + TON
TN = DIN + TON

nitrogen_pathways = {
    'NH4_Nitrification' : NH4_Nitrification,
    'NO3_Denit' : NO3_Denit,
    'ApUptakeFr_NH4': ApUptakeFr_NH4,
    'ApUptakeFr_NO3' : ApUptakeFr_NO3,
    'AbUptakeFr_NH4': AbUptakeFr_NH4,
}

print("dNH4dt", dNH4dt)
print("dNO3dt", dNO3dt)
print("dOrgNdt", dOrgNdt)

print("DIN", DIN)
print("TON", TON)
print("TKN", TKN)
print("TN", TN)

return DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt, nitrogen_pathways
