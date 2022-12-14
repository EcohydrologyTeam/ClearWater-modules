'''
=======================================================================================
Contaminant Simulation Module (CSM): Main Program
=======================================================================================

Developed by:
* Dr. Todd E. Steissberg (ERDC-EL)
* Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

This module computes the water quality of a single computational cell. The algorithms 
and structure of this program were adapted from the Fortran 95 version of this module, 
developed by:
* Dr. Billy E. Johnson (ERDC-EL)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

Version 1.0

Initial Version: June 15, 2021
Last Revision Date: June 15, 2021
'''

# ----------------------------------------------------------------------------------------------------------------------------------
# ZZ: a check needs to be included in GUI to avoid???
# all desorption rates can not be 0!
# Porosity cannot be 0!
# ka cannot be 0, otherwise KL will be 0, vv becomes infinity.
# pH cannot to 0 or 14, or else, CHH and COH will be zero.
# also MW cannot be 0!
# h2, depth, volume cannot be 0!
# I0pht, Lambdamax, KH, dt cannot be 0!
# kdeap, kdepom, kdepom2, kdep, kdep2 cannot be 0!
# ----------------------------------------------------------------------------------------------------------------------------------

from EnvironmentalSystems.ClearWater.wq_modules.python.src.water_quality_functions import TempCorrection
from EnvironmentalSystems.ClearWater.wq_modules.python.src.CSM.csm_globals import global_vals as globvals
import math


def csv():
    # Temperature correction parameters (47 = 4 + 10 + 9 + 2 + 1 + 12 + 9)

    # Ionization equilibrium constant for the formation
    # of singly charged anionic species (unitless)
    Ka1 = 0.0
    Ka1_tc = 0.0
    ka12_tc = 0.0

    # Ionization equilibrium constant for the formation
    # of doubly charged anionic species (unitless)
    Ka2 = 0.0
    Ka2_tc = 0.0
    ka22_tc = 0.0

    # Ionization equilibrium constant for the formation
    # of singly charged cationic species (unitless)
    Kb1 = 0.0
    Kb1_tc = 0.0
    kb12_tc = 0.0

    # Ionization equilibrium constant for the formation
    # of doubly charged cationic species (unitless)
    Kb2_tc = 0.0
    kb22_tc = 0.0

    # Adsorption rate of POM sorbed contaminant i in water column (L/µg/day)
    kadpom = 0.0
    kadpom_tc = 0.0

    # Adsorption rate of POM sorbed contaminant i in bed sediment (L/µg/day)
    kadpom2 = 0.0
    kadpom2_tc = 0.0

    # Desorption rate of POM sorbed contaminant i in water column (1/day)
    kdepom = 0.0
    kdepom_tc = 0.0

    # Desorption rate of POM sorbed contaminant i in bed sediment (1/day)
    kdepom2 = 0.0
    kdepom2_tc = 0.0

    # Adsorption rate of algae sorbed contaminant i in water column (L/µg/day)
    kadap = 0.0
    kadap_tc = 0.0

    # Desorption rate of algae sorbed contaminant i (1/day)
    kdeap = 0.0
    kdeap_tc = 0.0

    # Adsorption rate of solid j sorbed contaminant i in water column (L/µg/day)
    kadp = 0.0
    kadp_tc = 0.0

    # Adsorption rate of solid j sorbed contaminant i in bed sediment (L/µg/day)
    kadp2 = 0.0
    kadp2_tc = 0.0

    # Desorption rate of solid j sorbed contaminant i in water column (1/day)
    kdep = 0.0
    kdep_tc = 0.0

    # Desorption rate of solid j sorbed contaminant i in bed sediment (1/day)
    kdep2 = 0.0
    kdep2_tc = 0.0

    # n order decay rate of dissolved contaminant i in water column  (1/day)
    k1d = TempCorrection(0.0, 0.0)  # TODO Need to define

    # n order decay rate of dissolved contaminant i in sediment (1/day)
    k1d2 = TempCorrection(0.0, 0.0)  # TODO Need to define

    # n order decay rate of DOC sorbed contaminant i in water column (1/day)
    k1doc = TempCorrection(0.0, 0.0)  # TODO Need to define

    # n order decay rate of DOC sorbed contaminant i in sediment (1/day)
    k1doc2 = TempCorrection(0.0, 0.0)  # TODO Need to define

    # n order decay rate of algae sorbed contaminant i in water column (1/day)
    k1ap = TempCorrection(0.0, 0.0)  # TODO Need to define

    # n order decay rate of POM sorbed contaminant i in water column (1/day)
    k1pom = TempCorrection(0.0, 0.0)  # TODO Need to define

    # n order decay rate of POM sorbed contaminant i in sediment (1/day)
    k1pom2 = TempCorrection(0.0, 0.0)  # TODO Need to define

    # n order decay rate of solids sorbed contaminant i in water column (1/day)
    k1p = TempCorrection(0.0, 0.0)  # TODO Need to define

    # n order decay rate of solids sorbed contaminant i in bed sediment (1/day)
    k1p2 = TempCorrection(0.0, 0.0)  # TODO Need to define

    # Aquatic photolysis rate of dissolved contaminant i in water column (1/day)
    k1phtd = TempCorrection(0.0, 0.0)  # TODO Need to define

    # Aquatic photolysis rate of DOC sorbed contaminant i in water column (1/day)
    kphtdoc = TempCorrection(0.0, 0.0)  # TODO Need to define

    # Volatilization velocity of dissolved contaminant i in water column   (m/day)
    vv = TempCorrection(0.0, 0.0)  # TODO Need to define

    # Alkaline hydrolysis rate of dissolved contaminant i in water column (L/mol/day)
    khbd = 0.0
    khbd_tc = 0.0

    # Alkaline hydrolysis rate of dissolved contaminant i in bed sediment (L/mol/day)
    khbd2 = 0.0
    khbd2_tc = 0.0

    # Alkaline hydrolysis rate of DOC sorbed contaminant i in water column (L/mol/day)
    khbdoc_tc = 0.0
    khbdoc2 = 0.0

    # Alkaline hydrolysis rate of DOC sorbed contaminant i in bed sediment (L/mol/day)
    khbdoc2_tc = 0.0

    # Neutral hydrolysis rate of dissolved contaminant i in water column (1/day)
    khnd = 0.0
    khnd_tc = 0.0

    # Neutral hydrolysis rate of dissolved contaminant i in bed sediment (1/day)
    khnd2 = 0.0
    khnd2_tc = 0.0

    # Neutral hydrolysis rate of DOC sorbed contaminant i in water column (1/day)
    khndoc = 0.0
    khndoc_tc = 0.0

    # Neutral hydrolysis rate of DOC sorbed contaminant i in bed sediment (1/day)
    khndoc2 = 0.0
    khndoc2_tc = 0.0

    # Acid hydrolysis rate of dissolved contaminant i in water column (L/mol/day)
    khad = 0.0
    khad_tc = 0.0

    # Acid hydrolysis rate of dissolved contaminant i in bed sediment (L/mol/day)
    khad2 = 0.0
    khad2_tc = 0.0

    # Acid hydrolysis rate of DOC sorbed contaminant i in water column (L/mol/day)
    khadoc = 0.0
    khadoc_tc = 0.0

    # Acid hydrolysis rate of DOC sorbed contaminant i in bed sediment (L/mol/day)
    khadoc2 = 0.0
    khadoc2_tc = 0.0

    # 2nd reaction rate of dissolved contaminant i in water column (L/mg/day)
    kerd = 0.0
    kerd_tc = 0.0

    # 2nd reaction rate of dissolved contaminant i in sediment (L/mg/day)
    kerd2 = 0.0
    kerd2_tc = 0.0

    # 2nd reaction rate of DOC sorbed contaminant i in water column (L/mg/day)
    kerdoc = 0.0
    kerdoc_tc = 0.0

    # 2nd reaction rate of DOC sorbed contaminant i in sediment (L/mg/day)
    kerdoc2 = 0.0
    kerdoc2_tc = 0.0

    # 2nd reaction rate of algae sorbed contaminant i in water column (L/mg/day)
    kerap = 0.0
    kerap_tc = 0.0

    # 2nd reaction rate of POM sorbed contaminant i in water column (L/mg/day)
    kerpom = 0.0
    kerpom_tc = 0.0

    # 2nd reaction rate of POM sorbed contaminant i in sediment (L/mg/day)
    kerpom2 = 0.0
    kerpom2_tc = 0.0

    # 2nd reaction rate of solids sorbed contaminant i in water column (L/mg/day)
    kerp = 0.0
    kerp_tc = 0.0

    # 2nd reaction rate of solids sorbed contaminant i in bed sediment (L/mg/day)
    kerp2 = 0.0
    kerp2_tc = 0.0

    # Global parameters

    # 1D parameters (10)
    # Active Sediment layer thickness                                       (m)
    h2 = 0.0
    # Porosity or volume water per volume bed sediments                     (unitless)
    Por = 0.0
    # the average bioturbed depth in bed sediments                          (cm)
    z2 = 0.0
    # Biodiffusion coefficient representing particle diffusivity in the bed (cm2/d)
    Db = 0.0
    # Water-side benthic boundary layer mass transfer coefficient           (cm/d)
    Beta = 0.0
    # density of total bed sediments                                        (g/cm3)
    ps = 0.0
    # Settling velocity of algae in water column                            (m/d)
    vsap = 0.0
    # Settling velocity of POM in water column                              (m/d)
    vsom = 0.0
    # Coefficient to adjust light extinction coefficient                    (unitless)
    alpha = 0.0
    # Maximum relative error of Newton-Raphson or Bisection method          (unitless)
    res = 0.0

    # 2D or 3D parameters associated with specified contaminant (62 = 5 + 38 + 13 + 6)
    # Solubility of contaminant i                                           (mg/L)
    Csd = 0.0
    # Molecular weight of contaminant i                                     (g/mole)
    MW = 0.0
    # Molecular diffusivity of contaminant i                                (m2/s)
    Dm = 0.0
    # Sediment-water mass transfer velocity                                 (m/d)
    vm = 0.0
    # Henry's constant of contaminant i                                     (Pa m3/mol)
    KH = 0.0

    # Octanol-water partitioning coefficient                                             (unitless)
    Kow = 0.0
    # Dissolved organic carbon partition coefficient of contaminant i in water column    (L/kg)
    Kdoc = 0.0
    # Dissolved organic carbon partition coefficient of contaminant i in bed sediment    (L/kg)
    Kdoc2 = 0.0
    # Correction coefficient to compute Kdoc of contaminant i from Kow in water column   (unitless)
    adoc = 0.0
    # Correction coefficient to compute Kdoc2 of contaminant i from Kow in bed sediment  (unitless)
    adoc2 = 0.0
    # Particulate organic matter partition coefficient of contaminant i in water column  (L/kg)
    Kpom = 0.0
    # Particulate organic matter partition coefficient of contaminant i in bed sediment  (L/kg)
    Kpom2 = 0.0
    # Correction coefficient to compute Kpom of contaminant i from Kow in water column   (unitless)
    apom = 0.0
    # Correction coefficient to compute Kpom2 of contaminant i from Kow in bed sediment  (unitless)
    apom2 = 0.0
    # Langmuir adsorption constant of contaminant i on organic matter in water column    (L/mg)
    Klpom = 0.0
    # Langmuir adsorption constant of contaminant i on organic matter in bed sediment    (L/mg)
    Klpom2 = 0.0
    # Adsorption capacity of contaminant i on POM in water column                        (µg/g)
    qcpom = 0.0
    # Adsorption capacity of contaminant i on POM in bed sediment                        (µg/g)
    qcpom2 = 0.0
    # Freundlich adsorption constant of contaminant i on organic matter in water column  (L/kg)
    Kfpom = 0.0
    # Freundlich adsorption constant of contaminant i on organic matter in bed sediment  (L/kg)
    Kfpom2 = 0.0
    # Freundlich exponent of contaminant i on organic matter in water column             (unitless)
    bpom = 0.0
    # Freundlich exponent of contaminant i on organic matter in bed sediment             (unitless)
    bpom2 = 0.0
    # Algae partition coefficient of contaminant i in water column                       (L/kg)
    Kap = 0.0
    # Correction coefficient to compute Kap of contaminant i from Kow in water column    (unitless)
    aap = 0.0
    # Langmuir adsorption constant of contaminant i on algae in water column             (L/mg)
    Klap = 0.0
    # Adsoption capacity of contaminant i on algae                                       (µg/g)
    qcap = 0.0
    # Freundlich adsorption constant of contaminant i on algae                           (L/kg)
    Kfap = 0.0
    # Freundlich exponent of contaminant i on algae                                      (unitless)
    bap = 0.0
    # Solid partition coefficient of contaminant i in water column                       (L/kg)
    Kp = 0.0
    # Solid partition coefficient of contaminant i in bed sediment                       (L/kg)
    Kp2 = 0.0
    # Correction coefficient to compute Kp of contaminant i from Kow in water column     (unitless)
    ap = 0.0
    # Correction coefficient to compute Kp2 of contaminant i from Kow in bed sediment    (unitless)
    ap2 = 0.0
    # Langmuir adsorption constant of contaminant i on solid in water column             (L/mg)
    Klp = 0.0
    # Langmuir adsorption constant of contaminant i on solid in bed sediment             (L/mg)
    Klp2 = 0.0
    # Adsorption capacity of contaminant i on solid in water column                      (µg/g)
    qcp = 0.0
    # Adsorption capacity of contaminant i on solid in bed sediment                      (µg/g)
    qcp2 = 0.0
    # Freundlich adsorption constant of contaminant i on solid in water column           (L/kg)
    Kfp = 0.0
    # Freundlich adsorption constant of contaminant i on solid in bed sediment           (L/kg)
    Kfp2 = 0.0
    # Freundlich exponent of contaminant i on solid in water column                      (unitless)
    bp = 0.0
    # Freundlich exponent of contaminant i on solid in bed sediment                      (unitless)
    bp2 = 0.0
    # Arrhenius activation energy for adsorption rate                                    (kJ/mol)
    Eaad = 0.0
    # Arrhenius activation energy for desorption rate                                    (kJ/mol)
    Eade = 0.0
    # Reference temperature for which adsorption and desorption rate is reported         (C)
    Trade = 0.0

    # Reference temperature for which ionization constant is reported        (C)
    Trion = 0.0
    # Reaction enthalpy for the formation of anionic species 1               (kJ/mol)
    dHa1 = 0.0
    # Reaction enthalpy for the formation of anionic species 2               (kJ/mol)
    dHa2 = 0.0
    # Reaction enthalpy for the formation of cationic species 1              (kJ/mol)
    dHb1 = 0.0
    # Reaction enthalpy for the formation of cationic species 2              (kJ/mol)
    dHb2 = 0.0
    # Order of decay                                                         (unitless)
    nOrder = 0.0
    # Reference temperature for which hydrolysis rate is reported            (C)
    Trhyd = 0.0
    # Arrhenius activation energy for alkaline hydrolysis                    (kJ/mol)
    Eahb = 0.0
    # Arrhenius activation energy for neutral hydrolysis                     (kJ/mol)
    Eahn = 0.0
    # Arrhenius activation energy for acid hydrolysis                        (kJ/mol)
    Eaha = 0.0
    # Light intensity at which kpht is measured                              (W/m2)
    I0pht = 0.0
    # Reference temperature for which 2nd reaction rate is reported          (C)
    Trer = 0.0
    # Arrhenius activation energy for 2nd reaction rate                      (kJ/mol)
    Eaer = 0.0

    # Contaminant i to j mass yield coefficient from n order decay           (g/g)
    y1 = 0.0
    # Contaminant i to j mass yield coefficient from alkaline hydrolysis     (g/g)
    yhb = 0.0
    # Contaminant i to j mass yield coefficient from neutral hydrolysis      (g/g)
    yhn = 0.0
    # Contaminant i to j mass yield coefficient from acid hydrolysis         (g/g)
    yha = 0.0
    # Contaminant i to j mass yield coefficient from photolysis              (g/g)
    ypht = 0.0
    # Contaminant i to j mass yield coefficient from 2nd reaction            (g/g)
    yer = 0.0

    # Integer parameters (4)
    vm_option = 0
    kd_option = 0           # 1 measured partition coefficient; 2 input organic carbon partition coefficient; 3 compute organic carbon partition coefficient
    vv_option = 0           # 1 user defined; 2 compute
    Cd_solution_option = 0  # 1 Newton;       2 Bisection

    # Pathways
    C_Air_Deposition = 0.0
    C_Decay = 0.0
    C_Hydrolysis = 0.0
    C_Photolysis = 0.0
    C_Volatilization = 0.0
    C_2ndReaction = 0.0
    C_Transform_Decay = 0.0
    C_Transform_Hydrolysis = 0.0
    C_Transform_Photolysis = 0.0
    C_Transform_Reaction = 0.0
    C_C2_Settling = 0.0
    C_C2_Resuspension = 0.0
    C_C2_Transfer = 0.0

    C2_Decay = 0.0
    C2_Hydrolysis = 0.0
    C2_2ndReaction = 0.0
    C2_Transform_Decay = 0.0
    C2_Transform_Hydrolysis = 0.0
    C2_Transform_Reaction = 0.0
    C2_Settling = 0.0
    C2_Burial = 0.0
    C2_C_Resuspension = 0.0
    C2_C_Transfer = 0.0

    # Pathway index
    C_Air_Deposition_index = 0
    C_Decay_index = 0
    C_Hydrolysis_index = 0
    C_Photolysis_index = 0
    C_Volatilization_index = 0
    C_2ndReaction_index = 0
    C_Transform_Decay_index = 0
    C_Transform_Hydrolysis_index = 0
    C_Transform_Photolysis_index = 0
    C_Transform_Reaction_index = 0
    C_C2_Settling_index = 0
    C_C2_Resuspension_index = 0
    C_C2_Transfer_index = 0

    C2_Decay_index = 0
    C2_Hydrolysis_index = 0
    C2_2ndReaction_index = 0
    C2_Transform_Decay_index = 0
    C2_Transform_Hydrolysis_index = 0
    C2_Transform_Reaction_index = 0
    C2_Settling_index = 0
    C2_Burial_index = 0
    C2_C_Resuspension_index = 0
    C2_C_Transfer_index = 0
    

    # Local variables
    fion = 0.0
    fion2 = 0.0
    Cd_Species = 0.0
    Cdoc_Species = 0.0
    Cpom_Species = 0.0
    Cap_Species = 0.0
    Cd2_Species = 0.0
    Cdoc2_Species = 0.0
    Cpom2_Species = 0.0
    Cp_Species = 0.0
    Cp2_Species = 0.0

    CHH = 0.0             # concentration of H ion
    COH = 0.0             # concentration of OH ion
    gas_constant = 8.314  # universal gas constant (J/mole/K)
    Icpht = 0.0           # corrected solar radiation for photolysis
    TwaterK = 0.0         # water temperature (Kelvin)
    TsedK = 0.0           # sediment temperature (Kelvin)

    imax = 100
    i = 0
    j = 0
    k = 0
    id = 0

    IsWaterCell = False

    # Initialize all the input parameters (119 + 4) with default values
    # ionization(9 = 4 + 4 + 1)
    if globvals.use_Ionization_One:
        Ka1 = 1.0E-7
        Ka2 = 1.0E-7
        Kb1 = 1.0E-7
        Kb2 = 1.0E-7
        dHa1 = 25.0
        dHa2 = 25.0
        dHb1 = 25.0
        dHb2 = 25.0
        Trion = 25.0

    # Partition (35 = 1 + 4 + 10 + 10 + 10)
    if (globvals.use_DOCSorbed_One or ((globvals.use_Equilibrium_One(1) or globvals.use_Equilibrium_One(2)) and (globvals.nGS > 0 or globvals.use_POMSorbed_One or globvals.use_AlgaeSorbed_One))):
        Kow = 1.0E3

    # Linear equilibrium partition for DOC
    if globvals.use_DOCSorbed_One:
        Kdoc = 1.0E3
        adoc = 0.08
        if globvals.use_BedSediment:
            Kdoc2 = 1.0E3
            adoc2 = 0.08

    # Linear equilibrium partition for algae, POM, solids
    if globvals.use_Equilibrium_One(1):
        if globvals.use_AlgaeSorbed_One:
            Kap = 1.0E3
            aap = 0.2

        if globvals.use_POMSorbed_One:
            Kpom = 1.0E3
            apom = 0.2

        if globvals.nGS > 0:
            Kp = 1.0E3
            ap = 0.02

    if globvals.use_Equilibrium_One[1]:
        if globvals.use_POMSorbed_One:
            Kpom2 = 1.0E3
            apom2 = 0.2

        if globvals.nGS > 0:
            Kp2 = 1.0E3
            ap2 = 0.02

    if globvals.use_Langmuir_One[0]:
        if globvals.use_AlgaeSorbed_One:
            Klap = 1.0E3
            qcap = 1000.0

        if globvals.use_POMSorbed_One:
            Klpom = 1.0E3
            qcpom = 1000.0

        if globvals.nGS > 0:
            Klp = 1.0E3
            qcp = 1000.0

    if globvals.use_Langmuir_One[1]:
        if globvals.use_POMSorbed_One:
            Klpom2 = 1.0E3
            qcpom2 = 1000.0

        if globvals.nGS > 0:
            Klp2 = 1.0E3
            qcp2 = 1000.0

    # For algae, POM, solids
    if globvals.use_Freundlich_One[0]:
        if globvals.use_AlgaeSorbed_One:
            Kfap = 1.0E3
            bap = 1.0

        if globvals.use_POMSorbed_One:
            Kfpom = 1.0E3
            bpom = 1.0

        if globvals.nGS > 0:
            Kfp = 1.0E3
            bp = 1.0

    if globvals.use_Freundlich_One[1]:
        if globvals.use_POMSorbed_One:
            Kfpom2 = 1.0E3
            bpom2 = 1.0

        if globvals.nGS > 0:
            Kfp2 = 1.0E3
            bp2 = 1.0

    # Non-equilibrium (18 = 10 + 5 + 3 - 5?) for POM, algae, solids
    if globvals.use_NonEquilibrium_One[0]:
        if globvals.use_AlgaeSorbed_One:
            kadap = 1.0E-4
            kdeap = 0.1
            if not globvals.use_Langmuir_One[0]:
                qcap = 1000.0

        if globvals.use_POMSorbed_One:
            kadpom = 1.0E-4
            kdepom = 0.1

            if not globvals.use_Langmuir_One[0]:
                qcpom = 1000.0

        if globvals.nGS > 0:
            kadp = 1.0E-4
            kdep = 0.1

            if not globvals.use_Langmuir_One[0]:
                qcp = 1000.0

    if globvals.use_NonEquilibrium_One[1]:
        if globvals.use_POMSorbed_One:
            kadpom2 = 1.0E-4
            kdepom2 = 0.1
            if not globvals.use_Langmuir_One[1]:
                qcpom2 = 1000.0

        if globvals.nGS > 0:
            kadp2 = 1.0E-4
            kdep2 = 0.1
            if not globvals.use_Langmuir_One[1]:
                qcp2 = 1000.0

    if globvals.use_NonEquilibrium_One[0] or globvals.use_NonEquilibrium_One[1]:
        Eaad = 75.0
        Eade = 75.0
        Trade = 25.0

    # Decay (11 = 9 + 1 + 1)
    if globvals.use_nOrderDecay_One:
        k1d.rc20 = 0.1
        k1d.theta = 1.024
        if globvals.use_DOCSorbed_One:
            k1doc.rc20 = 0.1
            k1doc.theta = 1.024

        if globvals.use_AlgaeSorbed_One:
            k1ap.rc20 = 0.1
            k1ap.theta = 1.024

        if globvals.use_POMSorbed_One:
            k1pom.rc20 = 0.1
            k1pom.theta = 1.024

        if globvals.nGS > 0:
            k1p.rc20 = 0.1
            k1p.theta = 1.024

        if globvals.use_BedSediment:
            k1d2.rc20 = 0.1
            k1d2.theta = 1.024
            if globvals.use_DOCSorbed_One:
                k1doc2.rc20 = 0.1
                k1doc2.theta = 1.024

            if globvals.use_POMSorbed_One:
                k1pom2.rc20 = 0.1
                k1pom2.theta = 1.024

            if globvals.nGS > 0:
                k1p2.rc20 = 0.1
                k1p2.theta = 1.024

        nOrder = 1.0

        if globvals.nTransformProduct > 0:
            y1 = 1.0

    # Hydrolysis (19 = 12 + 4 + 3)
    if globvals.use_Hydrolysis_One:
        khbd = 1.0
        khnd = 1.0
        khad = 1.0

        if globvals.use_DOCSorbed_One:
            khbdoc = 1.0
            khndoc = 1.0
            khadoc = 1.0

        if globvals.use_BedSediment:
            khbd2 = 1.0
            khnd2 = 1.0
            khad2 = 1.0

            if globvals.use_DOCSorbed_One:
                khbdoc2 = 1.0
                khndoc2 = 1.0
                khadoc2 = 1.0

        Trhyd = 25.0
        Eahb = 75.0
        Eahn = 75.0
        Eaha = 75.0

        if globvals.nTransformProduct > 0:
            yhb = 1.0
            yhn = 1.0
            yha = 1.0

    # Photolysis (5 = 2 + 2 + 1)
    if globvals.use_Photolysis_One:
        kphtd.rc20 = 1.0
        kphtd.theta = 1.0
        if globvals.use_DOCSorbed_One:
            kphtdoc.rc20 = 1.0
            kphtdoc.theta = 1.0

        I0pht = 100.0
        alpha = 1.3

        if globvals.nTransformProduct > 0:
            ypht = 1.0

    # Volatilization (2)
    if globvals.use_Volatilization_One:
        vv.rc20 = 1.0
        vv.theta = 1.024
        KH = 0.1

    # 2nd reaction (12 = 9 + 2 + 1)
    if globvals.use_2ndReaction_One:
        kerd = 1.0
        if globvals.use_DOCSorbed_One:
            kerdoc = 1.0
        if globvals.use_AlgaeSorbed_One:
            kerap = 1.0
        if globvals.use_POMSorbed_One:
            kerpom = 1.0
        if globvals.nGS > 0:
            kerp = 1.0
        if globvals.use_BedSediment:
            kerd2 = 1.0
            if globvals.use_DOCSorbed_One:
                kerdoc2 = 1.0
            if globvals.use_POMSorbed_One:
                kerpom2 = 1.0
            if globvals.nGS > 0:
                kerp2 = 1.0

        Trer = 25.0
        Eaer = 75.0

        if globvals.nTransformProduct > 0:
            yer = 1.0

    # 1D parameters (9)
    if globvals.use_BedSediment:
        h2 = 0.1
        Por = 0.9
        z2 = 5.0
        Db = 0.03
        beta = 33.3
        ps = 2.7

    if globvals.use_AlgaeSorbed_One:
        vsap = 1.0

    if globvals.use_POMSorbed_One:
        vsom = 1.0

    if globvals.use_Langmuir_One[0] or globvals.use_Langmuir_One[1] or globvals.use_Freundlich_One[0] or globvals.use_Freundlich_One[1]:
        res = 0.001

    # 2D contaminant parameters (4)
    Csd = 100.0

    if globvals.use_BedSediment or globvals.use_Volatilization_One:
        MW = 100.0

    if globvals.use_BedSediment:
        Dm = 0.01
        vm = 1.0

    # Integer parameters (4)
    # 2D
    if globvals.use_BedSediment:
        vm_option = 1

    kd_option = 1

    if globvals.use_Langmuir_One[0] or globvals.use_Langmuir_One[1] or globvals.use_Freundlich_One[0] or globvals.use_Freundlich_One[1]:
        Cd_solution_option = 1

    if globvals.use_Volatilization_One:
        vv_option = 1
