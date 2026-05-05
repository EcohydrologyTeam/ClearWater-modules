# ClearWater Modules v3 NSM1 Phase 0: Per-Constituent Function-to-Method Mapping

**Status:** Gap analysis output (Phase 0.1)  
**Date:** 2026-05-04  
**Scope:** Function-to-method mapping for v1 NSM1 → v3 NSM1 port, organized by the 11 Process classes per the design spec Section 3.

---

## Overview

This document catalogs the v1 NSM1 kinetics functions (from `src/clearwater_modules/nsm1/processes.py`, ~290 functions across ~3,540 lines) and maps them to:

1. **v2 NSM1 methods** (where they exist; only 4 of 16 constituents are partially implemented in v2)
2. **v3 target Process.method** (proposed mapping to the 11 NSM1 Process classes)
3. **Kinetic role** (1-sentence description of what the function computes)
4. **Observable kinetic differences** between v1 and v2 code (bugs, hard-coded zeros, missing terms)

The 16 NSM1 constituents, grouped by Process class:

| Process Class | State Variables | v1 Status | v2 Status | v3 Target |
|---|---|---|---|---|
| `FloatingAlgae` | Ap | Fully implemented | Partial (bugs) | Extend |
| `BenthicAlgae` | Ab | Fully implemented | Partial (bugs) | Extend |
| `Nitrogen` | NH4, NO3, OrgN | Fully implemented | Partial (2 of 3) | Extend + Add |
| `Phosphorus` | TIP, OrgP | Fully implemented | Not in v2 | New |
| `Carbon` | POC, DOC, DIC | Fully implemented | Not in v2 | New |
| `POM` | POM | Fully implemented | Not in v2 | New |
| `CBOD` | CBOD (multi-group) | Fully implemented | Not in v2 | New |
| `DOX` | DOX | Fully implemented | Not in v2 | New |
| `Pathogen` | PX | Fully implemented | Not in v2 | New |
| `Alkalinity` | Alk | Declared but inactive | Not in v2 | New |
| `N2` | N2, TDG | Fully implemented | Not in v2 | New |

---

## 1. FloatingAlgae Process (Ap — Algal biomass, ug-Chla/L)

Floating phytoplankton, single-compartment.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `mu_max_tc` | 701-715 | Temperature-adjusted maximum growth rate (1/d) | — | — | FloatingAlgae.growth_rate_max_tc | None; kinetic formula same |
| `krp_tc` | 382-397 | Temperature-adjusted respiration rate (1/d) | — | — | FloatingAlgae.respiration_rate_tc | None |
| `kdp_tc` | 399-413 | Temperature-adjusted death rate (1/d) | — | — | FloatingAlgae.death_rate_tc | None |
| `FL` | 415-472 | Light limitation factor (Monod, dimensionless) | — | — | FloatingAlgae.light_limitation | None; v2 has equivalent limit_light() with multiple options |
| `FN` | 474-528 | Nitrogen limitation factor (Monod, dimensionless) | — | — | FloatingAlgae.nitrogen_limitation | None; computed via preference for NH4 vs NO3 |
| `FP` | 530-562 | Phosphorus limitation factor (Monod, dimensionless) | — | — | FloatingAlgae.phosphorus_limitation | None |
| `mu` | 564-604 | Growth rate under resource limitation (1/d) | — | — | FloatingAlgae.growth_rate | Multiplicative: μ = μ_max * FN * FP * FL; v2 has options (multiplicative, limiting nutrient, harmonic mean) |
| `ApGrowth` | 606-618 | Algal growth rate (ug-Chla/L/d) | — | — | FloatingAlgae.growth | Simple product μ * Ap |
| `ApRespiration` | 621-633 | Algal respiration rate (ug-Chla/L/d) | — | — | FloatingAlgae.respiration | Simple product krp_tc * Ap |
| `ApDeath` | 636-647 | Algal death rate (ug-Chla/L/d) | — | — | FloatingAlgae.death | Simple product kdp_tc * Ap |
| `ApSettling` | 650-664 | Algal settling loss (ug-Chla/L/d) | — | — | FloatingAlgae.settling | Product (vs/depth) * Ap; settling_velocity parameter |
| `dApdt` | 666-683 | Net rate of change of Ap (ug-Chla/L/d) | — | — | FloatingAlgae.change_algae | Sum of growth - respiration - death - settling |
| `Ap` | 685-700 | Final updated Ap state variable | — | — | FloatingAlgae.run (state integrator) | None |
| `rna` | 308-320 | Respiration coefficient for NH4 release (mg-N/ug-Chla) | — | — | FloatingAlgae.stoich_n_respiration | Constant ratio, used in ammonium coupling |
| `rpa` | 322-335 | Respiration coefficient for DIP release (mg-P/ug-Chla) | — | — | FloatingAlgae.stoich_p_respiration | Constant ratio |
| `rca` | 337-349 | Respiration coefficient for DIC release (mg-C/ug-Chla) | — | — | FloatingAlgae.stoich_c_respiration | Constant ratio |
| `rda` | 351-363 | Respiration coefficient for DOC production (mg-C/ug-Chla) | — | — | FloatingAlgae.stoich_c_doc | Constant ratio |
| `PN` | v1 constants.py | N:Chl-a ratio (mg-N/ug-Chla) | — | — | FloatingAlgae.DEFAULTS['pn'] | v2 floating_algae.py:43 has ratio_chla_nitrogen=7.2 |
| `ApUptakeFr_NH4` | 1206-1250 | Fraction of N uptake from NH4 vs NO3 (dimensionless, 0-1) | — | — | FloatingAlgae.nh4_uptake_fraction | v1 uses preference factor with Monod; v2 has threshold logic in nitrogen.py:239-247 |
| `NH4_ApRespiration` | 1472-1486 | NH4 source from algal respiration (mg-N/L/d) | ammonium_respiration() [FloatingAlgae:398-401] | floating_algae.py:398-401 | FloatingAlgae.ammonium_respiration() **IMPLEMENT** | **v2 returns 0 with TODO** — should compute rna * respiration |
| `NH4_ApGrowth` | 1488-1504 | NH4 sink from algal growth (mg-N/L/d) | ammonium_growth() [FloatingAlgae:403-405] | floating_algae.py:403-405 | FloatingAlgae.ammonium_growth() **IMPLEMENT** | **v2 returns 0 with TODO** — should compute rna * PN_calculated * growth |
| `DOC_algal_mortality` | 2565-2584 | DOC source from algal death (mg-C/L/d) | — | — | FloatingAlgae.doc_from_mortality() | Computed as rda * ApDeath |
| `POC_algal_mortality` | 2484-2503 | POC source from algal death (mg-C/L/d) | — | — | FloatingAlgae.poc_from_mortality() | Computed as (1-Fw) * rca * ApDeath, where Fw is fraction to water |
| `ApDeath_OrgN` | 1347-1360 | OrgN source from algal death (mg-N/L/d) | — | — | FloatingAlgae.orgn_from_mortality() | Computed as (1-Fw) * rna * ApDeath |
| `ApDeath_OrgP` | 1897-1912 | OrgP source from algal death (mg-P/L/d) | — | — | FloatingAlgae.orgp_from_mortality() | Computed as (1-Fw) * rpa * ApDeath |

### Key observations

1. **v2 FloatingAlgae has critical stub implementations:** `ammonium_respiration()` and `ammonium_growth()` both return 0 with TODO comments. This completely silences NH4 cycling through algal respiration and growth. **Must be implemented in v3** per v1 formulas at lines 1272-1276 and 1206-1220.

2. **Light limitation:** v1 uses v1's `FL` function (Beer-Lambert with self-shading). v2 has `limit_light()` with 3 options (half-saturation, Smith, Steele). v3 should adopt v2's flexibility while ensuring v1/v3 parity with v1's default (option 1 in v2).

3. **Nitrogen preference:** v1 line 1206 `ApUptakeFr_NH4` uses preference factor logic; v2 lines 239-247 have different threshold-based logic. **Observable difference:** v1 smooth transition; v2 threshold-based. v3 should adopt v1's preference factor for parity.

4. **Mortality routing:** v1 computes OrgN, OrgP, POC, DOC sources from mortality via dedicated functions (lines 1347, 1897, 2484, 2565). v2 declares `death_to_*` stub methods. **Must be implemented in v3** to route algal mortality to downstream constituents.

5. **v2 line 122 multiplicative integrator bug:** `algae = 0 + algae * rate * dt * 86400`. This is wrong on two counts: (a) multiplicative operator, (b) extra 86400 factor. v3 integrator fixes this to additive form.

---

## 2. BenthicAlgae Process (Ab — Benthic algal biomass, g-/m²)

Benthic algae, depth-integrated.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `mub_max_tc` | 701-715 | Temperature-adjusted max benthic growth (1/d) | — | — | BenthicAlgae.growth_rate_max_tc | None |
| `krb_tc` | 717-731 | Temperature-adjusted benthic respiration (1/d) | — | — | BenthicAlgae.respiration_rate_tc | None |
| `kdb_tc` | 733-747 | Temperature-adjusted benthic death (1/d) | — | — | BenthicAlgae.death_rate_tc | None |
| `FLb` | 803-862 | Light limitation for benthic algae (includes depth integration) | — | — | BenthicAlgae.light_limitation | v1 includes benthic-specific PAR/depth conversion; v2 absent |
| `FNb` | 864-917 | Nitrogen limitation for benthic (Monod) | — | — | BenthicAlgae.nitrogen_limitation | None |
| `FPb` | 919-951 | Phosphorus limitation for benthic (Monod) | — | — | BenthicAlgae.phosphorus_limitation | None |
| `FSb` | 953-982 | Sediment substrate availability (benthic-specific) | — | — | BenthicAlgae.sediment_limitation | Depth-integrated factor; absent in v2 |
| `mub` | 984-1023 | Growth rate under limitation (1/d, benthic) | — | — | BenthicAlgae.growth_rate | Multiplicative with sediment factor |
| `AbGrowth` | 1025-1038 | Benthic growth rate (g-/m²/d) | — | — | BenthicAlgae.growth | Product mub * Ab * fb (fraction bottom area) |
| `AbRespiration` | 1040-1052 | Benthic respiration (g-/m²/d) | — | — | BenthicAlgae.respiration | Product krb_tc * Ab * fb |
| `AbDeath` | 1054-1067 | Benthic death (g-/m²/d) | — | — | BenthicAlgae.death | Product kdb_tc * Ab * fb |
| `dAbdt` | 1069-1084 | Net rate of benthic algae change (g-/m²/d) | — | — | BenthicAlgae.change_algae | Sum growth - respiration - death - burial |
| `Ab` | 1086-1102 | Final updated Ab state | — | — | BenthicAlgae.run (state integrator) | None |
| `Chlb` | 1104-1120 | Chlorophyll-a from benthic algae (optional tracer) | — | — | BenthicAlgae.chlorophyll_diagnostic | Derived via stoichiometric ratio |
| `rnb` | 748-760 | Benthic respiration coefficient for NH4 (mg-N/g) | — | — | BenthicAlgae.stoich_n_respiration | Constant ratio |
| `rpb` | 762-774 | Benthic respiration coefficient for DIP (mg-P/g) | — | — | BenthicAlgae.stoich_p_respiration | Constant ratio |
| `rcb` | 776-788 | Benthic respiration coefficient for DIC (mg-C/g) | — | — | BenthicAlgae.stoich_c_respiration | Constant ratio |
| `rab` | 790-801 | Benthic death coefficient for OrgN (mg-N/g) | — | — | BenthicAlgae.stoich_n_death | Constant ratio |
| `AbUptakeFr_NH4` | 1263-1304 | NH4 uptake fraction for benthic (dimensionless) | — | — | BenthicAlgae.nh4_uptake_fraction | v1 depth-adjusted preference logic |
| `NH4_AbRespiration` | 1506-1525 | NH4 source from benthic respiration (mg-N/L/d) | ammonium_respiration() [BenthicAlgae] | benthic_algae.py | BenthicAlgae.ammonium_respiration() **IMPLEMENT** | **v2 equivalent expected but not verified** |
| `NH4_AbGrowth` | 1527-1547 | NH4 sink from benthic growth (mg-N/L/d) | ammonium_growth() [BenthicAlgae] | benthic_algae.py | BenthicAlgae.ammonium_growth() **IMPLEMENT** | **v2 equivalent expected** |
| `AbDeath_OrgN` | 1362-1380 | OrgN source from benthic death (mg-N/L/d) | — | — | BenthicAlgae.orgn_from_mortality() | Computed via benthic-specific ratio |
| `AbDeath_OrgP` | 1914-1935 | OrgP source from benthic death (mg-P/L/d) | — | — | BenthicAlgae.orgp_from_mortality() | Computed via benthic-specific ratio |
| `DOC_benthic_algae_mortality` | 2586-2612 | DOC from benthic death (mg-C/L/d) | — | — | BenthicAlgae.doc_from_mortality() | Includes depth integration (to L, not L/m²) |
| `POC_benthic_algae_mortality` | 2505-2530 | POC from benthic death (mg-C/L/d) | — | — | BenthicAlgae.poc_from_mortality() | Includes depth integration and Fb/Fw fractionation |

### Key observations

1. **Benthic depth integration:** v1 balances benthic state (g-/m²) with water-column updates (mg-N/L/d etc.) via `fb` (fraction bottom area) and depth normalization. v2 absent. **Critical for v3:** must thread depth through benthic methods to produce mg-/L/d rates consumable by water-column processes.

2. **Sediment limitation (FSb):** v1 includes sediment substrate availability (line 953). v2 absent. v3 should check if v1 really uses this or if it's historical artifact; if used, must implement.

3. **v2 benthic_algae.py is largely a stub.** Similar bug situation as FloatingAlgae: ammonium_respiration/growth are absent. Must be implemented in v3.

4. **Mortality fractionation:** v1 computes Fb/Fw-adjusted contributions (e.g., POC_benthic_algae_mortality at line 2505 uses `Fb * (1-Fw) * ...`). v2 absent. v3 must implement.

---

## 3. Nitrogen Process (NH4, NO3, OrgN)

Ammonium, nitrate, organic nitrogen. Extends existing v2 Nitrogen process.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `knit_tc` | 1122-1137 | Temperature-adjusted nitrification rate (1/d) | — | — | Nitrogen.nitrification_rate_tc | None |
| `rnh4_tc` | 1139-1154 | Temperature-adjusted OrgN to NH4 hydrolysis (1/d) | — | — | Nitrogen.orgn_hydrolysis_rate_tc | None |
| `vno3_tc` | 1156-1171 | Temperature-adjusted NO3 settling velocity (m/d) | — | — | Nitrogen.no3_settling_velocity_tc | None; v1 probably unused (NO3 doesn't settle) |
| `kon_tc` | 1173-1188 | Temperature-adjusted OrgN decay to NH4 (1/d) | — | — | Nitrogen.orgn_decay_rate_tc | None; duplicate of rnh4_tc conceptually |
| `kdnit_tc` | 1190-1204 | Temperature-adjusted denitrification rate (1/d) | — | — | Nitrogen.denitrification_rate_tc | None |
| `NitrificationInhibition` | 1422-1437 | Oxygen inhibition of nitrification (dimensionless 0-1) | — | — | Nitrogen.nitrification_inhibition | v1 uses exp(-KsOxdn * DOX); v2 line 365 similar |
| `NH4_Nitrification` | 1439-1455 | Nitrification rate (mg-N/L/d) | ammonium_nitrification() [Nitrogen:260-276] | nitrogen.py:260-276 | Nitrogen.nitrification | Rates computed identically; v2 has parameter wiring issues (line 191 hard-coded 1) |
| `NH4fromBed` | 1457-1470 | NH4 release from sediment (mg-N/L/d) | ammonium_from_bed() [Nitrogen:221-227] | nitrogen.py:221-227 | Nitrogen.nh4_from_sediment | None |
| `OrgN_NH4_Decay` | 1317-1331 | OrgN hydrolysis to NH4 (mg-N/L/d) | — | — | Nitrogen.orgn_hydrolysis | None |
| `OrgN_Settling` | 1333-1345 | OrgN settling loss (mg-N/L/d) | — | — | Nitrogen.orgn_settling | Product (vson_tc/depth) * OrgN |
| `ApDeath_OrgN` | 1347-1360 | OrgN contribution from Ap death | — | — | FloatingAlgae.orgn_from_mortality() | Routes from FloatingAlgae process |
| `AbDeath_OrgN` | 1362-1380 | OrgN contribution from Ab death | — | — | BenthicAlgae.orgn_from_mortality() | Routes from BenthicAlgae process |
| `dOrgNdt` | 1383-1403 | Net OrgN rate of change | — | — | Nitrogen.change_organic_nitrogen | New method; sum hydrolysis - settling + sources |
| `OrgN` | 1405-1420 | Updated OrgN state | — | — | Nitrogen.run (state integrator for OrgN) | New state; v1 computed but never returned |
| `NO3_Denit` | 1604-1638 | NO3 denitrification loss (mg-N/L/d) | nitrate_denitrification() [Nitrogen:291-313] | nitrogen.py:291-313 | Nitrogen.denitrification | **v2 lines 191 hard-coded 1 for half_saturation_oxygen — must wire parameter** |
| `NO3_BedDenit` | 1640-1655 | NO3 sediment denitrification (mg-N/L/d) | nitrate_bed_denitrification() [Nitrogen:315-328] | nitrogen.py:315-328 | Nitrogen.no3_from_sediment_denitrification | None |
| `NO3_ApGrowth` | 1657-1675 | NO3 uptake by Ap (mg-N/L/d) | nitrate_uptake_floating_algae() [Nitrogen:330-341] | nitrogen.py:330-341 | Nitrogen.no3_uptake_floating_algae | **v2 lines 204, 211 hard-coded 0 for algea_growth_rate — must read from registry** |
| `NO3_AbGrowth` | 1677-1698 | NO3 uptake by Ab (mg-N/L/d) | nitrate_uptake_benthic_algae() [Nitrogen:343-359] | nitrogen.py:343-359 | Nitrogen.no3_uptake_benthic_algae | **v2 lines 212, 213 hard-coded 0 — must read from registry** |
| `dNH4dt` | 1549-1585 | Net NH4 rate of change | change_ammonium() [Nitrogen:115-149] | nitrogen.py:115-149 | Nitrogen.change_ammonium (corrected) | **v2 line 101 multiplicative integrator: `ammonium = 0 + ammonium * rate * dt`** |
| `dNO3dt` | 1700-1730 | Net NO3 rate of change | change_nitrate() [Nitrogen:170-219] | nitrogen.py:170-219 | Nitrogen.change_nitrate (corrected) | **v2 line 112 multiplicative integrator: `nitrate = 0 + nitrate * rate * dt`** |
| `NH4` | 1587-1602 | Updated NH4 state | — | — | Nitrogen.run (state integrator for NH4) | None; v2 drops state update (no set_at_time) |
| `NO3` | 1732-1747 | Updated NO3 state | — | — | Nitrogen.run (state integrator for NO3) | None |

### Key observations

1. **Two critical v2 bugs in Nitrogen:**
   - **Lines 101, 112:** Multiplicative integrator `state = 0 + state * rate * dt` instead of additive `state = state + rate * dt`. This changes the kinetics from linear accumulation to exponential decay/growth. **MUST FIX in v3.**
   - **Lines 191, 204-212:** Hard-coded zeros for `half_saturation_oxygen=1` and `algea_growth_rate=0`. These silence oxygen inhibition on nitrification and algal NO3 uptake completely. **MUST WIRE PARAMETERS in v3.**

2. **NaN guards in v2:** Lines 147, 218, 250, 313 use `rate == np.nan` (always False per IEEE 754). v3 must replace with `.isnull()` check.

3. **OrgN is new in v3:** v1 computes `OrgN` (line 1405) but v2 Nitrogen doesn't have it. v3 adds as third state variable.

4. **State persistence in v2:** The integrator updates state variables (lines 101, 112) but never calls `set_at_time`, so changes are dropped. v3's pattern calls `set_at_time` after clip_negative check.

---

## 4. Phosphorus Process (TIP, OrgP)

Total inorganic phosphorus and organic phosphorus. New in v3.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `kop_tc` | 1833-1848 | Temperature-adjusted OrgP to TIP hydrolysis (1/d) | — | — | Phosphorus.orgp_hydrolysis_rate_tc | None |
| `rpo4_tc` | 1850-1864 | Temperature-adjusted OrgP settling velocity (m/d) | — | — | Phosphorus.orgp_settling_velocity_tc | None |
| `OrgP_DIP_decay` | 1866-1880 | OrgP hydrolysis to TIP (mg-P/L/d) | — | — | Phosphorus.orgp_hydrolysis | Product kon_tc * OrgP |
| `OrgP_Settling` | 1882-1895 | OrgP settling loss (mg-P/L/d) | — | — | Phosphorus.orgp_settling | Product (vsop_tc/depth) * OrgP; **v1 sentinel-999 default for vsop** |
| `ApDeath_OrgP` | 1897-1912 | OrgP source from Ap death (mg-P/L/d) | — | — | FloatingAlgae.orgp_from_mortality() | Routes from FloatingAlgae |
| `AbDeath_OrgP` | 1914-1935 | OrgP source from Ab death (mg-P/L/d) | — | — | BenthicAlgae.orgp_from_mortality() | Routes from BenthicAlgae |
| `dOrgPdt` | 1937-1957 | Net OrgP rate of change (mg-P/L/d) | — | — | Phosphorus.change_organic_phosphorus | Sum hydrolysis - settling + sources |
| `DIPfromBed` | 1959-1971 | DIP (TIP) release from sediment (mg-P/L/d) | — | — | Phosphorus.dip_from_sediment | Arrhenius-corrected parameterized flux |
| `fdp` | 290-306 | Fraction of TIP in dissolved phase (dimensionless 0-1) | — | — | Phosphorus.tip_dissolved_fraction (via shared utility) | v1 Solid-dependent; v2 hard-codes 0.5 at floating_algae.py:113 **TODO** |
| `TIP_Settling` | 1973-1988 | TIP (particulate phase) settling (mg-P/L/d) | — | — | Phosphorus.tip_settling | Product (vs/depth) * (1-fdp) * TIP; **v1 sentinel-999 default for vs** |
| `DIP_ApRespiration` | 1990-2003 | DIP source from Ap respiration (mg-P/L/d) | — | — | FloatingAlgae.dip_from_respiration() | Routes rpa * respiration |
| `DIP_ApGrowth` | 2005-2018 | DIP sink from Ap growth (mg-P/L/d) | — | — | FloatingAlgae.dip_uptake() | Routes via fdp fraction |
| `DIP_AbRespiration` | 2020-2037 | DIP source from Ab respiration (mg-P/L/d) | — | — | BenthicAlgae.dip_from_respiration() | Routes with depth integration |
| `DIP_AbGrowth` | 2039-2056 | DIP sink from Ab growth (mg-P/L/d) | — | — | BenthicAlgae.dip_uptake() | Routes with depth integration |
| `dTIPdt` | 2058-2093 | Net TIP rate of change (mg-P/L/d) | — | — | Phosphorus.change_tip | Sum settling - sources - sinks + sediment release |
| `TIP` | 2095-2109 | Updated TIP state | — | — | Phosphorus.run (state integrator for TIP) | None |
| `OrgP` | 2111-2124 | Updated OrgP state | — | — | Phosphorus.run (state integrator for OrgP) | None |

### Key observations

1. **fdp partitioning:** v1 line 290 computes `fdp` as Solid-dependent function. v2 hard-codes `phosphate_fraction_dissolved=0.5` at floating_algae.py:113 with TODO. v3 must implement `fdp` utility and wire to both Phosphorus and FloatingAlgae processes.

2. **Settling velocity sentinel defaults:** v1 constants.py sets `vs=999` and `vsop=999` (TIP and OrgP settling). v3 corrects to 0.1 m/d per spec Section 7.

3. **Algal uptake routing:** DIP_ApGrowth and DIP_AbGrowth (lines 2005, 2039) are separate functions. v3 routes through FloatingAlgae.dip_uptake() and BenthicAlgae.dip_uptake() methods reading from registry.

4. **Benthic depth integration:** OrgP and TIP contributions from benthic algae must include depth normalization (g-/m² → mg-P/L/d).

---

## 5. Carbon Process (POC, DOC, DIC)

Particulate, dissolved, and inorganic carbon pools with mineralization and reaeration.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `kpoc_tc` | 2439-2453 | Temperature-adjusted POC hydrolysis (1/d) | — | — | Carbon.poc_hydrolysis_rate_tc | None |
| `POC_hydrolysis` | 2455-2467 | POC to DOC conversion (mg-C/L/d) | — | — | Carbon.poc_hydrolysis | Product kpoc_tc * POC |
| `POC_settling` | 2469-2482 | POC settling loss (mg-C/L/d) | — | — | Carbon.poc_settling | Product (vspoc/depth) * POC |
| `POC_algal_mortality` | 2484-2503 | POC source from Ap death (mg-C/L/d) | — | — | FloatingAlgae.poc_from_mortality() | Computed as (1-Fw) * rca * ApDeath |
| `POC_benthic_algae_mortality` | 2505-2530 | POC source from Ab death (mg-C/L/d) | — | — | BenthicAlgae.poc_from_mortality() | Computed as Fb * (1-Fw) * rab * AbDeath with depth normalization |
| `dPOCdt` | 2532-2548 | Net POC rate of change (mg-C/L/d) | — | — | Carbon.change_poc | Sum hydrolysis - settling + sources |
| `POC` | 2550-2563 | Updated POC state | — | — | Carbon.run (state integrator for POC) | None |
| `kdoc_tc` | 2614-2627 | Temperature-adjusted DOC oxidation (1/d) | — | — | Carbon.doc_oxidation_rate_tc | None |
| `DOC_DIC_oxidation` | 2629-2649 | DOC to DIC oxidation (mg-C/L/d) | — | — | Carbon.doc_oxidation | Product kdoc_tc * Monod(DOX) * DOC; **monitored kinetic difference: v1 line 2639 uses `DOX/(DOX+KdocDOX)` Monod** |
| `DOC_algal_mortality` | 2565-2584 | DOC source from Ap death (mg-C/L/d) | — | — | FloatingAlgae.doc_from_mortality() | Computed as Fw * rda * ApDeath (fraction to water) |
| `DOC_benthic_algae_mortality` | 2586-2612 | DOC source from Ab death (mg-C/L/d) | — | — | BenthicAlgae.doc_from_mortality() | Computed as Fb * Fw * rda * AbDeath with depth normalization |
| `dDOCdt` | 2651-2669 | Net DOC rate of change (mg-C/L/d) | — | — | Carbon.change_doc | Sum POC hydrolysis - oxidation + sources |
| `DOC` | 2671-2685 | Updated DOC state | — | — | Carbon.run (state integrator for DOC) | None |
| `Henrys_k` | 2687-2696 | Henry's law constant for CO2 at temp (K atm / mol) | — | — | Carbon.henrys_k (or shared utility) | Constant per Weiss 1974 |
| `Atmospheric_CO2_reaeration` | 2698-2715 | CO2 reaeration from atmosphere (mg-C/L/d) | — | — | Carbon.dic_atmospheric_reaeration | Uses Henry's law and `ka_tc` reaeration rate |
| `DIC_algal_respiration` | 2717-2732 | DIC source from Ap respiration (mg-C/L/d) | — | — | FloatingAlgae.dic_from_respiration() | Routes rca * respiration |
| `DIC_algal_photosynthesis` | 2734-2749 | DIC sink from Ap photosynthesis (mg-C/L/d) | — | — | FloatingAlgae.dic_from_photosynthesis() | Routes rca * growth |
| `DIC_benthic_algae_respiration` | 2751-2770 | DIC source from Ab respiration (mg-C/L/d) | — | — | BenthicAlgae.dic_from_respiration() | Routes with depth integration |
| `DIC_benthic_algae_photosynthesis` | 2772-2791 | DIC sink from Ab photosynthesis (mg-C/L/d) | — | — | BenthicAlgae.dic_from_photosynthesis() | Routes with depth integration |
| `DIC_CBOD_oxidation` | 2793-2815 | DIC source from CBOD oxidation (mg-C/L/d) | — | — | Carbon.dic_from_cbod_oxidation (reads CBOD rate) | Routes stoichiometric C equivalent |
| `DIC_sed_release` | 2817-2832 | DIC source from sediment (mg-C/L/d) | — | — | Carbon.dic_from_sediment | Parameterized sediment release |
| `dDICdt` | 2834-2856 | Net DIC rate of change (mg-C/L/d) | — | — | Carbon.change_dic | Sum respiration + oxidation + sediment - photosynthesis - reaeration |
| `DIC` | 2858-2876 | Updated DIC state | — | — | Carbon.run (state integrator for DIC) | None |

### Key observations

1. **DOC oxidation kinetics:** v1 line 2639 uses Monod kinetics with DOX as limiting factor. **Kinetic difference:** must verify v1 formula `DOX/(DOX+KdocDOX)` is correct and implement in v3.

2. **Benthic depth integration:** POC, DOC, and DIC contributions from benthic algae include depth normalization (g-/m² → mg-C/L/d).

3. **Algal mortality fractionation:** v1 uses `Fw` (fraction to water) to split algal death between POC and DOC. v3 must route these via FloatingAlgae/BenthicAlgae methods.

4. **DIC reaeration:** v1 line 2698 uses Henry's law with atmosphere pCO2 and `ka_tc` reaeration. v3 must wire `ka_tc` from shared reaeration utility.

5. **CBOD coupling:** v1 line 2793 reads CBOD oxidation and produces DIC source. v3 reads `cbod_oxidation_rate` from registry.

---

## 6. POM Process (POM — Particulate Organic Matter, mg/L)

Refractory organic matter, settling and dissolution.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `kpom_tc` | 2185-2198 | Temperature-adjusted POM dissolution (1/d) | — | — | POM.dissolution_rate_tc | None |
| `POM_algal_settling` | 2200-2220 | POM source from Ap settling (mg/L/d) | — | — | FloatingAlgae.pom_from_settling() | Routes settling loss |
| `POM_dissolution` | 2222-2234 | POM dissolution to POC (mg/L/d) | — | — | POM.dissolution | Product kpom_tc * POM |
| `POM_POC_settling` | 2236-2255 | POM sink from POC settling (mg/L/d) | — | — | Carbon.poc_settling (reads into POM via registry) | Aliases POC settling |
| `POM_benthic_algae_mortality` | 2257-2279 | POM source from Ab death (mg/L/d) | — | — | BenthicAlgae.pom_from_mortality() | Routes with depth integration |
| `POM_burial` | 2281-2295 | POM burial loss (mg/L/d) | — | — | POM.burial | Product (vspom/depth) * POM; settling-like term |
| `dPOMdt` | 2297-2315 | Net POM rate of change (mg/L/d) | — | — | POM.change_pom | Sum dissolution - burial + sources |
| `POM` | 2317-2332 | Updated POM state | — | — | POM.run (state integrator for POM) | None |

### Key observations

1. **POM is largely independent:** depends only on dissolving to POC and receiving algal settling/death contributions. No inter-process couplings beyond algal mortality routing.

2. **Burial vs settling:** v1 line 2281 calls burial a settling term (`vspom/depth`), but conceptually it's sediment burial. v3 naming should clarify (likely just settling with different velocity parameter).

3. **POC settling coupling:** v1 line 2236 routes POC settling as both POC loss and POM source (same rate). v3 computes POC settling in Carbon.poc_settling() and writes to registry; POM.run() reads and applies.

---

## 7. CBOD Process (CBOD — Carbonaceous Biochemical Oxygen Demand, multi-group)

Labile organic matter groups (v1 uses up to 3 groups).

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `kbod_tc` | 2334-2350 | Temperature-adjusted CBOD hydrolysis (1/d) | — | — | CBOD.hydrolysis_rate_tc (per group) | None |
| `ksbod_tc` | 2352-2368 | Temperature-adjusted CBOD settling (m/d) | — | — | CBOD.settling_velocity_tc (per group) | None |
| `CBOD_oxidation` | 2370-2390 | CBOD aerobic decomposition (mg-C/L/d) | — | — | CBOD.oxidation (per group) | Product kbod_tc * Monod(DOX) * CBOD; **Monod kinetics with DOX dependency** |
| `CBOD_sedimentation` | 2392-2406 | CBOD settling loss (mg-C/L/d) | — | — | CBOD.settling (per group) | Product (ksbod_tc/depth) * CBOD |
| `dCBODdt` | 2408-2420 | Net CBOD rate of change (mg-C/L/d) | — | — | CBOD.change_cbod (per group) | Sum oxidation - settling |
| `CBOD` | 2422-2437 | Updated CBOD state (per group) | — | — | CBOD.run (state integrator for each group) | None |

### Key observations

1. **Multi-group structure:** v1 likely uses 1-3 groups; v3 should preserve loop-over-groups pattern. Each group has independent parameters and state.

2. **CBOD oxidation kinetics:** Uses Monod with DOX (similar to DOC oxidation). v1 line 2370 formula needs verification.

3. **Stoichiometric coupling:** CBOD oxidation produces DIC and consumes DOX. v3 writes `cbod_oxidation_rate` to registry (sum over groups) for DOX and Carbon processes.

4. **No v2 equivalent:** CBOD is absent in v2, so no parity issues.

---

## 8. DOX Process (DOX — Dissolved Oxygen, mg-O2/L)

The most highly coupled constituent, with photosynthesis/respiration, nitrification, organicmatter oxidation, reaeration, and SOD.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `pwv` | 2878-2888 | Partial pressure of water vapor (atm) at temp | — | — | DOX.water_vapor_pressure (or shared utility) | Temperature-dependent formula (Magnus) |
| `DOs_atm_alpha` | 2890-2899 | Solubility coefficient for O2 in water (mL/L/atm) | — | — | DOX.o2_solubility_coefficient (or shared utility) | APHA/QUAL2E formula with salinity/pressure correction |
| `DOX_sat` | 2901-2925 | Dissolved oxygen saturation (mg-O2/L) at temp/pressure/salinity | — | — | DOX.saturation (or shared utility) | Complex formula combining solubility, vapor pressure, barometric pressure |
| `Atm_O2_reaeration` | 2927-2940 | Atmospheric O2 reaeration source (mg-O2/L/d) | — | — | DOX.atmospheric_reaeration | Uses `ka_tc` reaeration rate and saturation deficit |
| `DOX_ApGrowth` | 2942-2960 | O2 production from Ap photosynthesis (mg-O2/L/d) | — | — | DOX.o2_from_photosynthesis_floating | Routes algal growth * stoich coeff; **NH4 vs NO3 fractionation** |
| `DOX_ApRespiration` | 2962-2978 | O2 consumption from Ap respiration (mg-O2/L/d) | — | — | DOX.o2_from_respiration_floating | Routes algal respiration * stoich coeff |
| `DOX_Nitrification` | 2980-3000 | O2 consumption from nitrification (mg-O2/L/d) | — | — | DOX.o2_from_nitrification (reads nitrification_rate) | Uses stoichiometric coefficient (typically 4.57 g-O2/g-N) |
| `DOX_DOC_oxidation` | 3002-3017 | O2 consumption from DOC oxidation (mg-O2/L/d) | — | — | DOX.o2_from_doc_oxidation (reads doc_oxidation_rate) | Uses stoichiometric coefficient |
| `DOX_CBOD_oxidation` | 3019-3030 | O2 consumption from CBOD oxidation (mg-O2/L/d) | — | — | DOX.o2_from_cbod_oxidation (reads cbod_oxidation_rate) | Uses stoichiometric coefficient; sums over CBOD groups |
| `DOX_AbGrowth` | 3032-3055 | O2 production from Ab photosynthesis (mg-O2/L/d) | — | — | DOX.o2_from_photosynthesis_benthic | Routes with depth integration; **NH4 vs NO3 fractionation** |
| `DOX_AbRespiration` | 3057-3079 | O2 consumption from Ab respiration (mg-O2/L/d) | — | — | DOX.o2_from_respiration_benthic | Routes with depth integration |
| `DOX_SOD` | 3081-3093 | Sediment oxygen demand (mg-O2/L/d) | — | — | DOX.sod | Uses `SOD_tc` shared utility; **v1 sentinel-999 default for SOD_theta causes blowup** |
| `dDOXdt` | 3095-3121 | Net DOX rate of change (mg-O2/L/d) | — | — | DOX.change_dox | Sum all sources and sinks |
| `DOX` | 3123-3139 | Updated DOX state | — | — | DOX.run (state integrator for DOX) | None |

### Key observations

1. **Photosynthesis/respiration coupling:** v1 lines 2942, 3032 handle algal O2 production with NH4/NO3 fractionation. When algae prefer NH4, O2 production is higher (redox stoichiometry). v3 must read `algal_nh4_uptake_fraction` from registry to correctly adjust stoichiometry.

2. **DO saturation calculation:** v1 lines 2901-2925 use complex formula with barometric pressure and salinity. v3 must port entire function to DOX process or shared utility.

3. **SOD sentinel-999 bug:** v1 constants.py sets `SOD_theta=999`, producing catastrophic exponential blowup at T > 20 °C (lines 3081-3093). **v3 corrects to 1.060 per spec Section 7.** Regression test required.

4. **Monod kinetics in oxidation sinks:** DOC and CBOD oxidation (lines 3002, 3019) use Monod with DOX. v3 must verify these formulas match v1.

5. **Multi-sink clipping risk:** DOX is subject to simultaneous sinks (SOD, nitrification, DOC oxidation, CBOD oxidation, respiration). Under large dt or stiff conditions, DOX can clip to zero. v3 may need semi-implicit treatment (Phase 5 per spec Section 3).

---

## 9. Pathogen Process (PX — Pathogenic organisms, counts/100mL or CFU/mL)

Natural decay, light-induced inactivation, settling.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `kdx_tc` | 3141-3156 | Temperature-adjusted pathogen decay (1/d) | — | — | Pathogen.decay_rate_tc | None |
| `PathogenDeath` | 3158-3170 | Natural pathogen mortality (counts/100mL/d) | — | — | Pathogen.natural_decay | Product kdx_tc * PX |
| `PathogenDecay` | 3172-3191 | Light-induced pathogen decay (counts/100mL/d) | — | — | Pathogen.light_decay | Exponential in light attenuation; uses `KEXT` from light utility |
| `PathogenSettling` | 3193-3207 | Pathogen settling loss (counts/100mL/d) | — | — | Pathogen.settling | Product (vs_px/depth) * PX |
| `dPXdt` | 3209-3225 | Net pathogen rate of change (counts/100mL/d) | — | — | Pathogen.change_px | Sum natural_decay + light_decay - settling |
| `PX` | 3227-3244 | Updated PX state | — | — | Pathogen.run (state integrator for PX) | None |

### Key observations

1. **Light-induced decay:** v1 line 3172 uses exponential dependence on light penetration. v3 must read `KEXT` (light extinction coefficient) from shared utility or registry.

2. **Simple process:** Pathogen is one of the simplest constituents with no coupling to other processes except reading light extinction.

3. **No v2 equivalent:** Pathogen is absent in v2.

---

## 10. Alkalinity Process (Alk — Alkalinity, mg-CaCO3/L or meq/L)

Tracer for acid-base chemistry, with nitrification/denitrification and algal growth/respiration coupling.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `Alk_denitrification` | 3246-3282 | Alkalinity production from denitrification (meq/L/d) | — | — | Alkalinity.alk_from_denitrification (reads denitrification_rate) | Uses stoichiometric coefficient (2 meq/mol NO3 reduced) |
| `Alk_nitrification` | 3284-3320 | Alkalinity consumption by nitrification (meq/L/d) | — | — | Alkalinity.alk_from_nitrification (reads nitrification_rate) | Uses stoichiometric coefficient (2 meq/mol N nitrified) |
| `Alk_algal_growth` | 3322-3343 | Alkalinity production from Ap growth (meq/L/d) | — | — | Alkalinity.alk_from_growth_floating | Routes via algal growth rate and N uptake fraction |
| `Alk_algal_respiration` | 3345-3362 | Alkalinity consumption from Ap respiration (meq/L/d) | — | — | Alkalinity.alk_from_respiration_floating | Routes via respiration rate |
| `Alk_benthic_algae_growth` | 3364-3389 | Alkalinity production from Ab growth (meq/L/d) | — | — | Alkalinity.alk_from_growth_benthic | Routes with depth integration |
| `Alk_benthic_algae_respiration` | 3391-3411 | Alkalinity consumption from Ab respiration (meq/L/d) | — | — | Alkalinity.alk_from_respiration_benthic | Routes with depth integration |
| `dAlkdt` | 3413-3433 | Net alkalinity rate of change (meq/L/d) | — | — | Alkalinity.change_alk | Sum all sources and sinks |
| `Alk` | 3435-3450 | Updated Alk state | — | — | Alkalinity.run (state integrator for Alk) | None |

### Key observations

1. **Alkalinity is declared but never updated in v1:** v1 computes `dAlkdt` and returns `Alk` state, but the model.py never applies these. v3 makes Alkalinity an active Process.

2. **Stoichiometric couplings:** Alkalinity reads `nitrification_rate` and `denitrification_rate` from registry, and `algal_growth_rate`, `algal_respiration_rate` from algal processes.

3. **pH solver deferred:** v3 1.0.0 Alkalinity is a tracer (integrates source/sink terms only). Full carbonate-pH solver comes in v3 1.1+ (NSM2 features).

4. **No v2 equivalent:** Alkalinity is absent in v2.

---

## 11. N2 / TDG Process (N2 — Dissolved nitrogen gas, mg-N2/L; TDG — Total dissolved gas as derived variable)

Atmospheric exchange via Henry's law, with denitrification source.

### v1 Functions

| v1 function | v1 lines | Kinetic role | v2 method | v2 lines | v3 target | Kinetic difference |
|---|---|---|---|---|---|---|
| `KHN2_tc` | 3452-3468 | Temperature-adjusted Henry's law constant for N2 (mol / (L·atm)) | — | — | N2.henrys_k_tc (or shared utility) | Temperature-dependent per Weiss 1970 |
| `N2sat` | 3470-3488 | N2 saturation concentration (mg-N2/L) | — | — | N2.saturation (or shared utility) | Uses Henry's law with atmospheric N2 partial pressure |
| `dN2dt` | 3490-3505 | Net N2 rate of change (mg-N2/L/d) | — | — | N2.change_n2 | Sum atmospheric reaeration - denitrification source |
| `N2` | 3507-3521 | Updated N2 state | — | — | N2.run (state integrator for N2) | None |
| `TDG` | 3523-3540 | Total dissolved gas as diagnostic (fraction saturation) | — | — | N2.tdg_diagnostic (derived variable) | TDG = N2/N2sat (or weighted with O2 if both present) |

### Key observations

1. **Atmospheric N2 exchange:** v1 line 3490 models N2 reaeration similar to DO, using `ka_tc * (N2sat - N2)`.

2. **Denitrification coupling:** N2 receives N2 source from denitrification. v3 reads `denitrification_rate` from Nitrogen process.

3. **TDG diagnostic:** v1 line 3523 computes TDG as derived variable (not state), useful for diagnostic of supersaturation.

4. **No v2 equivalent:** N2 is absent in v2.

---

## Shared Utilities and Physics Primitives

These functions appear in v1's `processes.py` but are not state-specific and should migrate to v3 utils:

| Utility | v1 lines | v3 destination | Status |
|---|---|---|---|
| `celsius_to_kelvin` | 9-10 | v2 utils.conversions (already exists) | Re-export from v2 |
| `arrhenius_correction` | 12-35 | v2 utils.conversions (already exists) | Re-export from v2 |
| `TwaterK` | 40-47 | v3 utils.temperature (wrapper, probably unused) | Deprecate or keep for compat |
| `kah_20` | 50-106 | **v3 utils.reaeration.py NEW** | 9 hydraulic options |
| `kah_tc` | 110-122 | **v3 utils.reaeration.py NEW** | Temperature-corrected hydraulic reaeration |
| `kaw_20` | 125-182 | **v3 utils.reaeration.py NEW** | 13 wind options |
| `kaw_tc` | 186-198 | **v3 utils.reaeration.py NEW** | Temperature-corrected wind reaeration |
| `ka_tc` | 202-214 | **v3 utils.reaeration.py NEW** | Combined reaeration rate |
| `SOD_tc` | 216-236 | **v3 utils.sediment.py NEW** | Temperature and DOX-corrected sediment oxygen demand |
| `L` | 238-272 | **v3 utils.light.py NEW** | Light extinction coefficient (Beer-Lambert) |
| `PAR` | 274-287 | **v3 utils.light.py NEW** | Photosynthetically active radiation |
| `fdp` | 290-306 | **v3 utils.partitioning.py NEW** | Fraction TIP dissolved phase |

---

## Cross-Constituent Dependencies and Function Call Graph

Functions that reference each other within processes.py (dependency order important for sequential reasoning):

### Floating Algae Dependency Chain
- `FL` → uses `L` (light extinction), `Ap` (state)
- `FN` → uses `ApUptakeFr_NH4` (preference)
- `FP` → independent
- `mu` → uses `mu_max_tc`, `FN`, `FP`, `FL`
- `ApGrowth` → uses `mu`
- `ApRespiration` → uses `krp_tc`
- `ApDeath` → uses `kdp_tc`
- `ApSettling` → independent of other functions
- `dApdt` → uses all four above
- `Ap` → integrates `dApdt`

### Nitrogen (Phased Aspect: Algae Coupling)
- `dNH4dt` → uses `NH4_Nitrification`, `NH4fromBed`, `NH4_ApRespiration`, `NH4_ApGrowth`, `NH4_AbRespiration`, `NH4_AbGrowth`
- `dNO3dt` → uses `NO3_Denit`, `NO3_BedDenit`, `NO3_ApGrowth`, `NO3_AbGrowth`
- `NO3_ApGrowth` → reads from Algae via preference factor
- `ApUptakeFr_NH4` → used by `NH4_ApGrowth` (fractionation logic)

### Carbon-DOX Coupling
- `dDOCdt` → uses `POC_hydrolysis`, `DOC_DIC_oxidation`, algal mortality sources
- `dDICdt` → uses `DOC_DIC_oxidation` (oxidation source), `DIC_algal_photosynthesis`, `DIC_algal_respiration`, `Atmospheric_CO2_reaeration`, `DIC_CBOD_oxidation`, `DIC_sed_release`
- `DOX_DOC_oxidation` → reads from `DOC_DIC_oxidation` (same rate, stoich difference)
- `dDOXdt` → uses `DOX_ApGrowth`, `DOX_ApRespiration`, `DOX_Nitrification`, `DOX_DOC_oxidation`, `DOX_CBOD_oxidation`, `DOX_AbGrowth`, `DOX_AbRespiration`, `DOX_SOD`, `Atm_O2_reaeration`

### Implications for v3 Dispatch Order
The spec defines dispatch order (Section 3.3):
1. FloatingAlgae (writes rates: growth, respiration, death, N uptake fraction)
2. BenthicAlgae (writes rates, benthic-specific)
3. Nitrogen (reads algal rates, writes nitrification, denitrification rates)
4. Phosphorus (reads algal P uptake)
5. Carbon (reads algal respiration/growth for DIC; reads nitrification for rates)
6. POM (reads algal settling from FloatingAlgae)
7. CBOD (independent, writes oxidation rate)
8. DOX (reads all upstream rates; computes saturation)
9. Pathogen (independent, reads light extinction)
10. Alkalinity (reads nitrification, denitrification, algal rates)
11. N2 (reads denitrification rate)

This order ensures all rate-variable producer→consumer edges are satisfied.

---

## Summary: Key Bugs and Differences Flagged

### Critical v2 Bugs (Must Fix in v3)

| Bug # | Location | Issue | v3 Fix |
|---|---|---|---|
| 1 | nitrogen.py:101 | Multiplicative integrator: `ammonium = 0 + ammonium * rate * dt` | Change to additive: `ammonium = ammonium + rate * dt` |
| 2 | nitrogen.py:112 | Multiplicative integrator for NO3 | Same fix |
| 3 | nitrogen.py:147, 218, 250, 313 | NaN guard `rate == np.nan` (always False) | Replace with `.isnull()` |
| 4 | floating_algae.py:122 | Multiplicative integrator + stray 86400: `algae = 0 + algae * rate * dt * 86400` | Remove both: `algae = algae + rate * dt` |
| 5 | nitrogen.py:191 | Hard-coded `half_saturation_oxygen=1` for nitrification inhibition | Wire parameter from config |
| 6 | nitrogen.py:204-212 | Hard-coded `algea_growth_rate=0` silences NO3 uptake | Read `algal_growth_rate` from registry |
| 7 | floating_algae.py:113 | Hard-coded `phosphate_fraction_dissolved=0.5` | Implement `fdp` utility and wire |
| 8 | floating_algae.py:398-401 | `ammonium_respiration()` returns 0 with TODO | Implement per v1 line 1272: `rna * respiration` |
| 9 | floating_algae.py:403-405 | `ammonium_growth()` returns 0 with TODO | Implement per v1 line 1206: `rna * PN_calculated * growth` |
| 10 | nitrogen.py, floating_algae.py | `set_at_time` never called; state updates dropped | Call `set_at_time` after integrator and clip check |

### Critical v1 Sentinel-999 Bugs (Correct in v3)

| Parameter | v1 default | v3 default | Risk |
|---|---|---|---|
| `SOD_theta` | 999 | 1.060 | 999^(T-20) blowup at T > 20 °C; catastrophic |
| `SOD_20` | 999 | 1.0 g-O2/m²/d | Wrong magnitude propagates |
| `vs` (TIP settling) | 999 | 0.1 m/d | Wrong magnitude propagates |
| `vsop` (OrgP settling) | 999 | 0.1 m/d | Wrong magnitude propagates |
| `kaw_20_user` | 999 | 0.0 (disabled) | Spurious reaeration unless explicitly set |
| `kah_20_user` | 999 | 0.0 (disabled) | Spurious reaeration unless explicitly set |

### Observable Kinetic Differences Between v1 and v2

1. **Nitrogen preference logic (lines 1206 vs 239):** v1 smooth preference factor; v2 threshold-based. Affects NH4 vs NO3 uptake fractionation.

2. **Light limitation (FL vs limit_light):** v1 single formula (Beer-Lambert); v2 offers 3 options. v3 adopts v2's flexibility with v1 default.

3. **Monod vs Limiting Nutrient (mu):** v1 multiplicative (μ = μ_max * FN * FP * FL); v2 offers limiting-nutrient and harmonic-mean options. v3 adopts v2 pattern with configurable option.

4. **DOC oxidation kinetics (line 2639):** v1 uses Monod(DOX); v3 must verify formula and implement.

5. **Benthic depth integration:** v1 carefully balances g-/m² state with mg-/L/d rates; v2 benthic process is largely a stub. v3 must implement full depth integration for benthic algae contributions to water-column processes.

6. **Algal mortality routing:** v1 routes to OrgN, OrgP, POC, DOC, POM via dedicated functions; v2 absent. v3 implements via FloatingAlgae.death_to_*() and BenthicAlgae.death_to_*() methods.

---

## Constituent Completeness Checklist

For each constituent, **green** = fully implemented in v1 and mapped to v3; **yellow** = v2 partial with known bugs; **red** = absent in v2.

| Constituent | v1 Status | v2 Status | v3 Target | Notes |
|---|---|---|---|---|
| **Ap** | Fully implemented | **Partial (bugs)** | Extend+fix | Bugs: integrator, ammonium_respiration/growth stubs, fdp hard-coded |
| **Ab** | Fully implemented | **Partial (stub)** | Extend+fix | Bugs: similar to Ap; benthic depth integration minimal |
| **NH4** | Fully implemented | **Partial (bugs)** | Extend+fix | Bugs: multiplicative integrator, hard-coded algae rates, NaN guards |
| **NO3** | Fully implemented | **Partial (bugs)** | Extend+fix | Same bugs as NH4 |
| **OrgN** | Fully implemented | **Absent** | New | Add as 3rd state in Nitrogen Process |
| **TIP** | Fully implemented | **Absent** | New | Implement with fdp partitioning |
| **OrgP** | Fully implemented | **Absent** | New | Pair with TIP in Phosphorus Process |
| **POC** | Fully implemented | **Absent** | New | Hydrolysis and settling; mortality source |
| **DOC** | Fully implemented | **Absent** | New | Oxidation with Monod(DOX); mortality source |
| **DIC** | Fully implemented | **Absent** | New | Reaeration, algal coupling, CBOD oxidation source |
| **POM** | Fully implemented | **Absent** | New | Dissolution, settling, burial |
| **CBOD** | Fully implemented | **Absent** | New | Multi-group; oxidation with Monod(DOX) |
| **DOX** | Fully implemented | **Absent** | New | Highly coupled; saturation, reaeration, SOD (sentinel bug) |
| **PX** | Fully implemented | **Absent** | New | Light-decay, settling; simple |
| **Alk** | Declared, not integrated | **Absent** | New | Nitrification/denitrification/algal coupling |
| **N2** | Fully implemented | **Absent** | New | Henry's law, denitrification source |

---

## Files Referenced

- **v1 NSM1:** `/src/clearwater_modules/nsm1/processes.py` (3,540 lines, ~290 functions)
- **v1 NSM1 constants:** `/src/clearwater_modules/nsm1/constants.py` (parameter defaults and sentinel-999 values)
- **v2 NSM1 Nitrogen:** `/src/clearwater_modules_v2/processes/nitrogen.py` (368 lines)
- **v2 NSM1 FloatingAlgae:** `/src/clearwater_modules_v2/processes/floating_algae.py` (406 lines)
- **v2 NSM1 BenthicAlgae:** `/src/clearwater_modules_v2/processes/benthic_algae.py` (expected similar structure)
- **v3 NSM1 base:** `/src/clearwater_modules_v3/processes/base.py` (Process class definition)
- **Design spec:** `/design/clearwater_modules_v3_nsm1_design_specification.md` (Sections 1-11, full context)

---

## End of Phase 0 Gap Analysis

**Next Steps:**
- Phase 1: Implement v3 utils (reaeration, sediment, light, partitioning, numerics)
- Phase 2: Fix 4 v2 processes and extend Nitrogen with OrgN
- Phases 3-6: Port remaining 12 constituents
- Phase 7: Validation, parity tests, end-to-end demo
- Phase 8: Documentation and review prep

**Contact:** Todd Steissberg (ERDC) with Claude

