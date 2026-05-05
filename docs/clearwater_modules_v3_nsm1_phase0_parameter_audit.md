# NSM1 v3 Phase 0.2 — Parameter Defaults Audit

**Phase 0 Task:** Complete inventory of v1 NSM1 parameter defaults, identifying any bad defaults beyond the documented sentinel-999 set.

**Source:** `src/clearwater_modules/nsm1/constants.py` (13 TypedDict groups, ~250 default entries).

**Date:** 2026-05-04  
**Scope:** Full parameter census with units inference from `processes.py` kinetic usage.

---

## Executive Summary: Newly Suspicious Parameters

Beyond the 6 known sentinel-999 items documented in the design spec (Section 7), this audit flagged **9 additional suspicious defaults** for Phase 1 review:

1. **rnh4_20=0, vno3_20=0** — Sediment release rates hardcoded to zero; feature is silently disabled unless overridden.
2. **rpo4_20=0** — Same as above for phosphorus; silently disabled.
3. **kdpo4=0.0** — Partitioning parameter set to 0; feature may be disabled.
4. **ksbod_20=0.0** — CBOD sedimentation rate zero; feature disabled (settled CBOD does not sink).
5. **apx=1, vx=1** — Pathogen parameters of 1 (unvalidated placeholders?).
6. **h2=0.1** — POM dissolution rate denominator; needs clarification on its physical role.
7. **vb=0.01** — Burial velocity (m/d); suspicious magnitude relative to settling velocities.
8. **dt=1, depth=1.5, flow=2, slope=2** — GlobalVars runtime defaults are toy placeholders, not sensible physical values.
9. **Alkalinity coefficients (r_alkXX)** — Computed as ratios; dimensionality check recommended.

---

## Complete Parameter Inventory by TypedDict Group

### Group 1: AlgaeStaticVariables (Floating Phytoplankton)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| AWd | 100 | mg-D/ug-Chla | Dry algal weight ratio; used in stoichiometry | OK |
| AWc | 40 | mg-C/ug-Chla | Algal carbon ratio | OK |
| AWn | 7.2 | mg-N/ug-Chla | Algal nitrogen ratio; N uptake & respiration coupling | OK |
| AWp | 1 | mg-P/ug-Chla | Algal phosphorus ratio | OK |
| AWa | 1000 | ug-Chla/ug-Chla | Chlorophyll per algal unit (self-ref); `rda = AWd/AWa` | suspicious-nonconst |
| KL | 10 | W/m² | Light limitation half-saturation (Smith/Steele models); reasonable photic range | OK |
| KsN | 0.04 | mg-N/L | N half-saturation; typical ~0.01–0.1 range | OK |
| KsP | 0.0012 | mg-P/L | P half-saturation; typical ~0.001–0.01 range | OK |
| mu_max_20 | 1 | 1/d | Max algal growth rate at 20 °C; ~0.5–2 d⁻¹ range | OK |
| kdp_20 | 0.15 | 1/d | Algal death rate at 20 °C | OK |
| krp_20 | 0.2 | 1/d | Algal respiration rate at 20 °C | OK |
| mu_max_theta | 1.047 | unitless | Arrhenius coeff for growth; standard van't Hoff value | OK |
| kdp_theta | 1.047 | unitless | Arrhenius coeff for death | OK |
| krp_theta | 1.047 | unitless | Arrhenius coeff for respiration | OK |
| vsap | 0.15 | m/d | Algal settling velocity; ~0.1–0.5 m/d typical | OK |
| growth_rate_option | 1 | 1=Mult, 2=Min, 3=Harmonic | Algal growth limitation method selector | OK |
| light_limitation_option | 1 | 1=HS, 2=Smith, 3=Steele | Light limitation formula selector | OK |

**Issues:** `AWa=1000` is used as denominator in `rda=AWd/AWa`; suspicious dimensionality (ug-Chla appears in both).

---

### Group 2: AlkalinityStaticVariables (Stoichiometric Ratios)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| r_alkaa | 14/(106×12×1000) ≈ 1.099e-4 | eq/mg-C | Alkalinity stoich: algal photosynthesis → Alk | OK |
| r_alkan | 18/(106×12×1000) ≈ 1.419e-4 | eq/mg-C | Alkalinity stoich: algal growth N uptake → Alk | OK |
| r_alkn | 2/(14×1000) ≈ 1.429e-4 | eq/mg-N | Alkalinity stoich: nitrification → Alk sink | OK |
| r_alkden | 4/(14×1000) ≈ 2.857e-4 | eq/mg-N | Alkalinity stoich: denitrification → Alk source | OK |
| r_alkba | 14/(106×12×1000) ≈ 1.099e-4 | eq/mg-C | Benthic algae photosynthesis → Alk | OK |
| r_alkbn | 18/(106×12×1000) ≈ 1.419e-4 | eq/mg-C | Benthic algae N uptake → Alk | OK |

**Issues:** All are computed ratios (stoichiometric formulas). Review dimensional consistency in processes.py alkalinity sources/sinks.

---

### Group 3: BalgaeStaticVariables (Benthic Algae)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| BWd | 100 | mg-D/g-D | Dry biomass weight | OK |
| BWc | 40 | mg-C/g-D | Carbon in benthic algae | OK |
| BWn | 7.2 | mg-N/g-D | Nitrogen in benthic algae | OK |
| BWp | 1 | mg-P/g-D | Phosphorus in benthic algae | OK |
| BWa | 3500 | g-D/m² | Benthic algae density at unit area; reaeration context | OK |
| KLb | 10 | W/m² | Benthic light half-sat; same as floating algae | OK |
| KsNb | 0.25 | mg-N/L | N half-sat for benthic algae; higher than floating (0.04) | OK |
| KsPb | 0.125 | mg-P/L | P half-sat; higher than floating (0.0012) | OK |
| Ksb | 10 | g-D/m² | Benthic algae half-sat for space limitation (FSb factor) | OK |
| mub_max_20 | 0.4 | 1/d | Max benthic algae growth; lower than floating (1) | OK |
| krb_20 | 0.2 | 1/d | Benthic algae respiration | OK |
| kdb_20 | 0.3 | 1/d | Benthic algae death | OK |
| mub_max_theta | 1.047 | unitless | Arrhenius coeff for growth | OK |
| krb_theta | 1.06 | unitless | Arrhenius coeff for respiration (1.06 vs 1.047 elsewhere) | OK |
| kdb_theta | 1.047 | unitless | Arrhenius coeff for death | OK |
| b_growth_rate_option | 1 | selector | Growth limitation method (same as floating) | OK |
| b_light_limitation_option | 1 | selector | Light limitation method | OK |
| Fw | 0.9 | unitless | Fraction of benthic photosynthesis as CO₂ release (0.1 as DOC) | OK |
| Fb | 0.9 | unitless | Fraction of benthic biomass released in death (0.1 as DOC) | OK |

**Issues:** `krb_theta=1.06` differs from the standard 1.047; documented in code or accidental?

---

### Group 4: NitrogenStaticVariables (Inorganic N Cycling)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| KNR | 0.6 | mg-O₂/L | Oxygen half-sat for nitrification inhibition (1 - exp(-KNR×DOX)) | OK |
| knit_20 | 0.1 | 1/d | Nitrification rate at 20 °C | OK |
| kon_20 | 0.1 | 1/d | Nitrification oxidation rate (ammonia oxidation) at 20 °C | OK |
| kdnit_20 | 0.002 | 1/d | Denitrification rate at 20 °C | OK |
| rnh4_20 | **0** | 1/d | Sediment NH₄ release rate at 20 °C | **suspicious-zero** |
| vno3_20 | **0** | 1/d | Sediment NO₃ denitrification rate at 20 °C | **suspicious-zero** |
| knit_theta | 1.083 | unitless | Arrhenius coeff for nitrification | OK |
| kon_theta | 1.074 | unitless | Arrhenius coeff for oxidation | OK |
| kdnit_theta | 1.08 | unitless | Arrhenius coeff for denitrification | OK |
| rnh4_theta | 1.047 | unitless | Arrhenius coeff for sediment NH₄ release | OK |
| vno3_theta | 1.045 | unitless | Arrhenius coeff for sediment NO₃ denit | OK |
| KsOxdn | 0.1 | mg-O₂/L | Oxygen half-sat for denitrification inhibition | OK |
| PN | 0.5 | unitless | Fraction of algal N uptake from NH₄ (vs NO₃) | OK |
| PNb | 0.5 | unitless | Fraction of benthic algal N uptake from NH₄ | OK |

**Issues:** `rnh4_20=0` and `vno3_20=0` disable sediment flux features; hardcoded to zero and not overrideable in v1 unless parameters are updated post-construction.

---

### Group 5: CarbonStaticVariables (POC, DOC, DIC)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| f_pocp | 0.9 | unitless | Fraction of algal death routed to POC (0.1 to DOC) | OK |
| kdoc_20 | 0.01 | 1/d | DOC decomposition rate at 20 °C | OK |
| kdoc_theta | 1.047 | unitless | Arrhenius coeff for DOC decay | OK |
| f_pocb | 0.9 | unitless | Fraction of benthic death routed to POC | OK |
| kpoc_20 | 0.005 | 1/d | POC decomposition rate at 20 °C | OK |
| kpoc_theta | 1.047 | unitless | Arrhenius coeff for POC decay | OK |
| KsOxmc | 1.0 | mg-O₂/L | Oxygen half-sat for mineralization (DOC & POC oxidation) | OK |
| pCO2 | 383.0 | ppm | Atmospheric CO₂ partial pressure; ~2024 value | OK |
| FCO2 | 0.2 | unitless | Fraction of DIC as CO₂ (vs HCO₃⁻, CO₃²⁻); pH-dependent; placeholder for NSM2 carbonate solver | suspicious-magnitude |
| roc | 32/12 ≈ 2.667 | mg-O₂/mg-C | Respiration stoich: O₂ consumed per C oxidized (ideal: 32/12) | OK |

**Issues:** `FCO2=0.2` is a rough average; carbonate speciation is NSM2 territory (v3 1.0.0 treats Alk as simple tracer, not full carbonate equilibrium).

---

### Group 6: CBODStaticVariables (Multi-Group CBOD)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| KsOxbod | 0.5 | mg-O₂/L | Oxygen half-sat for CBOD oxidation | OK |
| kbod_20 | 0.12 | 1/d | CBOD oxidation rate at 20 °C | OK |
| ksbod_20 | **0.0** | m/d | CBOD sedimentation (settling) rate at 20 °C | **suspicious-zero** |
| kbod_theta | 1.047 | unitless | Arrhenius coeff for CBOD decay | OK |
| ksbod_theta | 1.047 | unitless | Arrhenius coeff for CBOD settling | OK |

**Issues:** `ksbod_20=0.0` means CBOD never settles (sedimentation disabled). If CBOD is expected to settle, this is a bug. If CBOD is modeled as fully dissolved/suspended without settling, this is correct by design.

---

### Group 7: DOXStaticVariables (Dissolved Oxygen)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| ron | (2.0×32)/14 ≈ 4.571 | mg-O₂/mg-N | Stoich: O₂ consumed per N nitrified | OK |
| KsSOD | 1 | mg-O₂/L | Oxygen half-sat for sediment oxygen demand (SOD) | OK |

**Issues:** None identified; both are stoichiometric constants or Monod parameters.

---

### Group 8: N2StaticVariables (Dissolved Nitrogen Gas)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| (empty) | — | — | No parameters in v1 (N₂ atmospheric exchange via Henry's law, not parameterized) | OK |

---

### Group 9: POMStaticVariables (Particulate Organic Matter)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| kpom_20 | 0.1 | 1/d | POM dissolution rate at 20 °C | OK |
| h2 | **0.1** | m | Sediment burial/sedimentation depth denominator; unclear physical role | **unclear-role** |
| kpom_theta | 1.047 | unitless | Arrhenius coeff for POM dissolution | OK |

**Issues:** `h2=0.1` is used as a divisor in burial calculations (e.g., `vb * POM / h2`); needs clarification: is it a reference depth for sediment flux scaling, or something else?

---

### Group 10: PathogenStaticVariables (Indicator Pathogens)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| kdx_20 | 0.8 | 1/d | Pathogen decay (inactivation) rate at 20 °C | OK |
| kdx_theta | 1.07 | unitless | Arrhenius coeff for decay | OK |
| apx | **1** | unitless | Algal production of pathogen (stoich?); placeholder value | **suspicious-magnitude** |
| vx | **1** | m/d | Pathogen settling velocity | **suspicious-magnitude** |

**Issues:** `apx=1` and `vx=1` appear to be placeholder values; no literature basis documented in the code.

---

### Group 11: PhosphorusStaticVariables (Organic & Inorganic P)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| kop_20 | 0.1 | 1/d | Organic P decomposition rate at 20 °C | OK |
| rpo4_20 | **0** | 1/d | Sediment P release (DIP from bed) rate at 20 °C | **suspicious-zero** |
| kop_theta | 1.047 | unitless | Arrhenius coeff for OrgP decay | OK |
| rpo4_theta | 1.074 | unitless | Arrhenius coeff for sediment P release | OK |
| kdpo4 | **0.0** | L/kg | Partition coeff for DIP adsorption onto solids; feature disabled | **suspicious-zero** |

**Issues:** `rpo4_20=0` silently disables sediment P flux (like rnh4_20). `kdpo4=0.0` disables P partitioning.

---

### Group 12: GlobalParameters (Feature Flags)

All 16 parameters are boolean switches (True/False) enabling/disabling modules:
`use_NH4`, `use_NO3`, `use_OrgN`, `use_OrgP`, `use_TIP`, `use_SedFlux`, `use_POC`, `use_DOC`, `use_DOX`, `use_DIC`, `use_Algae`, `use_Balgae`, `use_N2`, `use_Pathogen`, `use_Alk`, `use_POM`.

| Parameter | Default | Comment | Sanity Flag |
|-----------|---------|---------|-------------|
| use_NH4...use_Alk | all True | All constituents enabled by default | OK |
| use_SedFlux | **False** | Sediment fluxes disabled by default (but see rnh4_20=0, etc., above) | OK |
| use_POM | True | POM enabled | OK |

**Issues:** `use_SedFlux=False` makes the zero values of `rnh4_20`, `vno3_20`, `rpo4_20` a non-issue *if* the flag properly gates all sediment flux calculations. Verify in processes.py that sediment fluxes are indeed gated by `use_SedFlux`.

---

### Group 13: GlobalVars (Runtime Scalars and Options)

| Parameter | Default | Units | Role | Sanity Flag |
|-----------|---------|-------|------|-------------|
| **Settling velocities** | | | | |
| vson | 0.01 | m/d | Organic N settling; very slow | OK |
| vsoc | 0.01 | m/d | POC settling; very slow | OK |
| vsop | **999** | m/d | Organic P settling | **sentinel-999** |
| vs | **999** | m/d | TIP settling (sediment flux feature) | **sentinel-999** |
| **Sediment oxygen demand** | | | | |
| SOD_20 | **999** | g-O₂/m²/d | SOD at 20 °C | **sentinel-999** |
| SOD_theta | **999** | unitless | Arrhenius coeff for SOD (catastrophic: 999^(T−20)) | **sentinel-999** |
| **Reaeration overrides** | | | | |
| kaw_20_user | **999** | m/d | User-override wind reaeration at 20 °C (only if option=1) | **sentinel-999** |
| kah_20_user | **999** | 1/d | User-override hydraulic reaeration at 20 °C (only if option=1) | **sentinel-999** |
| **Reaeration Arrhenius** | | | | |
| kaw_theta | 1.024 | unitless | Arrhenius coeff for wind reaeration | OK |
| kah_theta | 1.024 | unitless | Arrhenius coeff for hydraulic reaeration | OK |
| **Reaeration method selectors** | | | | |
| hydraulic_reaeration_option | 1 | 1–9 (9 formulae available) | Hydraulic reaeration calculation method | OK |
| wind_reaeration_option | 1 | 1–13 (13 formulae available) | Wind reaeration calculation method | OK |
| **Simulation and environmental inputs** | | | | |
| dt | **1** | d | Timestep (Global default; typically overridden by model) | **suspicious-nonconst** |
| depth | **1.5** | m | Water depth (toy value; spatially variable, overridden per cell) | **suspicious-nonconst** |
| TwaterC | 20 | °C | Water temperature (typical reference; overridden per cell/time) | OK |
| velocity | **1** | m/s | Velocity (toy value) | **suspicious-nonconst** |
| flow | **2** | m³/s | Flow (toy value) | **suspicious-nonconst** |
| topwidth | **1** | m | Top width (toy value) | **suspicious-nonconst** |
| slope | **2** | unitless | Slope (toy value; way too steep for realistic streams) | **suspicious-nonconst** |
| shear_velocity | **4** | m/s | Shear velocity (toy value) | **suspicious-nonconst** |
| pressure_mb | 2026.5 | hPa | Atmospheric pressure (roughly sea-level; ~1013 hPa actual) | **suspicious-magnitude** |
| wind_speed | **4** | m/s | Wind speed (toy value) | **suspicious-nonconst** |
| q_solar | **500** | W/m² (or day⁻¹?) | Solar radiation; units unclear (spec says "1/d" in comments) | **unclear-role** |
| **Light attenuation (Beer-Lambert)** | | | | |
| Solid | 1 | mg/L | Suspended solids (toy value) | **suspicious-nonconst** |
| lambda0 | 0.02 | 1/m | Background light extinction; reasonable for clear water | OK |
| lambda1 | 0.0088 | (1/m)/(ug-Chla/L) | Light extinction per unit chlorophyll; linear self-shading | OK |
| lambda2 | 0.054 | unitless | Light extinction from algae (non-linear term exponent 2/3); unitless coeff | OK |
| lambdas | 0.052 | L/(mg·m) | Light extinction from ISS (suspended solids); disabled in code (line 267: `# + lambdas * Solid`) | OK |
| lambdam | 0.0174 | L/(mg·m) | Light extinction from POM (suspended matter) | OK |
| Fr_PAR | 0.47 | unitless | Fraction of solar radiation in PAR (photosynthetically active radiation); ~47% typical | OK |

**Issues:**
- `vsop=999, vs=999, SOD_20=999, SOD_theta=999, kaw_20_user=999, kah_20_user=999` — **Known sentinel-999 bugs** (design spec Section 7).
- `dt, depth, velocity, flow, topwidth, slope, shear_velocity, wind_speed, Solid` — Toy placeholder values, not suitable for real simulations (overridden per cell/time).
- `pressure_mb=2026.5` — Unrealistic (~2× atmospheric); likely a data-entry error (should be ~1000–1050 hPa).
- `q_solar` units unclear ("1/d" in PAR function docstring; likely W/m² or MJ/m²/d).

---

## Summary: Sanity Flags and Recommendations

### Sentinel-999 Items (Known Bugs, Design Spec Section 7)

| Parameter | Issue | Phase 1 Fix |
|-----------|-------|------------|
| vsop | Multiplied into rate; wrong order of magnitude | Default to 0.1 m/d |
| vs | Same | Default to 0.1 m/d |
| SOD_20 | Wrong magnitude propagates | Default to 1.0 g-O₂/m²/d |
| SOD_theta | 999^(T-20) catastrophic Arrhenius blowup | Default to 1.060 |
| kaw_20_user | Spurious if set incorrectly | Default to 0.0 (disabled unless user opts in) |
| kah_20_user | Same | Default to 0.0 |

### Newly Flagged Items (Beyond Known Sentinel-999)

| Parameter | Group | Issue | Severity | Recommendation |
|-----------|-------|-------|----------|-----------------|
| rnh4_20 | Nitrogen | Hardcoded 0; sediment NH₄ flux disabled | Medium | Document as intentional (feature disabled by default); ensure `use_SedFlux` gates this correctly |
| vno3_20 | Nitrogen | Hardcoded 0; sediment NO₃ denit disabled | Medium | Same as above |
| rpo4_20 | Phosphorus | Hardcoded 0; sediment P flux disabled | Medium | Same as above |
| ksbod_20 | CBOD | Hardcoded 0; CBOD settling disabled | Medium | Clarify design intent: is CBOD modeled as fully suspended? If yes, OK; if no, change to ~0.001–0.005 m/d |
| kdpo4 | Phosphorus | Hardcoded 0; P partitioning feature disabled | Low | Document as intentional; NSM2 feature |
| h2 | POM | Unclear physical role (divisor in burial calc) | Low | Add docstring to constants.py and processes.py explaining h2's role |
| vx | Pathogen | Placeholder value (1 m/d); no literature basis | Low | Provide rationale or reference; plausible range ~0.01–1 m/d |
| apx | Pathogen | Placeholder value (1); dimensionality unclear | Low | Provide rationale or convert to a dimensionless fraction |
| pressure_mb | GlobalVars | Unrealistic (~2× atm); likely typo | High | Correct to ~1013 hPa (or make spatially/temporally variable) |

### Parameters with Toy/Placeholder Values (Expected for v1 Constants; Overridden at Runtime)

These are not bugs but design-as-intended defaults that get overridden by the model at runtime or per-cell:

- `dt, depth, velocity, flow, topwidth, slope, shear_velocity, wind_speed, Solid, q_solar`

**Action:** Verify in model initialization that these placeholders are always overridden before use in kinetics calculations.

---

## Detailed Review Notes for Phase 1

### 1. Sediment Flux Feature (`rnh4_20=0, vno3_20=0, rpo4_20=0`)

**Context:** The design spec (Section 14) states:
> Sediment-flux parameters (SOD_20, NH4fromBed, DIPfromBed, NO3_BedDenit, DIC_sed_release) are scalar globals in v3 1.0.0, applied uniformly to all cells, set in YAML (matching v1's pattern exactly). Per-cell spatially varying fluxes and dynamically computed fluxes both arrive in v3 1.1+ via the NSM2 sediment diagenesis Process.

**Finding:** In v1, sediment fluxes are controlled by `use_SedFlux=False` (default). The parameters `rnh4_20, vno3_20, rpo4_20` are zero, which means they're disabled even if `use_SedFlux` were True.

**Recommendation:** Verify that `processes.py` correctly gates sediment flux functions with `use_SedFlux`, and document that these three parameters are placeholders (to be replaced in NSM1 v3 with proper sediment flux objects in Phase 1).

### 2. CBOD Settling (`ksbod_20=0.0`)

**Context:** In processes.py, CBOD settling is calculated as:
```python
CBOD_sedimentation = CBOD * ksbod_tc  # then divided by h2 in burial calculation
```

**Finding:** With `ksbod_20=0`, CBOD never settles; the settling feature is disabled.

**Question:** Is this by design (CBOD modeled as fully suspended/dissolved without settling) or an oversight?

**Recommendation:** Code review with domain expert (LimnoTech) to confirm design intent. If CBOD should settle, change to ~0.001–0.005 m/d. If fully suspended is correct, document in v3 DEFAULTS.

### 3. Pathogen Parameters (`apx=1, vx=1`)

**Finding:** `apx` (algal production of pathogen) and `vx` (pathogen settling velocity) are both set to 1 with no units or literature basis in comments.

**Recommendation:** Add docstrings to constants.py clarifying units and providing a reference or rationale. Plausible ranges:
- `vx` (settling velocity): 0.01–1 m/d (similar to algae)
- `apx` (production stoich): needs clarity on dimensionality

### 4. Phosphorus Partitioning (`kdpo4=0.0`)

**Finding:** P adsorption onto solids is disabled (partition coeff = 0).

**Recommendation:** Document in constants.py as an NSM2 feature (when full sediment diagenesis is added). In v3 1.0.0, all TIP is treated as dissolved (no partitioning).

### 5. Light Attenuation (`lambdas` Unused)

**Finding:** In processes.py line 267, the `lambdas * Solid` term is commented out:
```python
L: xr.DataArray = lambda0  # + lambdas * Solid
```

**Recommendation:** Either uncomment (if intended) or remove the unused parameter from constants.py to reduce cognitive load. If kept, add a comment explaining why it's disabled.

### 6. POM Sedimentation Denominator (`h2`)

**Finding:** The parameter `h2=0.1` is used as a divisor in burial calculations (e.g., `vb * POM / h2`), but its physical meaning is not clearly documented.

**Recommendation:** Add a detailed docstring to constants.py:
```
h2: float = 0.1  # Reference depth (m) for sediment burial scaling. 
                 # Used in burial calculations as h2_reference = vb * state / h2.
                 # Dimensions: meters (same as water depth).
```

### 7. Atmospheric Pressure (`pressure_mb=2026.5`)

**Finding:** The value 2026.5 is approximately 2× atmospheric pressure (~1013 hPa at sea level). This is either a data-entry error or a non-standard choice.

**Recommendation:** High priority. Either:
- Correct to ~1000–1050 hPa (realistic sea-level range), OR
- Provide a reference for the unusual value, OR
- Make it spatially/temporally variable (loaded per simulation).

### 8. Solar Radiation (`q_solar`) Units

**Finding:** The parameter `q_solar=500` has unclear units (docstring in processes.py says "1/d", but likely means W/m² or MJ/m²/d).

**Recommendation:** Clarify units in constants.py docstring and standardize across v3 NSM1 (likely W/m² for consistency with PAR calculations).

---

## Conclusion

The parameter audit identified **9 newly suspicious defaults** beyond the 6 known sentinel-999 items. Most are either:
1. **Design-as-intended** (zero-valued features disabled by default, which is OK if gated by `use_SedFlux`),
2. **Placeholder toy values** (overridden at runtime), or
3. **NSM2-territory features** (partitioning, multi-pool sediment) deferred to v3 1.1+.

The highest-priority issues are:
- **`pressure_mb=2026.5`** — likely a typo; correct to ~1013 hPa.
- **`ksbod_20=0`** — confirm design intent with domain expert (settling disabled?).
- **Sediment flux gates** — verify `use_SedFlux=False` correctly disables all `rnh4_20`, `vno3_20`, `rpo4_20` terms.

All findings are documented in this audit for Phase 1 review and parameter-defaults-corrections.md (to be created in Phase 1).

