# ClearWater Modules v3 — NSM2 Features Implementation Plan (Multi-Group Algae + Remaining NSM2 Features)

**Status:** Draft, awaiting approval
**Author:** Todd Steissberg (ERDC), with Claude
**Date:** 2026-05-16
**Scope:** A single merged, phased implementation plan that carries forward **Step 1 — multi-group algae** (from `clearwater_modules_v3_multigroup_algae_design_specification.md`) and adds the previously-missing **Step 4 — remaining NSM2 features**: multi-pool organic matter, alkalinity/pH carbonate system, methane/sulfide, and silica. It also captures the **Step-3 sediment diagenesis design directive** (the recommended build approach and its extension points) and defines the Step-3↔Step-4 sediment-flux interface as a dependency boundary so Step 4 can ship on parameterized sediment fluxes without blocking on the Di Toro diagenesis module.

**Read this with two companions.** Umbrella conventions (env, branch, package architecture, integrator contract) live in `clearwater_modules_v3_architecture_specification.md`. The full Step-1 detail (component inventory, YAML schema, file-by-file change list, Step-1 open questions) lives in `clearwater_modules_v3_multigroup_algae_design_specification.md` and is **not duplicated** here — this document condenses Step 1 to its load-bearing conventions and focuses new material on Step 4.

> **River-first scope (binding).** ClearWater-Riverine and HEC-RAS-2D are depth-averaged and will remain so for the foreseeable future. Every algorithm in this plan is specified for a **vertically well-mixed water column**. Reservoir-stratification-only mechanisms (buoyancy-driven vertical migration, hypolimnetic anoxic accumulation, thermocline-gated fluxes, vertical nutrient gradients) are **out of scope**. The plan is nonetheless written so a future vertically-resolved transport driver can reuse every state variable, parameter, and kinetic term unchanged (Section 2.2).

---

## 1. Roadmap Placement

The four-step HAB / NSM2 roadmap (from the Step-1 spec, Section 1):

1. **Multi-group algae** — `algae_group` / `balgae_group` dimension and per-group parameterization. *Condensed in Section 3 of this plan.*
2. **HAB capabilities** — N-fixation, biomass floor + rewet seeding, toxin tracer, photo-inhibition, bell-curve T-response, P-luxury, source-attribution tracers, low-DO mortality. *Not covered here; Section 2.3 flags one Step-2-adjacent decision the river-scope question reopens.*
3. **Sediment diagenesis** — Di Toro multi-G (JNH4, JNO3, JCH4, JSO4, JH2S, JDIC, JDIP, SOD) as a v3-native process. *This plan now carries the Step-3 design directive (Section 4.E, author-approved 2026-05-16) — the recommended build approach and its extension points — together with the Step-3↔Step-4 interface boundary. Step 3 itself remains off this plan's critical path; a phased sketch is in Section 6.*
4. **Remaining NSM2 features** — multi-pool organic matter, alkalinity/pH, methane/sulfide, silica. **This is the new material (Section 4).**

This merged document exists because Step 4 had only roadmap placement and scattered deferral notes (gold-standard spec §8/§9; NSM1 design spec lines 78–79, 157, 356, 645–647), with no worked plan. Step 1 had a worked plan but in a separate spec. The two are now joined because Step 4 *consumes Step-1 machinery* — the multi-pool organic-matter pools replace the lumped `algal_orgn/orgp/poc/doc_from_mortality_rate` partition that Step 1's consumers `.sum("algae_group")`, and silica fills the `si_limitation_option` / `KsSi` / `AWsi` per-group hooks Step 1 deliberately reserved (Step-1 spec lines 102–104, 168–169, 431, 469–471).

### Authoritative sources and V&V baseline

The NSM2 Python modules on `main` (`src/clearwater_modules/nsm2/*.py`) are **empty scaffolding** (0-byte stubs); there is no prototype to port. The authoritative algorithmic references are, in priority order:

1. **NSM2 Fortran (HEC-RAS-WQ)** — the regression baseline. The V&V chain is **Fortran NSM2 → v3** (mirrors the NSM1 chain Fortran→v1→v3; v1 is not an oracle and there is no v1 NSM2).
2. **CE-QUAL-W2 v2026.02** — cross-check for the carbonate solver, silica (`ASI`/`AHSSI`), and multi-pool OM; local source path on file.
3. **QUAL2K/QUAL2Kw, WASP, CE-QUAL-ICM, Chapra (1997), Stumm & Morgan** — cross-checks and idea sources. Per project convention, a capability present in a cross-check but absent here is a *candidate*, not a deficiency; conformance is to NSM2 Fortran, not to QUAL2K/W2.

### 1.1 Merge gate and cutover sequencing (binding)

This plan (Step 1 and Step 4) must not merge into the integration target until the v3 1.0 cutover has run. The cutover sequence is the committed, decided plan in `clearwater_processes_library_plan.md` §6:

1. Close the open v3 NSM1 review CRITICALs (the v3-vs-v1 verification gate).
2. Retire v1 (`clearwater_modules`): delete the package; migrate parity tests to **frozen reference values**.
3. Merge `streaming` → upstream `EcohydrologyTeam/ClearWater-modules` `main` (the v3 1.0 integration).
4. In one discrete step at that point: rename the repository `ClearWater-modules` → `ClearWater-processes`; rename the package `clearwater_modules*` → `clearwater.processes` (a breaking import change, no redirect); set the APL identity.

That plan states its own execution gate: "nothing in §6 begins before step 1 (CRITICALs closed). The rename lands as part of the upstream-merge step, not piecemeal on the `streaming` fork beforehand."

**Binding rules for this plan:**

- **Develop on a separate branch off `streaming`** (descriptive kebab-case name, no `feature/` prefix, matching house style — e.g. `multigroup-algae`). Parallel development is fine; an unmerged branch is invisible to the v3-vs-v1 review, which runs on `streaming`.
- **Do not merge Step 1 or Step 4 into `streaming` or the integration target until after §6 step 3.** Merging earlier makes the algae and consumer code a moving target while it is under v3-vs-v1 review, and collides with the bit-identical baseline contract (next rule). This is the "piecemeal on the `streaming` fork beforehand" that the §6 gate forbids.
- **`nal=1` must reproduce the terminal gold-standard baseline exactly.** The coupled-demo parity check (`tests/v3/nsm1/baseline/check_baseline_parity.py`) compares against `baseline_coupled_trajectory_b51df71.nc` with `numpy.array_equal` — zero tolerance, no rtol or atol. The Step-1 spec's phrase "numerically identical within float tolerance" is weaker than the contract the test actually enforces; the binding bar is **bit-identical, with no new baseline**. Achieve it by squeezing the length-1 group dimension at the registry boundary and preserving floating-point operand order in the rewritten algae rate methods. Re-baselining a structural refactor is **not** the intended path: per `tests/v3/nsm1/baseline/README.md`, operand-reordering / broadcast-shift is exactly the regression class this contract exists to catch, and re-baseline is reserved for deliberate kinetics changes or reviewed default-parameter changes. A signed-off re-baseline (separate commit, hash in filename, prior baselines retained) is a documented last resort only — used only if a written operand-order analysis proves bit-identity is unattainable.
- **The parity-test harness code must change.** The hard-coded two-dimensional `(n_substeps+1, n_cells)` trajectory buffer and `check_baseline_parity.py` must squeeze the group dimension back to the exact prior shape. This is reviewable test infrastructure in the branch; it must squeeze, never loosen tolerance.
- **Two porting boundaries cross the cutover.** A long-lived pre-cutover branch incurs both: (i) §6 step 2 retires v1 and converts parity to frozen reference values — any `nal=1` "v1 parity" contract must be re-pointed to those frozen references; (ii) §6 step 4 renames the package — `clearwater_modules_v3` / `clearwater_modules*` → `clearwater.processes` across every touched module and the baseline harness (`check_baseline_parity.py`, `capture_baseline_trajectory.py`, the demo builder). Lower-risk option: start the branch after §6 step 4, since the Step-1 estimate (~9–12 days) is short relative to the porting cost, and multi-group is a gold-standard-spec §8 deferral — explicitly post-gate, not a merge blocker.

**Authoritative override.** The multi-group spec and earlier text in this plan anchor `nal=1` parity to "v1 parity" and to the `clearwater_modules_v3` package path. Both referents are invalidated by §6 steps 2 and 4. Where they conflict, this section governs: the parity contract targets the frozen reference values and the post-rename `clearwater.processes` path.

---

## 2. Scoping Principles (Binding)

### 2.1 Reaction module is stratification-agnostic; transport owns vertical structure

The contract: **the reaction kernel computes per-cell source/sink rates given that cell's state; it never assumes how cells are arranged vertically.** Vertical structure (layering, stratification, hypolimnetic transport) is the transport driver's responsibility (ClearWater-Riverine today; a future reservoir driver later). Consequences:

- All Step-4 kinetics are local thermodynamic/kinetic functions of the cell's own concentrations and temperature. None reads a "layer below me" or a "thermocline" state.
- Settling and deposition are expressed as a **per-cell downward flux** (`vs · C`) returned to the framework, not as an internal layer-to-layer transfer. A depth-averaged 2D river treats it as a bed-deposition sink; a future layered driver routes the same flux between layers. The kinetic code is identical in both.
- "Cell" is the unit of computation. Nothing collapses an array to a single layer inside the kernel. If a future driver introduces `(cell, layer)`, the kernel broadcasts over it exactly as Step 1 broadcasts over `algae_group`.

### 2.2 Reservoir-readiness ("don't paint into corners")

| Keep (reservoir-ready by construction) | Exclude now (stratification-only; not river-relevant) |
|---|---|
| Carbonate-system inputs include **ionic strength** (from salinity / specific conductance), so the solver is valid for freshwater rivers *and* future brackish/estuarine/hypolimnetic water without re-derivation | Vertical migration / buoyancy-driven repositioning of biomass between transported layers |
| CH₄/H₂S **water-column kinetics + sediment source flux + degassing** — identical equations in a well-mixed river reach or a future anoxic hypolimnion cell | Hypolimnetic accumulation logic, oxycline-gated CH₄/H₂S build-up (a transport/stratification outcome, not a kinetic) |
| Settling/deposition as per-cell flux (Section 2.1) | Internal vertical re-suspension/settling cascades across layers |
| Per-cell parameters that *could* later vary by `(cell, layer)` | Thermocline-dependent or depth-class parameter switching |
| Nutrient/Si/OM kinetics that are vertically uniform within a river cell | Vertical nutrient gradients (a stratified-reservoir feature) |

The discipline: **never hard-code "well-mixed"** in a way a 1DV/3D driver would have to tear out. The kernel is well-mixed *because the transport grid is*, not because the kinetics assume it.

### 2.3 Flagged divergence — buoyant cyanobacteria and vertical light position

The Step-1 spec treats buoyancy as one blanket deferral: `is_buoyant` is "buoyancy/migration flag (step 2; lake/3D only)" (Step-1 spec line 100; Appendix A line 433). The river-scope question shows that framing is too coarse. There are **two separable mechanisms**:

- **(a) Vertical light-position sub-model — river-admissible.** Buoyant cyanobacteria concentrate near the surface (scum) and therefore experience a *higher mean irradiance* than a vertically-uniform population. This can be represented as a within-cell, sub-grid **light-weighting**: compute the group's effective light-limitation factor assuming it occupies the top fraction `f_surf` of the depth, integrating the Beer–Lambert profile over that sub-layer only. It touches **only the local light-limitation term** — a non-conservative, cell-local quantity. It introduces **no transported vertical dimension and no vertical nutrient gradient**, which is physically defensible in a turbulent river: solutes mix vertically far faster than positively-buoyant cells re-stratify, so treating nutrients as vertically uniform while light is depth-weighted is internally consistent.
- **(b) Buoyancy-driven vertical migration — reservoir-only, excluded.** Moving biomass mass between transported layers in response to light/nutrient history. This is transport, requires vertical resolution, and is correctly deferred.

**Proposed disposition:** keep (b) excluded; reclassify (a) as a **river-admissible optional per-group light sub-model**, default off, delivered in **Step 2 (HAB)**, *not* Step 4. It is recorded here because the river-scope question raised it and it must not be lost to the Step-1 blanket deferral. **Open question Q1** asks the author to confirm this split and the Step-2 placement. (Answering the literal question: yes — depth-weighted light with vertically-uniform nutrients is a sound approximation for a depth-averaged river; it is the *migration*, not the light weighting, that needs stratification.)

---

## 3. Step 1 — Multi-Group Algae (Condensed)

Full detail in `clearwater_modules_v3_multigroup_algae_design_specification.md`. The load-bearing conventions Step 4 depends on:

- **Group dimension.** `nal` floating groups, `nalb` benthic groups (globals in `parameters/global_parameters.py`). State variables `algae_floating` → `(cell, algae_group)`, `benthic_algae` → `(cell, balgae_group)`. Default `nal=nalb=1` is numerically identical to current single-group v3 (within float tolerance).
- **Per-group parameters.** `GROUP_DEFAULTS` template in `parameters/{algae,balgae}.py`; YAML accepts a flat dict (broadcast to length 1) or a list of length `nal`. `DEFAULTS = GROUP_DEFAULTS` alias preserved for migration.
- **Consumer collapse.** Downstream processes `.sum("algae_group")` before applying their own kinetics. **This is the seam Step 4 modifies:** the mortality partition currently summed (`algal_orgn/orgp/poc/doc_from_mortality_rate`) is repartitioned into the multi-pool OM pools (Section 4.A).
- **Reserved Step-4 hooks (already in `GROUP_DEFAULTS`).** `si_limitation_option` (1=unlimited; Step 4 enables 2=Monod), `KsSi`, `AWsi`. Step 4 wires these — no YAML migration for users who configured them under Step 1.
- **Phased convention & estimate.** Step 1 = Phases 0–7, ~9–12 working days, with `nal=1` parity as the regression contract. Carried into the merged schedule (Section 6) unchanged.
- **Step-1 open questions** (default-block broadcast, audit-bug sequencing, consumer-sum helper, BenthicAlgae release cadence, label representation, per-group-flag vs. W2 `MIGRATE_GROUP`) remain in the Step-1 spec and are **not** re-litigated here.

---

## 4. Step 4 — Remaining NSM2 Features (New)

Each subsection gives: NSM2/W2 reference, new state variables, river-valid kinetics, parameters, the sediment-flux dependency, and the linkage to Steps 1–3. All kinetics use the v3 Forward-Euler integrator and Arrhenius temperature correction `θ^(T−20)` consistent with NSM1 v3.

### 4.A Multi-pool organic matter

**Reference:** NSM2 Fortran organic-matter modules; CE-QUAL-ICM / WASP lineage. Replaces NSM1's lumped `OrgN`, `OrgP`, `POC`, `DOC`, single `POM`.

**New state variables (per cell):**

| Pool | Symbol | Replaces |
|---|---|---|
| Refractory particulate organic N | `RPON` | part of NSM1 `OrgN` |
| Labile particulate organic N | `LPON` | part of NSM1 `OrgN` |
| Dissolved organic N | `DON` | part of NSM1 `OrgN` |
| Refractory/labile particulate organic P | `RPOP`, `LPOP` | part of NSM1 `OrgP` |
| Dissolved organic P | `DOP` | part of NSM1 `OrgP` |
| Refractory/labile particulate organic C | `RPOC`, `LPOC` | NSM1 `POC` |
| Refractory/labile dissolved organic C | `RDOC`, `LDOC` | NSM1 `DOC` |

`POM` (particulate organic matter mass) is retained as the mass-balance bookkeeping pool.

**River-valid kinetics** (first-order, T-corrected):

- **Hydrolysis (particulate → dissolved):** `RPON → DON`, `LPON → DON` (labile faster); analogous for P and C. Rate `k_hyd · θ^(T−20) · pool`.
- **Mineralization (dissolved → inorganic):** `DON → NH4`, `DOP → DIP`, `LDOC/RDOC → DIC` with O₂ consumption (couples to DOX). Optionally DO-limited via a Monod `DO/(KsOxmn+DO)` switch (NSM2 option).
- **Settling/deposition:** all particulate pools settle with per-group/per-pool `vs` as a **per-cell downward flux** (Section 2.1). In the depth-averaged river this is a bed-deposition sink; the flux value is exactly the input the Step-3 diagenesis G-classes will later consume (Section 4.E).
- **Algal/benthic mortality routing:** the Step-1 mortality partition is repartitioned by per-group stoichiometric fractions (`f_rpon`, `f_lpon`, `f_don`, and P/C analogues) into the new pools, replacing the lumped `.sum("algae_group")` targets. Fractions are per-algae-group (cyano vs. diatom detritus differ in lability).

**Conservation requirement:** total N, P, C conserved across algae → multi-pool OM → DIN/DIP/DIC. Tier-1 closed-system test (Section 7).

**River/reservoir note:** purely concentration-and-temperature driven; identical in a river reach or a future reservoir layer. No stratification dependence.

### 4.B Alkalinity / pH carbonate system

**Reference:** CE-QUAL-W2 carbonate routine; Chapra (1997) Ch. 22; Stumm & Morgan. NSM1 v3 already carries `Alk` as a **tracer** with source/sink terms (nitrification consumption, denitrification production, algal growth/respiration coupling — NSM1 design spec line 645). Step 4 **adds the equilibrium solver on top of the existing source/sink terms; it does not replace them.**

**New diagnostic outputs (not transported state):** `pH`, `[CO2*]` (free CO₂ + H₂CO₃), `[HCO3⁻]`, `[CO3²⁻]`, un-ionized ammonia fraction `f_NH3` and `NH3` concentration, `pCO2`.

**Solver (cell-local, vertically-agnostic):**

1. Inputs: `Alk`, `DIC` (from the Carbon process), `T`, **ionic strength `I`** (from salinity / specific conductance — kept as an explicit input for reservoir/estuary readiness, Section 2.2). Freshwater rivers: low `I`; use freshwater-appropriate apparent constants (Plummer & Busenberg / Millero) for `K1(T,I)`, `K2(T,I)`, `Kw(T,I)`, Henry's `KH(T,I)`.
2. Solve carbonate alkalinity balance `Alk ≈ [HCO3⁻] + 2[CO3²⁻] + [OH⁻] − [H⁺]` for `[H⁺]` (bracketed Newton; robust over pH 4–11).
3. Speciate: `[CO2*] = α0·DIC`, `[HCO3⁻] = α1·DIC`, `[CO3²⁻] = α2·DIC`.
4. **CO₂ atmospheric exchange:** `k_CO2 = k_O2 · (Sc_CO2/Sc_O2)^(−1/2)` — reuses the existing NSM1 v3 reaeration menu (river-appropriate piston velocity); flux `k_CO2·([CO2*]_sat − [CO2*])` updates DIC.
5. **Un-ionized ammonia:** `f_NH3 = 1/(1 + 10^(pKa(T) − pH))`; outputs `NH3` for toxicity reporting (and optional NH₃ volatilization, default off in rivers).
6. **CO₂ growth limitation:** `f_CO2 = [CO2*] / (KsCO2 + [CO2*])`, wired into the **Step-1 multi-group growth-limitation** chain (multiplicative/min/harmonic, per `growth_rate_option`). This closes Step-1 non-goal #7 ("pH-driven CO₂ limitation … deferred to v3.x when alkalinity/pH lands") — the more defensible bloom-collapse mechanism.

**River/reservoir note:** the solver is a per-cell thermodynamic calculation — no vertical dependence. Keeping `I` as an input is the single reservoir-readiness hook; everything else is freshwater-river default.

### 4.C Methane and sulfide

**Reference:** NSM2 Fortran CH₄/H₂S; CE-QUAL-W2; Di Toro (sediment source). 

**New state variables (per cell):** `CH4`, `H2S`. (`SO4` optional; in rivers sulfate is typically abundant and parameterized rather than tracked — **open question Q4**.)

**River-valid kinetics:**

- **Sources:** sediment release `J_CH4`, `J_H2S`. **Parameterized now** (`CH4_sed_release`, `H2S_sed_release` — scalar global or per-cell, exactly matching the NSM1 v3 parameterized-bed-flux pattern, design spec line 646); upgraded to mechanistic Di Toro flux when Step 3 lands (Section 4.E) — *no water-column state or kinetic change at that swap.* Optional in-situ anaerobic mineralization of `LDOC` under low DO (simple electron-acceptor switch; default off for rivers — **open question Q4**).
- **Oxidation (O₂ sink, DIC source):** `CH4 + 2 O2 → CO2 + 2 H2O`; `H2S + 2 O2 → SO4²⁻ + 2 H⁺`. First-order in substrate, DO-limited (`DO/(Ks+DO)`). O₂ demand routes into the **DOX** process; CH₄ oxidation carbon routes into **DIC** (and thus 4.B).
- **Degassing:** CH₄ ebullition/volatilization and H₂S volatilization as `k_gas·(C − C_sat)` with river piston velocity (Schmidt scaling from reaeration, as in 4.B step 4).

**River/reservoir note (important):** the *kinetics* above are identical in a well-mixed river reach and a future anoxic hypolimnion cell — they are reservoir-ready by construction. What is **excluded** is any logic that *accumulates* CH₄/H₂S because of stratification (oxycline-gated hypolimnetic build-up). That accumulation is an emergent transport outcome in a layered driver, not a kinetic; this module never encodes it. In a depth-averaged river, CH₄/H₂S arise from the sediment source + low-DO reaches and are consumed by oxidation/degassing — the physically correct river behavior.

### 4.D Silica

**Reference:** NSM2 Fortran silica; CE-QUAL-W2 `ASI(JA)` / `AHSSI(JA)`. **Fills the Step-1 reserved diatom hooks.**

**New state variables (per cell):** `DSi` (dissolved available silica), `PBSi` (particulate biogenic silica / opal).

**River-valid kinetics:**

- **Diatom uptake:** `−rSi · μ · A` for Si-limited groups, where `rSi` ← Step-1 `AWsi[g]` (Si:Chl-a or Si:C stoichiometry). Removes `DSi`.
- **Si growth limitation:** for groups with Step-1 `si_limitation_option[g] = 2`, `f_Si = DSi/(KsSi[g] + DSi)` enters the multi-group growth-limitation chain alongside N, P, light, and (4.B) CO₂. Groups with `si_limitation_option=1` are Si-unlimited (current behavior) — backward compatible.
- **Return:** diatom mortality/respiration → `PBSi` (and a small direct `DSi` fraction).
- **Dissolution:** `PBSi → DSi` at `k_diss·θ^(T−20)·PBSi`.
- **Settling:** `PBSi` settles with `vs_Si` as a per-cell flux (Section 2.1); bed sink now, Step-3 Si return later.
- **Sediment Si flux:** parameterized `Si_sed_release` now; mechanistic via Step 3 later (Section 4.E).

**Linkage:** this is the feature that closes the Step-1↔Step-4 loop — Step 1 reserved `AWsi`, `KsSi`, `si_limitation_option` per group precisely so this lands with no YAML migration. River-relevant (diatom/cyano succession via Si:N:P); no stratification dependence.

### 4.E Step-3 sediment diagenesis — design directive and the Step-4 interface boundary

#### Step-4 interface boundary (the dependency contract)

Step 4's particulate-OM deposition, CH₄/H₂S source, and Si return all interact with the bed. To keep Step 4 **independent of Step 3**, this plan fixes a stable flux interface:

- **Now (Step 4 ships on this):** sediment fluxes are **parameterized** — scalar-global or per-cell input values (`SOD_20`, `NH4fromBed`, `DIPfromBed`, `NO3_BedDenit`, `DIC_sed_release`, plus new `CH4_sed_release`, `H2S_sed_release`, `Si_sed_release`). This is exactly the NSM1 v3 1.0.0 pattern (design spec lines 79, 646) — no new mechanism.
- **Later (Step 3):** the Di Toro multi-G model computes those same flux symbols from the deposited particulate-OM G-classes. **The water-column state variables and kinetics in 4.A–4.D do not change** when this swap happens — only the *provenance* of the flux values. The deposition flux defined per pool in 4.A is precisely the reactivity-class input Step 3 consumes.
- **`kdpo4` / TIP sorption** (gold-standard §8: "needs multi-class solids + NSM2 sediment-flux coupling"): inorganic-P sorption to suspended solids is river-relevant (turbid rivers) but needs a solids/TSS state variable. **Scoped out of this plan** as a separate solids-coupled item; flagged so it is not silently bundled — **open question Q5**.

#### Step-3 design directive (author-approved 2026-05-16)

From the W2 source review (`Diagenesis Sediment Flux Model 05.f90` and companion modules) and the NSM2 validation record. Sediment diagenesis across the lineage is **three tiers**, not "NSM2 vs W2":

1. **NSM1** (v3 1.0.0 now): prescribed constant bed fluxes. No diagenesis.
2. **NSM2**: the canonical **Di Toro / Cerco / SedFlux** two-layer, three-reactivity-class model. Not crude — intercompared against SedFlux and CE-QUAL-ICM under Chesapeake-Bay forcing in **ERDC/EL TR-16-11** (Johnson & Zhang 2016), with good agreement. TR-16-11 is the authoritative NSM2 sediment-diagenesis verification document; it is a code intercomparison, not a field validation of sediment state variables.
3. **CE-QUAL-W2 v2026.02**: the same Di Toro core **plus added mechanisms** — methane/H₂S/NH₃/CO₂ ebullition (`Diagenesis Bubbles Code 01.f90`), a fresh-deposit transient surface layer (`Diagenesis FFT Layer 01.f90`; "FFT" here means fluff, not Fourier), bed consolidation/erosion, optional iron/manganese, dynamic active-layer depth.

The limiting assumptions are mostly **shared by both NSM2 and W2**: two-layer collapse; quasi-steady pore water (W2 relaxes over elapsed simulation time, not the coupling step); constant reactivity-class stoichiometry; linear instantaneous sorption; an either/or methane-versus-sulfide switch on a single sulfate threshold; constant burial velocity; an empirical oxygen-dependent phosphate trap standing in for iron-redox. Two capabilities are **missing from both**: benthic-stress / particle-mixing hysteresis (the dominant control on nutrient recycling after repeated low-oxygen events — the largest gap for eutrophic stratified reservoirs) and silica diagenesis (no silica flux anywhere in W2; this is the bed side of the Step-4 silica cycle). Copying W2 therefore does not avoid the design dead-end.

**Recommended Step-3 build approach:**

1. **Core model:** a ClearWater-native Di Toro two-layer, three-reactivity-class model — not a port of NSM2 or W2, not a new formulation, and not a depth-resolved sediment column (the last is research scope, unnecessary for rivers). It is the accepted standard and the right level of detail for vertically-averaged rivers.
2. **The flux interface is the central contract.** The sediment process delivers the fixed flux set (ammonia, nitrate, methane, sulfide/sulfate, phosphate, inorganic carbon, silica, SOD) and consumes the settling organic-matter flux split by reactivity class plus overlying-water concentrations and temperature. Any provider meeting this interface is interchangeable: prescribed constants (now), time series or spatial maps, the Di Toro model, or an extended version. Water-column kinetics never change when the provider is swapped.
3. **Named extension points, not a single block of code.** Write the state schema and per-step solve so each of these is a localized later addition rather than a rewrite — none built initially: a benthic-stress state and its effect on particle mixing; a biogenic-silica reactivity pool and silica flux (designed jointly with the Step-4 silica feature so Step 4 does not bake in a silica dead-end); a methane gas store and bubble-release pathway; a fresh-deposit transient surface layer and bed consolidation/erosion coupling; an iron/pH-coupled phosphate path to replace the simple oxygen-dependent trap.
4. **True transient integration via the operator split.** Both legacy models relax toward quasi-steady; W2 integrates over elapsed simulation time, a known weakness. ClearWater already separates reaction from transport with a controlled sub-step — advance the slow diagenesis pools as genuine differential equations over that sub-step while keeping the fast pore-water partition implicit as Di Toro does. This removes the one limiting assumption neither legacy model fixes and is what makes the implementation defensible as more than a re-port.
5. **Fail gracefully.** Keep the iterative SOD solve as Di Toro specifies, but do not copy W2's hard stop on a singular solver matrix.
6. **Validate in three layers:** reproduce the TR-16-11 intercomparison (against SedFlux and CE-QUAL-ICM) as the primary regression anchor; cross-check against W2 v2026.02 on the carbon/nitrogen/phosphorus/sulfur fluxes, knowing W2's documented deviations so differences are explained, not chased (W2 is a cross-check, not the oracle); add the numerical-correctness tests the published record lacks — the analytical sediment-methane-at-saturation case, a manufactured-solution test for the transient solver, and closed-system mass conservation from deposition through diagenesis to flux.
7. **Supersede, do not patch, the existing draft.** A v3/v2 "Draft sediment flux implementation" (commit `303d285`, 2025-07-28) carries a transcription error (`0.0061` versus the Fortran `0.0432`; `clearwater_modules_v3_review_findings.md:94`). Treat it as a non-authoritative sketch and rebuild under the interface so that error class does not survive.

**Defects not to inherit from W2:** `SD_KappaH2Sp1` declared but never assigned (zeroes marine particulate-sulfide oxidation); `Lin_Sys` hard-stops on a singular matrix; an older-W2 `CSODmax` exponent error (current W2 is fixed — relevant only if reading older W2); iron/manganese matrix divide-by-near-zero (treat the W2 Fe/Mn path as the least mature part).

**Two judgment calls for the author** — proposed defaults in **Q7** (iron/manganese/pH-coupled phosphorus scope) and **Q8** (methane bubble-release ordering).

**Off the critical path.** None of the above is built for Step 4 or v3 1.0. Step 4 ships on the parameterized interface above; the Di Toro model lands as Step 3 after the §6 cutover gate (Section 1.1) and replaces the flux provider with no change to water-column code. A compact Step-3 phased sketch is in Section 6.

---

## 5. Process Ordering

The NSM1 v3 `ComputeKinetics` order must be updated for NSM2 (NSM1 design spec line 157: "SedFlux moves earlier; methane-sulfide goes between carbon and DOX; alkalinity moves to the end"). Because v3 declares process order in YAML, this is a configuration change, not a code change. Step-4 order of record:

```
SedFlux (parameterized or Step-3)  →  Algae (multi-group, Step 1)  →  Benthic algae
  →  Multi-pool OM (4.A)  →  Nitrogen  →  Phosphorus  →  Silica (4.D)
  →  Carbon/DIC  →  Methane–Sulfide (4.C)  →  DOX  →  Alkalinity/pH solver (4.B, last)
```

Rationale: SedFlux first so bed fluxes are available to all consumers; Si after N/P (parallel nutrient); CH₄/H₂S between Carbon and DOX so their O₂ demand and DIC source are correctly sequenced; the carbonate/pH solver **last** so it sees the fully updated `Alk` and `DIC` for the timestep.

---

## 6. Merged Phased Implementation Plan

Step 1 phases are carried from the multi-group spec unchanged (summarized); Step 4 phases are new. Estimates assume Claude does the coding, in the multi-group spec's accounting convention.

### Step 1 (Phases 0–7) — Multi-group algae

Per `clearwater_modules_v3_multigroup_algae_design_specification.md` §6: Phase 0 gap analysis → Phase 1 parameter/config plumbing → Phase 2 FloatingAlgae per-group kinetics → Phase 3 downstream-consumer refactor → Phase 4 BenthicAlgae mirror → Phase 5 boundary conditions/output → Phase 6 tests/validation → Phase 7 docs/review. **`nal=1` parity is the regression contract.** *~9–12 working days.*

### Step 4 (Phases S4-0 … S4-6) — Remaining NSM2 features

**Prerequisite:** Step 1 merged (Step 4 modifies the Step-1 consumer-collapse seam).

#### Phase S4-0 — State-variable migration plan & gap analysis (1 day)
Catalog every site reading the lumped `OrgN/OrgP/POC/DOC/POM` and the Step-1 mortality-partition seam. Define the backward-compatible migration: a `multipool_om: false` config flag keeps NSM1 lumped behavior (single OrgN/OrgP/POC/DOC) so existing NSM1 applications are unaffected; `true` enables the seven-pool set. **Deliverable:** `design/nsm2_multipool_om_gap_analysis.md` + migration/back-compat decision.

#### Phase S4-1 — Multi-pool organic matter (3–4 days)
New pools, hydrolysis/mineralization/settling kinetics, per-algae-group mortality repartition replacing the Step-1 `.sum("algae_group")` lumped targets. Wire DOC oxidation O₂ demand into DOX, DIC source into Carbon. **Deliverable:** seven-pool OM with N/P/C closed-system conservation passing; `multipool_om:false` reproduces NSM1 lumped trajectories.

#### Phase S4-2 — Silica (1–2 days)
`DSi`, `PBSi` state; diatom uptake/return/dissolution/settling; wire `f_Si` into the Step-1 multi-group growth-limitation chain via the reserved `AWsi`/`KsSi`/`si_limitation_option` hooks. **Deliverable:** Si cycle with diatom-group limitation; `si_limitation_option=1` reproduces Step-1 behavior exactly.

#### Phase S4-3 — Alkalinity/pH carbonate solver (3–4 days)
Equilibrium solver on top of the existing NSM1 `Alk` source/sink terms; `pH`, `[CO2*]`, speciation, `f_NH3`, CO₂ atmospheric exchange via the reaeration menu; `f_CO2` into the Step-1 growth-limitation chain. Ionic-strength input plumbed (reservoir-ready). **Deliverable:** carbonate solver validated against analytical equilibria and W2 cross-check; CO₂ limitation toggles cleanly.

#### Phase S4-4 — Methane & sulfide (2–3 days)
`CH4`, `H2S` state; parameterized sediment source (Section 4.E interface); DO-limited oxidation routed to DOX and DIC; degassing. **Deliverable:** CH₄/H₂S kinetics with O₂/DIC coupling; reservoir-ready kinetics, no stratification logic.

#### Phase S4-5 — Integration, reordering, conservation & regression (2–3 days)
Apply the Section-5 process order in YAML; full-system mass-conservation (N, P, C, Si) and an MMS test (closes the gold-standard MMS item for the NSM2 set); **Fortran NSM2 → v3 regression** on a representative river reach. **Deliverable:** integrated NSM2 feature set; conservation + MMS + Fortran-regression green.

#### Phase S4-6 — Documentation, migration, review prep (1 day)
YAML schema for new state variables and `multipool_om` flag; migration notes (NSM1 lumped → NSM2 multi-pool, including how to split lumped OrgN into RPON/LPON/DON); the post-hoc-pH worked example promoted to a live solver (supersedes NSM1 design spec line 645 stopgap); Appendix-A mapping extended (Section 8 here). **Deliverable:** Step-4 NSM2 features ready for LimnoTech review and a v3.x release.

**Step-4 estimated wall-clock: 13–19 working days.** Combined Step 1 + Step 4: **~22–31 working days**, sequenced (Step 4 after Step 1 merge).

### Step 3 (sediment diagenesis) — off critical path, phased sketch

Not on this plan's critical path; the design directive is in Section 4.E. It can proceed in parallel with or after Step 4, and when it lands it swaps the parameterized flux provider for the Di Toro model with no 4.A–4.D code change. Indicative phases:

- **S3-0 — Interface & schema (1–2 days):** lock the flux interface from Section 4.E; define the state schema with the named extension points reserved (benthic-stress, biogenic-Si pool, CH₄ gas store, transient-vs-fixed active layer).
- **S3-1 — Di Toro core (5–7 days):** two-layer, three-reactivity-class diagenesis; transient pool integration over the operator-split sub-step; iterative SOD solve with graceful failure; supersede the `303d285` draft.
- **S3-2 — Provider swap (1 day):** replace the parameterized flux provider with the Di Toro provider; confirm 4.A–4.D water-column code is unchanged.
- **S3-3 — Validation (3–4 days):** TR-16-11 intercomparison (SedFlux / CE-QUAL-ICM) as the regression anchor; W2 cross-check; analytical methane-at-saturation, manufactured-solution, and closed-system conservation tests.
- **S3-4 — Extension points landed as needed (scoped per application):** benthic-stress hysteresis and silica diagenesis are the priority extensions for eutrophic-reservoir and diatom applications; ebullition and the fresh-deposit layer follow per Q8.

**Step-3 indicative wall-clock for the core path (S3-0…S3-3): ~10–14 working days**, excluding the optional extensions.

---

## 7. Testing and Validation

Extends the Step-1 validation tiers (multi-group spec §4.2). Per project V&V convention, **NSM2 Fortran is the regression baseline**; W2/QUAL2K/WASP/Chapra are cross-checks.

- **Tier 1 — Conservation:** closed-system N, P, C, Si mass balance across algae → multi-pool OM → DIN/DIP/DIC/DSi. Zero leakage to float tolerance.
- **Tier 2 — Backward-compat parity:** `multipool_om=false`, `si_limitation_option=1`, carbonate solver off → byte/float-equivalent to current NSM1 v3. The non-regression contract for existing applications.
- **Tier 3 — Analytical:** carbonate solver vs. closed-form equilibria across pH 4–11 and a temperature/ionic-strength grid; un-ionized ammonia fraction vs. tabulated `pKa(T)`.
- **Tier 4 — MMS:** method-of-manufactured-solutions for the coupled OM/Si/carbonate/CH₄–H₂S kinetics (closes the gold-standard MMS item for the NSM2 set).
- **Tier 5 — Fortran regression:** NSM2 Fortran → v3 on a representative depth-averaged river reach (the authoritative baseline).
- **Tier 6 — Cross-checks (idea source, not gate):** carbonate solver vs. CE-QUAL-W2 and Chapra worked examples; silica vs. W2 `ASI`/`AHSSI`; multi-pool OM vs. WASP/CE-QUAL-ICM lineage. Divergence is investigated and explained, not auto-failed.

---

## 8. Risks and Mitigations

| Risk | Prob. | Impact | Mitigation |
|---|---|---|---|
| Carbonate solver non-convergence at pH extremes / very low ionic strength | Medium | Medium | Bracketed Newton with pH 4–11 bounds + fallback bisection; Tier-3 analytical grid includes low-`I` freshwater corners |
| Multi-pool OM migration breaks existing NSM1 applications | Medium | High | `multipool_om=false` default = exact NSM1 lumped behavior; Tier-2 parity is a merge gate |
| Step-4 work blocks on Step 3 (sediment diagenesis) | Medium | High | Section 4.E interface: parameterized flux now (NSM1 pattern), mechanistic later with no water-column change — Step 4 never waits on Step 3 |
| Scope creep: stratified-reservoir mechanisms pulled into river kinetics | Medium | Medium | Section 2 scoping table is a binding review checklist; buoyancy-migration explicitly excluded (Q1); reviewer rejects any "layer-aware" kinetic |
| CH₄/H₂S O₂-demand double-counts with CBOD/SOD | Medium | Medium | Explicit O₂-budget accounting test in Tier 1; CH₄/H₂S oxidation is a *distinct* O₂ sink from CBOD and parameterized SOD — verified non-overlapping |
| `f_CO2`/`f_Si` change calibrated NSM1 algae behavior | Medium | Medium | Both default off (`si_limitation_option=1`, CO₂ limitation opt-in); Tier-2 parity proves no change when off |
| Fortran NSM2 reference ambiguity (no v1 NSM2, stubs only) | Medium | Medium | Treat Fortran NSM2 as sole baseline; where Fortran is itself ambiguous, document against W2/Chapra and flag for author decision (cross-check ≠ conformance) |
| Step-3 built as a fixed simplified model, locking out reservoir capability later | Medium | High | Section 4.E directive: Di Toro core behind the stable flux interface; named extension points for benthic-stress hysteresis, silica diagenesis, ebullition, fresh-deposit layer, transient solver; TR-16-11 regression anchor |
| Existing v3/v2 draft sediment-flux code (`303d285`) carries a transcription error | High | Medium | Supersede under the new interface, do not patch (Section 4.E item 7); the draft is non-authoritative |
| Plan merged before the §6 cutover gate corrupts the v3-vs-v1 review and the bit-identical baseline | Medium | High | Section 1.1 binding merge gate: unmerged branch off `streaming`, no merge until after §6 step 3, bit-identical `nal=1` with no re-baseline; re-baseline only as documented last resort |

---

## 9. Open Questions (Proposed Defaults — Confirm or Override)

Following the Step-1 spec's "proposed default + please confirm" house style.

1. **Buoyant-cyanobacteria light position (Section 2.3).** Split buoyancy into (a) river-admissible vertical-light-position sub-model and (b) excluded vertical migration; deliver (a) as an optional per-group light sub-model in **Step 2**, default off. **Proposed default: yes, split; (a) → Step 2, (b) excluded.** Confirms the river-scope answer and prevents the Step-1 blanket deferral from burying (a).
2. **Multi-pool OM back-compat switch.** `multipool_om` config flag, default `false` (NSM1 lumped behavior preserved). **Proposed default: yes, default false.** Existing NSM1 applications must be untouched until a user opts in.
3. **CO₂ growth limitation default.** Off by default (opt-in via a per-run flag), on when alkalinity/pH solver is enabled and the user requests it. **Proposed default: solver computes `[CO2*]`/`pH` always when enabled; `f_CO2` enters growth only when explicitly opted in.** Avoids silently changing calibrated algae kinetics.
4. **Sulfate (`SO4`) and in-situ anaerobic mineralization.** Track `SO4` as state and model in-situ anaerobic `LDOC` mineralization, or parameterize sediment CH₄/H₂S only? **Proposed default: parameterized sediment source only; no `SO4` state, no in-situ anaerobic pathway in Step 4** (rivers are typically sulfate-replete and oxic; the anaerobic cascade is a Step-3/diagenesis concern). Revisit if a river application needs it.
5. **`kdpo4` / TIP sorption scope.** Inorganic-P sorption to suspended solids needs a solids/TSS state. **Proposed default: out of scope for this plan**, tracked as a separate solids-coupled item (gold-standard §8). Flag, don't silently bundle.
6. **State-variable naming.** Adopt NSM2 Fortran names (`RPON`, `LPON`, `DON`, …, `DSi`, `PBSi`) or W2 names? **Proposed default: NSM2 Fortran names** (this is an NSM2 port; matches the regression baseline), with a W2-name mapping table (Section 8 / Appendix). Consistent with Step-1's choice to mirror NSM2 structure.
7. **Step-3 iron/manganese/pH-coupled phosphorus scope.** Ship the simple oxygen-dependent phosphate trap (sufficient for rivers; what W2 mostly relies on), or build the mechanistic iron path (matters for seasonally anoxic reservoirs)? **Proposed default: ship the simple trap; reserve the interface and an extension point for the mechanistic iron path.** Confirm the scoping.
8. **Step-3 methane bubble-release ordering.** Bubble release is the dominant carbon export in organic-rich shallow sediments; the full model (fracture-mechanics release plus turbulence feedback, as in W2) is a large addition. **Proposed default: ship aqueous methane flux first, with the gas store and ebullition pathway as a reserved extension — not in the first Step-3 delivery.** Confirm whether ebullition must instead be in the first Step-3 delivery.

---

## 10. Appendix A — State-Variable & Parameter Mapping (extends Step-1 Appendix A)

| v3 (Step 4) | NSM2 (Fortran) | CE-QUAL-W2 | Notes |
|---|---|---|---|
| `RPON`,`LPON`,`DON` | RPON, LPON, DON | LPOM-N / RPOM-N / DOM-N analogues | replaces lumped NSM1 `OrgN` |
| `RPOP`,`LPOP`,`DOP` | RPOP, LPOP, DOP | — | replaces lumped NSM1 `OrgP` |
| `RPOC`,`LPOC`,`RDOC`,`LDOC` | RPOC, LPOC, RDOC, LDOC | LDOM/RDOM-C | replaces NSM1 `POC`/`DOC` |
| `DSi`,`PBSi` | available Si, biogenic Si | `DSI`, `PSI` | wires Step-1 `AWsi`/`KsSi`/`si_limitation_option` |
| `CH4`,`H2S` | CH4, H2S | `CH4`, `H2S` | parameterized sediment source (4.E) |
| `pH`,`[CO2*]`,`f_NH3` | derived | W2 carbonate routine | diagnostic, not transported |
| `multipool_om` (flag) | (implicit; always multi-pool) | (implicit) | back-compat switch; default false |
| `CH4_sed_release`,`H2S_sed_release`,`Si_sed_release` | Di Toro `JCH4`/`JH2S`/`JSi` | sediment compartment | parameterized now; Step-3 mechanistic later |
| `KsCO2`,`KsSi[g]`,`AWsi[g]` | `KsCO2`, `KsSip(i,r)`, Si stoich | `AHSSI(JA)`, `ASI(JA)` | growth-limitation half-sats / stoichiometry |
| `J_NH4`,`J_NO3`,`J_CH4`,`J_SO4`,`J_H2S`,`J_PO4`,`J_DIC`,`J_Si`,`SOD` | Di Toro flux symbols | W2 `SedimentFlux` outputs | Step-3↔Step-4 stable flux interface; provider-swappable (4.E) |
| benthic-stress state; biogenic-Si reactivity pool; CH₄ gas store | (absent in NSM2) | (W2: `SD_BEN_STR*` declared-unused; no silica flux; `Bubbles`/`FFT` modules) | Step-3 extension points — reserved, not built in v3 1.x (4.E) |
| TR-16-11 intercomparison | NSM2 ≈ SedFlux ≈ CE-QUAL-ICM | — | authoritative Step-3 regression anchor (intercomparison, not field validation) |

---

## 11. Approval Criteria

Approved when the author has reviewed and accepted:

1. River-first scope and the reaction-module/transport split (Section 2.1–2.2) as binding.
2. The Section 1.1 merge gate and cutover sequencing as binding.
3. The buoyancy split and its Step-2 placement (Section 2.3 / Q1).
4. Step-4 feature designs 4.A–4.D and the Step-3↔Step-4 interface boundary (Section 4.E).
5. The Step-3 sediment diagenesis design directive (Section 4.E) — author-approved 2026-05-16; reflected here for the record.
6. The Section-5 process reordering.
7. The merged phased plan and estimates, including the Step-3 phased sketch (Section 6).
8. The validation tiers with NSM2 Fortran as the regression baseline (Section 7).
9. Risks (Section 8) and resolutions to the eight open questions (Section 9).

Once approved, Step 4 implementation begins at Phase S4-0, after Step 1 is merged.
