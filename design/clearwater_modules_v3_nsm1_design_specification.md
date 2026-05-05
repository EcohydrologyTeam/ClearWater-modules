# ClearWater Modules v3 — NSM1 Design Specification

**Status:** Draft for review
**Author:** Todd Steissberg (ERDC), with Claude
**Date:** 2026-05-04
**Read this with the umbrella spec.** Quick start, env setup, branch conventions, package architecture, integrator-pattern contract, retirement plan, and umbrella risks are documented in `clearwater_modules_v3_architecture_specification.md`. TSM-specific design lives in `clearwater_modules_v3_tsm_design_specification.md`. This document covers NSM1-specific design only.

**Scope:** NSM1 module within `clearwater_modules_v3`.

---

## 1. Background and Motivation

NSM1 (Nutrient Simulation Module, version 1) is a coupled water-quality kinetics package covering 16 constituents: phytoplankton, benthic algae, ammonium, nitrate, organic nitrogen, dissolved nitrogen gas, total inorganic phosphorus, organic phosphorus, particulate organic carbon, dissolved organic carbon, dissolved inorganic carbon, particulate organic matter, multi-group CBOD, dissolved oxygen, pathogen, and alkalinity.

The current state of NSM1 across the two parallel codebases is markedly asymmetric:

### NSM1 v1 (`clearwater_modules`)

- **All 16 constituents implemented** with working kinetics (~3,540 lines in `processes.py`, ~290 process functions)
- **~430 passing regression tests** with expected-value tables in `tests/NSM Manual Calcs/`
- Reaeration menu (9 hydraulic options + 13 wind options) in `shared/processes.py`
- Beer-Lambert light extinction (with ISS, POC, Chl-a contributions) and PAR computation
- Parameter `TypedDict` library (~250 default entries across 13 groups)
- Shared physics primitives library (`shared/processes.py`)
- Optional numba `@njit` decorators on inner kernels
- Frozen since 2024-09-24 except for hotstart kwargs (Todd, 2026-05-01)
- **Known bug:** Sentinel-`999` defaults in `nsm1/constants.py` for `vsop`, `vs`, `SOD_20`, `SOD_theta`, `kaw_20_user`, `kah_20_user`. `SOD_theta=999` produces catastrophic Arrhenius blowup (`999^(T-20)`) at T > 20 °C. Documented in `docs/NSM1_DOX_rate_bug_investigation.md` (commit `b5f2c8b`, 2026-02-20).

### NSM1 v2 (`clearwater_modules_v2`)

- **4 of 16 constituents partially implemented** (NH4, NO3 in `nitrogen.py`; Ap in `floating_algae.py`; Ab in `benthic_algae.py`)
- **12 constituents absent:** OrgN, TIP-as-process, OrgP, POC, DOC, DIC, POM, CBOD, DOX, PX, Alk, N2
- Built on the v2 framework patterns established by TSM v2: `Process` subclass, `ProcessFactory.register`, per-process `time_step`, `registry.get_at_time`/`set_at_time`
- First v2 NSM1 commit: 2025-11-11 by Paul Tomasula (4 months after v2 framework scaffolding)
- Last v2 NSM1 commit on `memory-refactor-pytestUpdate`: 2026-01-29 by Paul Tomasula ("Single Registry Instance"); v2 NSM1 development paused
- **Critical bugs in the 4 implemented constituents:**
  1. Multiplicative integrator at `nitrogen.py:101` and `nitrogen.py:115`: `ammonium = 0 + ammonium * rate * dt` instead of `ammonium + rate * dt`
  2. Multiplicative integrator + stray `* 86400` at `floating_algae.py:122`
  3. NaN guards at four locations using `rate == np.nan` (always False per IEEE 754, so the guard never fires)
  4. Hard-coded `half_saturation_oxygen=1` at `nitrogen.py:197` with `# TODO: need argument`
  5. Hard-coded `algea_growth_rate=0` at `nitrogen.py:211, 217` (silently disables algal nitrate uptake)
  6. Death rate placeholder at `nitrogen.py:60` with `# TODO: this should come from floating algae process`
  7. `floating_algae.ammonium_respiration()` returns 0 with `# TODO: implement`
  8. `floating_algae.ammonium_growth()` returns 0 with `# TODO: implement`
  9. Hard-coded `phosphate_fraction_dissolved=0.5` at `floating_algae.py:113` with `# TODO`
- Author has fixed 4 of these on the `streaming` branch (`set_at_time` persistence, `time_step` typo, 2 of 4 NaN guards) but those fixes are not on `memory-refactor-pytestUpdate`

### What v3 NSM1 must accomplish

Unlike TSM v3, which was a merge of two divergent full implementations, NSM1 v3 has a fundamentally different shape:

1. **Foundation work:** complete and correct the 4 partial v2 NSM1 constituents (fix the 9+ known bugs, wire up the parameters that were left as `# TODO` placeholders)
2. **Major port:** translate the 12 missing constituents from v1's function-style processes into v3's `Process` class pattern, preserving the working kinetics
3. **Shared infrastructure:** port v1's reaeration menu, light extinction, fdp partitioning, and other shared physics primitives into v3
4. **Parameter library:** establish a v3 pattern for the ~250 NSM1 parameters that v1 organizes via `TypedDicts`
5. **Bug correction during port:** replace v1's sentinel-`999` defaults with sensible values; do not carry the bug forward
6. **Orchestration improvements:** v3 NSM1 inherits the v3 `Model` (kernel optimization, wet-mask, hotstart) being designed for TSM v3

---

## 2. Goals and Non-Goals

### Goals

1. All 16 NSM1 constituents implemented in v3 framework (v2 `Process` class pattern), runnable end-to-end coupled with Riverine on a representative test case (Sumwere Creek with nutrients, or equivalent).
2. The 4 partial v2 NSM1 implementations completed and corrected: integrator bug fixed, parameters wired, algal-uptake stubs implemented, NaN guards corrected.
3. The 12 missing constituents ported from v1 with kinetics-formula equivalence (numerical outputs match v1 within floating-point tolerance for the regression test suite).
4. Sentinel-`999` defaults from v1's `nsm1/constants.py` corrected during the port; documented as fixed.
5. Reaeration menu (9 hydraulic + 13 wind options), Beer-Lambert light extinction, fdp partitioning, and other shared physics primitives ported into v3.
6. Test coverage at least matching v1's NSM1 calculation suite (~11 per-constituent test files in `tests/NSM Manual Calcs/`), adapted to v3's `Process`-based dispatch.
7. v2/v3 parity tests for the 4 constituents that exist in both, demonstrating that v3 reproduces v2's intended behavior (with v2 bugs fixed) and v1's behavior under matched conditions.
8. YAML configuration backward compatibility: existing v2 NSM1 YAML configs run on v3; new YAML keys for the 12 newly-ported constituents follow the v2 pattern.
9. NSM1 v3 inherits TSM v3's orchestration improvements (kernel optimization, wet-mask, hotstart) without redundant per-process implementation.

### Non-Goals

1. NSM2 features are explicitly out of scope. Multi-pool organic matter (RPON/LPON/DON, RPOP/LPOP/DOP, RPOC/LPOC/RDOC/LDOC), alkalinity solver with pH chemistry, methane and sulfide kinetics, and silica cycle are deferred to subsequent v3 releases (1.1, 1.2, ...).
2. Sediment diagenesis (Di Toro multi-G) is out of scope for v3 1.0.0. v3 1.0.0 uses v1's parameterized sediment fluxes (`NH4fromBed`, `DIPfromBed`, `NO3_BedDenit`, `DIC_sed_release`) consistent with v1 NSM1's actual implementation. Full sediment diagenesis comes with NSM2 features in v3 1.1+.
3. v3 NSM1 will not introduce a new architectural pattern beyond what TSM v3 establishes. The v2 `Process` class pattern is the architectural baseline.
4. v3 NSM1 will not break v2's public API for the 4 constituents that exist in both. Process names, configuration keys, and registry names remain compatible.
5. Performance optimization beyond inheriting TSM v3's `Model` improvements is out of scope for the initial release. If profiling reveals NSM1-specific bottlenecks, they are addressed in a subsequent release.

---

## 3. Architectural Approach

NSM1 v3 sits inside the v3 package alongside TSM v3, using the same thin-overlay strategy and the same v3 `Model` orchestration. This section covers only NSM1-specific architectural decisions; for general v3 architecture see the architecture spec (`clearwater_modules_v3_architecture_specification.md`).

### Process organization

NSM1 v3 follows the v2 NSM1 pattern of grouping related state variables into a single `Process` class rather than one Process per state variable. The mapping is:

| Process class | State variables managed | Source |
|---|---|---|
| `Nitrogen` | NH4, NO3, OrgN | Extends existing v2 `Nitrogen` |
| `Phosphorus` | TIP, OrgP | New in v3 |
| `Carbon` | POC, DOC, DIC | New in v3 |
| `POM` | POM | New in v3 |
| `CBOD` | CBOD (multi-group) | New in v3 |
| `DOX` | DOX | New in v3 |
| `Pathogen` | PX | New in v3 |
| `Alkalinity` | Alk | New in v3 (declared but inactive in v1) |
| `N2` | N2, TDG | New in v3 |
| `FloatingAlgae` | Ap | Extends existing v2 `FloatingAlgae` |
| `BenthicAlgae` | Ab | Extends existing v2 `BenthicAlgae` |

Eleven Process classes total. Each subclasses `clearwater_modules_v2.processes.base.Process` (re-exported from v3).

### Inter-process coupling

NSM1 has substantial inter-process coupling: algae produce O2 (consumed by DOX), nitrification consumes O2 (DOX) and produces alkalinity sink, denitrification produces N2 and consumes NO3, organic matter mineralization consumes O2 and produces DIC, and so on. The coupling pattern in v3:

1. **Earlier-running processes write rate contributions to the registry as named variables.** For example, `FloatingAlgae.run` writes `algal_growth_rate`, `algal_respiration_rate`, `algal_death_rate`, and `algal_nitrate_uptake_rate` to the registry.
2. **Later-running processes read these rates and apply them to their own state.** For example, `DOX.run` reads `algal_growth_rate`, multiplies by stoichiometric coefficients, and adds the photosynthesis O2 source to its DOX rate equation.
3. **The YAML config declares processes in the correct dispatch order.** v3's `Model` honors that order; any process declared earlier in the YAML runs earlier in each timestep.

This is the v2 pattern, generalized. The 4 existing v2 NSM1 processes already use this pattern (e.g., `Nitrogen.change_ammonium` reads `floating_algae_process.ammonium_respiration()` and `ammonium_growth()`). v3 makes the pattern uniform across all 11 processes by routing through the registry rather than direct cross-process method calls.

#### Within-timestep update semantics — Jacobi for state, Gauss-Seidel for rate variables

The Registry distinguishes two value classes with different lifecycles:

- **State variables** (Ap, Ab, NH4, NO3, OrgN, TIP, OrgP, POC, DOC, DIC, POM, CBOD groups, DOX, PX, Alk, N2): time-indexed. Within a step processing `t_current → t_current + dt`, every `Process.run` reads state via `get_at_time(t=t_current)` and writes via `set_at_time(t=t_current + dt)`. State reads always see pre-update (time-`n`) values regardless of dispatch order — Jacobi semantics. Two dispatch orderings that respect rate-variable producer→consumer edges produce identical state evolution.
- **Rate variables** (`algal_growth_rate`, `algal_respiration_rate`, `algal_death_rate`, `algal_nh4_uptake_fraction`, `nitrification_rate`, `denitrification_rate`, `doc_dic_oxidation_rate`, `cbod_oxidation_rate`, the sediment-flux variables): step-scoped, not time-indexed. Cleared at the start of each step. Producers write via `set_rate_variable(name, value)`; consumers read via `get_rate_variable(name)`. Reading a rate variable that has not been written in the current step raises an error — Gauss-Seidel semantics with strict producer-precedes-consumer enforcement, which catches dispatch-order bugs immediately.

The Model.run loop is therefore:

```
for each step (t_current → t_current + dt):
    registry.clear_rate_variables()
    for process in dispatch_order:
        process.run(t_current, dt)   # state reads at t_current; rate vars within-step
    registry.advance_time(t_current + dt)
```

NSM1 Processes never read state via a "latest available time" accessor; reads always pass `t=t_current` explicitly. The naming convention from Appendix A (`_rate` / `_fraction` suffixes) doubles as the lifecycle marker: state variables don't carry those suffixes; rate variables do.

This convention matches v1 NSM1's xarray-DAG behavior (compute all dynamic quantities in dependency order using time-`n` state, then apply state updates), which tightens v1/v3 parity tolerances. v2's accidental sequential-read pattern is closer to pure Gauss-Seidel; v2/v3 parity tests must constrain conditions where the state-read difference is moot, or accept small differences attributable to v2's looser contract.

### Process dispatch order

The NSM1 dispatch order, declared in the YAML and honored by v3 `Model`:

1. `FloatingAlgae` — computes algal growth, respiration, death, nutrient uptake rates
2. `BenthicAlgae` — same for benthic algae
3. `Nitrogen` — consumes algal N uptake; computes nitrification (consumes DOX), denitrification (consumes NO3), OrgN hydrolysis
4. `Phosphorus` — consumes algal P uptake; computes OrgP hydrolysis, TIP partitioning and settling
5. `Carbon` — POC and DOC hydrolysis/mineralization, DIC reaeration with atmosphere, photosynthesis/respiration coupling
6. `POM` — settling, dissolution
7. `CBOD` — multi-group oxidation (consumes DOX), sedimentation
8. `DOX` — reaeration, photosynthesis (from algae rates), respiration (from algae rates), nitrification (from Nitrogen rates), DOC oxidation, CBOD oxidation, SOD
9. `Pathogen` — natural and light-induced decay, settling
10. `Alkalinity` — nitrification consumption, denitrification production, algal growth/respiration coupling
11. `N2` — atmospheric exchange via Henry's law; consumes denitrification

This order matches v1 NSM1's `ComputeKinetics` sequence. NSM2 features (when added in v3 1.1+) will require reordering: SedFlux moves earlier; methane-sulfide goes between carbon and DOX; alkalinity moves to the end. The v3 framework handles ordering via the YAML declaration, so future reordering is a config change.

Note that under the Jacobi-state semantics defined above, dispatch order is required only to respect rate-variable producer→consumer dependencies — state variable evolution is order-independent. Two orderings that both respect the rate-DAG produce identical state outputs. The order above reflects the rate-DAG; alternative orderings that satisfy the DAG (e.g., swapping Pathogen and Alkalinity, which share no rate variables) are equivalent.

### Integrator pattern

The v2 NSM1 multiplicative-integrator bug is a symptom of an unclear contract for what `Process.run` should do. v3 establishes the contract explicitly:

1. Each `Process.run` reads its state variables from the registry at the current time.
2. Each `Process.run` computes net rate of change for each state variable, with units of `[state] / second` (additively combining sources and sinks).
3. Each `Process.run` applies the rate to the state via Forward Euler: `state_new = state_old + rate * self.time_step.total_seconds()`.
4. Each `Process.run` calls the shared `clip_negative_state(state_new, name, cell_mask, diagnostics)` utility (in `clearwater_modules_v3/utils/numerics.py`), which clips negatives to zero, emits a structured log entry per affected cell (with rate-limiting for high-volume events), and increments a per-Process per-variable counter on the run-level `diagnostics` object. Clipping is a safety net for off-design conditions (unphysical parameters, dt too large, bad initial conditions); under reasonable closed-system inputs it should never fire — see Section 9 Tier 1.
5. Each `Process.run` writes the (post-clip) updated state back to the registry via `set_at_time`.

Notes on the clipping contract:
- Clip target is exactly 0, not a small epsilon. Monod ratios `C/(C+K)` are well-defined at `C=0`. If a kinetics formula would divide directly by a clipped state, the formula is reformulated rather than the clip target adjusted.
- DOX is the most clip-prone constituent (multiple first-order sinks: SOD + nitrification + CBOD oxidation + DOC oxidation). Phase 5 evaluates whether DOX should opt into a semi-implicit treatment of its first-order sinks (`DOX_new = (DOX_old + sources*dt) / (1 + sum_of_first_order_sink_coefs*dt)`), which is positive by construction. This is a per-Process opt-in, not a framework-wide change.
- Adaptive substepping when a clip would fire is deferred to v3 1.1+ if profiling shows it is needed.

This is the corrected version of what v2 NSM1 attempted. The integrator is the same shape across all 11 processes; only the rate computation differs.

### Parameter library

v1 NSM1's `constants.py` provides 13 `TypedDict` groups (algae, alkalinity, balgae, carbon, CBOD, DOX, nitrogen, POM, N2, phosphorus, pathogen, global_parameters, global_vars) with ~250 entries total. v2's pattern is per-Process `__init__` arguments.

v3 NSM1 adopts a hybrid pattern:

1. Each `Process` class accepts a `parameters: dict` argument in its YAML config block.
2. The `parameters` dict is unpacked into the `Process.__init__` via `**config['parameters']`.
3. Inside the `Process`, parameters are stored as `self.<name>` attributes for use in kinetic methods.
4. Default parameter values are class-level `DEFAULTS: dict[str, float]` attributes, merged with user-provided parameters at construction time.

This preserves v2's per-Process YAML pattern while organizing parameters in groups that match v1's `TypedDict` structure for migration ease.

Example YAML block:

```yaml
- nitrogen:
    time_step: '30s'
    parameters:
      knit_20: 0.10
      knit_theta: 1.083
      kdnit_20: 0.10
      kdnit_theta: 1.045
      KsOxdn: 0.1
      vno3_20: 0.0
      vno3_theta: 1.08
      KNR: 0.6
      kon_20: 0.1
      kon_theta: 1.074
      vson_20: 0.1
      vson_theta: 1.024
      PN: 0.5
      NH4fromBed: 0.0
      NO3_BedDenit: 0.0
      use_NH4: True
      use_NO3: True
      use_OrgN: True
```

Default values are taken from v1's `nitrogen.py` `Nitrogen` `TypedDict` defaults; sentinel-`999` defaults are corrected to physically reasonable values during the port (see Section 7).

---

## 4. Component Inventory

### 4.1 Existing v2 NSM1 components to fix and extend

| Component | Action | Notes |
|---|---|---|
| `Nitrogen` (NH4, NO3) | Fix integrator bugs (lines 101, 115); fix `time_step_frequency` typo; replace `rate == np.nan` with `.isnull()` (4 places); wire up `half_saturation_oxygen` parameter; remove hard-coded `algea_growth_rate=0` and route through algal processes; **extend to add OrgN** | The integrator fix changes simulated values |
| `FloatingAlgae` (Ap) | Fix integrator bug (line 122) including stray `* 86400`; implement `ammonium_respiration()` and `ammonium_growth()` (currently return 0); replace hard-coded `phosphate_fraction_dissolved=0.5` with TIP partitioning function; wire up algal-mortality routing to OrgN/OrgP/POC/DOC pools | The implementations of `ammonium_respiration` and `ammonium_growth` are essential — without them, algal N coupling silently doesn't happen |
| `BenthicAlgae` (Ab) | Same integrator and NaN guard fixes as FloatingAlgae; implement equivalent `ammonium_respiration` and `ammonium_growth` for benthic algae; add benthic-algae burial and Fw/Fb mortality fractionation routing | Currently a stub with limited coverage of v1's benthic algae kinetics |

### 4.2 Missing constituents to port from v1

The 12 constituents to port. Each becomes a v3 `Process` class. Order in the table reflects suggested implementation order (simpler/independent first).

| Constituent | v1 source | Complexity | Dependencies |
|---|---|---|---|
| Pathogen (PX) | `nsm1/processes.py:3141-3227` | Small | Independent (uses light extinction) |
| Organic Nitrogen (extension to existing Nitrogen) | `nsm1/processes.py:1173-1405` | Small-Medium | Algal mortality routing |
| Organic Phosphorus + TIP (`Phosphorus`) | `nsm1/processes.py:1833-2168` | Medium | Algal uptake; fdp partitioning |
| POM | `nsm1/processes.py:2185-2317` | Medium | Algal mortality, settling |
| CBOD | `nsm1/processes.py:2334-2422` | Medium | DOX (for oxidation) |
| Pathogen, OrgN, P, POM, CBOD subtotal | — | — | — |
| Carbon (POC, DOC, DIC) | `nsm1/processes.py:2439-2858` | Large | Algal photosynthesis/respiration; CO2 reaeration; sediment release |
| DOX | `nsm1/processes.py:2876-3123` | Large | Almost everything else (algal photosynthesis/respiration, nitrification, DOC oxidation, CBOD oxidation, SOD) |
| Alkalinity | `nsm1/processes.py:3246-3435` | Medium | Nitrification, denitrification, algal growth/respiration |
| N2/TDG | `nsm1/processes.py:3452-3540` | Small-Medium | Denitrification (consumes NO3, produces N2); Henry's law for atmospheric exchange |

### 4.3 Shared physics primitives to port

v1's `shared/processes.py` contains physics primitives used by multiple NSM1 processes. These need v3 equivalents:

| Primitive | v1 location | v3 destination |
|---|---|---|
| `arrhenius_correction` | `shared/processes.py:16` | Already in `clearwater_modules_v2/utils/conversions.py` (re-exported in v3) |
| `celsius_to_kelvin`, `kelvin_to_celsius` | `shared/processes.py:6,11` | Already in v2 utils |
| `kah_20` (9 hydraulic reaeration options) | `shared/processes.py:65` | New: `clearwater_modules_v3/utils/reaeration.py` |
| `kaw_20` (13 wind reaeration options) | `shared/processes.py:101` | Same file |
| `ka_tc` (temperature-corrected combined reaeration) | `shared/processes.py:165` | Same file |
| `SOD_tc` (temperature-corrected sediment oxygen demand) | `shared/processes.py:239` | New: `clearwater_modules_v3/utils/sediment.py` |
| `L` (Beer-Lambert light extinction with ISS, POC, Chl-a contributions) | `shared/processes.py:202` | New: `clearwater_modules_v3/utils/light.py` |
| `PAR` (photosynthetically active radiation derivation) | `shared/processes.py:222` | Same file |
| `fdp` (TIP solid-dissolved partitioning) | `shared/processes.py:257` | New: `clearwater_modules_v3/utils/partitioning.py` |

These primitives are stateless functions that take their inputs as `ArrayLike` and return scalars or arrays. They are shared utilities used across multiple `Process` classes.

### 4.4 Parameter library

v3 NSM1 parameters are organized into the following groups, mapping to v1's `TypedDict` structure:

| Group | v1 TypedDict | v3 Process consumer |
|---|---|---|
| algae | `Algae` | `FloatingAlgae` |
| balgae | `Balgae` | `BenthicAlgae` |
| nitrogen | `Nitrogen` | `Nitrogen` |
| phosphorus | `Phosphorus` | `Phosphorus` |
| carbon | `Carbon` | `Carbon` |
| POM | `POM` | `POM` |
| CBOD | `CBOD` | `CBOD` |
| DOX | `DOX` | `DOX` |
| pathogen | `Pathogen` | `Pathogen` |
| alkalinity | `Alkalinity` | `Alkalinity` |
| N2 | `N2` | `N2` |
| global_parameters | `GlobalParameters` | All processes (passed via dependency injection or registry) |
| global_vars | `GlobalVars` | Initial conditions for state variables |

Parameters that need correction during port (see Section 7 for the canonical list of seven critical default-value corrections, including the Phase 0 finding `pressure_mb=2026.5 → 1013.25 hPa`).

The corrections are documented in v3's `parameter_defaults_corrections.md` (a new doc) so the rationale is on record. The DOX bug investigation already documented in `NSM1_DOX_rate_bug_investigation.md` is referenced.

#### Other defaults under Phase 1 audit

The Phase 0 parameter audit (`docs/clearwater_modules_v3_nsm1_phase0_parameter_audit.md`) flagged 8 additional defaults that do not block the port but warrant Phase 1 review and documentation:

- `rnh4_20=0`, `vno3_20=0`, `rpo4_20=0` — sediment release silently disabled; verify gated by `use_SedFlux`
- `kdpo4=0.0` — TIP partitioning feature disabled; NSM2 territory
- `ksbod_20=0.0` — CBOD never settles; confirm with LimnoTech whether intentional
- `apx=1`, `vx=1` — pathogen placeholders without literature basis; document units and source
- `h2=0.1` — POM dissolution depth, unclear physical role; add docstring
- `vb=0.01` — burial velocity magnitude; document with reference
- `q_solar=500` units — code docstring says `1/d` but actually W/m²; standardize
- `lambdas` — light extinction parameter defined in constants but disabled in code; remove or activate

Phase 1 dispositions each item in `parameter_defaults_corrections.md` alongside the seven critical corrections.

---

## 5. Per-Constituent Design Notes

This section captures non-obvious design decisions for each constituent. Routine ports (kinetics formula matches v1) are not enumerated here.

### Nitrogen (extending existing v2 Nitrogen with OrgN)

- The existing v2 `Nitrogen.change_ammonium` and `change_nitrate` methods are corrected, not rewritten. The body of each method already enumerates the right source/sink terms; only the integrator and parameter wiring need fixing.
- OrgN is added as a third state variable. `change_organic_nitrogen` is added as a method, including hydrolysis to NH4, settling, and contributions from algal mortality.
- The algal mortality routing reads `floating_algae_process.death_to_orgn()` and `benthic_algae_process.death_to_orgn()` (new methods on those processes that return the rate of N transferred from algal death to OrgN).

### FloatingAlgae and BenthicAlgae (extensions)

- `ammonium_respiration()` returns `rna * AlgalRespiration` per v1 line 1273
- `ammonium_growth()` returns `rna * PN_calculated * AlgalGrowth` where `PN_calculated` uses v1's preference function (line 1206)
- `death_to_orgn()` and `death_to_orgp()` return per-element mortality routing for downstream processes
- Light limitation uses v3's `Beer-Lambert L(...)` from the new shared utilities

### Phosphorus (TIP, OrgP)

- TIP partitioning uses `fdp` shared utility
- Settling velocity for inorganic P: `vs/depth * (1-fdp) * TIP` matches v1
- OrgP follows same pattern as OrgN: hydrolysis to TIP, settling, algal mortality contribution

### Carbon (POC, DOC, DIC)

- POC and DOC mineralization with Monod DO attenuation per v1
- DIC reaeration with atmosphere uses Henry's law: `0.923 * ka_tc * (KH * pCO2/1e6 - Fco2 * DIC)` (v1 line 2742)
- Algal photosynthesis/respiration coupling: DIC sink from photosynthesis, source from respiration
- Sediment release: `JDIC/depth/12000` if SedFlux enabled, else parameterized constant

### DOX

- `O2sat` calculation uses APHA / QUAL2E formulation with pressure and salinity correction per v1 lines 2935-2956
- Atmospheric reaeration: `ka_tc * (O2sat - DOX)` where `ka_tc` comes from the v3 reaeration menu utility
- Photosynthesis/respiration source: reads algal rates from FloatingAlgae and BenthicAlgae, applies stoichiometric coefficients with NH4 vs NO3 uptake fractionation
- Nitrification O2 sink: reads `nitrification_rate` from Nitrogen process
- DOC oxidation O2 sink: reads `doc_dic_oxidation_rate` from Carbon process
- CBOD oxidation O2 sink: reads `cbod_oxidation_rate` (sum over groups) from CBOD process
- SOD: uses `SOD_tc` shared utility (with the corrected `SOD_theta` parameter)

### Pathogen (PX)

- Direct port of v1 logic: natural death + light-induced death + settling
- Light-induced death uses v3's `KEXT` from the light extinction shared utility

### Alkalinity (newly active in v3)

- v1 declares `Alk` state variable but never updates it
- v3 implements the alkalinity ODE per v1 lines 3246-3435: nitrification consumption, denitrification production, algal growth/respiration coupling
- Note: full carbonate-pH solver is NSM2 territory; v3 1.0.0 just integrates Alk as a tracer with the source/sink terms above. pH solver comes in v3 1.1+ as part of NSM2 features integration.

### N2 / TDG

- N2 saturation from Henry's law: `N2sat = 2.8e4 * KHN2_tc * 0.79 * (pressure_atm - p_wv)`
- Atmospheric exchange: `1.034 * ka_tc * (N2sat - N2)`
- Denitrification source: reads `denitrification_rate` from Nitrogen process
- TDG derived variable: `TDG = N2/N2sat` (or weighted with O2 if both are present)

---

## 6. Bug Fix Inventory for Existing v2 Constituents

Consolidated list of bugs that v3 fixes in the 4 existing v2 NSM1 constituents, with source files and line numbers:

| # | Bug | Location | Status on `streaming` | Status on `memory-refactor-pytestUpdate` | v3 fix |
|---|---|---|---|---|---|
| 1 | Multiplicative integrator for NH4: `ammonium = 0 + ammonium * rate * dt` | `nitrogen.py:101` | BROKEN | BROKEN | Replace with `ammonium = ammonium + rate * dt` |
| 2 | Multiplicative integrator for NO3: same pattern | `nitrogen.py:115` | BROKEN | BROKEN | Same fix |
| 3 | `time_step_frequency` typo on line 115 | `nitrogen.py:115` | FIXED | BROKEN | Already fixed on streaming; carry forward |
| 4 | Multiplicative integrator + stray `* 86400` for Ap | `floating_algae.py:122` | BROKEN | BROKEN | Replace with `algae = algae + rate * dt` (no extra 86400) |
| 5 | NaN guard `rate == np.nan` (always False) | `nitrogen.py:147` | FIXED | BROKEN | Use `.isnull()` |
| 6 | NaN guard same | `nitrogen.py:218` | FIXED | BROKEN | Use `.isnull()` |
| 7 | NaN guard same | `nitrogen.py:250` | BROKEN | BROKEN | Use `.isnull()` |
| 8 | NaN guard same | `nitrogen.py:313` | BROKEN | BROKEN | Use `.isnull()` |
| 9 | Hard-coded `half_saturation_oxygen=1` | `nitrogen.py:197` | TODO | TODO | Wire up parameter |
| 10 | Hard-coded `algea_growth_rate=0` for nitrate uptake | `nitrogen.py:211` | TODO | TODO | Read from FloatingAlgae process |
| 11 | Hard-coded `algea_growth_rate=0` for benthic uptake | `nitrogen.py:217` | TODO | TODO | Read from BenthicAlgae process |
| 12 | Death rate placeholder | `nitrogen.py:60` | TODO | TODO | Wire up to floating_algae_process |
| 13 | `ammonium_respiration()` returns 0 | `floating_algae.py:401` | TODO | TODO | Implement per v1 line 1273 |
| 14 | `ammonium_growth()` returns 0 | `floating_algae.py:407` | TODO | TODO | Implement per v1 line 1206 |
| 15 | `phosphate_fraction_dissolved=0.5` hard-coded | `floating_algae.py:113` | TODO | TODO | Replace with `fdp` from shared utility |
| 16 | `set_at_time` persistence (state computed but dropped) | `nitrogen.py:run`, `floating_algae.py:run` | FIXED | BROKEN | Already fixed on streaming; carry forward |

The fixes from `streaming` (#3, #5, #6, #16) need to be brought forward into v3 alongside the new fixes.

---

## 7. Critical Default-Value Corrections

v1's `nsm1/constants.py` contains 7 parameter defaults that v3 corrects at the port: 6 sentinel-`999` values and 1 magnitude error in `pressure_mb` (identified during Phase 0 parameter audit, 2026-05-04).

| Parameter | v1 default | Risk | v3 default | Source |
|---|---|---|---|---|
| `vsop` (OrgP settling velocity) | 999 | Multiplied into rate; wrong magnitude propagates | `0.1` m/d | Typical 0.01–1.0 m/d |
| `vs` (TIP settling velocity) | 999 | Same | `0.1` m/d | Typical 0.01–1.0 m/d |
| `SOD_20` (SOD at 20 °C) | 999 | Wrong magnitude propagates | `1.0` g-O2/m²/d | Reasonable for moderate organic loading |
| `SOD_theta` (Arrhenius theta for SOD) | 999 | `999^(T-20)` — catastrophic blowup at T > 20 °C | `1.060` | Chapra 1997 standard |
| `kaw_20_user` (user-override wind reaeration at 20 °C) | 999 | Used only when option set; spurious if set incorrectly | `0.0` m/d | Disabled unless user opts in |
| `kah_20_user` (user-override hydraulic reaeration at 20 °C) | 999 | Same | `0.0` m/d | Disabled unless user opts in |
| `pressure_mb` (atmospheric pressure) | 2026.5 | ~2× sea-level pressure; biases `O2sat`, `N2sat`, atmospheric reaeration | `1013.25` hPa | Standard sea-level pressure (ISO 2533) |

Each fix is documented in `clearwater_modules_v3/parameter_defaults_corrections.md` (new doc) with the rationale. Regression tests added: `test_sod_theta_no_blowup.py` confirms that with default parameters, a 30 °C simulation produces finite SOD values; `test_pressure_correction.py` confirms `O2sat` and `N2sat` track standard tabulated values at the corrected pressure.

Phase 0 also flagged 8 lower-priority parameter findings (sediment-flux gating, `ksbod_20=0`, pathogen placeholders, units/docstring gaps) that do not require spec amendments; see `docs/clearwater_modules_v3_nsm1_gap_analysis.md` Section 3.2 for the full list and Phase 1 actions.

---

## 8. Migration Strategy

### From v2 NSM1 to v3 NSM1

For users who have v2 NSM1 code or YAML configs:

| v2 reference | v3 equivalent |
|---|---|
| `from clearwater_modules_v2.processes.nitrogen import Nitrogen` | `from clearwater_modules_v3.processes.nitrogen import Nitrogen` |
| `from clearwater_modules_v2.processes.floating_algae import FloatingAlgae` | `from clearwater_modules_v3.processes.floating_algae import FloatingAlgae` |
| `from clearwater_modules_v2.processes.benthic_algae import BenthicAlgae` | `from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae` |
| YAML: `- nitrogen: time_step: '30s'` | YAML: same, plus `parameters:` block now required for full functionality |

Behavioral differences after migration:
- v3 NSM1 actually computes algal nitrate uptake (v2 silently set this to zero); expect lower steady-state nitrate concentrations in any simulation with active algae
- v3 NSM1 actually persists state changes (v2 dropped them); expect non-zero state evolution where v2 produced near-constant outputs
- v3 NSM1 applies the additive integrator instead of the multiplicative one; numerical outputs differ substantially

### From v1 NSM1 to v3 NSM1

For users who have v1 NSM1 code (the larger migration):

| v1 reference | v3 equivalent |
|---|---|
| `from clearwater_modules.nsm1.model import NutrientBudget` | Use `clearwater_modules_v3.config.init_from_file` with a YAML that lists all NSM1 processes |
| Code-defined parameters: `algae_parameters={...}, ...` | YAML `parameters:` blocks per process |
| Direct invocation: `nsm.increment_timestep()` | YAML-driven `model.run()` |
| Hotstart: `NutrientBudget(hotstart_dataset=ds, hotstart_timestep=0, ...)` | YAML `hotstart:` block with `dataset_path` and `timestep` |
| State variable access: `nsm.dataset['NH4']` | `model._Model__registry.get('ammonium')` (or v3-provided accessor TBD) |

A `migration_v1_to_v3_nsm1.md` doc accompanies v3 with worked examples.

### Retirement timeline (consistent with overall v3 plan)

NSM1 v3 1.0.0 ships by 2026-05-31 (target). At that point:
- v1 NSM1 (`clearwater_modules.nsm1`) is frozen
- v2 NSM1 (`clearwater_modules_v2.processes.nitrogen` etc.) is frozen
- All NSM1 development goes to v3

v3 1.1.0 (post-2026-05-31) removes v1 and v2 NSM1 source from the tree alongside the broader v2 retirement.

---

## 9. Testing and Validation

### Test infrastructure

v3 NSM1 inherits the test directory structure and adds:

- `tests/v3/nsm1/test_<constituent>_calculations_v3.py` for each of 11 processes — uses v1's expected-value tables in `tests/NSM Manual Calcs/` (which are framework-independent), exercising v3's `Process.run` flow
- `tests/v3/nsm1/test_v1_v3_parity.py` — A/B comparison: instantiate v1 NSM1 (`NutrientBudget`) and v3 NSM1 with matched parameters, run both for a few timesteps, assert state-variable evolution matches within tolerance for each constituent
- `tests/v3/nsm1/test_v2_v3_parity_existing.py` — for the 4 constituents that exist in both v2 and v3, run v2 (with bugs intact) and v3 (with bugs fixed) on a controlled case where the fixes don't matter (e.g., zero algae growth → algal coupling differences vanish), assert outputs match
- `tests/v3/nsm1/test_coupled_nsm1_riverine_v3.py` — programmatic end-to-end test of a coupled NSM1+Riverine simulation, asserting mass conservation and reasonable evolution
- `tests/v3/nsm1/test_sod_theta_no_blowup.py` — regression test for the sentinel-`999` fix
- `tests/v3/nsm1/test_validation_tier1_conservation.py` — closed-system mass balance per the no-reference validation strategy: total N, P, C, O₂-equivalents conserved to roundoff when boundary fluxes are zero

### Validation tiers

Per the validation strategy (`Validation_strategy_no_reference.md`):

- **Tier 1 (conservation):** mass balance for N, P, C, O₂-equivalents, alkalinity. Closed-system test: no boundaries, no settling, all source/sink pairs internally balanced. Assert totals constant to roundoff AND `model.diagnostics.clip_events == 0`. A clip event under closed-system test conditions indicates either unphysical parameters or a malformed test case — the test fails in that situation, which is the correct diagnostic. **Mandatory before merge.**
- **Tier 2 (analytical limits):** each kinetic function reduces to known closed-form solution under simplified conditions (first-order decay, Streeter-Phelps DO sag, Monod growth at constant nutrient, Henry's-law equilibration for DIC and N2). **Mandatory before merge.**
- **Tier 3 (steady-state algebra):** constant-inflow steady-state matches algebraic balance. **Recommended.**
- **Tier 5 (sensitivity):** parameter sweeps confirm physical direction of response. **Recommended.**

Tiers 4 (MMS) and 6 (published case-study comparison) are deferred to NSM2 features integration in v3 1.1+.

---

## 10. Performance Targets

The current v2 NSM1 doesn't run end-to-end (4 of 16 constituents), so there is no v2 baseline to measure against directly. Reference points:

- **v1 NSM1** with the 418× kernel optimization: ~3 ms per timestep on a 5-cell mesh (per the kernel-optimization commit benchmark). For Sumwere Creek (~600 cells), that scales to roughly 1-2 seconds per timestep.
- **v2 TSM** coupled with Riverine on Sumwere Creek (the verified baseline): 89 seconds for 4,320 timesteps (~20 ms/step including transport). NSM1 is more complex than TSM (~16 constituents vs 1), so per-step time will be higher.

v3 NSM1 targets:

- **Must:** complete a 4,320-timestep coupled NSM1+Riverine simulation on Sumwere Creek within 30 minutes (≤ 415 ms/step).
- **Should:** complete the same simulation within 10 minutes (≤ 138 ms/step).
- **Aspirational:** within 3 minutes (≤ 41 ms/step), competitive with v1's optimized performance.

If the "must" target is not met, performance optimization is added in v3 1.0.1.

---

## 11. Phased Implementation Plan

### Phase 0 — Gap analysis and parameter audit (1 day)

- Per-constituent diff: v1 NSM1 source vs v2 NSM1 source (where overlap exists). Catalog which v1 functions become which v3 methods.
- Parameter audit: verify v1's default parameter values for each of the 13 TypedDict groups. Identify any other bad defaults beyond the documented sentinel-`999` set.
- Test fixture audit: enumerate v1's NSM1 calculation tests and the expected-value tables they use.

**Deliverable:** `docs/clearwater_modules_v3_nsm1_gap_analysis.md` with the per-constituent and parameter inventory.

### Phase 1 — Foundation: shared physics primitives + parameter library (1-2 days)

- Implement `clearwater_modules_v3/utils/reaeration.py` with `kah_20`, `kaw_20`, `ka_tc`
- Implement `clearwater_modules_v3/utils/sediment.py` with `SOD_tc`
- Implement `clearwater_modules_v3/utils/light.py` with `L`, `PAR`
- Implement `clearwater_modules_v3/utils/partitioning.py` with `fdp`
- Establish parameter-library pattern: `Process.DEFAULTS` class-level dicts merged with YAML overrides
- Apply sentinel-`999` corrections in DEFAULTS

**Deliverable:** v3 utils package with shared primitives; parameter handling demonstrated on existing v2 Nitrogen process.

### Phase 2 — Fix 4 existing v2 NSM1 processes (2-3 days)

- Apply all 16 bug fixes from Section 6 to `Nitrogen`, `FloatingAlgae`, `BenthicAlgae`
- Implement `ammonium_respiration`, `ammonium_growth` on FloatingAlgae and BenthicAlgae
- Wire up algal mortality routing (`death_to_orgn`, `death_to_orgp`, `death_to_poc`, `death_to_doc` methods)
- Add OrgN as a third state variable in Nitrogen
- Replace `phosphate_fraction_dissolved=0.5` with `fdp` shared utility
- Tier 1 conservation test for N (closed system) must pass

**Deliverable:** existing 4 v2 processes corrected and extended; OrgN now part of Nitrogen process; Tier 1 N conservation passes.

### Phase 3 — Port simple constituents (2-3 days)

- Pathogen (small)
- POM (medium)
- CBOD (medium)
- N2/TDG (small-medium)

Each follows the pattern: read v1 process, translate kinetics into a v3 `Process` subclass, implement `run` with the standard integrator pattern, wire registry inputs/outputs, add Tier 1 conservation test.

**Deliverable:** 4 new Process classes; YAML can declare and run a simulation including these; Tier 1 conservation passes for each.

### Phase 4 — Port Phosphorus (1-2 days)

- TIP and OrgP as a single Phosphorus Process
- TIP partitioning via `fdp` utility
- Algal P uptake routing
- Tier 1 P conservation test

**Deliverable:** Phosphorus Process; YAML extended; Tier 1 P passes.

### Phase 5 — Port Carbon and DOX (3-4 days, the largest single chunk)

- Carbon Process (POC, DOC, DIC) with all source/sink couplings
- DOX Process with all sink couplings (algal photosynthesis/respiration, nitrification, DOC oxidation, CBOD oxidation, SOD, atmospheric reaeration, DIC reaeration coupling)
- Tier 1 C and O2-equivalent conservation tests
- Streeter-Phelps analytical-limit test for DOX

**Deliverable:** Carbon and DOX Processes; coupled simulation runs end-to-end with all 11 processes; Tier 1 and Tier 2 tests pass.

### Phase 6 — Port Alkalinity (1 day)

- Alk Process with nitrification, denitrification, algal growth/respiration source/sink terms
- Note: full carbonate solver is NSM2 territory (v3 1.1+); v3 1.0.0 just integrates Alk as a tracer
- Tier 1 alkalinity equivalence test

**Deliverable:** Alkalinity Process; Tier 1 alk passes.

### Phase 7 — Tests, validation, end-to-end demo (2-3 days)

- Port v1 calculation tests for all 11 processes (Tier 0 — value matching against v1 expected outputs)
- v1/v3 parity test
- v2/v3 parity test for existing 4 constituents
- Coupled NSM1+Riverine demo notebook (extending the TSM+Riverine demo to also exercise nutrients and algae)
- Performance benchmark against the "must" target

**Deliverable:** Full test suite passing; coupled demo runs and produces sensible nutrient/biomass timeseries.

### Phase 8 — Documentation and review prep (1 day)

- v3 NSM1 README
- Migration notes (v1→v3, v2→v3)
- Parameter defaults corrections doc
- LimnoTech review materials

**Deliverable:** v3 NSM1 1.0.0 ready for LimnoTech sign-off.

**Total estimated wall-clock with Claude doing the coding: 14-20 working days.** Fits within the 18-day window from May 14 (after TSM v3 ships) through May 31 if reviews are tight and no major surprises surface.

---

## 12. Coordination with LimnoTech

### What to discuss with Anthony (concurrent with v3 NSM1 work)

1. The v2 NSM1 partial implementation's bug status (the 16 items in Section 6) and whether any are LimnoTech's intent that should be preserved differently than my fix proposes.
2. The decision to make NSM1 v3 use the same v3 `Model` as TSM v3 (so kernel optimization, wet-mask, hotstart apply to NSM1 too).
3. The integrator-pattern contract (additive Forward Euler, units-per-second rates) as the v3 standard.
4. The parameter library pattern (`Process.DEFAULTS` + YAML override merge) for sign-off as the v3 standard.
5. The dispatch-order question: are any LimnoTech demonstrations relying on a different process order than the one in Section 3?

### What can proceed without LimnoTech input

- Phase 0 (gap analysis): pure analysis
- Phase 1 (shared utilities): physics primitives are independent of framework discussion
- Phase 2 items 1-4 (integrator fixes, NaN guards, persistence): pure bug fixes with no design ambiguity
- Phase 3 items (simple constituent ports): standard kinetics, no design ambiguity
- All test work

Items that genuinely benefit from LimnoTech input:

- Algal-uptake stub implementation (Phase 2 items): `ammonium_respiration` and `ammonium_growth` semantics — confirm LimnoTech intends the v1 formulations
- Parameter library pattern (Phase 1)
- Coupled-process method conventions (`death_to_orgn` etc.)

---

## 13. Risks and Mitigations

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| 14-20 day estimate slips beyond May 31 | Medium-High | Medium | Tight phase gates; if Phase 5 (Carbon+DOX, the largest) overruns, defer Pathogen+N2 to v3 1.0.1 to ship core nutrient cycle on time |
| v1/v3 parity test reveals systematic differences (beyond floating-point noise) | Medium | Medium-High | Build the parity test early in Phase 2 (alongside Nitrogen extension); each new constituent in Phase 3-6 is added to the parity test as it lands |
| Inter-process coupling via registry rates is awkward or slow | Medium | Medium | If coupling is awkward, fall back to the v2 pattern of direct `process_other.method()` calls; if slow, profile and optimize |
| Sentinel-`999` correction values turn out to differ from what users expect | Low | Low | Document each correction with rationale; users can override via YAML; corrections are discoverable in `parameter_defaults_corrections.md` |
| Phase 5 DOX port is harder than estimated due to coupling complexity | Medium | Medium | Start Phase 5 by building the DOX skeleton with all couplings stubbed out, then fill in formulas one at a time |
| Tier 1 conservation tests fail in subtle ways | Medium | High | Build the test harness once in Phase 1; run after every constituent addition; failures surface immediately, not at the end |
| Performance target ("must" of 30 min for 4,320 steps) not met | Low-Medium | Medium | If exceeded, port the kernel optimization concept from TSM v3 to NSM1 dispatch; profile to identify hot loops; numba decoration available as a fallback |
| Hotstart and wet-mask design from TSM v3 doesn't extend cleanly to NSM1's multi-process state | Low-Medium | Medium | NSM1 v3 design is built on top of TSM v3's `Model`; if the design surfaces issues in NSM1, revise the v3 `Model` design with NSM1 in mind |
| LimnoTech objects to v3 NSM1 design after substantial work has been done | Low (with early communication) | High | Send the v2 NSM1 bug list to Anthony alongside the v3 design spec; their reply informs the design before code is written |

---

## 14. Design Decisions

All design questions surfaced during spec review are resolved. Three of the LimnoTech-input items (Alkalinity scope, sediment-flux structure, single-compartment algae) are tentatively decided as listed below; LimnoTech confirmation closes them. None block the start of Phase 0 or Phase 1.

**Resolved:**
- OrgN is added as a third state variable on the existing `Nitrogen` Process (extension, not a separate `OrganicNitrogen` Process), consistent with the topical-grouping rule applied to `Phosphorus` (TIP+OrgP) and `Carbon` (POC+DOC+DIC). Decided 2026-05-04.
- Inter-process rate variables in the registry follow the names listed in Appendix A (e.g., `algal_growth_rate`, `algal_respiration_rate`, `algal_death_rate`, `algal_nh4_uptake_fraction`, `nitrification_rate`, `denitrification_rate`, `doc_dic_oxidation_rate`, `cbod_oxidation_rate`, `sod_rate`, `nh4_from_bed`, `dip_from_bed`, `no3_from_bed_denit`, `dic_from_bed`). Convention: snake_case, suffix `_rate` for time-derivative quantities, `_fraction` for dimensionless ratios, source-named prefixes for sediment fluxes. Decided 2026-05-04.
- Negative-state handling uses the **clip-with-log** contract: every `Process.run` calls `clip_negative_state(...)` after the Forward Euler step; clips are counted on `model.diagnostics`; Tier 1 closed-system tests assert `clip_events == 0`. Clip target is exactly 0 (not epsilon). Adaptive substepping is deferred to v3 1.1+; semi-implicit sink treatment is a Phase 5 per-Process opt-in (most likely for DOX). Decided 2026-05-04.
- Within-timestep update semantics: **Jacobi for state, Gauss-Seidel for rate variables.** State reads always pass `t=t_current` (pre-update values, order-independent); rate variables are step-scoped, written by producers and read by consumers within the same step, cleared at start of each step, and reading an unset rate variable raises an error. Naming convention from Appendix A (`_rate` / `_fraction` suffixes) marks the lifecycle class. Decided 2026-05-04.
- Package location: v3 NSM1 lives in `src/clearwater_modules_v3/` alongside v3 TSM. The existing layout (`config/`, `processes/`, `utils/`, `model.py`, `README.md`, `__init__.py`) is the foundation; Phase 1 adds `utils/{reaeration,sediment,light,partitioning,numerics}.py` and `parameters/<group>.py`; Phase 2+ extends `processes/` with the NSM1 Process classes. Decided 2026-05-04.
- v3 NSM1 ships a working default `modules.yml` that exercises all 11 Process classes end-to-end. Built and validated as part of Phase 7 alongside the coupled NSM1+Riverine demo notebook; lives at `src/clearwater_modules_v3/config/nsm1_default.yml` (or equivalent path under `config/`). Decided 2026-05-04.
- Alkalinity in v3 1.0.0 is a **simple tracer**: state variable with source/sink terms (nitrification consumption, denitrification production, algal growth/respiration coupling) integrated by Forward Euler. No carbonate equilibrium, no pH solver. Full pH chemistry (carbonate equilibrium, NH3/NH4+ partitioning, free-CO2 fraction, carbonate speciation) is NSM2 territory in v3 1.1+. v3 1.0.0 documentation includes a worked example of post-hoc pH computation from `Alk`, `DIC`, `T`, salinity output trajectories for users who need a quick number. Decided 2026-05-04 (pending LimnoTech confirmation that no current applications require pH from day one).
- Sediment-flux parameters (`SOD_20`, `NH4fromBed`, `DIPfromBed`, `NO3_BedDenit`, `DIC_sed_release`) are **scalar globals** in v3 1.0.0, applied uniformly to all cells, set in YAML (matching v1's pattern exactly). Per-cell spatially varying fluxes and dynamically computed fluxes both arrive in v3 1.1+ via the NSM2 sediment diagenesis Process. Decided 2026-05-04 (pending LimnoTech confirmation that no current applications use spatially varying SOD or bed fluxes).
- Phytoplankton and benthic algae in v3 1.0.0 remain **single-compartment** (one `Ap` and one `Ab`), matching v1/v2. NSM2's multi-group capability lands as new Process classes (e.g., `PhytoplanktonGroups`) added *alongside* the single-compartment `FloatingAlgae`/`BenthicAlgae`, not as in-place extension. The YAML-driven Process registration framework supports this additive extension natively, so framework extensibility is not a blocker for 1.0.0. Decided 2026-05-04 (pending LimnoTech confirmation that no near-term applications require multi-group before v3 1.1+).

---

## 15. Approval Criteria

This specification is complete enough to proceed if the reviewer agrees that:

1. The motivation for v3 NSM1 (Section 1) accurately characterizes the asymmetric state of v1 and v2 NSM1.
2. The 11-process organization (Section 3) is the right granularity (Process per topic, with multiple state variables per Process where physically grouped).
3. The integrator-pattern contract (Section 3, additive Forward Euler) is correct.
4. The component inventory (Section 4) accurately captures what v3 NSM1 will contain.
5. The bug fix inventory (Section 6) is correct and complete.
6. The sentinel-`999` corrections (Section 7) are reasonable defaults.
7. The phased plan (Section 11) is the right sequence and the May 31 target is achievable with the listed mitigations.
8. The risks (Section 13) are correctly identified and the mitigations are reasonable.
9. The open questions (Section 14) are the right ones to surface for discussion.

If any section is wrong, incomplete, or misframed, mark up this document directly and the spec will be revised before any code is written.

---

## Appendix A: v3 NSM1 process registry-coupling cheatsheet

Variables written to the registry by each Process for downstream consumers:

| Producer | Registry variable | Consumed by |
|---|---|---|
| FloatingAlgae | `algal_growth_rate` | DOX (photosynthesis source), Carbon (DIC sink), Nitrogen (NH4/NO3 sink), Phosphorus (TIP sink), Alkalinity |
| FloatingAlgae | `algal_respiration_rate` | DOX (sink), Carbon (DIC source), Nitrogen (NH4 source) |
| FloatingAlgae | `algal_death_rate` | OrgN, OrgP, POC, DOC, POM (mortality routing) |
| FloatingAlgae | `algal_nh4_uptake_fraction` (PN_calculated) | Nitrogen (split between NH4 and NO3 uptake) |
| BenthicAlgae | analogous to FloatingAlgae | Same downstream consumers, with depth-integration adjustment |
| Nitrogen | `nitrification_rate` | DOX (sink), Alkalinity (source) |
| Nitrogen | `denitrification_rate` | Alkalinity (source), N2 (source) |
| Carbon | `doc_dic_oxidation_rate` | DOX (sink) |
| CBOD | `cbod_oxidation_rate` (sum over groups) | DOX (sink) |
| Sediment (parameterized in v3 1.0.0) | `sod_rate`, `nh4_from_bed`, `dip_from_bed`, `no3_from_bed_denit`, `dic_from_bed` | DOX, Nitrogen, Phosphorus, Carbon |

These names are the v3 NSM1 registry convention (decided 2026-05-04). Phase 0 verifies completeness against the per-constituent kinetics audit; additions during Phase 0 follow the same naming rules (snake_case, `_rate` for time derivatives, `_fraction` for dimensionless ratios, source-named prefixes for sediment fluxes).

---

## Appendix B: Parameter library file structure (proposal)

```
src/clearwater_modules_v3/
├── parameters/
│   ├── __init__.py
│   ├── algae.py              # FloatingAlgae DEFAULTS dict
│   ├── balgae.py             # BenthicAlgae DEFAULTS dict
│   ├── nitrogen.py           # Nitrogen DEFAULTS
│   ├── phosphorus.py         # Phosphorus DEFAULTS
│   ├── carbon.py             # Carbon DEFAULTS
│   ├── pom.py                # POM DEFAULTS
│   ├── cbod.py               # CBOD DEFAULTS
│   ├── dox.py                # DOX DEFAULTS
│   ├── pathogen.py           # Pathogen DEFAULTS
│   ├── alkalinity.py         # Alkalinity DEFAULTS
│   ├── n2.py                 # N2 DEFAULTS
│   └── global.py             # Global/shared parameters
└── parameter_defaults_corrections.md   # rationale for each sentinel-999 fix and any other corrections
```

Each `parameters/*.py` defines a single `DEFAULTS: dict[str, float]` exported by name. Each Process imports its DEFAULTS at construction time and merges with user-provided values:

```python
from clearwater_modules_v3.parameters.nitrogen import DEFAULTS as NITROGEN_DEFAULTS

class Nitrogen(Process):
    def __init__(self, parameters: dict, time_step: timedelta = ..., ...):
        merged = {**NITROGEN_DEFAULTS, **parameters}
        for k, v in merged.items():
            setattr(self, k, v)
        Process.__init__(self, time_step)
```

This keeps parameter values discoverable, documented (each `DEFAULTS` dict can have inline comments), and easy to override per simulation via YAML.
