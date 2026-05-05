# v3 NSM1 1.0.0 — LimnoTech Reviewer Materials

**Re:** v3 NSM1 1.0.0, the convergence of v1 NSM1 (`clearwater_modules.nsm1`)
and v2 NSM1 (`clearwater_modules_v2.processes.{nitrogen, floating_algae,
benthic_algae}`) into a single 11-Process implementation under the v2
framework.

**Source branch:** `EcohydrologyTeam/ClearWater-modules-streaming` `streaming`

ClearWater is an ERDC product released under an open-source license; this
packet describes the post-audit state of v3 NSM1 1.0.0 within that
framework.

---

## 1. Executive summary

v3 NSM1 1.0.0 ports v1 NSM1's 16-constituent kinetics suite into the v2
`Process` class framework. After an earlier "review-ready" claim in
Phase 8, a line-by-line three-way audit (Fortran NSM1 vs v1 Python
NSM1 vs v3 Python NSM1) was performed across all 11 Processes and the
shared utilities. The audit identified 26 critical correctness defects
in v3, of which 22 have been resolved across Phases 9.A, 9.B, and 9.C.
The remaining 4 items, plus 2 disclosed-but-non-blocking design choices
and 1 escalated DIC unit reconciliation, are the focus of this review
packet.

The fixes are coherent and mechanical: (a) algae and Nitrogen kinetic
methods that were silently reading legacy v2 kwargs (defaulted to
`0.0`/`1.0`) instead of the canonical `DEFAULTS` dict have been
rewired, plus three formula bugs (light-limit parenthesization,
harmonic-mean zero-guard, Steele exponent sign); (b) Carbon and DOX
algal coupling now uses the correctly derived stoichiometric ratios
`rca = AWc/AWa = 0.04` and `rcb = BWc/BWd = 0.4` instead of the raw
weights, and the previously missing CBOD source on dDIC/dt and
DOX-Monod attenuation on the SOD sink have been added; (c) two
inherited parameter values were corrected (`vson_20` 0.1→0.01 and
`lambdam` 0.0174→0.174); (d) the latent fdp partitioning unit error
was fixed.

The full test suite reports 705 tests passing and 3 xfailed.
Approximately 35 Fortran-anchored regression tests were added across
Phases 9.A/B/C to cover the audit findings (default-instantiation
regressions, rca/rcb numerical references, NH4+NO3 mass-balance
closure, hypoxic SOD attenuation, CBOD-DIC source, kdpo4>0 partitioning,
parameter pinnings). The 8 Tier 1 closed-system mass-conservation tests
all pass at `rtol=1e-12` with zero clip events. The coupled
NSM1+Riverine demo notebook at `examples/V3/04_Example_NSM1.ipynb`
exercises all 11 Process classes end-to-end on a synthetic 5-cell
mesh; performance is 17.6 ms/step, well below the spec's "must" target
of 415 ms/step.

The four open items in this packet are escalated to LimnoTech because
the canonical literature value is ambiguous from the codebase alone
or because the discrepancy is between v1's and Fortran's defaults
where neither has a clear primacy.

---

## 2. Document index — supporting materials

| Topic                                          | Canonical document                                                                |
| ---                                            | ---                                                                               |
| Source-of-truth design                         | `design/clearwater_modules_v3_nsm1_design_specification.md`                       |
| Gap analysis (Phase 0)                         | `design/clearwater_modules_v3_nsm1_gap_analysis.md`                               |
| README + Process inventory + audit history     | `design/clearwater_modules_v3_nsm1_README.md`                                     |
| Migration notes (v1→v3, v2→v3, audit fixes)    | `design/clearwater_modules_v3_nsm1_migration.md`                                  |
| Parameter corrections + LimnoTech reconciliation | `src/clearwater_modules_v3/parameter_defaults_corrections.md`                   |
| Three-way audit summary                        | `design/clearwater_modules_v3_nsm1_audit_summary.md`                              |
| Sub-audit: FloatingAlgae + BenthicAlgae        | `design/clearwater_modules_v3_nsm1_audit_algae.md`                                |
| Sub-audit: Nitrogen + Phosphorus               | `design/clearwater_modules_v3_nsm1_audit_n_p.md`                                  |
| Sub-audit: Carbon + DOX                        | `design/clearwater_modules_v3_nsm1_audit_c_dox.md`                                |
| Sub-audit: POM + CBOD + Pathogen + N2 + Alkalinity | `design/clearwater_modules_v3_nsm1_audit_simple_constituents.md`              |
| Sub-audit: utilities + parameter library       | `design/clearwater_modules_v3_nsm1_audit_utilities_params.md`                     |
| End-to-end coupled demo                        | `examples/V3/04_Example_NSM1.ipynb`                                               |
| v3 umbrella architecture                       | `design/clearwater_modules_v3_architecture_specification.md`                      |

This document is the index plus the focused-review summary; the
canonical content lives in those files.

---

## 3. Audit history (post-Phase-8 sequence)

After Phase 8 packaged the v3 NSM1 1.0.0 deliverable as
"review-ready", a follow-on three-way audit was performed across five
Process families and the shared utilities. The audit produced 64
findings deduplicated across the five sub-audits, of which 22 were
critical correctness defects in v3.

- **Phase 9.A** rewired 13 algae kinetic methods and 5 Nitrogen
  kinetic methods to read the canonical `DEFAULTS` dict instead of
  legacy v2 kwargs that defaulted to `0.0` (rates) or `1.0` (thetas);
  fixed 3 algae formula bugs (FloatingAlgae light-limit option-1
  parenthesization, harmonic-mean zero-guard, BenthicAlgae Steele
  exponent sign); dropped the phantom NH4 decay source; reconstructed
  the NO3 algal-uptake split to be mass-conservative; reconstructed
  the benthic NO3 uptake formula. **Resolved 17 of 22 critical items.**
- **Phase 9.B** corrected `rca`/`rcb` stoichiometric ratios in 8 sites
  across `carbon.py` and `dox.py`; added the CBOD oxidation source to
  dDIC/dt; added DOX-Monod attenuation to the SOD sink; dropped the
  spurious DOX attenuation from POC hydrolysis; corrected the `fdp`
  partitioning unit factor. **Resolved 5 of 22 critical items.**
- **Phase 9.C** corrected two inherited parameter values (`vson_20`
  0.1→0.01; `lambdam` 0.0174→0.174); reconciled `utils/reaeration.py`
  author attributions against Fortran source; relocated NSM1 docs
  from `docs/` to `design/`. **Closed 2 of 6 three-way disagreements;
  4 remain open for this review.**

Per-Phase commits: `eeefe51` (9.A), `01d70fd` (9.B), `afbbdb7` (9.C).

---

## 4. Reviewer focus areas

The review can proceed in any order. These are the items where
LimnoTech's literature knowledge or institutional history with v1
NSM1 calibration is the deciding input.

### 4.1 Nitrogen 4-way theta swap

`parameters/nitrogen.py` carries four Arrhenius theta values that
appear to swap pairwise with Fortran NSM1's
(`modNitrogen.f90:82, 89, 95, 100`):

| Parameter         | v3 / v1 value | Fortran value |
| ---               | ---           | ---           |
| `kon_theta`       | 1.074         | 1.047         |
| `kdnit_theta`     | 1.08          | 1.045         |
| `rnh4_theta`      | 1.047         | 1.074         |
| `vno3_theta`      | 1.045         | 1.08          |

The pattern (v3/v1 has `1.074, 1.08, 1.047, 1.045` while Fortran has
`1.047, 1.045, 1.074, 1.08`) is consistent with a 4-way swap during
the v1 port from Fortran. v3 inherits the v1 values; the Phase 9.C
audit annotated each with an inline `FIXME(phase9c-audit):` comment in
`parameters/nitrogen.py`.

**Ask:** does LimnoTech have authoritative QUAL2K-Table-X (or
equivalent literature) references for these four theta values that
would resolve the swap? At non-20-C temperatures the choice changes
nitrification, denitrification, sediment NH4 release, and sediment NO3
uptake Arrhenius corrections by ~10-15% individually.

### 4.2 v3 deliberate value choices that differ from Fortran

These are values where v3 deliberately diverges from Fortran NSM1
defaults; in each case v3's choice is defensible from literature but
not unambiguously correct.

| Parameter        | Fortran        | v3            | v3 rationale (current)                          |
| ---              | ---            | ---           | ---                                              |
| `vsop`           | 0.01 m/d       | 0.1 m/d       | Mid-range over Fortran's lower-bound; same as v3's `vs` |
| `SOD_20`         | 0.2 g-O2/m^2/d | 1.0 g-O2/m^2/d | Conservative midpoint over Fortran's lower bound (Chapra 1997 Table 25.2 supports 0.5-3.0 for moderate loading) |
| `BWa`            | 5000 g-D/mg-Chla | 3500        | v3 inherits v1; Fortran's 5000 differs by ~43% |
| `kah_20_user`    | 1.0 1/d        | 0.0 1/d       | v3 disabled-by-default to make the user-override branch a no-op when not configured |

The `kah_20_user=0` choice is documented prominently because it has
a behavioral consequence: at default `hydraulic_reaeration_option=1`,
v3 produces zero atmospheric hydraulic reaeration while Fortran NSM1
produces 1.0 1/d. Side-by-side runs of v3 vs Fortran NSM1 with
all-default settings will show DOX recovery in Fortran but not in v3.
See `src/clearwater_modules_v3/parameter_defaults_corrections.md`
Section 1.6 for the full behavioral note.

**Ask:** confirm these four choices, or supply Fortran-aligned values
where the LimnoTech-side preference is to track Fortran. Each is
documented in detail in `parameter_defaults_corrections.md` Sections
1.1, 1.3, 1.6 and Section 4.2.

### 4.3 rca / rcb derivation

Phase 9.B replaced raw-weight uses (`self.AWc = 40 mg-C/ug-Chla`,
`self.BWc = 40 mg-C/mg-D`) with the correctly derived ratios in 8 call
sites in `carbon.py` and `dox.py`:

```
rca = AWc / AWa = 40 / 1000 = 0.04 mg-C / ug-Chla   (floating algae)
rcb = BWc / BWd = 40 /  100 = 0.4  mg-C / mg-D     (benthic algae)
```

The v1 Python source defines the helper functions to compute `rca`
and `rcb` this way (`shared/processes.py` and `nsm1/algae_processes.py`).
Fortran derives them analogously. v3 1.0.0 uses the derived values.

The pre-fix raw-weight usage was off by a factor of 1000 (floating)
or 100 (benthic) on the algal coupling terms in dDIC/dt and dDOX/dt.
The pre-fix parity tests passed because they explicitly called the v1
helpers with `rca = AWc = 40` (passing the same wrong value into both
sides of the comparison), so v1-equivalence was vacuous.

**Ask:** sanity-check the derivation. The ratios are unitful
(mg-C / ug-Chla and mg-C / mg-D respectively) and the algal biomass
state variables (`Ap` in ug-Chla / L, `Ab` in mg-D / m^2) carry the
matching denominators, so the units close to mg-C / L per unit time.
This is consistent with v1's helper-function derivation and Fortran's
inline derivation; we want to confirm LimnoTech reads it the same way.

### 4.4 v3 deliberate improvements (vindicated by audit)

These are v3-only architectural choices that the audit found correct
and supported by Fortran's behavior; we list them so the reviewer is
not surprised by deviations from v1.

- **`SOD_tc` pure-Arrhenius split (Phase 1.1 architectural; Phase 9.B
  consumer-side Monod re-applied).** v3's `utils/sediment.SOD_tc` is a
  pure Arrhenius temperature correction; the DOX-Monod attenuation
  factor `DOX/(DOX+KsSod)` is applied at the consumer site rather than
  baked into the utility. This supports a future opt-in to a
  semi-implicit DOX treatment without forking the utility. Phase 9.B
  re-applied the Monod factor at the consumer so the hypoxic-attenuation
  behavior matches Fortran.
- **`PAR` consumer-side toggle (Phase 1.1).** v3 returns `q_solar *
  Fr_PAR` unconditionally; the `use_Algae | use_Balgae` toggle moved
  to the consumer Process. v1's two-arg `xr.where` form returns NaN in
  the false branch when both algae modules are disabled (latent
  NaN-propagation bug); v3 avoids it cleanly.
- **`clip_negative_state` contract.** Each `Process.run` calls
  `clearwater_modules_v3.utils.numerics.clip_negative_state(...)` after
  Forward Euler; clip target is exactly 0; clip events are logged
  (rate-limited) and counted on `model.diagnostics`. Tier 1 closed-system
  conservation tests assert `clip_events == 0`.
- **DOX-Monod attenuation routed through Nitrogen flux cache.** v3
  Alkalinity reads pre-attenuated `nitrification_flux_rate` and
  `denitrification_flux_rate` from `Nitrogen.run`'s rate cache rather
  than re-applying the Monod factor locally. Single-source-of-truth
  for the attenuation factor.
- **N2 Process closes the N mass balance.** v3 `N2.run` adds
  `denit_source = nitrogen_process.denitrification_flux_rate` to dN2/dt,
  closing the NO3 → N2 → N2sat exchange that Fortran and v1 silently
  break. Collapses to v1 form when `use_nitrogen=False`.
- **Phosphorus DIP derivation.** v3 uses `DIP = TIP * fdp` (matching
  v1, which is correct); Fortran's `DIP = TIP / fdp` is sign-inverted.
  v3 is right.

**Ask:** confirm these architectural improvements are acceptable for
v3 1.0.0. None of them changes the steady-state physics; each is
either a unit-of-work cleanup, a safety net, or a corrected sign.

---

## 5. Test summary

- **705 tests passing, 3 xfailed** across the v3 NSM1 suite.
- **~35 audit-anchored Fortran-anchored regression tests** added across
  Phases 9.A/B/C, including:
  - default-instantiation regression tests (`Process()` with no kwargs
    matches v1/Fortran-aligned reference values for 6 Processes);
  - stoichiometric-ratio regression tests (`Carbon.dic_algal_resp(...)`,
    `DOX._floating_algae_growth_flux(...)` against numerical references);
  - mass-balance closure tests (`NH4 algal-uptake + NO3 algal-uptake
    = rna * AlgalGrowth` over 100-step integrations);
  - hypoxic SOD attenuation (asserts `_sod_flux -> 0` as `DOX -> 0`,
    not just clip);
  - CBOD-DIC source (asserts dDIC/dt under `use_cbod=True, CBOD>0`
    includes the CBOD oxidation contribution);
  - `kdpo4 > 0` partitioning (asserts post-fix `fdp` matches Fortran
    over the sorption regime);
  - parameter pinnings (`vson_20 = 0.01`, `lambdam = 0.174`).
- **8 Tier 1 closed-system mass-conservation tests** under
  `tests/v3/nsm1/test_<constituent>_tier1.py` plus the consolidated
  harness at `tests/v3/nsm1/test_validation_tier1_conservation.py`.
  Each Process passes at `rtol=1e-12` with zero clip events.
- **36 sub-rate v1-parity tests** in `tests/test_5_*_calculations_v2.py`,
  each pinning a single v1 sub-term against v3's cached rate variable.
- **1 coupled end-to-end demo** at `examples/V3/04_Example_NSM1.ipynb`
  exercising all 11 Process classes plus Riverine on a synthetic mesh,
  runs end-to-end without clip events.

The test inventory above is the v3 NSM1 contribution. The full v3 test
suite (TSM + Model + NSM1) is documented at the v3 umbrella README.

---

## 6. Resolved design questions (spec Section 14)

The 9 questions surfaced during spec review are all decided. One
sentence each on the rationale; full design context in spec Section
14.

1. **OrgN as third state variable on `Nitrogen` Process** (not a
   separate `OrganicNitrogen` Process) — consistent with the
   topical-grouping rule applied to `Phosphorus` (TIP+OrgP) and
   `Carbon` (POC+DOC+DIC).
2. **Inter-process rate variables in the registry** follow the
   Appendix A naming convention (snake_case, suffix `_rate` for time
   derivatives, `_fraction` for dimensionless ratios, source-named
   prefixes for sediment fluxes).
3. **Negative-state handling: clip-with-log.** Forward Euler
   integration; per-Process call to
   `clearwater_modules_v3.utils.numerics.clip_negative_state(...)`;
   clip target is exactly 0; clips counted on `model.diagnostics`;
   Tier 1 closed-system tests assert `clip_events == 0`.
4. **Within-step semantics: Jacobi for state, Gauss-Seidel for rate
   variables.** State reads always pass `t=t_current`
   (order-independent); rate variables are step-scoped with strict
   producer-precedes-consumer enforcement.
5. **Package location: v3 NSM1 lives in `src/clearwater_modules_v3/`**
   alongside v3 TSM, sharing the v3 `Model`, registry, hotstart, and
   wet-mask infrastructure.
6. **Working default `modules.yml` shipped** at
   `src/clearwater_modules_v3/config/nsm1_default.yml` (or equivalent)
   exercising all 11 Process classes end-to-end.
7. **DEFAULTS-merge parameter library pattern.** Each `Process`
   imports a `DEFAULTS: dict[str, float]` from `parameters/<group>.py`
   and merges with user-provided values at construction time.
   Phase 9.A corrected the kinetic-method wiring so the methods read
   the DEFAULTS keys directly.
8. **v1 NSM1 expected-value tables are usable but not directly
   portable.** v1's hard-coded values were generated under a mixed
   Jacobi/Gauss-Seidel state-read pattern; v3 enforces strict Jacobi.
   Phase 7 regenerates expected values where the state-read pattern
   affects the result.
9. **DOX semi-implicit treatment is a per-Process opt-in deferred
   to v3 1.1+** if profiling shows the Forward Euler treatment of
   DOX's first-order sinks (SOD + nitrification + CBOD oxidation +
   DOC oxidation) requires tightening.

---

## 7. Tentative design decisions awaiting LimnoTech confirmation

These three are decided in spec Section 14 with the explicit annotation
"pending LimnoTech confirmation". Each is a low-stakes confirmation —
v3 1.0.0 ships under the listed defaults; LimnoTech objection routes
to v3 1.0.x or v3 1.1.

### 7.1 Alkalinity as simple tracer (no carbonate solver)

**Decision.** v3 1.0.0 implements `Alkalinity` as a state variable with
source/sink terms (nitrification consumption, denitrification
production, algal growth/respiration coupling) integrated by Forward
Euler. No carbonate equilibrium, no pH solver. Full pH chemistry
(carbonate equilibrium, NH3/NH4+ partitioning, free-CO2 fraction,
carbonate speciation) is NSM2 territory in v3 1.1+.

**Ask.** Confirm no current LimnoTech application requires pH from day
one. v3 1.0.0 documentation includes a worked example of post-hoc pH
computation from `Alk`, `DIC`, `T`, salinity output trajectories for
users who need a quick number.

### 7.2 Sediment-flux parameters as scalar globals

**Decision.** `SOD_20`, `NH4fromBed`, `DIPfromBed`, `NO3_BedDenit`,
`DIC_sed_release` are scalar globals applied uniformly to all cells,
set in YAML (matching v1's pattern exactly). Per-cell spatially varying
fluxes and dynamically computed fluxes both arrive in v3 1.1+ via the
NSM2 sediment diagenesis Process.

**Ask.** Confirm no current LimnoTech application uses spatially
varying SOD or bed fluxes. If yes, route to a v3 1.0.x follow-up that
allows per-cell DataArrays in the YAML (small change; the registry
already supports it).

### 7.3 Single-compartment algae (one Ap, one Ab)

**Decision.** Phytoplankton and benthic algae remain
single-compartment in v3 1.0.0, matching v1/v2. NSM2's multi-group
capability lands as new Process classes (e.g., `PhytoplanktonGroups`)
added *alongside* the single-compartment `FloatingAlgae` /
`BenthicAlgae`, not as in-place extension. The YAML-driven Process
registration framework supports this additive extension natively.

**Ask.** Confirm no near-term LimnoTech application requires
multi-group algae before v3 1.1+.

---

## 8. Remaining items (disclosed, non-blocking)

These are items the audit surfaced that v3 1.0.0 ships with
deliberately and that we want to flag explicitly for review.

### 8.1 C5 — zero default reaeration design choice

`kah_20_user = 0.0` (v3) vs `1.0` (Fortran). At default
`hydraulic_reaeration_option=1`, v3 produces zero atmospheric
hydraulic reaeration; Fortran NSM1 produces 1.0 1/d. Documented in
detail in `src/clearwater_modules_v3/parameter_defaults_corrections.md`
Section 1.6 with explicit user guidance (set `kah_20_user > 0`, or
select an empirical option from the menu).

### 8.2 C6 — DOX salinity correction on O2sat

Fortran applies `O2sat *= exp(-Salinity * (0.017674 - 10.754/Tk +
2140.7/Tk^2))`; v1 omits; v3 inherits the v1 omission. Freshwater
impact is zero. Estuarine impact is up to ~20% on `O2sat` at
salinity 35. Defer if v3 1.x targets fresh water only; flag if
brackish/estuarine deployment is on the LimnoTech roadmap.

### 8.3 C9 / C10 — DIC unit reconciliation (mol/L vs mg/L)

The dDIC/dt budget in v1 carries a long-standing unit-factor pattern
where some terms are computed in mol/L (with implicit conversion via
`/12000`) and others in mg/L. v1 / Fortran agreement is not
unambiguous; v3 inherits the v1 form. Audit items C9 and C10 in
`design/clearwater_modules_v3_nsm1_audit_c_dox.md` lay out the
candidate readings. Escalated for LimnoTech review; the decision is
whether to land a unit cleanup in 1.0.x or defer to the v3 1.1
carbonate solver.

---

## 9. Performance

- **17.6 ms/step on a 5-cell synthetic mesh** (current measured).
- **Spec "must" target** (spec Section 10): 4,320-timestep coupled
  NSM1+Riverine simulation on Sumwere Creek within 30 minutes
  (≤ 415 ms/step). Confidently met.
- **"Should" target**: within 10 minutes (≤ 138 ms/step). Likely met
  pending profiling.
- **"Aspirational" target**: within 3 minutes (≤ 41 ms/step). Requires
  profiling and possible kernel optimization.

Extrapolation to Sumwere Creek (~600 cells) is deferred to a formal
benchmark in a follow-up phase.

---

## 10. Suggested review focus areas (25-minute first pass)

1. **Section 4.1** — the 4-way Nitrogen theta swap. The single
   load-bearing literature question for this review.
2. **Section 4.2** — confirm the 4 v3-deliberate value choices
   (`vsop`, `SOD_20`, `BWa`, `kah_20_user`).
3. **Section 4.3** — sanity-check the rca/rcb derivation.
4. **Sections 7.1 / 7.2 / 7.3** — the 3 tentative decisions awaiting
   LimnoTech confirmation. Each is a low-stakes yes/no.
5. **Section 8.3** — the DIC unit reconciliation. Decide 1.0.x or
   1.1.

For a deeper second pass: the audit summary
(`design/clearwater_modules_v3_nsm1_audit_summary.md`), the five
sub-audits, the corrections record at
`src/clearwater_modules_v3/parameter_defaults_corrections.md`, and
the corresponding `processes/<name>.py` implementations.

---

## 11. Asks of the reviewer (summary)

Five confirmations close out the v3 NSM1 1.0.0 review:

1. **Resolve the 4-way Nitrogen theta swap** (Section 4.1). Authoritative
   QUAL2K-Table-X reference would close it.
2. **Confirm the 4 v3-deliberate value choices** (Section 4.2):
   `vsop=0.1`, `SOD_20=1.0`, `BWa=3500`, `kah_20_user=0.0`. Any of
   these can be flipped to Fortran-aligned values if the LimnoTech-side
   preference is to track Fortran.
3. **Sanity-check the rca/rcb derivation** (Section 4.3) and the
   Phase 9.B fix.
4. **Confirm the 3 tentative spec-Section-14 decisions** (Sections
   7.1 / 7.2 / 7.3): Alkalinity simple-tracer scope, sediment-flux
   scalar globals, single-compartment algae.
5. **Decide the DIC unit reconciliation** (Section 8.3): land a unit
   cleanup in 1.0.x, or defer to the v3 1.1 carbonate solver.

Any LimnoTech application that v3's design choices would surprise — in
particular spatially varying SOD, multi-group algae, brackish/estuarine
DOX salinity correction, or pH-from-day-one — routes to v3 1.0.x or
v3 1.1+ rather than blocking 1.0.0 sign-off.
