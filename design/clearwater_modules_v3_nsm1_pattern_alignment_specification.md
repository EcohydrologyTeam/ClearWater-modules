# ClearWater Modules v3 NSM1 — Pattern Alignment Specification

**Status:** Approved for implementation (2026-05-13)
**Author:** Todd Steissberg (ERDC)
**Date:** 2026-05-13
**Scope:** Retro-apply the v3 TSM `Temperature.run()` design patterns to every v3 NSM1 `Process`, including all tracers. Restores conformance with `clearwater_modules_v3_nsm1_design_specification.md` Section 14 (registry rate-variable convention, clip-with-log contract, Jacobi/GS substep semantics) and harmonises all 11 Processes against a single canonical run-method shape.

**Read this with:**

- `clearwater_modules_v3_nsm1_design_specification.md` — the 1.0.0 spec; Section 14 design decisions and Appendix A registry-coupling cheatsheet are load-bearing inputs to this work.
- `clearwater_modules_v3_architecture_specification.md` — umbrella v3 architecture; chunking, hotstart, and wet-mask contracts.
- `clearwater_modules_v3_tsm_design_specification.md` — TSM is the canonical exemplar referenced throughout this spec.
- `src/clearwater_modules_v3/processes/temperature.py` — the reference implementation that defines patterns A–J below.

---

## 1. Motivation

The v3 NSM1 1.0.0 audit packet (`clearwater_modules_v3_nsm1_limnotech_review.md`) shipped with the kinetics-correctness items resolved across Phases 9.A–9.C. A subsequent structural-conformance pass against the v3 TSM exemplar identified that the eleven NSM1 `Process` classes implement the v2 LimnoTech *shape* of `run()` (read forcings, compute rate, integrate, persist) but do not uniformly carry the v3 TSM polish that the NSM1 1.0.0 design spec assumed:

- The Appendix A registry rate-variable convention exists in code as `self.<rate>` instance attributes but is not exposed to the registry except in `N2.run` (one of eleven).
- The clip-with-log contract is implemented inconsistently — five Processes apply it unconditionally; four wrap it in `isinstance(... DataArray) and self.diagnostics is not None` guards with `xr.where(state < 0, 0, ...)` fallbacks.
- The fused `(delta, components)` helper that v3 TSM uses (`_temperature_change_with_factors`) — which is the single point of truth for both the integrator update and the per-component diagnostics — has no NSM1 analogue; sub-fluxes are composed inline in `run` and component diagnostics are not exposed.
- Auxiliary inconsistencies across modules (NaN-scrub helper used in 6/11; clip-import location varies; missing-input idioms differ; one module re-invokes a kinetic helper redundantly; another caches in the kinetic helper instead of `run`).

The eleven Processes are coherent and pass 705 tests, but the surface area they expose for calibration, downstream coupling, validation, and post-processing is not what the design spec called for. This specification defines the retrofit.

---

## 2. Goals and Non-Goals

### Goals

1. Every NSM1 `Process` exposes a fused `(delta, components)` helper analogous to TSM's `_temperature_change_with_factors`, returning the integrator-ready delta(s) and a `components: dict` of all sub-rate/sub-flux diagnostics consumed by the integrator step.
2. Every NSM1 `Process` writes every Appendix A rate variable it produces to the registry **opportunistically**: cached on `self.<name>` for sibling consumption, and written via `set_at_time(name, time, value)` **only when the variable is pre-registered**.
3. Every NSM1 `Process` uses unconditional clip-with-log via `clip_negative_state(...)`; the `isinstance/diagnostics-not-None` guard branches are removed.
4. Auxiliary inconsistencies (NaN-scrub style, import location, sibling discovery in `init_process`, cache-write location, redundant kinetic-helper invocations) are resolved with a single canonical convention applied uniformly.
5. **ZERO REGRESSIONS** — the kinetics outputs are bit-identical for any run that does not pre-register the new diagnostic variables. This is non-negotiable, gated at every phase commit, and the load-bearing invariant of this work. See §12 for the full contract.
6. The Appendix A registry-coupling cheatsheet in `clearwater_modules_v3_nsm1_design_specification.md` is the authoritative names list, and every name in Appendix A maps to exactly one (Process, cached attribute, opportunistic-write target).

### Non-Goals

1. No kinetics-formula changes. This is a structural refactor; numerical outputs are unchanged for unchanged registry configurations.
2. No new parameters, no new `DEFAULTS` keys, no parameter renames.
3. No retroactive changes to v1 or v2 NSM1 source.
4. No changes to the Forward Euler integrator contract or the Jacobi/GS substep semantics.
5. No introduction of NSM2 features (multi-pool organic matter, carbonate solver, methane-sulfide, silica, sediment diagenesis).
6. No changes to the `Model`-layer wet-mask, hotstart, or chunking — those are the architecture spec's concern, not this one.

---

## 3. The Canonical Pattern (A–J)

The TSM exemplar (`temperature.py:410-583`) is the canonical implementation. Each pattern below has a one-line statement, the TSM reference, and the NSM1 generalisation that every Process must satisfy after this work.

### Pattern A — Read all state/forcing at top via `get_at_time`

- **Statement:** `run` reads each state variable and each forcing variable from the registry once, near the top of the method, into a local name. No mid-method `get_at_time` calls.
- **TSM:** `temperature.py:417-443`.
- **NSM1 rule:** A single contiguous block of `get_at_time` calls at the top of `run`, one per name in `self.variables` plus any optional inputs. Optional inputs use the `if "<name>" in registry: ... else: <fallback>` idiom (see pattern G note for the consistent membership-check form).

### Pattern B — Fused `(delta, components)` helper

- **Statement:** `run` delegates to a single helper `_change_with_components(...)` (per-Process convention; TSM uses `_temperature_change_with_factors`) that returns a tuple `(delta_state_1, ..., delta_state_N, components: dict)`. The helper computes every sub-rate/sub-flux exactly once and returns them in `components` for caching and opportunistic registry exposure.
- **TSM:** `temperature.py:1186-1289` (`_temperature_change_with_factors`).
- **NSM1 rule:** `run` calls `self._change_with_components(...)` exactly once. The integrator step uses `delta_state_*` directly; sibling-process caches (pattern F) and opportunistic registry writes (pattern G) read from `components`. No sub-rate is recomputed.
- **NSM1 naming convention:** helper name is `_change_with_components` for all multi-state processes; for the rare single-state Process (`Pathogen`), the helper name is `_rate_with_components` since the integrator already uses a per-second rate. Both forms return `(integrator_argument..., components: dict)`.

### Pattern C — Forward Euler additive update

- **Statement:** Integration is `state_new = state_old + delta` (delta-form) or `state_new = state_old + rate * dt_days` (rate-form). No multiplicative `state_old * rate * dt` forms anywhere.
- **TSM:** `temperature.py:484`.
- **NSM1 rule:** Whichever form (delta or rate-times-dt) is in use today is preserved. This is C/contract, already satisfied everywhere post Phase 2.B.

### Pattern D — Unconditional clip-with-log

- **Statement:** Every state update is followed by `state_new = clip_negative_state(state_new, name, self.diagnostics, step=step_index)`. No conditional fallbacks.
- **TSM:** Implicit at the `run` level — TSM's water_temperature is energy-balanced and naturally bounded; an explicit clip is not used but the Tier 1 conservation test asserts the equivalent invariant.
- **NSM1 rule:** Every state-variable assignment is wrapped in `clip_negative_state`. The `isinstance(... xr.DataArray) and self.diagnostics is not None` branches in Alkalinity (`alkalinity.py:485-490`), DOX (`dox.py:730-735`), N2 (`n2.py:367-370`), and Nitrogen's `_clip` wrapper (`nitrogen.py:371-376`) are removed. `clip_negative_state` is hardened to no-op gracefully on `None` diagnostics and on scalar/ndarray inputs so the guard is unnecessary at every call site.
- **Out-of-scope clarification:** `Pathogen` already uses unconditional `clip_negative_state`; no change.

### Pattern E — `set_at_time` for primary outputs

- **Statement:** Every state variable updated by the Process is persisted via `registry.set_at_time(name, time, state_new)`.
- **TSM:** `temperature.py:486, 535-537`.
- **NSM1 rule:** Already in place for every Process post Phase 2.B (Bug #16 fix). This pattern is verified, not changed.

### Pattern F — Cache step-scoped rates on `self.<name>`

- **Statement:** Every sub-rate/sub-flux that any sibling Process consumes is assigned to a `self.<name>` attribute inside `run` (not inside a kinetic helper). The cache is step-scoped: overwritten each substep, no persistence across substeps.
- **TSM:** `temperature.py:544-550`.
- **NSM1 rule:** Cache writes live in `run`, immediately after the call to `_change_with_components`. Cache names match Appendix A exactly. `POM` is the one current divergence (`pom.py:280` writes inside `rate()`); the cache assignment moves to `run`.

### Pattern G — Opportunistic registry writes

- **Statement:** After cache writes, iterate the Process's registry-exposed names and write each one with `set_at_time` **only if** the name is currently registered. Users who pre-register a name receive it on their output; users who don't pay no cost.
- **TSM:** `temperature.py:551-563`.
- **NSM1 rule:** Each Process declares a class-level `REGISTRY_DIAGNOSTICS: tuple[str, ...]` listing the Appendix A names it owns. `run` ends with:
  ```python
  for name in self.REGISTRY_DIAGNOSTICS:
      if name in registry:
          registry.set_at_time(name, time, getattr(self, name))
  ```
  This is the single mechanical change that makes Appendix A real on the output side.

### Pattern H — No per-process wet-mask

- **Statement:** No per-process `xr.where(volume > 0, ...)` guard. Wet/dry is the `Model.__apply_wet_mask` layer's responsibility.
- **TSM:** `temperature.py:466-485` (commentary documenting the removal).
- **NSM1 rule:** Audit confirms all 11 Processes are already clean. Verified, not changed.

### Pattern I — Magnitudes-only sub-functions

- **Statement:** Kinetic sub-functions return positive magnitudes (or signed-by-physical-gradient values where physics dictates, e.g., `T_air − T_water`); signs are applied at composition time inside `_change_with_components`.
- **TSM:** `temperature.py:585-607` (the energy-balance section header).
- **NSM1 rule:** Audit confirms all 11 Processes are already aligned. Verified, not changed.

### Pattern J — `init_process` for sibling discovery + diagnostics capture

- **Statement:** `init_process(model, registry)` probes for optional sibling Processes via `model.has_process("<Class>")` and caches handles via `model.get_process("<Class>")`. Sets `self.use_<sibling> = bool(...)`. Captures `self.diagnostics = getattr(model, "diagnostics", None)`.
- **NSM1 rule:** Carbon, DOX, Nitrogen, Phosphorus, FloatingAlgae, POM, N2, and Pathogen already follow this. `BenthicAlgae` hard-codes `self.use_* = True` at `benthic_algae.py:222-224`; this must be replaced with `has_process`/`get_process` discovery to match `FloatingAlgae`'s `init_process`. `Alkalinity` and `CBOD` capture `self.diagnostics` but do not discover siblings; they should adopt the same `has_process` / `get_process` pattern for any Process whose rates they consume (the audit table in §5 enumerates).

### Auxiliary conventions (apply uniformly post-refactor)

- **NaN scrub:** `sanitize_rate(value)` from `utils.numerics` is the canonical helper (used by Alkalinity, Carbon, DOX, N2, Pathogen, Phosphorus today). Applied once at the *components-dict* level inside `_change_with_components`, not piecemeal at each sub-flux. FloatingAlgae, POM, BenthicAlgae, CBOD adopt it; CBOD's inline `xr.where(...isnull(), 0, ...)` is removed.
- **Import location for `clip_negative_state`:** module-level at top of file. FloatingAlgae (`floating_algae.py:422`) and POM (`pom.py:219`) currently use local-in-function imports; lift to module level.
- **Missing-input idiom:** `if "<name>" in registry:` followed by `get_at_time` else a `xr.zeros_like(<reference>)` fallback. CBOD's `try/except KeyError` block (`cbod.py:206-216`) is replaced with the membership-check form. Pathogen's `_get_optional` helper-with-warn-once-latch is also acceptable for variables where the warning is useful; document the choice per Process in code.
- **Step-index argument to `clip_negative_state`:** the substep counter from the Process or `model.current_step`; the placeholder `step=0` currently used in `Nitrogen._clip` is replaced with the real value (resolution in Phase 0; depends on how Model exposes the step counter — see §10 Open Questions).

---

## 4. Per-Process Component Inventory

This table is the load-bearing artifact of this spec. Each row defines (a) which Appendix A names a Process owns, (b) what its `_change_with_components` returns, and (c) which sibling caches it reads. Rows are ordered by phase (heavy integrators first, tracers last) to match §6.

| Process | State variables | Owned registry diagnostics (Appendix A) | `_change_with_components` returns | Sibling cache reads |
|---|---|---|---|---|
| **Carbon** | POC, DOC, DIC | `poc_hydrolysis_rate`, `doc_dic_oxidation_rate`, `dic_atm_exchange_rate`, `dic_sed_release_rate`, `carbon_algal_resp_rate`, `carbon_balgae_resp_rate`, `carbon_algal_photo_rate`, `carbon_balgae_photo_rate`, `carbon_cbod_oxidation_rate` | `(d_poc, d_doc, d_dic, components: dict)` | `floating_algae.algal_growth_rate`, `floating_algae.algal_respiration_rate`, `benthic_algae.balgae_growth_rate`, `benthic_algae.balgae_respiration_rate`, `cbod.cbod_oxidation_rate` |
| **DOX** | DOX | `dox_sat`, `atm_reaeration_rate`, `dox_nitrification_rate`, `dox_sod_rate`, `dox_doc_oxidation_rate`, `dox_cbod_oxidation_rate`, `dox_algal_photo_rate`, `dox_algal_resp_rate`, `dox_balgae_photo_rate`, `dox_balgae_resp_rate` | `(delta_dox, components: dict)` | `nitrogen.nitrification_flux_rate`, `carbon.doc_dic_oxidation_rate`, `cbod.cbod_oxidation_rate`, `floating_algae.algal_growth_rate / algal_respiration_rate`, `benthic_algae.balgae_growth_rate / balgae_respiration_rate`, `floating_algae.algal_nh4_uptake_fraction`, `benthic_algae.balgae_nh4_uptake_fraction` |
| **Nitrogen** | NH4, NO3, OrgN | `nitrification_flux_rate`, `denitrification_flux_rate`, `nh4_from_bed`, `no3_from_bed_denit`, `orgn_hydrolysis_rate`, `orgn_settling_rate`, `nh4_algal_growth_rate`, `no3_algal_growth_rate`, `nh4_algal_resp_rate`, `nh4_balgae_resp_rate` | `(d_nh4, d_no3, d_orgn, components: dict)` | `floating_algae.algal_growth_rate / algal_respiration_rate / algal_orgn_from_mortality_rate / algal_nh4_uptake_fraction`, `benthic_algae.balgae_*` |
| **FloatingAlgae** | Ap (chlorophyll) | `algal_growth_rate`, `algal_respiration_rate`, `algal_death_rate`, `algal_settling_rate`, `algal_orgn_from_mortality_rate`, `algal_orgp_from_mortality_rate`, `algal_poc_from_mortality_rate`, `algal_doc_from_mortality_rate`, `algal_pom_from_settling_rate`, `algal_nh4_uptake_fraction`, `algal_light_limitation`, `algal_nutrient_limitation_n`, `algal_nutrient_limitation_p` | `(d_ap, components: dict)` | `nitrogen.ammonium / nitrate` via registry; `phosphorus.tip` via registry; no sibling cache reads |
| **BenthicAlgae** | Ab (dry mass per area) | `balgae_growth_rate`, `balgae_respiration_rate`, `balgae_death_rate`, `balgae_orgn_from_mortality_rate`, `balgae_orgp_from_mortality_rate`, `balgae_poc_from_mortality_rate`, `balgae_doc_from_mortality_rate`, `balgae_nh4_uptake_fraction`, `balgae_light_limitation`, `balgae_nutrient_limitation_n`, `balgae_nutrient_limitation_p` | `(d_ab, components: dict)` | same as FloatingAlgae |
| **Phosphorus** | TIP, OrgP | `orgp_hydrolysis_rate`, `orgp_settling_rate`, `tip_settling_rate`, `dip_from_bed`, `orgp_algal_mortality_rate`, `tip_algal_growth_rate`, `tip_balgae_growth_rate` | `(d_tip, d_orgp, components: dict)` | `floating_algae.algal_growth_rate / algal_orgp_from_mortality_rate`, `benthic_algae.balgae_*` |
| **POM** | POM | `pom_hydrolysis_rate`, `pom_settling_rate`, `pom_algal_mortality_rate`, `pom_balgae_mortality_rate` | `(d_pom, components: dict)` | `floating_algae.algal_pom_from_settling_rate`, `benthic_algae.balgae_pom_from_mortality_rate` (if applicable) |
| **CBOD** | CBOD (per group) | `cbod_oxidation_rate`, `cbod_settling_rate` (per group; sums optional) | `(d_cbod, components: dict)` | none |
| **N2** | N2 | `n2_atm_exchange_rate`, `n2_sat`, `total_dissolved_gas` (existing); add `n2_denit_source_rate` | `(delta_n2, components: dict)` | `nitrogen.denitrification_flux_rate` |
| **Pathogen** | PX | `pathogen_natural_death_rate`, `pathogen_light_death_rate`, `pathogen_settling_rate` | `(rate_px, components: dict)` (rate-form integrator) | none |
| **Alkalinity** | Alk | `alk_nitrification_sink_rate`, `alk_denitrification_source_rate`, `alk_algal_growth_rate`, `alk_algal_respiration_rate`, `alk_balgae_growth_rate`, `alk_balgae_respiration_rate` | `(delta_alk, components: dict)` | `nitrogen.nitrification_flux_rate / denitrification_flux_rate`, `floating_algae.algal_growth_rate / algal_respiration_rate / algal_nh4_uptake_fraction`, `benthic_algae.balgae_*` |

**Naming convention reminder** (Appendix A, design spec §14 resolution):

- snake_case throughout.
- Suffix `_rate` for time-derivative quantities (mg/L/d, ug-Chla/L/d, etc.).
- Suffix `_fraction` for dimensionless ratios.
- Source-named prefixes for sediment fluxes (e.g., `nh4_from_bed`).
- Process-name prefix on diagnostics that disambiguate consumer attribution (`dox_nitrification_rate` is "the contribution of nitrification to DOX", computed inside DOX; `nitrification_flux_rate` is the bare nitrification flux computed inside Nitrogen).

**Cross-reference rule:** any name in the table above that does not appear in design spec Appendix A is a candidate Appendix A addition. Phase 0 (§6) produces the diff and amends Appendix A.

---

## 5. Per-Process Refactor Plan

Each Process gets the same structural treatment; the table below records the per-Process delta against the current code.

| Process | A | B | D | F | G | J | Aux fixes |
|---|---|---|---|---|---|---|---|
| **Carbon** | clean | **add** `_change_with_components` | clean | clean | **add** loop over `REGISTRY_DIAGNOSTICS` | clean | none |
| **DOX** | clean | **add** | **harmonise** (drop guard) | clean | **add** | clean | none |
| **Nitrogen** | clean | **add** | **harmonise** (replace `_clip` with direct `clip_negative_state`) | clean | **add** | clean | resolve step-index placeholder |
| **FloatingAlgae** | clean | **add** | clean | clean | **add** | clean | lift `clip_negative_state` import to module level; adopt `sanitize_rate` |
| **BenthicAlgae** | clean | **add** | clean | **move cache writes from `_cache_benthic_mortality_rates` to `run`; drop the redundant second `rate_death` invocation** | **add** | **fix sibling discovery** (replace hardcoded `use_* = True`) | adopt `sanitize_rate`; document `init_process` change |
| **Phosphorus** | clean | **add** | clean | clean | **add** | clean | none |
| **POM** | clean | **add** | clean | **move cache write from `rate()` to `run`** | **add** | clean | lift `clip_negative_state` import; adopt `sanitize_rate` |
| **CBOD** | clean | **add** | clean | clean | **add** | **add** sibling discovery for DOX (currently `try/except`) | replace `try/except KeyError` with `if "<name>" in registry`; replace inline `xr.where(isnull, 0, ...)` with `sanitize_rate` |
| **N2** | clean | **add** | **harmonise** (drop guard) | clean | clean (already present for `total_dissolved_gas`); extend loop to cover full Appendix A set | clean | none |
| **Pathogen** | clean | **add** `_rate_with_components` | clean | n/a (no siblings consume) | **add** (writes the three death/settling rates) | clean | none |
| **Alkalinity** | clean | **add** | **harmonise** (drop guard) | clean | **add** | **add** sibling discovery for Nitrogen, FloatingAlgae, BenthicAlgae | drop the locally-defined `sanitize_rate` if it shadows utils; use the canonical one |

Patterns C (additive integrator), E (`set_at_time` persistence), H (no wet-mask), I (magnitudes-only) require no changes — they are verified as conforming by the audit and re-verified by Phase 5 regression tests.

---

## 6. Phased Implementation Plan

The phasing is heavy integrators first (Carbon, DOX, Nitrogen) — those are the ones whose component diagnostics calibration users will actually subscribe to — then algae and Phosphorus, then tracers. Each phase has a deliverable, a test gate, and a "no-cost-when-unused" regression check.

### Phase 0 — Spec acceptance, prep, baseline capture, Appendix A diff (1 day)

- Reviewer accepts this spec.
- **Baseline capture** (load-bearing for the zero-regression contract; see §12):
  - Run the full v3 test suite on the pre-refactor branch; save the verbose pytest output, exact numeric assertion values, and per-test wall-clock to `tests/v3/nsm1/baseline_<commit>.txt` and `tests/v3/nsm1/baseline_<commit>.json`.
  - Run the coupled NSM1+Riverine demo for 4,320 steps with a fixed RNG seed; save the full state-variable trajectories (every substep, every cell, every state variable) as a NetCDF dataset at `tests/v3/nsm1/baseline_coupled_trajectory_<commit>.nc`. This is the *gold reference* against which every subsequent phase's outputs are compared bit-for-bit.
  - Run each Process's Tier 1 closed-system conservation test; record the achieved tolerance and clip-event count for each. Save as `tests/v3/nsm1/baseline_tier1_<commit>.json`.
  - Capture `pip freeze` / `pixi list` to pin the exact dependency versions used for the baseline. Any phase that changes a dependency invalidates the baseline and must re-capture.
- Diff §4's component inventory against `clearwater_modules_v3_nsm1_design_specification.md` Appendix A; amend Appendix A with the new names introduced here (most are new; the Appendix A list was scoped to inter-process *coupling* names, not the broader diagnostic set).
- Resolve the open question on `clip_negative_state` step-index source (§10).
- Resolve the open question on `clip_negative_state` graceful-no-op behaviour when diagnostics is None (§10).
- Decide: utility helpers (`sanitize_rate`, `clip_negative_state`) finalised at this point; no churn in later phases.

**Deliverable:** spec marked approved; baseline artifacts committed under `tests/v3/nsm1/baseline_*`; Appendix A amended; `utils/numerics` finalised (clip-with-log no-op contract + step-index handling).

### Phase 1 — Mechanical alignment pass across all 11 Processes (1 day)

Applies the §3 auxiliary conventions and the §5 column-D harmonisation, without touching kinetics:

- Harmonise pattern D: replace `_clip` wrapper / inline-guard variants with unconditional `clip_negative_state` (Alkalinity, DOX, N2, Nitrogen).
- Lift `clip_negative_state` imports to module level (FloatingAlgae, POM).
- Adopt `sanitize_rate` uniformly; remove CBOD's inline NaN scrub.
- Fix BenthicAlgae `init_process` sibling discovery.
- Fix BenthicAlgae redundant `rate_death` re-invocation.
- Move POM cache write from `rate()` to `run()`.
- Replace CBOD `try/except KeyError` DOX read with `if "<name>" in registry` form.
- Add `init_process` sibling discovery to Alkalinity and CBOD.

**Test gate (§12 zero-regression contract enforced):**
- Full test suite passes (`pytest tests/v3 -x`). No skipped tests beyond the pre-refactor xfailed set.
- `test_baseline_parity.py` (added Phase 0) passes: every state variable in the 4,320-step coupled trajectory matches the baseline NetCDF *bit-identically* (`xr.testing.assert_identical`).
- Per-test numeric assertions match the baseline JSON within ULP-level tolerance (`rtol=0, atol=0` where the pre-refactor tolerance was already `0`; pre-refactor tolerance otherwise).
- Tier 1 conservation: tolerance achieved and clip-event count match the baseline JSON exactly.
- BenthicAlgae sibling-discovery change: review every test under `tests/v3/nsm1/test_benthic_*.py` for setup assumptions; if any test relied on the old hardcoded behaviour, the test is the regression — investigate, do not "fix" the test until the regression is understood.

**Deliverable:** all 11 Processes pass pattern D + J + auxiliary; commit gated by every clause of the test gate.

### Phase 2 — Carbon `_change_with_components` (1.5 days)

- Add `Carbon._change_with_components(...)` returning `(d_poc, d_doc, d_dic, components)`.
- Add `Carbon.REGISTRY_DIAGNOSTICS` class attribute.
- Rewire `Carbon.run` to call the helper, cache the components, and run the opportunistic-write loop.
- **Refactor rule for this phase (applies to every Process refactor through Phase 9):** the legacy inline composition must move into `_change_with_components` *without rearranging operand order, replacing intermediate variables, or "simplifying" expressions*. The diff must be code motion, not code rewrite. Any algebraic simplification, even one that is mathematically identical, requires a separate commit and a re-run of the full §12 contract.

**Test gate (§12 zero-regression contract enforced):**
- Full test suite passes.
- `test_baseline_parity.py` passes bit-identically against the Phase 0 baseline NetCDF, with no Carbon `REGISTRY_DIAGNOSTICS` names pre-registered.
- **`test_carbon_helper_vs_inline.py`** (new, deleted after Phase 10): the legacy inline composition is preserved as `_change_legacy_inline` for the duration of Phases 2–9 and the test asserts `_change_with_components` returns bit-identical `(d_poc, d_doc, d_dic)` for a parametrised matrix of state/forcing inputs covering: zero state; uniform state; randomised state with fixed seed; edge cases (zero DOX, zero temperature, very high temperature, very thin depth). Tolerance: `rtol=0, atol=0`.
- **`test_carbon_registry_diagnostics.py`** (new, retained): when a `REGISTRY_DIAGNOSTICS` name is pre-registered, it is written each substep with finite values; when none are pre-registered, the state-variable trajectory is bit-identical to the no-diagnostics-subscribed baseline.
- Tier 1 carbon conservation: tolerance and clip-event count match baseline exactly.

**Deliverable:** Carbon conformant; helper-vs-inline parity verified; zero-change-when-not-used invariant verified at the substep level.

### Phases 3–9 — Per-Process refactor (each phase follows the Phase 2 template)

Every per-Process phase below follows the **identical structure and test gate** defined in Phase 2:

1. Add `_change_with_components` (or `_rate_with_components` for Pathogen) returning `(delta(s)..., components)`.
2. Add `REGISTRY_DIAGNOSTICS` class attribute.
3. Rewire `run` to call the helper, cache the components, and run the opportunistic-write loop.
4. **Code-motion-only refactor rule** (Phase 2 statement applies verbatim).
5. **§12 zero-regression test gate** applies (full suite + baseline parity NetCDF + `test_<process>_helper_vs_inline.py` with `rtol=0, atol=0` + diagnostics-not-subscribed bit-identical state-trajectory + Tier 1 conservation tolerance/clip-count match).

### Phase 3 — DOX `_change_with_components` (1.5 days)

DOX is the second-heaviest integrator (eight sub-fluxes) and has the most calibration-interesting per-component diagnostics. The helper-vs-inline parity matrix in `test_dox_helper_vs_inline.py` must also cover the DOX-Monod attenuation regime (`DOX → 0`) and the hypoxic regime where SOD attenuation kicks in.

**Deliverable:** DOX conformant; Phase 2 test gate satisfied; hypoxic-attenuation parity verified specifically.

### Phase 4 — Nitrogen `_change_with_components` (2 days)

Three state variables (NH4, NO3, OrgN), and Nitrogen owns the nitrification/denitrification flux caches that DOX and Alkalinity consume. The helper must produce those caches in `components` while preserving the existing `self.nitrification_flux_rate / denitrification_flux_rate` attribute names exactly (no rename, no semantic shift; sibling reads continue to work without changes). The parity test matrix must cover: `use_OrgN={True,False}`, `use_SedFlux={True,False}` (False only — True is `NotImplementedError`), and the NH4/NO3 mass-balance closure case.

**Deliverable:** Nitrogen conformant; Phase 2 test gate satisfied; DOX and Alkalinity verified to read identical values from the preserved sibling-cache attributes (Phase 4 commit triggers a re-run of the Phase 3 DOX baseline parity).

### Phase 5 — FloatingAlgae and BenthicAlgae (2.5 days, 1.25 each)

- FloatingAlgae: integrate the eight existing mortality-routing caches (`algal_*_from_mortality_rate`) into the components dict; the helper must preserve those attribute names exactly.
- BenthicAlgae: consolidate the `_cache_benthic_mortality_rates` helper into the components-dict path (one call, no duplicate `rate_death`). The duplicate-elimination must produce bit-identical mortality rates — the duplicate was a code defect, not a semantic.

**Deliverable:** both algae Processes conformant; Phase 2 test gate satisfied; downstream Nitrogen/Phosphorus/Carbon parity re-verified (algae caches are widely consumed; the Phase 5 commit triggers re-run of Phase 2/3/4 baselines).

### Phase 6 — Phosphorus (1 day)

Existing `tip_settling_rate`, `orgp_settling_rate`, `orgp_to_tip_hydrolysis_rate` caches feed into the components dict; attribute names preserved.

**Deliverable:** Phosphorus conformant; Phase 2 test gate satisfied.

### Phase 7 — POM and CBOD (1.5 days total)

- POM: cache write moves from `rate()` to `run()` (Phase 1 already moved it, but Phase 7 wraps the consolidated rate-and-cache flow into `_change_with_components`).
- CBOD: per-group diagnostics aggregated; per-group bit-identicality required.

**Deliverable:** POM and CBOD conformant; Phase 2 test gate satisfied; per-CBOD-group output verified.

### Phase 8 — N2 and Pathogen (1 day total)

- N2: `_change_with_components`; extend the existing opportunistic-write loop to cover the full Appendix A set (N2 already does this for `total_dissolved_gas` — *extend*, don't *replace* the loop).
- Pathogen: `_rate_with_components` returning `(rate_px, components)`; per-component death/settling rates exposed.

**Deliverable:** N2 and Pathogen conformant; Phase 2 test gate satisfied.

### Phase 9 — Alkalinity (1 day)

`_change_with_components` returning `(delta_alk, components)`. Six nitrification/denitrification/algal-coupling rates exposed.

**Deliverable:** Alkalinity conformant; Phase 2 test gate satisfied.

### Phase 10 — End-to-end conformance + parity + perf check; legacy-inline removal (1.5 days)

- Conformance test: a single test (`tests/v3/nsm1/test_pattern_conformance.py`) iterates every Process class and asserts: (a) `_change_with_components` (or `_rate_with_components`) exists; (b) `REGISTRY_DIAGNOSTICS` exists and is a non-empty tuple; (c) `init_process` captures `self.diagnostics`; (d) clip-with-log call sites match a fixed regex (no `isinstance/diagnostics-not-None` guards remain).
- **Final baseline parity run** (the load-bearing check):
  - Replay the Phase 0 baseline scenario exactly (same RNG seed, same dependency versions, same YAML config, same chunk size).
  - Verify the produced NetCDF trajectory equals the Phase 0 baseline NetCDF *bit-identically* across all 11 state variables × all cells × all 4,320 substeps.
  - Verify per-test numeric assertion values in the JSON baseline match exactly.
  - Verify Tier 1 conservation tolerance and clip-event counts match exactly.
- **Diagnostics-subscription smoke run:** the same demo, with all Appendix A names pre-registered, produces a populated diagnostics output dataset; assert finite values, expected shapes, and that the state-variable subset of the output is *still* bit-identical to the no-subscription baseline.
- Perf check: per-step wall time within 5% of the pre-refactor baseline on the 5-cell synthetic mesh (current 17.6 ms/step), measured no-subscription. With full subscription, wall time within 15% (the higher overhead reflects the `set_at_time` calls).
- **Legacy-inline-cleanup commit:** once the final baseline parity passes, the `_change_legacy_inline` shadow methods added in Phases 2–9 are deleted in a single commit, and the `test_<process>_helper_vs_inline.py` tests are deleted alongside them. This commit re-runs the full §12 contract one final time.

**Deliverable:** v3 NSM1 1.0.0 + pattern-alignment release candidate; full test suite green; bit-identical parity verified end-to-end; perf within tolerance; legacy-inline scaffolding removed.

### Phase 11 — Docs, migration notes, LimnoTech materials (0.5 day)

- Update `clearwater_modules_v3_nsm1_README.md` with the registry-diagnostics convention.
- Update `clearwater_modules_v3_nsm1_design_specification.md` Section 14 Q10 closure to point at this spec as the structural realisation.
- Update `clearwater_modules_v3_nsm1_migration.md` with a "subscribing to diagnostics" worked example.
- Update or refresh `clearwater_modules_v3_nsm1_limnotech_review.md` (decision on §11 versioning — see §10 Open Questions).

**Deliverable:** documentation aligned; review packet refreshed.

**Total wall-clock with Claude doing the coding: 11–13 working days.** Mostly mechanical per-Process passes (Phases 2–9) on a single canonical template.

---

## 7. Testing and Validation

### Per-phase tests (added in the phase that introduces the pattern for the Process)

- `test_<process>_components.py` — fused helper returns the same integrator delta(s) as the legacy inline composition (matrix of state/forcing inputs).
- `test_<process>_registry_diagnostics.py` — pre-registering a `REGISTRY_DIAGNOSTICS` name causes it to be written each substep; not pre-registering yields bit-identical state-variable outputs.
- `test_<process>_clip_unconditional.py` — `clip_negative_state` is called unconditionally; verify with a `unittest.mock.patch` of `clip_negative_state` that it is invoked once per state variable per substep regardless of `self.diagnostics`.

### Cross-Process tests (added in Phase 10)

- `test_pattern_conformance.py` — structural assertions over every Process class.
- `test_appendix_a_completeness.py` — for every name in design spec Appendix A, exactly one Process declares it in `REGISTRY_DIAGNOSTICS`.
- `test_coupled_demo_parity.py` — bit-identical state-variable outputs pre-vs-post refactor over a 4,320-step coupled NSM1+Riverine run, no diagnostics pre-registered.
- `test_coupled_demo_diagnostics.py` — full Appendix A pre-registered; verify diagnostics dataset has finite values and reasonable dimensional shape.

### Tier 1 conservation tests

The 8 existing Tier 1 closed-system mass-conservation tests at `tests/v3/nsm1/test_<constituent>_tier1.py` must continue to pass at `rtol=1e-12` with zero clip events after each phase. This is the most aggressive parity check the suite carries; phases 2–9 are individually gated on it.

---

## 8. Performance Targets

The refactor introduces zero new work on the integrator hot path when no `REGISTRY_DIAGNOSTICS` names are pre-registered. The opportunistic-write loop is `n_diagnostics × O(1)` membership checks; on the 5-cell mesh this is below measurement noise.

When all Appendix A names are pre-registered, the cost is `n_diagnostics × n_cells × set_at_time` — bounded by the number of names per Process (max ~13 for FloatingAlgae) and the registry write cost. Phase 10 perf check budgets ≤5% wall-clock overhead on the 5-cell baseline (≤18.5 ms/step from the current 17.6 ms).

If the perf check fails, the most likely culprit is `set_at_time` overhead per name; mitigation is to batch the writes into a single `update` call per Process. This is a Phase 10 contingency, not a Phase 0 risk.

---

## 9. Risks and Mitigations

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| `_change_with_components` rewrite introduces subtle kinetics changes (e.g., evaluation order, broadcast semantics) | Medium | High | §11 zero-regression contract enforced at every phase commit: bit-identical state trajectory + helper-vs-inline parity with `rtol=0, atol=0` + code-motion-only refactor discipline. Failures roll back the phase; no silent tolerance relaxation. |
| `clip_negative_state` graceful-no-op on `None` diagnostics breaks Tier 1 tests by absorbing what was previously a hard failure | Low-Medium | Medium | Tier 1 tests run after Phase 0 hardening; if any newly passes that previously failed, investigate. |
| BenthicAlgae `init_process` sibling-discovery fix uncovers tests that assumed `use_floating_algae = True` unconditionally | Low | Low-Medium | Audit `tests/v3/nsm1/test_benthic_*.py` in Phase 1 for setup assumptions; update fixtures if needed. |
| Step-index argument to `clip_negative_state` cannot be resolved cleanly from `Process.run` without a Model dependency | Medium | Low | Phase 0 resolution: introduce a `self._substep_counter` increment in `Process.run` base class or pass `step` via `model` on `init_process`; the existing `step=0` placeholder is correct enough that diagnostics rate-limiting still works in the worst case. |
| Appendix A diff turns out larger than expected and the diagnostic naming reveals inconsistencies | Low | Low | Phase 0 produces the diff before any code change; if naming churn is needed, do it once at Phase 0, not mid-stream. |
| LimnoTech objects to the registry-diagnostics convention's surface area or naming | Low-Medium | Medium | Spec is sent for review at Phase 0; objections route to Phase 0 amendment before any code is written. The convention is the resolution of design spec §14 Q2/Q10, so objections would be a re-litigation of that decision, not a new question. |
| Performance regression beyond the 5% budget when all diagnostics are subscribed | Low | Low-Medium | Phase 10 contingency: batch `set_at_time` writes per Process. Falls outside the parity invariant, so the workaround is purely internal. |

---

## 10. Open Questions

1. **Step-index argument to `clip_negative_state`.** Current `Nitrogen._clip` passes `step=0` as a placeholder. The diagnostics rate-limiter relies on the step counter to deduplicate clip-event logs. Options: (a) introduce a `self._substep_counter` incremented in a base-class `Process.run` wrapper; (b) capture `self._step_ref = model.current_step_ref` in `init_process` and read it; (c) accept the placeholder and document the rate-limiter degraded behaviour. **Recommended: (b)**, since the Model already maintains a step counter and `init_process` is the natural binding point. Resolve in Phase 0.

2. **`clip_negative_state` graceful no-op when diagnostics is `None`.** Today the `isinstance/diagnostics-not-None` guards exist precisely because `clip_negative_state` was assumed to require diagnostics. The utility must support `None`-diagnostics with a fast pure-`xr.where` path so call sites can drop the guards. Resolve in Phase 0; verify against Tier 1 tests.

3. **`Process` versioning for this work.** Two viable framings:
   - **NSM1 1.1.0** — semver-minor; this is a structural alignment plus diagnostics-surface addition, no kinetics behaviour change. Diagnostics subscription is additive.
   - **NSM1 1.0.1** — semver-patch; the diagnostics-subscription surface is invisible to users who don't opt in, and bit-identical outputs when not subscribed make the patch reading defensible.
   **Recommended: 1.0.1 with the understanding that 1.1.0 is reserved for first-kinetics-additive work** (e.g., NSM2 multi-pool organic matter). Decide before Phase 11.

4. **LimnoTech review-packet reopening.** Two viable framings:
   - **Reopen 1.0.0 packet** — refresh `clearwater_modules_v3_nsm1_limnotech_review.md` with the structural-alignment section before sign-off.
   - **Ship as 1.0.1 follow-up** — 1.0.0 review proceeds with its current scope; this work lands as a separate 1.0.1 review.
   **Recommended: ship as 1.0.1 follow-up** — the 1.0.0 packet is content-correct and the structural items don't load-bear on the four open kinetics-correctness asks (4-way theta swap, value choices, rca/rcb derivation, DIC units). Decide at Phase 11.

5. **Pathogen integrator form.** Pathogen currently uses a rate-form integrator (`state_new = state + rate * dt_days`). Its `_rate_with_components` is the natural fit, but uniformity with the delta-form would name it `_change_with_components` returning `(delta, components)` instead. Choice: (a) per-Process naming reflects integrator form (chosen above); (b) uniform name across all Processes, with the rate-vs-delta distinction internal. **Recommended: (a)** for explicit-is-better; the convention shows up in the signature, not just the name.

---

## 11. Zero-Regression Contract (load-bearing)

This contract is the single most important section of this specification. Every phase commits under its terms; any phase commit that does not satisfy it is rejected and rolled back before the next phase begins. The contract has six clauses.

### 11.1 The Phase 0 baseline is the gold reference

The Phase 0 baseline-capture artifacts (`baseline_<commit>.txt`, `baseline_<commit>.json`, `baseline_coupled_trajectory_<commit>.nc`, `baseline_tier1_<commit>.json`) are the authoritative pre-refactor reference. They are committed to the repository under `tests/v3/nsm1/` and version-pinned by commit hash. They are never overwritten — if a later baseline re-capture becomes necessary (e.g., dependency upgrade), it is committed as a separate file with the new hash, and the comparison-target is updated explicitly in the test suite, not silently.

### 11.2 Bit-identical state-trajectory invariant

For every phase commit, the state-variable trajectory of the 4,320-step coupled run (no `REGISTRY_DIAGNOSTICS` names pre-registered) must equal the Phase 0 baseline NetCDF *bit-identically*. The comparison uses `xr.testing.assert_identical` (which requires equal dtype, equal coordinates, equal attrs, and `numpy.array_equal` on the values). Tolerance: `rtol=0, atol=0`. No exceptions for "tiny" differences. A single-bit difference in a single cell at a single substep fails the phase.

Rationale: bit-identicality is the only invariant that catches every kinetics regression. Approximate-equality tolerances (`rtol=1e-12`, etc.) accumulate over 4,320 substeps and can mask a regression that produces a `~1e-14` per-step error — which is exactly the size of error that operand-reordering, intermediate-variable elimination, or Numpy/xarray broadcast-rule shifts can introduce.

### 11.3 Per-Process helper-vs-inline parity

Phases 2–9 each retain a `_change_legacy_inline(...)` shadow method that holds the pre-refactor inline composition verbatim. The phase's `test_<process>_helper_vs_inline.py` asserts the new `_change_with_components` and the legacy `_change_legacy_inline` produce bit-identical outputs (`rtol=0, atol=0`) across a parametrised input matrix:

- Zero state (`state = 0` for all state variables).
- Uniform state (all cells, all variables at typical mid-range values).
- Randomised state, fixed RNG seed (state and forcing sampled from realistic ranges).
- Edge cases per process (zero DOX, zero depth, very high temperature, near-zero algae, etc.; the per-Process Phase docstring enumerates).

The shadow methods are removed in Phase 10 only after the final end-to-end baseline parity passes.

### 11.4 Full test suite, no skips, every phase

Every phase commit runs `pytest tests/v3 -x` and the exit code must be zero. The number of `xfailed` tests must equal the Phase 0 baseline count exactly (currently 3). A test that becomes xpassed (xfailed → passed) is investigated, not silently accepted. A test that becomes xfailed (passed → xfailed) is a regression and the phase is rolled back.

The Tier 1 closed-system conservation tests (`tests/v3/nsm1/test_<constituent>_tier1.py` + `tests/v3/nsm1/test_validation_tier1_conservation.py`) must achieve the same tolerance and the same clip-event count as the Phase 0 baseline JSON. `rtol=1e-12` with zero clip events is the documented contract; if a phase changes either, it is a regression.

### 11.5 Rollback contract

If any phase fails any clause of §12.2 / §12.3 / §12.4, the phase commit is rolled back via `git revert` (not `git reset --hard`, so the rollback is visible in the history). The phase is then re-planned, not retried verbatim. The next phase does not begin until the failing phase is either re-landed cleanly or formally deferred.

Rolling back never silently relaxes the contract. If a phase truly cannot satisfy bit-identicality (e.g., a Numpy version change in a dependency, an unavoidable broadcast-rule shift), the failure is documented, the user is notified, and a deliberate decision is made — never absorbed into the test tolerances without explicit sign-off.

### 11.6 Refactor-discipline rules (apply within each phase)

The following code-level rules are non-negotiable inside each `_change_with_components` refactor:

1. **Code motion, not code rewrite.** Lines move from `run` (or wherever they live today) into the helper. Operand order is preserved. Intermediate variables are preserved with their original names where they cross the helper boundary.
2. **No algebraic simplification.** Even mathematically identical rearrangements (`a * b + a * c` → `a * (b + c)`) are forbidden in this refactor. They land in a separate commit, after the phase's parity gate passes, and require a re-run of the §12 contract for the affected Process.
3. **No reordering of `xr.where` branches.** `xr.where(cond, A, B)` is not converted to `xr.where(~cond, B, A)`, even though it is algebraically identical, because Numpy/xarray broadcast-rule differences in dtype promotion can produce different floating-point outputs in the false-cell padding.
4. **No replacement of `np.minimum`/`np.maximum` with `xr.where` or `.clip(...)`, or vice-versa.** The exact ufunc the original code used is preserved.
5. **No replacement of in-method constants with class-level constants** (or vice versa) inside the refactor commit; cosmetic refactors land in a separate commit.
6. **No changes to dtype.** If the legacy code computes in `float64`, the helper computes in `float64`. If a `float32` cast was implicit somewhere, it remains implicit at the same call site.
7. **The `components` dict is populated with the values used by the integrator** — not recomputed downstream. The dict is the single source of truth for both the integrator-update arithmetic and the registry-diagnostics export. Recomputing a sub-flux for the dict (vs. caching it from the integrator path) risks subtle drift and is forbidden.

These rules are checked by code review at every phase commit. A diff that violates any rule is rejected at review, regardless of whether the test gate passes.

---

## 12. Approval Criteria

This specification is complete enough to proceed if the reviewer agrees that:

1. The TSM pattern (§3, A–J) is the right canonical pattern for NSM1.
2. The per-Process component inventory (§4) covers the right Appendix A diagnostics for each Process; the additions to Appendix A are reasonable.
3. The per-Process refactor delta (§5) matches the current code state.
4. The phased plan (§6) sequences the heavy integrators first and the tracers last, with bit-identical parity gates between phases.
5. The testing strategy (§7) catches both kinetics regressions and structural divergences.
6. The performance target (§8) is acceptable.
7. The risks (§9) are correctly identified.
8. The five open questions (§10) are the right ones to surface before code begins.
9. **The zero-regression contract (§11) is the right load-bearing safety net** — bit-identical state trajectories, per-phase rollback, code-motion-only refactor discipline.

If any section is wrong, incomplete, or misframed, mark up this document directly and the spec will be revised before any code is written.

---

## Appendix A — Cross-reference to NSM1 1.0.0 design-spec Section 14 resolutions

This work realises the structural side of the following Section 14 resolutions that were left as code TODOs in 1.0.0:

- **Q2 — Inter-process rate variables in the registry.** Spec §14 resolved naming convention; this spec resolves the *exposure mechanism* (opportunistic registry writes from a `REGISTRY_DIAGNOSTICS` tuple per Process).
- **Q3 — Negative-state handling: clip-with-log.** Spec §14 resolved the contract; this spec realises the uniform application (pattern D harmonisation across all Processes).
- **Q4 — Jacobi for state, Gauss-Seidel for rate variables.** Spec §14 resolved the semantics; this spec realises the producer-precedes-consumer enforcement by making rate caches step-scoped via `_change_with_components` (the cache is rebuilt every substep, never persisting across substeps).
- **Q10 — Step-scoped GS-rate plumbing.** Spec §14 deferred the registry-side plumbing as a "Phase 2.A.1 follow-up"; this spec is that follow-up.

The kinetics-correctness side of the 1.0.0 review (4-way theta swap, value choices, rca/rcb derivation, DIC units) is unaffected by this work and remains the canonical content of `clearwater_modules_v3_nsm1_limnotech_review.md`.
