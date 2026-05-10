# ClearWater Modules v3 — Multi-Group Algae Design Specification

**Status:** Draft, awaiting approval
**Author:** Todd Steissberg (ERDC), with Claude
**Date:** 2026-05-09
**Scope:** Refactor `FloatingAlgae` and `BenthicAlgae` in `clearwater_modules_v3` to support a configurable number of algal functional groups. Step 1 of the four-step HAB / NSM2 roadmap.

**Read this with the umbrella spec.** Quick start, env setup, branch conventions, package architecture, integrator-pattern contract, and umbrella risks live in `clearwater_modules_v3_architecture_specification.md`. This document covers only the multi-group algae refactor.

---

## 1. Background and Roadmap Placement

This work is **step 1 of a four-step sequence** that brings ClearWater toward feature parity with W2 and NSM2 for HAB-relevant studies on rivers, with a focus on the ClearWater-Riverine + HEC-RAS-2D coupling:

1. **Multi-group algae** — this spec.
2. **HAB capabilities** — N-fixation, biomass floor + rewet seeding, toxin tracer, photo-inhibition, bell-curve T-response, P-luxury, source-attribution tracers, low-DO mortality (documented as calibration proxy).
3. **Sediment diagenesis** — Di Toro multi-G (JNH4, JNO3, JCH4, JSO4, JH2S, JDIC, JDIP, SOD) as a v3-native process module.
4. **Remaining NSM2 features** — multi-pool organic matter, alkalinity/pH, methane/sulfide, silica.

Step 1 introduces the structural mechanism (an `algae_group` dimension and per-group parameterization). Steps 2–4 layer features on top of this structure as per-group parameters and as new state variables. This spec does **not** add HAB-specific kinetics, sediment diagenesis fluxes, alkalinity/pH coupling, multi-pool OM, or silica limitation. Those are explicitly out of scope.

### Reference implementations surveyed

| Source | Group control | State per group | Cyano-specific path | Si handling |
|---|---|---|---|---|
| **CE-QUAL-W2 v2026.02** | `NAL` (integer; per-group arrays via `(K,I,JA)` indexing) | Algae as concentration; per-group stoichiometry, settling, T-optima, light response | `CRIT_TIN(JA)` per-group N-fixation trigger; `MIGRATE_GROUP(MIGI)` opts groups into buoyancy migration | `ASI(JA)` per-group Si stoichiometry; `AHSSI(JA)` per-group Si half-sat |
| **NSM2 (HEC-RAS-WQ Fortran)** | `nAp` (integer; dynamically allocated `Ap(nAp)`; 1-based loop) | Single biomass state in µg-Chl-a/L; derived stoichiometry per group | None — N-fixation simulated by setting `KsNp(i,r) = 0` (no explicit flag) | `Si_limitation_option(i,r)` per group (1=unlimited, 2=Monod); `KsSip(i,r)` per group |
| **NSM1 v3 (current ClearWater)** | None — single generic group | `floating_algae` (µg-Chl-a/L); `benthic_algae` (g-D/m²) | None | None |

The design below mirrors the **NSM2 structural pattern** (per-group config block, derived stoichiometry, per-group selectors for limitation models) while adopting the **W2 naming convention** (`NAL` → `nal` as the group count) for downstream readability when users move between models. The HAB-specific features that W2 and the planned step-2 work add (N-fixation, buoyancy, toxin) are introduced as per-group fields in step 2; their parameter slots are reserved here but defaulted to "off."

### Current state in v3

`FloatingAlgae` and `BenthicAlgae` are single-group, with parameter dicts at `clearwater_modules_v3/parameters/{algae,balgae}.py` and process implementations at `clearwater_modules_v3/processes/{floating_algae,benthic_algae}.py`. Several documented wiring bugs from the NSM1 audit (`design/clearwater_modules_v3_nsm1_audit_algae.md`) remain open: light-option-1 parenthesis (F5), Steele exponent sign (B6), growth-option-3 P logic (F14 — partially fixed in `floating_algae.py:545–556`). These are tracked separately in the NSM1 audit and the in-progress Phase 2.A NSM1 rewrite. **This refactor does not pre-empt those fixes**, but it must not regress them either; the multi-group rate methods inherit whichever single-group formulation lands at merge time. See Section 7 for the bug-fix sequencing decision.

---

## 2. Goals and Non-Goals

### Goals

1. **Configurable group count.** Users specify `nal` (number of floating algae groups) and `nalb` (number of benthic algae groups) in the YAML config. Default `nal = 1`, `nalb = 1` produces behavior numerically identical to the current single-group implementation (within floating-point tolerance).
2. **Per-group parameters.** Every parameter currently in `parameters/algae.py` and `parameters/balgae.py` becomes per-group. The YAML `parameters` block accepts a list of per-group dicts, with omitted entries falling back to the v3 default.
3. **Group-aware kinetics.** Rate methods (`rate_growth`, `rate_death`, `rate_respiration`, `rate_settling`, `limit_*`) broadcast over the group dimension. Downstream consumers (Nitrogen, Phosphorus, Carbon, DOX, POM) sum group contributions before applying their own kinetics.
4. **xarray-native group dimension.** State variables gain an `algae_group` dim (for floating) and `balgae_group` dim (for benthic). Existing v2/v3 registry names (`algae_floating`, `benthic_algae`) are preserved; only the shape changes.
5. **Reserved hooks for step 2.** Per-group fields for N-fixation flag, T-optimum, photo-inhibition selector, buoyancy flag, and toxin yield are present in the schema with default values that disable the feature. Step 2 wires them into rate equations; step 1 reserves the slots.
6. **Test coverage.** A new `tests/v3/test_multigroup_algae_v3.py` exercises 1-group (parity with current), 2-group (cyano + non-cyano), and 3-group (W2-style: cyano, greens, diatoms) configurations with regression assertions against single-group baseline outputs.
7. **Documentation.** Migration notes for v2/v3 single-group users; YAML schema documentation; mapping tables to W2 `NAL` and NSM2 `nAp` parameter names.

### Non-Goals

1. **No HAB-specific kinetics in step 1.** N-fixation logic, biomass floor, toxin tracer, photo-inhibition, T-bell-curve, P-luxury (Droop), source-attribution tracers, and low-DO mortality are step 2.
2. **No silica.** Silica limitation, Si state variables, and diatom-specific Si stoichiometry are deferred to the step-4 NSM2 feature work. Step 1 reserves the parameter slot but does not enable it.
3. **No per-cell-per-group spatial heterogeneity in parameters.** Parameters in step 1 are per-group, scalar across cells. Per-cell parameter overrides (e.g., site-specific cyano kinetics) remain a v3.x feature.
4. **No legacy v1 changes.** v1 stays single-group and is on its existing deprecation track per the architecture spec.
5. **No retirement of v2.** v2 stays single-group and frozen per architecture spec Section 5.
6. **No new integrator pattern.** The Forward-Euler integrator from Phase 2.A is retained; only the array shape changes.
7. **No alkalinity/pH coupling.** pH-driven CO₂ limitation (the more defensible bloom-collapse mechanism) is deferred to v3.x when alkalinity/pH lands.

---

## 3. Component Inventory

### 3.1 New / changed parameter files

#### `clearwater_modules_v3/parameters/algae.py`

The flat `DEFAULTS` dict becomes a **per-group default template** plus a global section:

```python
GROUP_DEFAULTS: dict[str, float | int | bool] = {
    # Stoichiometry (mirrors NSM2 modAlgae.f90:15-46)
    'AWd': 100.0, 'AWc': 40.0, 'AWn': 7.2, 'AWp': 1.0, 'AWa': 1000.0,

    # Kinetics (single-group v3 values; published consensus midpoints)
    'mu_max_20': 2.0, 'mu_max_theta': 1.047,
    'kdp_20': 0.05,  'kdp_theta':  1.047,
    'krp_20': 0.10,  'krp_theta':  1.047,
    'vsap':   0.15,

    # Limitation half-sats and constants
    'KL':  10.0, 'KsN': 0.04, 'KsP': 0.0012,

    # Selectors (per group)
    'growth_rate_option':     1,   # 1=Multiplicative, 2=Min, 3=Harmonic
    'light_limitation_option': 1,  # 1=Half-Sat, 2=Smith, 3=Steele

    # NH4 vs. NO3 preference (was implicit; promoted to per-group, default v1 behavior)
    'PN': 0.5,

    # Step-2 reserved hooks (defaults disable the feature; wired in step 2)
    'is_n_fixer': False,            # cyano N-fixation flag
    'crit_tin':   0.0,              # mg-N/L; TIN threshold for N-fixation onset
    'toxin_yield': 0.0,              # ug-toxin / ug-Chla; intracellular toxin
    't_opt':      20.0,              # °C; optimal temperature (placeholder; bell-curve in step 2)
    'kt1':        0.0,               # 1/°C^2; below-optimum slope (step 2)
    'kt2':        0.0,               # 1/°C^2; above-optimum slope (step 2)
    'photoinhibit_option': 0,        # 0=off, 1=Steele-style, 2=Platt-Jassby (step 2)
    'is_buoyant':   False,            # buoyancy/migration flag (step 2; lake/3D only)
    # Si reserved for step 4
    'si_limitation_option': 1,       # 1=unlimited (step 4 enables option 2=Monod)
    'KsSi':                 0.03,    # mg-Si/L; ignored when si_limitation_option=1
    'AWsi':                 0.0,     # mg-Si / ug-Chla; ignored when si_limitation_option=1
}

# Optional per-group label hints (purely cosmetic; no behavioral effect)
GROUP_LABEL_DEFAULTS: dict[int, str] = {
    0: "phytoplankton",   # used when nal=1
}

# Also retained for back-compat with v2-style flat-dict imports during migration:
DEFAULTS = GROUP_DEFAULTS
```

**Rationale**: keeping `DEFAULTS` aliased to `GROUP_DEFAULTS` lets the existing single-group code and tests continue to work during phased rollout. The class-level merge in `FloatingAlgae.__init__` (currently lines 213–221) gains a small wrapper that builds a per-group list of merged dicts.

#### `clearwater_modules_v3/parameters/balgae.py`

Identical pattern. `BWd`, `BWc`, `BWn`, `BWp`, `BWa`, `KLb`, `KsNb`, `KsPb`, `Ksb`, `mub_max_*`, `krb_*`, `kdb_*`, `Fw`, `Fb` all become per-group. The existing comment block on the WASP7-canonical `BWa = 1000` correction (Phase 9.E) is preserved.

#### `clearwater_modules_v3/parameters/global_parameters.py` (new entries)

```python
GLOBAL_DEFAULTS_ADDITIONS = {
    'nal':  1,   # number of floating algae groups
    'nalb': 1,   # number of benthic algae groups
}
```

These live in `global_parameters.py` rather than `algae.py` because the orchestrator and the registry need them at scaffold time, before the algae process is instantiated.

### 3.2 YAML schema additions

**Single-group (v2/v3 backward-compatible):**
```yaml
processes:
  floating_algae:
    parameters:
      mu_max_20: 2.0
      KsN: 0.04
      # ... existing flat dict
```
Behavior: `nal` defaults to 1; the flat parameters dict is wrapped to a 1-element list. Output identical to current v3.

**Multi-group (new):**
```yaml
global:
  nal: 3
  nalb: 1

processes:
  floating_algae:
    parameters:
      - label: cyanobacteria      # group 0; optional label for logs/output
        mu_max_20: 1.8
        KsN: 0.02                 # lower N half-sat -> wins at low DIN
        is_n_fixer: true          # step-2 hook (no-op in step 1)
        crit_tin: 0.05            # step-2 hook (no-op in step 1)
        t_opt: 28.0               # warm optimum (step-2 hook)
        toxin_yield: 0.1          # step-2 hook
      - label: greens             # group 1
        mu_max_20: 2.5
        KsN: 0.05
      - label: diatoms            # group 2
        mu_max_20: 2.0
        KsN: 0.04
        AWsi: 0.5                 # step-4 hook (no-op in step 1)
        si_limitation_option: 2   # step-4 hook (no-op in step 1)
```

**Schema rule**: when `parameters` is a list, its length must equal `nal` (validated at config load). When `parameters` is a dict, `nal` must be 1 (or the dict is implicitly broadcast to all groups, which is the W2 default-block convention — see open question Q1).

### 3.3 Process refactor: `processes/floating_algae.py`

The class gains a `nal` attribute and an internal list of per-group merged parameter dicts (`self.group_params: list[dict]`). Per-group attributes are stored as numpy arrays of length `nal` keyed by parameter name, e.g., `self.mu_max_20: np.ndarray` with shape `(nal,)`. xarray broadcasting handles cell × group automatically.

| Current | After refactor |
|---|---|
| `self.AWn: float` | `self.AWn: np.ndarray` shape `(nal,)` |
| `algae = registry.get_at_time("algae_floating", time)` returns `(cell,)` DataArray | Returns `(cell, algae_group)` DataArray |
| `rate_growth(...)` returns `(cell,)` | Returns `(cell, algae_group)` — xarray broadcasts cleanly |
| `algal_growth_rate` cache stores `(cell,)` | Stores `(cell, algae_group)` |
| Mortality routing: `rna * ap_death` returns `(cell,)` | Returns `(cell, algae_group)`; downstream consumers `.sum("algae_group")` before adding to their own state |

The new pattern in `_cache_mortality_rates` (currently lines 431–473):

```python
def _cache_mortality_rates(self, algae, water_temperature):
    ap_death = self.rate_death(algae, water_temperature)  # (cell, algae_group)

    rna = self.AWn / self.AWa  # (algae_group,)
    rpa = self.AWp / self.AWa
    rca = self.AWc / self.AWa

    self.algal_death_rate = ap_death                          # per-group
    self.algal_orgn_from_mortality_rate = (rna * ap_death)    # per-group
    self.algal_orgp_from_mortality_rate = (rpa * ap_death)
    self.algal_poc_from_mortality_rate  = (self.f_pocp * rca * ap_death)
    self.algal_doc_from_mortality_rate  = ((1 - self.f_pocp) * rca * ap_death)

    self.algal_pom_from_settling_rate = (
        self.vsap * algae * (self.AWd / self.AWa) / self.h2
    )  # per-group
```

Downstream consumers (`Nitrogen`, `Phosphorus`, `Carbon`, `POM`, `DOX`) sum over groups when building their water-column rate:

```python
# in Nitrogen.run, replacing the current scalar read:
algal_orgn = self.floating_algae.algal_orgn_from_mortality_rate.sum("algae_group")
```

The `.sum("algae_group")` collapse is the single behavioral change forced on every consumer. **Open question Q3** asks whether this should be encapsulated in a helper to avoid scattering it across consumers.

### 3.4 BenthicAlgae refactor: `processes/benthic_algae.py`

Same pattern with `nalb` and `balgae_group` dim. Benthic-specific group concept maps cleanly to functional types: diatom biofilm (low T-optimum, Si-limited in step 4), green-algae mats, cyanobacterial mats (N-fixing in step 2). The `balgae_group` dim is independent of `algae_group`.

### 3.5 Registry / state-variable changes

| Variable | Before | After |
|---|---|---|
| `algae_floating` | `(cell,)` | `(cell, algae_group)` |
| `benthic_algae` | `(cell,)` | `(cell, balgae_group)` |
| Boundary conditions / inflow concentrations | `(time, cell)` | `(time, cell, algae_group)` |

Boundary conditions need user-supplied per-group inflow time series. **For the default `nal = 1` case, the existing single-tributary-Chla time series is broadcast to the 1-element group dim transparently.** Multi-group configs require per-group inflow time series; documentation will guide users on splitting an aggregate Chla measurement into functional-group fractions (a non-trivial but established practice — W2 user manual gives the same guidance).

### 3.6 Output schema

Output xarray Datasets gain the `algae_group` and `balgae_group` dims. The optional `label` field from the YAML config is written as a coordinate so downstream analysis can refer to groups by name (`ds.sel(algae_group="cyanobacteria")`) rather than index.

### 3.7 Files affected (initial estimate)

| File | Change type |
|---|---|
| `parameters/algae.py` | Rename `DEFAULTS` → `GROUP_DEFAULTS`; add reserved-hook keys; alias `DEFAULTS = GROUP_DEFAULTS` |
| `parameters/balgae.py` | Same pattern |
| `parameters/global_parameters.py` | Add `nal`, `nalb` defaults |
| `processes/floating_algae.py` | Per-group attribute storage; broadcast over group dim |
| `processes/benthic_algae.py` | Same |
| `processes/nitrogen.py` | Add `.sum("algae_group")` at consumer sites |
| `processes/phosphorus.py` | Same |
| `processes/carbon.py` | Same |
| `processes/cbod.py` | Same |
| `processes/dox.py` | Same |
| `processes/pom.py` | Same |
| `config/init.py` | Validate `parameters` list-vs-dict shape; build per-group merged param list |
| `model.py` | Plumb `nal`/`nalb` from global config to processes at scaffold time |
| `tests/v3/test_multigroup_algae_v3.py` | New |
| `tests/v3/test_v2_v3_parity.py` (existing) | Add 1-group parity assertion for multigroup-refactored code |

---

## 4. Testing and Validation

### 4.1 Test infrastructure

- **`tests/v3/test_multigroup_algae_v3.py`** — three configurations exercised:
  - `nal=1` parity: outputs match pre-refactor single-group baseline within float tolerance
  - `nal=2` (cyano + non-cyano with identical kinetics): output sums to 2× the 1-group baseline; per-group state tracks correctly
  - `nal=3` (W2-style functional types with distinct kinetics): regression test against hand-computed expected values
- **`tests/v3/test_multigroup_balgae_v3.py`** — same three configurations for benthic algae
- **`tests/v3/test_multigroup_consumers_v3.py`** — verifies that `Nitrogen`, `Phosphorus`, `Carbon`, `DOX`, `POM`, `CBOD` correctly sum group contributions before applying their own kinetics. A 2-group config with identical kinetics must produce DIN/DIP/DO trajectories identical to a 1-group config with twice the biomass IC.
- **`tests/v3/test_multigroup_yaml_schema.py`** — config-load validation: list length must equal `nal`; missing keys fall back to defaults; unknown keys warn-and-ignore (matching current v3 behavior in `floating_algae.py:223–230`).

### 4.2 Validation tiers

- **Tier 1 (conservation):** mass conservation of N, P, C across the algae→OM→DIN cycle in a 2-group closed-system test. Total biomass change in one group must produce the right OrgN/OrgP/POC/DOC partition with no leakage.
- **Tier 2 (single-group parity):** with `nal=1`, every regression test for the existing single-group implementation passes unchanged.
- **Tier 3 (group additivity):** two groups with identical parameters must produce the same downstream consumer effects as one group with doubled IC.
- **Tier 5 (group differentiation):** with `nal=2` and contrasting kinetics (low-`KsN` cyano vs. high-`KsN` greens), low-DIN regimes must produce cyano-dominated steady state.

### 4.3 Regression suite

The existing v3 NSM1 regression suite (cf. `design/clearwater_modules_v3_nsm1_design_specification.md`) runs unchanged against the `nal=1` default. Any divergence is treated as a refactor bug.

---

## 5. Performance Considerations

xarray broadcasting over a length-`nal` dim adds memory proportional to `nal` and roughly `nal × O(per-cell cost)` compute. For typical HAB studies with `nal ∈ {2, 3, 4}` the slowdown should be near-linear. Targets:

- **Must:** `nal=1` runtime within 1.10× of pre-refactor single-group v3 baseline (a 10% overhead for the broadcasting machinery is acceptable).
- **Should:** `nal=3` runtime within 3.5× of `nal=1` baseline.
- **Aspirational:** `nal=3` runtime within 3.0× of `nal=1` baseline (linear scaling with no broadcasting overhead).

Profile after Phase 2; if `nal=1` overhead exceeds 10%, investigate xarray dim-squeeze paths or special-case the 1-group code path in the hot kernel.

---

## 6. Phased Implementation Plan

Following the TSM-spec convention.

### Phase 0 — Gap analysis (½ day)

Catalog every site in `processes/{floating_algae,benthic_algae}.py` and the seven downstream consumer modules where a single-group assumption is encoded. Produce a Markdown table at `design/multigroup_algae_gap_analysis.md`. Each row: file:line, current shape assumption, refactor disposition (`broadcast`, `sum("algae_group")`, `iterate`, `keep`).

**Deliverable:** the gap-analysis table.

### Phase 1 — Parameter and config plumbing (1 day)

- Rename `DEFAULTS` → `GROUP_DEFAULTS` in `algae.py`, `balgae.py`; add reserved-hook keys; preserve `DEFAULTS` alias.
- Add `nal`, `nalb` to `global_parameters.py`.
- Update `config/init.py` to accept list-of-dicts under `parameters`; build per-group merged param lists; validate length matches `nal`/`nalb`.
- No process-side changes yet; existing single-group code continues to work via the alias.

**Deliverable:** YAML accepts both flat-dict (length 1 broadcast) and list-of-dicts; v3 still runs with `nal=1` exactly as before.

### Phase 2 — FloatingAlgae per-group kinetics (2–3 days)

- Refactor `FloatingAlgae` per-group attribute storage (numpy arrays of length `nal` for scalars; broadcasting handled by xarray).
- Update rate methods (`rate`, `rate_growth`, `rate_death`, `rate_respiration`, `rate_settling`, `limit_*`, `_cache_mortality_rates`) to broadcast over the `algae_group` dim.
- Update state-variable shape: `algae_floating` carries `algae_group` dim.
- Run existing regression suite with `nal=1`; assert byte-identical (or floating-point equivalent) outputs.

**Deliverable:** `FloatingAlgae` works at any `nal`; single-group regression passes.

### Phase 3 — Downstream consumer refactor (2 days)

- Add `.sum("algae_group")` (or helper, see open question Q3) at every consumer site.
- Update consumer rate caches and registry writes to keep their existing shapes (i.e., the algae-group dim is collapsed at the boundary; nitrogen/phosphorus/carbon/DOX state remains `(cell,)`).
- Test with the 2-group identical-kinetics regression: consumer trajectories must match 1-group-doubled-IC.

**Deliverable:** Multi-group floating algae produces correct consumer-side fluxes.

### Phase 4 — BenthicAlgae mirror refactor (1–2 days)

- Apply Phases 1–3 to `BenthicAlgae` with `nalb` / `balgae_group` dim.
- Same regression and additivity tests.

**Deliverable:** BenthicAlgae multi-group ready.

### Phase 5 — Boundary conditions and output (1 day)

- Update boundary-condition / inflow time-series readers to accept per-group time series; broadcast 1-group time series to `nal=1` configs.
- Add `label` coordinate to output Datasets when YAML supplies group labels.

**Deliverable:** End-to-end multi-group run on a representative test case.

### Phase 6 — Tests and validation (2 days)

- Implement the four test files in Section 4.1.
- Run all four validation tiers in Section 4.2.

**Deliverable:** Test suite green; A/B comparison vs. pre-refactor v3 single-group baseline.

### Phase 7 — Documentation and review prep (½ day)

- Update v3 README with multi-group YAML schema examples.
- Migration notes for v2/v3 single-group users.
- Mapping table to W2 `NAL` and NSM2 `nAp` parameter names for users transitioning between models.
- Prepare materials for LimnoTech review.

**Deliverable:** Multi-group algae ready for v3.1.0 release.

**Total estimated wall-clock with Claude doing the coding: 9–12 working days.**

---

## 7. Risks and Mitigations

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| `nal=1` parity test fails after refactor (broadcast machinery introduces float drift) | Medium | Low–Medium | Use exact-shape preservation when `nal=1`; if drift unavoidable, document tolerance bound and assert in test |
| Downstream consumer (e.g., DOX, CBOD) has hidden single-group assumption that escapes Phase 0 audit | Medium | Medium | Use Tier-3 group-additivity test as catch-all: any consumer that fails it has a hidden assumption |
| Broadcasting overhead > 10% in `nal=1` case | Low–Medium | Medium | Profile early in Phase 2; squeeze the dim in single-group hot path if needed |
| Open NSM1 audit bugs (F5, F14, B6) interact with refactor; fixing during refactor inflates scope | Medium | Medium | Sequencing decision in Q2 below: refactor first against current (buggy) kinetics, then audit fixes land separately on the multi-group structure |
| Boundary-condition handlers in `clearwater_riverine` (downstream consumer of these state shapes) break when `algae_group` dim appears | Medium | Medium | Coordinate with riverine repo: confirm xarray dim handling in their advection/transport code before Phase 5; if breakage, broadcast at the riverine boundary as a temporary shim |
| Multi-group config files diverge stylistically across studies (no convention) | Low | Low | Ship 3 reference YAML configs (1-group, 2-group cyano-vs-non-cyano, 3-group functional types) as templates |

### Coordination

- **LimnoTech review needed for:** YAML schema (list-vs-dict ambiguity), boundary-condition shape changes, downstream-consumer `.sum("algae_group")` pattern (vs. helper encapsulation, Q3).
- **ClearWater-Riverine coordination needed for:** advection/transport handling of the new `algae_group` dim. If they treat it as just-another-constituent-dim, no change. If they have any per-state-variable special-casing, that needs review.
- **Items that can proceed without external input:** Phase 0 (gap analysis), Phase 1 (parameter plumbing internal to v3), Phase 2 (process-internal refactor with regression backstop), Phase 6 tests.

---

## 8. Open Questions

1. **Default-block broadcast vs. strict length match.** When YAML supplies a single `parameters:` dict with `nal > 1`, do we (a) broadcast the dict to all groups (W2 default-block convention; ergonomic for "all groups identical except one"), or (b) raise a config error requiring explicit list of length `nal`? **Proposed default: (a) broadcast**, with a config-load warning so users see what's happening. This matches W2's user mental model and reduces friction for "all groups same except cyano" experiments.

2. **NSM1 audit bug fixes — sequenced before, during, or after the refactor?** Three options: (i) land audit fixes first on single-group v3, then refactor; (ii) bundle fixes into the refactor PRs; (iii) refactor first, fixes after on the multi-group structure. **Proposed default: (iii) refactor first, fixes after.** Rationale: the refactor is structural and reversible; the bug fixes are mechanistic and need their own validation. Mixing them complicates regression analysis. The Phase 2 single-group parity test will pin down "no behavioral change vs. current v3," which is exactly the contract we want.

3. **Consumer-side group sum: inline `.sum("algae_group")` or helper utility?** Inline is simpler but scatters the dim name across modules. A helper (e.g., `floating_algae.consumer_rate("algal_orgn_from_mortality_rate")`) hides the dim and centralizes the pattern. **Proposed default: inline `.sum("algae_group")`** for step 1 — it's transparent and easy to audit; if it becomes painful in step 2 (HAB consumers add more rate caches), refactor to a helper at that point.

4. **BenthicAlgae and FloatingAlgae: same release or staggered?** Staggered would let FloatingAlgae go to v3.1.0 sooner, with BenthicAlgae following in v3.1.1. Same release ships them together as a coherent feature. **Proposed default: same release.** They share enough refactor pattern that staggering doesn't save much, and downstream modules (Nitrogen, Phosphorus, Carbon) need to handle both anyway.

5. **Group label representation.** Numeric (`0..nal-1`) only, or numeric + optional `label` field? The label is purely cosmetic — for log messages, output coordinates, error messages — but materially improves usability for HAB studies where users will routinely think in terms of "cyano group" vs. "diatom group." **Proposed default: numeric + optional label**, with the label stored as an xarray coordinate on the `algae_group` dim. No behavioral effect; pure ergonomics.

6. **W2-style `MIGRATE_GROUP` indirection vs. per-group flag.** W2 has a separate `MIGRATE_GROUP` array indexing which algae groups have buoyancy migration; NSM2 has none (no migration). For step 2's HAB scope, we add `is_n_fixer`, `is_buoyant`, `photoinhibit_option` as per-group flags directly. **Proposed default: per-group flags** (NSM2-style), not a separate index array. Flags are simpler and clearer in YAML. This also extends naturally to the source-attribution tracers (per-group flag for "track this group from source X").

---

## 9. Approval Criteria

This spec is approved when the author has reviewed and accepted:

1. Roadmap placement (Section 1) — multi-group algae is step 1; HAB / sediment diagenesis / NSM2 are step 2/3/4.
2. Goals and non-goals (Section 2) — explicitly excludes HAB kinetics, silica, alkalinity/pH, sediment diagenesis.
3. Component inventory (Section 3) — `nal`/`nalb` global, per-group `GROUP_DEFAULTS`, list-of-dicts YAML schema, xarray `algae_group`/`balgae_group` dims, `.sum("algae_group")` at consumer sites.
4. Testing plan (Section 4) — three-tier: 1-group parity, 2-group additivity, 3-group functional differentiation.
5. Performance targets (Section 5) — `nal=1` overhead < 10%, `nal=3` linear-ish scaling.
6. Phased plan (Section 6) — Phases 0 through 7, ~9–12 working days.
7. Risks (Section 7) and proposed mitigations.
8. Resolutions for the six open questions (Section 8). Defaults are proposed; please confirm or override.

Once approved, implementation begins with Phase 0.

---

## Appendix A: Parameter Mapping — v3 ↔ NSM2 ↔ W2

| v3 multi-group | NSM2 (Fortran) | W2 (Fortran) | Description |
|---|---|---|---|
| `nal` | `nAp` | `NAL` | Number of floating algae groups |
| `nalb` | (not present; benthic single-group) | `NEP` (epiphyton) | Number of benthic algae groups |
| `mu_max_20[g]` | `mu_max(i,r)` | `AG(JA)` | Max growth rate per group |
| `KsN[g]` | `KsNp(i,r)` | `AHSN(JA)` | N half-saturation per group |
| `KsP[g]` | `KsPp(i,r)` | `AHSP(JA)` | P half-saturation per group |
| `vsap[g]` | `vsap(i,r)` | `AS(JA)` | Settling velocity per group |
| `KL[g]` | `KL(i,r)` | `ASAT(JA)` | Light half-saturation per group |
| `AWn[g]` | `AWn(i,r)` | `AN(JA)` × stoichiometry | N stoichiometry per group |
| `is_n_fixer[g]` (step 2) | (set `KsNp=0` per group) | `CRIT_TIN(JA)` | N-fixation flag |
| `t_opt[g]`, `kt1[g]`, `kt2[g]` (step 2) | `T0p(i,r)`, `ktp1(i,r)`, `ktp2(i,r)` | `AT1..AT4(JA)` | Temperature response per group |
| `light_limitation_option[g]` | `light_limitation_option(i,r)` | (single global option) | Light model selector per group |
| `growth_rate_option[g]` | (multiplicative only) | (multiplicative only) | Growth model selector per group |
| `si_limitation_option[g]` (step 4) | `Si_limitation_option(i,r)` | per-group `ASI(JA)` | Si limitation per group |
| `toxin_yield[g]` (step 2) | (not present) | `CTP(J,JA) × CTB(J,JA)` | Toxin yield per group |
| `is_buoyant[g]` (step 2; lake/3D only) | (not present) | `MIGRATE_GROUP(MIGI)` | Buoyancy migration flag |

---

## Appendix B: Reference YAML — three-group HAB-ready config

```yaml
global:
  nal: 3
  nalb: 1

processes:
  floating_algae:
    parameters:
      - label: cyanobacteria
        mu_max_20: 1.8
        mu_max_theta: 1.07
        KsN: 0.02
        KsP: 0.0008
        vsap: 0.05
        AWn: 7.2
        AWp: 1.0
        # step-2 hooks (no-op until step 2 lands)
        is_n_fixer: true
        crit_tin: 0.05
        t_opt: 28.0
        toxin_yield: 0.10
      - label: greens
        mu_max_20: 2.5
        KsN: 0.05
        vsap: 0.20
      - label: diatoms
        mu_max_20: 2.0
        KsN: 0.04
        vsap: 0.30
        # step-4 hooks (no-op until step 4 lands)
        AWsi: 0.5
        si_limitation_option: 2
        KsSi: 0.03

  benthic_algae:
    parameters:
      - label: epiphyton
        # uses balgae GROUP_DEFAULTS
```

In step 1, the step-2 and step-4 fields are accepted by the schema and stored on the process instance but have no effect on kinetics. Step 2 wires `is_n_fixer`, `crit_tin`, `t_opt`, `toxin_yield`, etc., into rate equations. Step 4 wires `AWsi`, `si_limitation_option`, `KsSi`. This forward-compatibility lets users start configuring HAB-ready experiments now without YAML migration when the step-2 release lands.
