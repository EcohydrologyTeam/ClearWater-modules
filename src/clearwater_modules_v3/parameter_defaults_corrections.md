# v3 NSM1 Parameter Defaults: Corrections and Audit Findings

This document records every default-value correction applied to v3 NSM1's
parameter library at the port from v1 (`src/clearwater_modules/nsm1/constants.py`)
and the lower-priority audit findings that were left as-is for v3 1.0.0 with
inline `FIXME(phase1-audit):` comments in the relevant `parameters/<group>.py`
module.

The canonical reference for the seven critical corrections is
`design/clearwater_modules_v3_nsm1_design_specification.md` Section 7. The
broader Phase 0 inventory and rationale live in
`docs/clearwater_modules_v3_nsm1_phase0_parameter_audit.md`.

---

## Section 1: Critical default-value corrections (7 items, applied at port)

These are applied directly in the v3 `DEFAULTS` dicts. Each inline comment in
the relevant `parameters/<group>.py` module records the v1 original.

### 1.1 `vsop` — organic-P settling velocity
- **Module:** `parameters/phosphorus.py`
- **v1 default:** `999` m/d
- **v3 default:** `0.1` m/d
- **Rationale:** `vsop` is multiplied directly into the OrgP loss rate
  (`vsop * OrgP / depth`). At 999 m/d the water column would be entirely
  evacuated of organic phosphorus on every timestep. Typical values for
  organic-P settling sit in the 0.01-1.0 m/d range; 0.1 m/d is a defensible
  midpoint consistent with literature values for medium-sized organic
  particles.

### 1.2 `vs` — TIP settling velocity
- **Module:** `parameters/phosphorus.py`
- **v1 default:** `999` m/d
- **v3 default:** `0.1` m/d
- **Rationale:** Same magnitude error as `vsop`. `vs` controls the settling
  loss of total inorganic phosphorus adsorbed onto suspended solids; physically
  reasonable values fall in the 0.01-1.0 m/d range. 0.1 m/d is the same
  midpoint chosen for `vsop` and is consistent with Chapra (1997) recommended
  ranges for inorganic-P partitioning onto fine sediments.

### 1.3 `SOD_20` — sediment oxygen demand at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` g-O2/m^2/d
- **v3 default:** `1.0` g-O2/m^2/d
- **Rationale:** Realistic SOD values for moderate organic loading sit in the
  0.5-3.0 g-O2/m^2/d range (Chapra 1997, Table 25.2). 1.0 g-O2/m^2/d is a
  defensible default for an unspecified moderate-loading site; users with
  field measurements should override per simulation. The v1 value of 999
  drives DOX immediately negative on any wet-bed timestep.

### 1.4 `SOD_theta` — Arrhenius coefficient for SOD
- **Module:** `parameters/dox.py`
- **v1 default:** `999` (unitless)
- **v3 default:** `1.060` (unitless)
- **Rationale:** Arrhenius-style temperature corrections take the form
  `theta^(T-20)`. With `theta=999` and any temperature above 20 C, the
  corrected SOD blows up by orders of magnitude per degree (catastrophic).
  Chapra (1997) recommends `theta=1.060` for SOD as the standard literature
  value, equivalent to a Q10 of approximately 1.79. This was Phase 0's most
  severe finding.

### 1.5 `kaw_20_user` — user-override wind reaeration coefficient at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` m/d
- **v3 default:** `0.0` m/d
- **Rationale:** This parameter is consulted only when `wind_reaeration_option`
  selects the user-override branch (option 1). Setting the default to zero
  means the user-override branch is *off* by default; if a user opts into the
  override path without supplying their own value they will see no reaeration
  rather than a runaway 999 m/d.

### 1.6 `kah_20_user` — user-override hydraulic reaeration coefficient at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` 1/d
- **v3 default:** `0.0` 1/d
- **Rationale:** Same as `kaw_20_user`. Defaulting to zero makes the override
  branch a no-op when not configured, instead of a numerical disaster.

### 1.7 `pressure_mb` — atmospheric pressure
- **Module:** `parameters/global_parameters.py`
- **v1 default:** `2026.5` hPa
- **v3 default:** `1013.25` hPa
- **Rationale:** Standard sea-level atmospheric pressure is `1013.25` hPa
  (ISO 2533). The v1 value is approximately 2x correct; this biases every
  downstream computation that depends on it, notably `O2sat` (DO saturation
  via Henry's law) and `N2sat` (N2 saturation), as well as any atmospheric
  reaeration formula that uses partial pressure. This was identified as a
  Phase 0 finding (2026-05-04) and added to the canonical correction list.

---

## Section 2: Lower-priority audit findings under review (8 items)

These eight items were flagged during Phase 0 as suspicious or unclear but
left at v1 values for the v3 1.0.0 port. Each is marked with a
`FIXME(phase1-audit):` inline comment in the relevant parameter module.
Disposition for each is described below.

### 2.1 `rnh4_20=0`, `vno3_20=0`, `rpo4_20=0` — sediment-flux release rates disabled
- **Modules:** `parameters/nitrogen.py` (`rnh4_20`, `vno3_20`),
  `parameters/phosphorus.py` (`rpo4_20`).
- **Issue:** All three sediment-release/uptake rates are zero, which silently
  disables sediment-flux contributions to NH4, NO3, and PO4 budgets. v1 also
  defaults `use_SedFlux=False`, which should make this moot, but Phase 0
  could not confirm that every code path consuming these parameters is gated
  by `use_SedFlux`.
- **Disposition (v3 1.0.0):** Kept at zero, consistent with v1. Phase 1.3
  Process implementation must verify that `use_SedFlux=False` gates all
  consumers of these parameters; until that audit is complete, leaving the
  values at zero preserves v1 behavior. If `use_SedFlux=True` is set without
  also providing site-specific rates, results will silently exclude sediment
  fluxes.

### 2.2 `kdpo4=0.0` — TIP partitioning coefficient
- **Module:** `parameters/phosphorus.py`
- **Issue:** Phosphorus partitioning between dissolved (DIP) and particulate
  (sorbed) phases is governed by `kdpo4` (L/kg). At zero, no DIP adsorbs onto
  suspended solids; the TIP partitioning feature is effectively disabled.
- **Disposition (v3 1.0.0):** Kept at zero. Full DIP-solid partitioning is
  NSM2 territory; v3 NSM1 1.0.0 maintains the v1 simplification. Documented
  for clarity.

### 2.3 `ksbod_20=0.0` — CBOD settling rate
- **Module:** `parameters/cbod.py`
- **Issue:** CBOD never settles in v1 because the settling rate is hardcoded
  to zero. If CBOD is conceptually fully dissolved this is correct; if CBOD
  represents a mixture including particulate fractions, this is a bug.
- **Disposition (v3 1.0.0):** Kept at zero pending clarification from
  LimnoTech on the intended interpretation of CBOD groups in NSM1.

### 2.4 `apx=1`, `vx=1` — pathogen placeholder values
- **Module:** `parameters/pathogen.py`
- **Issue:** Both pathogen production/coupling (`apx`) and settling velocity
  (`vx`) are set to 1.0 with no documented literature basis. Likely
  placeholder values inherited from v1.
- **Disposition (v3 1.0.0):** Kept as-is; updating defaults to literature-
  backed values is deferred to a future audit pending site-specific
  pathogen-tracer studies.

### 2.5 `h2=0.1` — POM dissolution depth denominator
- **Module:** `parameters/pom.py`
- **Issue:** `h2` is used as a divisor in burial/sedimentation terms (e.g.,
  `vb * POM / h2`); its physical role is unclear, and the units (meters) do
  not obviously line up with the surrounding terms. Phase 0 flagged it as
  needing clarification.
- **Disposition (v3 1.0.0):** Kept as-is; flagged with `FIXME(phase1-audit):`.
  Full clarification is part of the Phase 1.3 process implementation.

### 2.6 `vb=0.01` — burial velocity magnitude
- **Module:** `parameters/global_vars.py`
- **Issue:** Burial velocity of 0.01 m/d has no documented validation; some
  sediment-flux literature uses values one to two orders of magnitude
  different.
- **Disposition (v3 1.0.0):** Kept as-is; flagged with `FIXME(phase1-audit):`.
  Without tied-to-site sediment data, no defensible alternative default
  exists.

### 2.7 `q_solar=500` units mismatch
- **Module:** `parameters/global_vars.py`
- **Issue:** v1's docstring documents `q_solar` as having units of `1/d`, but
  the value 500 and its consumption pattern (Beer-Lambert PAR computation)
  are consistent with W/m^2. The docstring is wrong, not the value.
- **Disposition (v3 1.0.0):** Value kept at 500; inline comment in
  `global_vars.py` documents the correct units (W/m^2). Docstring will be
  corrected in the Process docstring during Phase 1.3.

### 2.8 `lambdas` — light extinction parameter defined but not used
- **Module:** `parameters/global_vars.py`
- **Issue:** `lambdas=0.052` represents the suspended-solids contribution to
  Beer-Lambert light extinction. v1's code path computing the effective
  extinction coefficient comments out this term: `lambdas * Solid` is
  defined but not added to the sum. Either the term should be enabled or
  the parameter should be removed.
- **Disposition (v3 1.0.0):** Parameter retained for backward compatibility
  with any user YAML that already sets it; flagged with `FIXME(phase1-audit):`
  in `global_vars.py`. Phase 1.3 light/extinction utility implementation
  should decide whether to re-enable the term (and document the choice
  here).
