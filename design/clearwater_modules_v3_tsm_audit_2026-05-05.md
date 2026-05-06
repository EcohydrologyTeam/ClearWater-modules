# v3 TSM Four-Way Audit Findings

**Date:** 2026-05-05
**Last updated:** 2026-05-05 (resolution status added)
**Method:** Three parallel reviewer agents partitioned by physics domain.
**Sources audited:**
- **v3 TSM:** `src/clearwater_modules_v3/processes/temperature.py`
- **v1 TSM:** `src/clearwater_modules/tsm/processes.py`
- **Fortran-A (HEC-RAS-WQ):** `HEC-RAS-WQ/RAS-1D-WQ/Kinetics Libraries/{TemperatureEnergyBudget,TemperatureEquilibrium}/Source files/modTemperature.f90`
- **Fortran-B (WQM1D):** `ClearWater/WQM1D/TEMP/Source files/modHeatFlux.f90`

**Framing:** the legacy Fortran code was not adequately validated and likely contains flaws. v3 has been actively improving on both v1 and Fortran. Deviation from Fortran does not imply v3 is wrong; the audit is about reasoning which version is physically correct.

## Resolution status (added 2026-05-05)

All findings against v3 and all open questions have been resolved on the
`streaming` branch. The table below maps each item to its commit hash; full
details are inline at each item below.

| ID | Severity | Status | Commit | Title |
|---|---|---|---|---|
| F2 | MAJOR | Resolved | `808facb` | depth-ramp/rate-cap break water-sediment energy conservation |
| F-Richardson-doc | MINOR | Resolved | `808facb` | Richardson docstring sign-convention rewording |
| F-flux-sediment-doc | MINOR | Resolved | `808facb` | `/86400` docstring wording fix |
| F-sign-convention | MINOR | Resolved | `9d8ebfc` | flux methods refactored to magnitudes-only |
| Q1 (F2 fix) | — | Option A | `808facb` | symmetric ramp+cap on sediment side |
| Q2 (TeqC) | — | Resolved | `a00065a` | vectorized Newton-Raphson, opt-in registry diagnostic |
| Q3 (per-component fluxes) | — | Resolved | `a00065a` | seven `q_*` outputs + `flux_components` dict |
| Q4 (wind-function provenance) | — | Resolved | `a8874e7` | Edinger Brady & Geyer (1974) cited |
| Q5 (273.16 → 273.15) | — | Resolved | `a8874e7` | v3 SI Kelvin offset; v3 NSM1 modules already correct |
| Q6 (v1 wind defaults) | — | Resolved | `a8874e7` | `wind_a/b/c = 0.3 / 1.5 / 3.0` defaults |
| Q7 (sign-convention) | — | Resolved | `9d8ebfc` | refactored to v1/Fortran magnitudes-only |

Test status after all fixes: 248 passed, 1 xfailed on the v3 suite (was 199 at audit time; +49 new regression tests).

**Fortran-side defects** (1 CRITICAL in both Fortrans, 1 MAJOR in Fortran-B) are tracked in §2 below as informational findings against the audited Fortran reference snapshots. They are out of scope for v3 Python work, and the upstream Fortran sources are not maintained as part of this project. v3 Python already has the correct forms of both calculations.

---

## Executive summary

Three parallel agents covering energy fluxes (Agent A), thermodynamics + atmospheric stability (Agent B), and sediment/integration/constants (Agent C) collectively examined every formula, coefficient, and sign convention in the four sources. The four implementations agree on the structural form of every audited quantity. **Most divergences between v3 and Fortran are v3 improvements** (corrected unit bugs, defensive guards, NaN handling). **Two genuine v3 defects were identified**, one MAJOR and one MINOR, both in code I added recently.

**Findings tally against v3 (the code under audit):**
- 0 CRITICAL
- 1 MAJOR (F2: depth-ramp/rate-cap break sediment-water energy conservation in shallow cells)
- 3 MINOR (Richardson sign-convention docstring; flux_sediment `/86400` docstring wording; mixed sign-convention bookkeeping in flux returns)

**Findings tally against the legacy Fortran sources:**
- 1 CRITICAL (both Fortrans): Lv polynomial fed Kelvin instead of Celsius — ~26 % underestimate at 20 °C
- 1 MAJOR (Fortran-B only): Richardson out-of-bounds constants `12.3` / `0.03` returned with the clamp commented out — non-physical discontinuity
- Several MINOR (273.16 vs 273.15 offset; missing mixing-ratio guard; missing degenerate-`sediment_thickness` guard; missing NaN propagation through Richardson)

**v3 verified-correct items:** all energy-flux formulas, all thermodynamic polynomials (Brutsaert e_sat, latent-heat polynomial coefficients in correct Celsius form, water density, water cp, air density, density_air_sat, Richardson formula and stability function with correct cutoffs/exponents, wind function structure), all constants, sediment heat flux formula, dynamic T_sed evolution algebra (matches Fortran exactly), water-T per-substep update equivalence to v1 and both Fortrans, thin-water guards algebraically equivalent to v1.

---

## 1. v3 defects identified by the audit

### F2 (MAJOR, Agent C) — **RESOLVED in `808facb`**. Depth ramp + rate cap break water-sediment energy conservation in shallow cells.

**Resolution:** Option A applied. `temperature_change` was refactored to
delegate to a new `_temperature_change_with_factors` helper that returns
`(delta, ramp, clip_ratio)`. `Temperature.run` now multiplies the
sediment-side delta by `ramp * clip_ratio` so the guards apply
symmetrically. Defaults make `ramp = 1` and `clip_ratio = 1` in typical
deep cells; the fix is a no-op outside the guard regions. Five new
regression tests pin the contract:
`test_F2_helper_returns_correct_ramp_when_active`,
`test_F2_helper_returns_ramp_one_when_disabled_or_deep`,
`test_F2_helper_returns_clip_ratio_when_cap_fires`,
`test_F2_water_sediment_conservation_under_depth_ramp`,
`test_F2_sediment_delta_scales_by_ramp_and_clip_ratio_in_run`.

**Location:** `src/clearwater_modules_v3/processes/temperature.py`. The depth ramp and rate cap are applied to the water-side per-substep delta T inside `temperature_change` (lines 564-585), but `sediment_temperature_change` (lines 617-627) applies neither.

**Mechanism.** When `depth < q_net_depth_ramp_ref` (e.g., a 0.05 m cell with the default 0.3 m reference depth, ramp = 1/6 ≈ 0.167):

- Water column receives `q_sediment · ramp` of energy per substep.
- Sediment loses energy at the **unramped** rate: `dE_sediment = pb · Cps · h2 · ΔT_s = pb · Cps · α / (0.5 · h2) · (T_w − T_s) · dt / 86400`.
- The pair-cancellation derived in the C10 fix derivation (`design/clearwater_modules_v3_review_findings.md` C10 closure) **fails by a factor of `(1 − ramp)`**.

In the example above, sediment relaxes toward water 6× faster than water absorbs the corresponding energy. T_sed converges toward T_water at a non-physical rate while T_water lags. The result: a one-way energy sink at the wet/dry margin proportional to the ramp deficit and the persistent T_sed − T_water gradient.

**Same defect for the rate cap.** When `dTdt_max_per_hour` clips the water-side delta but the sediment-side delta is uncapped, the cancellation breaks proportionally to the clip ratio.

**Real-world impact:**
- For deep riverine cells (depth ≫ 0.3 m), `ramp = 1` and the cap rarely fires; the defect is dormant.
- For shallow / floodplain / drying cells, the bias is non-trivial and accumulates monotonically. This is exactly v3's target use case (mission-critical riverine simulation).
- For Sumwere Creek (where I observed the wet/dry margin artifact), the depth-ramp-active cells coincide with the same wet/dry margin cells.

**How this got introduced.** The C10 fix (dynamic sediment-T evolution restored from Fortran) was added in commit `200a26e` (sediment heat-exchange parity with canonical Fortran TSM). The energy-conservation derivation in that commit's test (`tests/v3/test_tsm_sediment_v3.py::test_water_sediment_energy_conservation_per_substep`) verifies the unguarded-path cancellation. The thin-water guards predated the dynamic-sediment-T work; nobody re-verified the conservation contract under the guards. This is exactly the kind of invariant breakage that an end-to-end MMS test would have caught — and which Phase R-5 deferred as `MMS test for energy conservation` to v3.x NSM1 work.

**Two recommended fixes (your choice):**

**Option A — Symmetric application (preferred for energy conservation):**
Apply the same `ramp` factor to the sediment-side delta inside `sediment_temperature_change`. This preserves the per-cell water-sediment cancellation invariant. The rate cap should similarly apply proportionally to both sides (clip the `(T_water − T_sed)` driving force when the water-side cap fires).

**Option B — Document the trade-off:**
Acknowledge that the thin-water guards intentionally relax the energy-conservation invariant in shallow cells in exchange for stability. Add a docstring at `temperature_change` noting this and add a regression test that quantifies the bias.

**My recommendation: Option A.** The depth ramp's purpose is to attenuate the surface flux on shallow cells where the explicit-Euler step would be unstable; physically, that means slowing the heat transfer at both interfaces (water-air via the surface flux, water-sediment via the sediment flux). Applying the ramp symmetrically to both `flux_net` and `sediment_temperature_change` matches that physical intent and preserves conservation. The rate cap is more nuanced (it caps the resulting ΔT, not a flux), so cap-symmetry needs care — applying the same proportional clip to ΔT_s as to ΔT_w preserves cancellation.

### F-Richardson-doc (MINOR, Agent B) — **RESOLVED in `808facb`**. Richardson docstring claims "removed −1 factor" but the sign is actually carried by `GRAVITY = -9.806`.

**Resolution:** Docstring at `temperature.py:796-800` rewritten to make the
GRAVITY sign convention explicit: "v3 stores GRAVITY = -9.806 m/s²; the
formula uses it directly. v1 has the same convention. Fortran-A and
Fortran-B store gravity = +9.806 and apply an explicit `-gravity`. Do
not 'normalize' GRAVITY to +9.806 without re-auditing every consumer."

**Location:** `src/clearwater_modules_v3/processes/temperature.py:796-800`.

The v3 docstring asserts the leading `−1` was deleted "per Jason Rutyna's investigation" and that v1 had no such factor. The factual situation is:
- v3 and v1 both store `GRAVITY = -9.806` and use it directly: `richardson_number = GRAVITY * (rho_air − rho_air_sat) * 2 / (rho_air * U²)`.
- Fortran-A and Fortran-B store `gravity = +9.806` (positive) and apply an explicit leading `−` in the formula: `RichardsonNo = -gravity * (rho_air − rho_air_sat) * 2 / (rho_air * U²)`.
- Algebraically identical. The "removed −1" framing is misleading because it implies the sign was a bug-prone explicit factor, when in fact the sign is carried by the constant.

**Risk.** A future contributor "normalizing" `GRAVITY` to `+9.806` (the SI convention) without auditing the sign-bookkeeping would silently flip every Richardson regime in v3.

**Fix.** Reword the comment at lines 796-800 to make the sign convention explicit: "v3 stores `GRAVITY = -9.806 m/s²`; the Richardson formula uses it directly. v1 has the same convention. Fortran-A and Fortran-B store `gravity = +9.806` and apply an explicit `-gravity`. The two conventions produce algebraically identical Richardson numbers; do not 'normalize' `GRAVITY` to `+9.806` without re-auditing every consumer." Or, alternatively, flip `GRAVITY` to `+9.806` in `clearwater_modules_v3/utils/constants.py` and reintroduce the explicit `-` in the Richardson formula, matching Fortran convention. Either resolution is fine; the documentation defect is the priority.

### F-flux-sediment-doc (MINOR, Agent C) — **RESOLVED in `808facb`**. `flux_sediment` docstring mis-states what `/86400` converts.

**Location:** `src/clearwater_modules_v3/processes/temperature.py:402-404`.

Current text: *"The `/ 86400` converts the product of diffusivity (m²/s) and bulk thermal capacity into the per-substep flux units expected by the energy balance."*

This is wrong:
- `alphas` (sediment diffusivity) is in **m²/day** (correctly stated in the constructor docstring at lines 113-119).
- `/86400` converts day to seconds, giving the flux a `1/s` time component matching W/m² output.
- The flux is already in W/m²; "per-substep flux units" is misleading because no per-substep multiplication is happening here.

v1's wording at `tsm/processes.py:414` is correct: *"86400 converts the sediment thermal diffusivity from units of m²/d to m²/s"*. **Fix:** rewrite the v3 docstring to match v1's wording.

### F-sign-convention (MINOR, Agent A) — **RESOLVED in `9d8ebfc`**. Mixed pre-negate vs compose-time-negate in v3 flux returns.

**Resolution:** Refactored to v1 / Fortran magnitudes-only convention.
`flux_upwelling_longwave` now returns the positive Stefan-Boltzmann
magnitude; `flux_latent_heat` returns the unnegated form (positive in
the evaporative regime, negative for condensation). `flux_net`
composition is now `+ sensible + solar + sediment + LW_down - LW_up
- latent`, matching v1 `tsm/processes.py:q_net` and Fortran-A
`modTemperature.f90:257`. Math is identical end-to-end. New file
`tests/v3/test_tsm_sign_convention_v3.py` (7 tests) pins the contract.

**Location:** `src/clearwater_modules_v3/processes/temperature.py`. Sign convention is inconsistent across the flux methods:

- `flux_upwelling_longwave` (lines 308-315): pre-negated inside the function (returns negative for outgoing flux).
- `flux_latent_heat` (lines 349-369): pre-negated inside the function.
- `flux_sensible`, `flux_atmospheric_longwave`, `flux_sediment`, solar input: signs come from argument structure (e.g., `T_air − T_water` for sensible, `T_sed − T_water` for sediment).

`flux_net` (line 478) sums these: `sensible + solar + sediment + atmospheric + upwelling + latent`. Mathematically correct, but a reader auditing the signs has to remember which functions are pre-negated.

v1 and both Fortran sources use a **consistent** convention: all flux functions return magnitudes; signs are applied at composition time (`q_net = q_sens + q_solar + q_sed + q_LW_down − q_LW_up − q_latent`).

**Fix:** either (a) add a one-time docstring at the top of the energy-balance section documenting which functions are pre-negated, or (b) refactor to match v1/Fortran convention (return magnitudes; apply signs in `flux_net`). Option (a) is cheaper and lower-risk.

---

## 2. Defects in the audited Fortran reference snapshots (informational)

These are not v3 issues — v3 has already corrected them. The Fortran files audited on 2026-05-05 are reference snapshots only; they are not maintained as part of this project, and their current upstream state is unknown. Recorded here in case anyone tries to use those snapshots as a baseline in the future and needs to know which numerical conventions they should NOT inherit.

### Fortran-A and Fortran-B: Lv polynomial fed Kelvin instead of Celsius (CRITICAL)

**Location:** Fortran-A `modTemperature.f90:393-397`; Fortran-B `modHeatFlux.f90:399-403`.

```fortran
mf_latent_heat_vaporization = 2499999 - 2385.74 * TwaterK    ! Kelvin input
```

The polynomial `2,499,999 − 2385.74 · T` is calibrated for Celsius (`Lv ≈ 2.50 MJ/kg @ 0 °C`, `2.45 MJ/kg @ 20 °C`). Feeding Kelvin produces nonsensical values: at 20 °C (`T_K = 293.15`), Fortran computes `Lv ≈ 2,499,999 − 2385.74 · 293.15 = 1.80 × 10⁶ J/kg` — about **27 % below the physical value**. This systematically underestimates evaporative cooling and biases simulated water temperature warm.

v3 (`temperature.py:654`) and v1 streaming-branch (`processes.py:213-228`) both apply the polynomial in Celsius. **v3 is correct; both Fortrans are wrong.**

### Fortran-B only: Richardson out-of-bounds constants (MAJOR)

**Location:** Fortran-B `modHeatFlux.f90:353-391`.

The Richardson clamp to `[−1, 2]` (present in Fortran-A and the Python implementations) is commented out in Fortran-B, replaced by a fallthrough returning constant values `12.3` and `0.03` for `Rn < −1` and `Rn ≥ 2` respectively. These constants approximately match the function evaluated at the bounds (`(1 − 22·(−1))^0.80 ≈ 13.5`; `(1 + 34·2)^(−0.80) ≈ 0.0335`) but are non-physical: they introduce a discontinuity between the in-range branch and the constants. Whether `12.3` / `0.03` were calibration values matching the function at the threshold or debug placeholders is unclear from the source.

v3 (with the clamp-then-evaluate structure inherited from v1 and matching Fortran-A) is correct.

### Several minor Fortran omissions

- **No `e_air > P_air` guard in `mixing_ratio_air`** (Fortran-A, Fortran-B, v1). v3's guard at `temperature.py:691-729` defends against forcing-data errors that produce negative denominators. v3 is the only source with this protection.
- **No NaN propagation through `richardson_function`** (all three other sources). v3 forces `Ri_function = NaN` when `Ri_number = NaN` at lines 864-868 so missing meteorology forcing produces a visible defect rather than a finite-but-wrong stability multiplier.
- **No degenerate-`sediment_thickness ≤ 0` guard** in either Fortran source or v1. v3's guard at `temperature.py:417-435` and `610-627` returns 0 rather than producing inf/NaN that would poison the entire field.
- **`+273.16` vs `+273.15` Kelvin offset** in all four sources. The triple point of water is 273.16 K; the SI absolute-temperature offset for 0 °C is 273.15 K. The 0.01 K bias propagates into every Kelvin-evaluated quantity (`e_sat`, longwave T⁴ terms, air densities). Relative error ~3.4 × 10⁻⁵ at 293 K, well below other model uncertainties, but worth documenting for v4.
- **Fortran wind-function defaults `1.0/1.0/1.0`** are likely placeholder values intended to be overridden by application-supplied parameter files. v1 uses `0.3/1.5/3.0`; v3 has no default and forces explicit constructor args (a defensible safety posture).

---

## 3. Open questions for the human — **ALL RESOLVED 2026-05-05**

Audit-surfaced questions, with the user's resolution and shipping commit:

1. **F2 fix decision** — **Option A** (symmetric ramp+cap on sediment side).
   Shipped in `808facb`.
2. **TeqC equilibrium-temperature exposure** — **Implement.** Vectorized
   Newton-Raphson, opt-in registry diagnostic. Shipped in `a00065a` as
   `Temperature.equilibrium_temperature`.
3. **Per-component flux outputs** — **Implement.** Seven Fortran-A pathway
   outputs (`q_sensible`, `q_latent`, `q_longwave_up`, `q_longwave_down`,
   `q_solar`, `q_sediment`, `q_net`) cached on the process and
   opportunistically written to the registry when pre-registered. Shipped
   in `a00065a` as `Temperature.flux_components`.
4. **Wind-function provenance** — **Recover and cite.** Edinger, J.E.,
   D.K. Brady, and J.C. Geyer (1974), *Heat exchange and transport in the
   environment*, Report 14, Cooling Water Discharge Research Project
   (RP-49), Electric Power Research Institute, Palo Alto, CA, 125 pp.
   Citation added to `Temperature.__init__` and `wind_function` docstrings
   in `a8874e7`.
5. **`+273.16` Kelvin offset** — **Switch to 273.15** (SI canonical
   offset; 273.16 is the triple point, the wrong reference for 0 °C).
   v3 `utils/conversions.py` now defines `celsius_to_kelvin` locally with
   `KELVIN_OFFSET = 273.15` rather than re-exporting v2's 273.16. v3
   NSM1 modules already used literal 273.15 and are unchanged. v2 is
   untouched (still 273.16 internally for v1-utility-parity tests).
   Shipped in `a8874e7`.
6. **v3 wind defaults** — **Adopt v1's `0.3 / 1.5 / 3.0`.** `Temperature`
   constructor no longer requires explicit wind args; v1-aligned defaults
   apply unless overridden. Shipped in `a8874e7`.
7. **Sign-convention bookkeeping** — **Refactor to magnitudes-only.**
   Shipped in `9d8ebfc` (see F-sign-convention above).

---

## 4. Action plan — **status as of 2026-05-05**

### Immediate (shipped in this session):
- ✅ **F2 (MAJOR):** `808facb` — symmetric ramp+cap on sediment side, 5 regression tests.
- ✅ **F-Richardson-doc (MINOR):** `808facb` — `temperature.py` Richardson docstring rewritten.
- ✅ **F-flux-sediment-doc (MINOR):** `808facb` — `/86400` wording corrected.
- ✅ **F-sign-convention (MINOR):** `9d8ebfc` — flux methods refactored to v1/Fortran magnitudes-only convention; new `tests/v3/test_tsm_sign_convention_v3.py` (7 tests).
- ✅ **Q1 (F2 fix decision):** Option A applied.
- ✅ **Q2 (TeqC):** `a00065a` — `Temperature.equilibrium_temperature` Newton-Raphson.
- ✅ **Q3 (per-component fluxes):** `a00065a` — `Temperature.flux_components` returns dict; cached on instance and opportunistically written to registry.
- ✅ **Q4 (wind-function provenance):** `a8874e7` — Edinger Brady & Geyer (1974) cited.
- ✅ **Q5 (273.16 → 273.15):** `a8874e7` — v3 SI Kelvin offset.
- ✅ **Q6 (v1 wind defaults):** `a8874e7` — `Temperature` defaults `0.3 / 1.5 / 3.0`.
- ✅ **Q7 (sign-convention refactor):** `9d8ebfc` (same as F-sign-convention).

### v3.x or v4:
- MMS energy-conservation test (was deferred at Phase R-5; an MMS test would have caught F2).
- LimnoTech reconciliation of the 273.16 vs 273.15 convention if v2 parity is renegotiated (currently v2 retains 273.16 internally).

### Not in v3 scope:
- **Fortran-side fixes** (Lv-Kelvin in both Fortrans, Fortran-B Richardson clamp) — out of scope for the v3 Python port. The audited Fortran files are reference snapshots, not maintained as part of this project, and may not represent current upstream state. v3 Python already has the correct forms; §2 below records the legacy defects for context only.

---

## 5. Per-agent reports

The full per-agent reports (~6,000 words combined) are available in the conversation transcript. This document consolidates and prioritizes their findings.

- **Agent A (energy fluxes):** zero findings against v3 (after self-correction); identified the Fortran Lv-Kelvin defect and the per-component-flux-output design question.
- **Agent B (thermodynamics + atmospheric stability):** zero CRITICAL/MAJOR against v3; identified the Richardson docstring issue and the Fortran-B Richardson-clamp issue.
- **Agent C (sediment + integration + constants):** identified F2 (the only MAJOR finding against v3) and the flux_sediment docstring inaccuracy.

---

**End of audit findings.** v3's energy-budget kernel is in good shape across the board. The single MAJOR defect (F2) is a real regression worth fixing, but it's contained, well-localized, and has a clear fix path (Option A, ~10-20 lines of code plus a regression test).
