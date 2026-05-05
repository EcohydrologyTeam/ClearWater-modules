# v3 TSM Four-Way Audit Findings

**Date:** 2026-05-05
**Method:** Three parallel reviewer agents partitioned by physics domain.
**Sources audited:**
- **v3 TSM:** `src/clearwater_modules_v3/processes/temperature.py`
- **v1 TSM:** `src/clearwater_modules/tsm/processes.py`
- **Fortran-A (HEC-RAS-WQ):** `HEC-RAS-WQ/RAS-1D-WQ/Kinetics Libraries/{TemperatureEnergyBudget,TemperatureEquilibrium}/Source files/modTemperature.f90`
- **Fortran-B (WQM1D):** `ClearWater/WQM1D/TEMP/Source files/modHeatFlux.f90`

**Framing:** the legacy Fortran code was not adequately validated and likely contains flaws. v3 has been actively improving on both v1 and Fortran. Deviation from Fortran does not imply v3 is wrong; the audit is about reasoning which version is physically correct.

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

### F2 (MAJOR, Agent C). Depth ramp + rate cap break water-sediment energy conservation in shallow cells.

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

### F-Richardson-doc (MINOR, Agent B). Richardson docstring claims "removed −1 factor" but the sign is actually carried by `GRAVITY = -9.806`.

**Location:** `src/clearwater_modules_v3/processes/temperature.py:796-800`.

The v3 docstring asserts the leading `−1` was deleted "per Jason Rutyna's investigation" and that v1 had no such factor. The factual situation is:
- v3 and v1 both store `GRAVITY = -9.806` and use it directly: `richardson_number = GRAVITY * (rho_air − rho_air_sat) * 2 / (rho_air * U²)`.
- Fortran-A and Fortran-B store `gravity = +9.806` (positive) and apply an explicit leading `−` in the formula: `RichardsonNo = -gravity * (rho_air − rho_air_sat) * 2 / (rho_air * U²)`.
- Algebraically identical. The "removed −1" framing is misleading because it implies the sign was a bug-prone explicit factor, when in fact the sign is carried by the constant.

**Risk.** A future contributor "normalizing" `GRAVITY` to `+9.806` (the SI convention) without auditing the sign-bookkeeping would silently flip every Richardson regime in v3.

**Fix.** Reword the comment at lines 796-800 to make the sign convention explicit: "v3 stores `GRAVITY = -9.806 m/s²`; the Richardson formula uses it directly. v1 has the same convention. Fortran-A and Fortran-B store `gravity = +9.806` and apply an explicit `-gravity`. The two conventions produce algebraically identical Richardson numbers; do not 'normalize' `GRAVITY` to `+9.806` without re-auditing every consumer." Or, alternatively, flip `GRAVITY` to `+9.806` in `clearwater_modules_v3/utils/constants.py` and reintroduce the explicit `-` in the Richardson formula, matching Fortran convention. Either resolution is fine; the documentation defect is the priority.

### F-flux-sediment-doc (MINOR, Agent C). `flux_sediment` docstring mis-states what `/86400` converts.

**Location:** `src/clearwater_modules_v3/processes/temperature.py:402-404`.

Current text: *"The `/ 86400` converts the product of diffusivity (m²/s) and bulk thermal capacity into the per-substep flux units expected by the energy balance."*

This is wrong:
- `alphas` (sediment diffusivity) is in **m²/day** (correctly stated in the constructor docstring at lines 113-119).
- `/86400` converts day to seconds, giving the flux a `1/s` time component matching W/m² output.
- The flux is already in W/m²; "per-substep flux units" is misleading because no per-substep multiplication is happening here.

v1's wording at `tsm/processes.py:414` is correct: *"86400 converts the sediment thermal diffusivity from units of m²/d to m²/s"*. **Fix:** rewrite the v3 docstring to match v1's wording.

### F-sign-convention (MINOR, Agent A). Mixed pre-negate vs compose-time-negate in v3 flux returns.

**Location:** `src/clearwater_modules_v3/processes/temperature.py`. Sign convention is inconsistent across the flux methods:

- `flux_upwelling_longwave` (lines 308-315): pre-negated inside the function (returns negative for outgoing flux).
- `flux_latent_heat` (lines 349-369): pre-negated inside the function.
- `flux_sensible`, `flux_atmospheric_longwave`, `flux_sediment`, solar input: signs come from argument structure (e.g., `T_air − T_water` for sensible, `T_sed − T_water` for sediment).

`flux_net` (line 478) sums these: `sensible + solar + sediment + atmospheric + upwelling + latent`. Mathematically correct, but a reader auditing the signs has to remember which functions are pre-negated.

v1 and both Fortran sources use a **consistent** convention: all flux functions return magnitudes; signs are applied at composition time (`q_net = q_sens + q_solar + q_sed + q_LW_down − q_LW_up − q_latent`).

**Fix:** either (a) add a one-time docstring at the top of the energy-balance section documenting which functions are pre-negated, or (b) refactor to match v1/Fortran convention (return magnitudes; apply signs in `flux_net`). Option (a) is cheaper and lower-risk.

---

## 2. Defects in the legacy Fortran code (informational)

These are not v3 issues — v3 has already corrected them. Recorded here as a reviewer-facing artifact when v3 ships back to LimnoTech / ERDC, in case the Fortran modules are used as a baseline for any future validation.

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

## 3. Open questions for the human

Audit-surfaced questions that benefit from your input:

1. **F2 fix decision.** Option A (apply ramp + cap symmetrically to sediment-side delta, preserving energy conservation) or Option B (document the trade-off and ship the asymmetric version as-is)? Recommend Option A.
2. **TeqC equilibrium-temperature exposure.** Fortran-A implements a Newton-Raphson loop computing equilibrium temperature (`modTemperature.f90:209-262`). v3 does not expose `TeqC`. Was this an intentional simplification, or is `TeqC` needed as a future output?
3. **Per-component flux outputs.** Should v3 expose `q_latent`, `q_sensible`, `q_longwave_up`, `q_longwave_down`, `q_sediment` as registry diagnostics matching Fortran-A's `TempPathwayOutput`? Useful for calibration and validation; v3 currently bundles them into `flux_net` only.
4. **Wind-function provenance.** None of the four sources cites the calibration reference for `wind_a`, `wind_b`, `wind_c`, or the `1e6` divisor. Edinger Brady & Geyer (1974) is the most likely source. Recommend recovering and citing.
5. **`+273.16` Kelvin offset.** Confirmed intentional ("for testing consistency with v1") or inherited typo? Schedule a `273.15` migration for v4 with a single-pass parity test?
6. **v3 wind defaults.** Adopt v1's `0.3 / 1.5 / 3.0` defaults (matching v1 calibration) or maintain "explicit args required" for safety? This is a usability/design call.
7. **Sign-convention bookkeeping.** Refactor v3 to match v1/Fortran's "magnitudes-only, signs-at-composition" convention, or document the v3 mixed convention?

---

## 4. Recommended action plan

### Immediate (before next sponsor demo or LimnoTech PR):
- **Fix F2 (MAJOR):** apply depth ramp symmetrically to `sediment_temperature_change`. Add a regression test extending `test_water_sediment_energy_conservation_per_substep` to the guarded-path case.
- **Fix F-Richardson-doc (MINOR):** reword `temperature.py:796-800` to make the GRAVITY-sign convention explicit.
- **Fix F-flux-sediment-doc (MINOR):** rewrite `temperature.py:402-404` to match v1's correct wording.

### v3 1.0.x:
- F-sign-convention (MINOR): add docstring documenting which flux methods are pre-negated, or refactor to magnitudes-only convention.
- Open question 1 (TeqC) — design decision; default to "ship without TeqC" unless a specific use case appears.
- Open question 3 (per-component flux diagnostics) — small surface-area addition; useful for future calibration.

### v3.x or v4:
- 273.16 → 273.15 alongside v1-parity-test decommissioning (already deferred per Phase R-4).
- Wind-function provenance documentation.
- MMS energy-conservation test (was deferred at Phase R-5; an MMS test would have caught F2).

### Not in v3 scope:
- Fortran-side fixes (Lv-Kelvin, Fortran-B Richardson clamp): record in the LimnoTech reviewer materials and offer to share v3's corrections back. Not v3 work.

---

## 5. Per-agent reports

The full per-agent reports (~6,000 words combined) are available in the conversation transcript. This document consolidates and prioritizes their findings.

- **Agent A (energy fluxes):** zero findings against v3 (after self-correction); identified the Fortran Lv-Kelvin defect and the per-component-flux-output design question.
- **Agent B (thermodynamics + atmospheric stability):** zero CRITICAL/MAJOR against v3; identified the Richardson docstring issue and the Fortran-B Richardson-clamp issue.
- **Agent C (sediment + integration + constants):** identified F2 (the only MAJOR finding against v3) and the flux_sediment docstring inaccuracy.

---

**End of audit findings.** v3's energy-budget kernel is in good shape across the board. The single MAJOR defect (F2) is a real regression worth fixing, but it's contained, well-localized, and has a clear fix path (Option A, ~10-20 lines of code plus a regression test).
