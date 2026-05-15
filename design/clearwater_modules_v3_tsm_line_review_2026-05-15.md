# v3 TSM Line-Level Code Review

**Date:** 2026-05-15
**Reviewer:** Todd Steissberg (ERDC)
**Scope:** `src/clearwater_modules_v3/processes/temperature.py` (1,651 lines)
reviewed line-by-line against the v1 reference
(`src/clearwater_modules/tsm/processes.py`),
`src/clearwater_modules/tsm/constants.py`,
`src/clearwater_modules_v3/utils/constants.py`, and
`src/clearwater_modules_v3/utils/conversions.py`.
**Companion documents:** `clearwater_modules_v3_tsm_design_specification.md`,
`clearwater_modules_v3_tsm_gap_analysis.md`,
`clearwater_modules_v3_tsm_audit_2026-05-05.md`,
`clearwater_modules_v3_tsm_wind_function_specification.md`.
**Note:** This review was performed independently of (and concurrently
with) the parallel NSM1 v3 review. The shared orchestration layer
(`model.py`) is in NSM1 review scope; this document covers TSM only,
with the one exception noted in §5 (chunking TODO resolution spot-check).

---

## Verdict

The xarray refactoring is **complete**, the algorithms **match v1**
except where v3 deliberately improves on them, the documentation is
**thorough and accurate**, and **no live TODO or unresolved-issue
comment remains** in `temperature.py`. Every fixed defect carries a
comment stating it was corrected and citing the fix
(`F2 fix (audit 2026-05-05)`, `m2 / M3 / m6 / m9 fix
(review-findings 2026-05-04)`, etc.). The single `TODO` token in the
file (line 1580) is a *historical* note recording that the Richardson
`-1` factor "was resolved per Jason Rutyna's January 2026 diff
investigation" — correctly worded as resolved, not outstanding.

Core v3 TSM test files pass: **88 passed** across
`test_tsm_sign_convention_v3`, `test_tsm_sediment_v3`,
`test_tsm_latent_heat_v3`, `test_tsm_stability_ramp_v3`,
`test_5_tsm_calculations_v3`, `test_tsm_wind_function_spec_v3`. The 32
warnings emitted are the *intentional* `wind_c = 3.0` UserWarning from
the v1-parity stability-ramp ports — expected behavior, not a defect.

---

## 1. Algorithm equivalence vs. v1 (all verified)

| Quantity | v1 ↔ v3 | Notes |
|---|---|---|
| `latent_heat_vaporization` | Equivalent | v3 takes Celsius directly; v1 converts K→C first. Same polynomial (`2,499,999 − 2385.74·T_C`), same result. v1 latent-heat unit fix correctly carried. |
| `saturation_vapor_pressure` | Identical | Horner form; Brutsaert (1982) coefficients `A0…A6` match `tsm/constants.py` exactly. |
| `richardson_number` + stability function | Identical | Formula identical (`GRAVITY = −9.806` carries the sign in both v1 and v3; no explicit leading `−1`). Bounds `[−1, 2]`, regime cutoffs (`±0.01`), exponents (`0.80`, `−0.80`), coefficients (`22`, `34`) all match. v3 adds NaN propagation (improvement, M3). |
| `wind_function` | Reduces to v1 | Edinger/Brady/Geyer form identical; v3's shelter + log-law height correction is a no-op under default constructor values. |
| `flux_sediment`, `flux_sensible`, `flux_latent_heat`, `flux_atmospheric_longwave`, `flux_upwelling_longwave`, `flux_net` | Identical | `/0.5/h2`, `/86400`, and the `+sensible +solar +sediment +LW_down −LW_up −latent` composition all match v1 `q_net` and Fortran-A. |
| `temperature_change` | Algebraically identical | Equivalent to v1 `dTdt_water_c`; the `86400·dt` factor correctly relocated from `q_net` (v1) to `temperature_change` (v3) per gap-analysis N7. Depth ramp and rate cap match v1 with the same disable semantics (`q_net_depth_ramp_ref=0.0`, `dTdt_max_per_hour=+inf`). |
| `water_density`, `water_specific_heat`, `density_air`, `density_air_sat`, `mixing_ratio_air` | Identical + guards | Formulas identical; v3 adds defensive guards (degenerate-thickness, denominator ≤ 0, NaN propagation) that only change behavior on already-pathological inputs. |

### Independently re-derived invariants

- **Water–sediment energy conservation.** Confirmed the implicit
  active-layer convention `C_sed = ρ_b·c_ps·h2` (temperature gradient
  taken over the half-thickness `0.5·h2`, thermal storage over the
  full `h2`). The F2 fix is mathematically sound:
  `delta_clipped = delta_unclipped · clip_ratio` *exactly* (the clip
  is a uniform proportional scaling of the linear flux sum), so
  applying the same `ramp · clip_ratio` to the sediment-side delta
  preserves per-substep `dE_water + dE_sediment = 0` under both the
  depth ramp and the rate cap.
- **Equilibrium-temperature Newton-Raphson.** Derivative signs are
  correct (`d_qnet_dT = −d_upwelling −d_latent +d_sensible
  +d_sediment`, strictly negative for realistic forcing). The
  Brutsaert-polynomial derivative `d_esat_dT` is the correct
  termwise derivative. Kelvin-based derivatives with a Celsius Newton
  step are valid (the offset is constant, so `d/dT_C ≡ d/dT_K`).

---

## 2. Audit findings — independently confirmed resolved in code

The 2026-05-05 audit claimed all findings resolved. Verified against
the current source (not trusting the audit document):

| Finding | Status in code | Evidence |
|---|---|---|
| F2 (MAJOR) — sediment energy conservation under guards | Resolved | `_temperature_change_with_factors` returns `(delta, ramp, clip_ratio, components)`; `run` applies `ramp_factor * cap_clip_ratio` to the sediment delta (lines 504–524). |
| F-Richardson-doc (MINOR) | Resolved | Lines 1569–1582 make the `GRAVITY = −9.806` sign convention explicit and warn against "normalizing" to `+9.806`. |
| F-flux-sediment-doc (MINOR) | Resolved | `flux_sediment` docstring (lines 708–718) correctly states `/86400` converts `sediment_diffusivity` from m²/day to m²/s. |
| F-sign-convention (MINOR) | Resolved | Magnitudes-only convention block (lines 585–607); `flux_upwelling_longwave` and `flux_latent_heat` return unnegated magnitudes; signs applied in `flux_components`. |
| Q1–Q5, Q7 | Resolved | Option A (F2); `equilibrium_temperature` Newton-Raphson; `flux_components` 7-output dict; Edinger/Brady/Geyer (1974) citation; `KELVIN_OFFSET = 273.15`; magnitudes-only refactor. |

Referenced fix commits all exist: `808facb`, `9d8ebfc`, `a00065a`,
`a8874e7`, `200a26e`.

---

## 3. Most significant intentional divergence from v1

**`wind_c` default is `2.0`, not v1's `3.0`.** This is the single
largest numerical difference from v1 and is **deliberate**:

- Extensively justified in the constructor docstring (CE-QUAL-W2 manual
  default `CFW = 2.0`; QUAL2K Brady-Graves-Geyer default `2.0`; all
  seven W2 example case studies use `CFW = 2.0`).
- Backed by `clearwater_modules_v3_tsm_wind_function_specification.md`.
- Validated: the Willamette River Santiam-Salem case study results
  looked good.
- `c = 3.0` remains accepted at the upper bound for back-compat with a
  `UserWarning`.

Reviewers and downstream users should be aware that v3 ≠ v1 here by
design. The v1-parity test ports that pass `wind_c=3.0` explicitly
still pass (and correctly emit the warning).

---

## 4. Stale documentation (not code)

The audit findings document
(`clearwater_modules_v3_tsm_audit_2026-05-05.md`, Q6 at lines 31 and
222–224) and the gap analysis still state the v3 wind defaults as
`0.3 / 1.5 / 3.0`. The code moved the `wind_c` default to `2.0` in
later commits (`201a364` → `42f381b` → `7efa946`). These are
historical findings documents frozen at their dates, so this is
expected — but if the `3.0` figure is cited in
sponsor/LimnoTech-facing materials it is now incorrect. Recommend a
one-line "superseded by the wind-function specification (wind_c = 2.0)"
annotation in those documents where they state the default.

---

## 5. Genuinely deferred items (all documented, none blocking)

1. **`np.select` dask incompatibility** in `water_specific_heat`
   (Gemini review finding 4). No dask dispatch via xarray's ufunc
   registry; materializes the chunk eagerly. Inert under the current
   in-memory `Model`; the limitation and the migration path
   (`xr.where` chain or `xr.apply_ufunc(..., dask="allowed")`) are
   documented in-code. v1 used `np.select` here too — not a
   regression.
2. **MMS energy-conservation test** deferred to v3.x/v4 (audit §4). An
   end-to-end MMS test would have caught F2 before it shipped; this is
   the highest-value deferred test.
3. **273.15-vs-273.16 LimnoTech reconciliation** if v2 numerical
   parity is renegotiated (v2 retains 273.16 internally).
4. **Wind `a`/`b` calibration study** — "tracked as future work" in
   the `__init__` docstring. The `a = 0.3, b = 1.5` magnitude
   coefficients are inherited from v1 and are unit-coupled to the
   `/1e6` normalization and to `c`; a focused calibration against
   observations is recommended but not done.

Spot-checked the shared `model.py`: the four v2 `__process_loop_chunked`
TODOs (gap-analysis O7) are genuinely resolved with a documented
resolution and the C7 integer-step-index fix — not left as TODOs.

---

## 6. Minor consistency note (optional, not a bug)

`_temperature_change_with_factors` computes `volume / surface_area`
inside an `xr.where` without an `np.errstate` guard, so dry/zero-area
cells can emit a benign `RuntimeWarning` before being masked to `0.0`.
`richardson_number` suppresses the analogous divide-by-zero warning
via `np.errstate`; v1's `dTdt_water_c` also did not suppress it here.
Harmonizing is optional and purely cosmetic — the masked result is
numerically correct either way.

---

## Conclusion

Nothing in this review blocks the current state of v3 TSM. The xarray
refactoring and v2-severance are complete for TSM (no v1/v2 imports
remain in `temperature.py`). The F2 conservation fix and the
magnitudes-only sign refactor are genuine improvements over v1,
correctly implemented and test-covered. The principal action items are
documentation hygiene (§4) and the deferred MMS test (§5.2), not code
defects.
