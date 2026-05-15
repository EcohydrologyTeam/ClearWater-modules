# v3 NSM1 CBOD and POM -- Source-Code and Science-Correctness Review

**Review date:** 2026-05-15
**Reviewer:** Claude Code (scientific software review agent, claude-sonnet-4-6)
**Branch:** `streaming`
**Commit:** 54f2b12
**Scope:** CBOD and POM processes only

| File | LOC | Role |
|---|---|---|
| `src/clearwater_modules_v3/processes/cbod.py` | 328 | CBOD process class |
| `src/clearwater_modules_v3/processes/pom.py` | 487 | POM process class |
| `src/clearwater_modules_v3/parameters/cbod.py` | 51 | CBOD parameter defaults |
| `src/clearwater_modules_v3/parameters/pom.py` | 28 | POM parameter defaults |

**V1 reference read:** `src/clearwater_modules/nsm1/processes.py` lines 2185--2329 (POM) and 2334--2434 (CBOD).

**Authoritative docs read:**
- `design/clearwater_modules_v3_nsm1_audit_simple_constituents.md`
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` (targeted sections via grep)

---

## 1. Summary Verdict

Both the CBOD and POM process classes are scientifically sound and algorithmically faithful to v1 under the conditions that are expected to occur in practice (default `ksbod_20 = 0`, `vb = 6.85e-6 m/d`). The xarray refactor is structurally complete: no Python-level cell loops, no scalar-only logic in hot paths, no `== np.nan` antipatterns, and no raw `if` tests on DataArray values inside the kernel. Forward Euler integration and `clip_negative_state` are applied consistently. The CBOD `ksbod_tc / depth` form is a documented, intentional architectural deviation from v1 that is silent under defaults and carries a known recalibration requirement for nonzero `ksbod_20` users.

The review finds **zero Critical findings**, **one Major finding**, **two Minor findings**, and **three Observations**. The single Major finding is a misleading diagnostic name (`pom_settling_rate` actually stores burial, not water-column settling) that exposes incorrect output if any downstream code or user analysis consumes that registry slot by name. The two Minor findings concern (1) a legacy `rate()` method in `pom.py` that duplicates logic without being guarded as deprecated, and (2) the `use_DOX` flag not being configurable via the DEFAULTS dict or `parameters` constructor argument. Both stale FIXME markers required by the scope mandate (one in `parameters/cbod.py` and one in `parameters/pom.py`) were confirmed as resolved and correctly labeled "FIXME cleared."

---

## 2. Findings

### Critical

None.

### Major

**F1 -- `pom_settling_rate` diagnostic name is misleading and maps to the wrong physical process**

- **File:line:** `src/clearwater_modules_v3/processes/pom.py:149`, `pom.py:355`, `pom.py:400`
- **Category:** Documentation-to-code fidelity / interface contract
- **Description:** The `REGISTRY_DIAGNOSTICS` tuple at line 149 declares `"pom_settling_rate"` as a named diagnostic. Inside `_change_with_components`, the variable `rate_burial` is computed as `self.vb * pom / self.h2` (line 355) -- the sediment burial velocity term -- and is then stored in `components["pom_settling_rate"]` at line 400. This is not a settling flux; it is burial removal of bed-sediment POM out of the `h2` layer. Water-column POC settling into the POM compartment (`rate_poc_settling`) is a separate, correctly named intermediate but is not exposed in `REGISTRY_DIAGNOSTICS` at all.

  Any downstream code, diagnostic consumer, or post-processing script that reads `pom_settling_rate` from the registry will receive the burial rate and interpret it as a settling rate, with incorrect physical attribution. The dual-path legacy `rate()` method (lines 410--487) makes the same assignment via the same `components` dict pattern even though it predates the refactor, so the error is consistent across both code paths.

- **Consequence:** Incorrect physical attribution in all diagnostic outputs that reference `pom_settling_rate`. Budget closure checks that split POM loss into "settling" and "burial" components will produce wrong partitioning. The error is not silent: consumers of the registry slot will read a value, but its physical meaning will be wrong.
- **Recommendation:** Rename `"pom_settling_rate"` to `"pom_burial_rate"` in `REGISTRY_DIAGNOSTICS` (line 149), in the `components` dict key at line 400, and in the parallel assignment inside `rate()` at the corresponding line (approximately line 444 region -- `rate_burial` returned as `"pom_settling_rate"` in the legacy path also). If the POC-to-POM settling flux is of diagnostic interest, add `"pom_poc_settling_rate"` as a separate entry and populate it from `rate_poc_settling`. Update any consumer that reads `pom_settling_rate` by that name.

---

### Minor

**F2 -- Legacy `rate()` method in `pom.py` duplicates `_change_with_components` without a deprecation guard**

- **File:line:** `src/clearwater_modules_v3/processes/pom.py:410--487`
- **Category:** Code quality / maintainability
- **Description:** The `rate()` method is a 78-line near-verbatim copy of the logic in `_change_with_components`. The class docstring (Phase 7 note at lines 325--329) acknowledges the duplication and states that `rate()` "retains its prior behaviour for back-compat with external tests that call `pom.rate(...)` directly." However, the method carries no `@deprecated` decorator, no `DeprecationWarning` call, and no `# noqa: deprecated` annotation. The duplicate also writes `self.pom_doc_source_rate` as a side effect (line 439), meaning that a caller who invokes `rate()` directly (outside the normal `run` path) will silently overwrite the cached DOC source rate. There is no guarantee the two implementations remain synchronized through future refactors.
- **Consequence:** Any future change to the rate formula (e.g., adding a new source or correcting a coefficient) must be applied in two places. If the change is made to `_change_with_components` but not to `rate()`, unit tests that call `rate()` directly will silently test the stale formula. The side-effect write to `pom_doc_source_rate` in `rate()` is particularly hazardous because it will silently shadow the canonical cache with a value computed on different inputs than those seen by `_change_with_components`.
- **Recommendation:** Add a `DeprecationWarning` to `rate()` at its first line (e.g., `warnings.warn("POM.rate() is deprecated; use POM._change_with_components() or POM.run() instead", DeprecationWarning, stacklevel=2)`) and, in a follow-on cleanup, remove the method once the one known external test that uses it is migrated to the `run` path. At minimum, remove the `self.pom_doc_source_rate` side-effect write from `rate()` so the canonical cache can only be set by the `run` path.

**F3 -- `use_DOX` flag is not user-configurable through the parameter DEFAULTS dict or constructor argument**

- **File:line:** `src/clearwater_modules_v3/processes/cbod.py:160--166`
- **Category:** Interface / usability
- **Description:** The `use_DOX` flag is hardcoded as `self.use_DOX = True` inside `__init__` and is not a key in `CBOD_DEFAULTS`. A user who passes `parameters={"use_DOX": False}` will trigger the unknown-key warning at line 142--146 and the value will be silently ignored (the `merged` dict is built from `DEFAULTS` keys only, so the user override never reaches `setattr`). The only way to disable Monod DO attenuation is to subclass `CBOD` or monkey-patch the instance after construction. In contrast, v1 exposes `use_DOX` as a per-column DataArray flag (v1 `processes.py:2384--2386`), giving it full per-cell override capability.
- **Consequence:** Users cannot disable the Monod DO attenuation via the standard configuration pathway. This is a supported use case (pure first-order CBOD decay in data-limited contexts where DOX is not modeled), and the v1 API supports it. The silent parameter-ignore behavior is particularly hazardous because no error is raised.
- **Recommendation:** Add `'use_DOX': True` to `CBOD_DEFAULTS` in `parameters/cbod.py`. The `__init__` assignment `self.use_DOX = True` is then superseded by the DEFAULTS-merge loop (lines 149--150) and user overrides will propagate correctly. No logic change to `_change_with_components` is required.

---

### Observations

**O1 -- CBOD `ksbod_tc / depth` dimensional form deviates from v1/Fortran (documented, dormant under defaults)**

- **File:line:** `src/clearwater_modules_v3/processes/cbod.py:310`; `parameters/cbod.py:9--43`
- **Description:** v3 implements CBOD settling as `ksbod_tc / depth * cbod`, treating `ksbod_20` as a settling velocity (m/d). v1 (`processes.py:2392--2404`) and Fortran (`modCBOD.f90:114`) implement it as `ksbod_tc * cbod`, treating `ksbod_20` as a first-order rate constant (1/d). The parameter module docstring (lines 30--42 of `parameters/cbod.py`) documents this deviation and the related `ksbod_theta = 1.047` vs Fortran `1.024` mismatch. Under the default `ksbod_20 = 0`, both forms are identically zero. The deviation is also recorded in `parameter_defaults_corrections.md` Section 3.5 and pinned in `tests/test_5_cbod_calculations_v2.py:173--200`.

  This observation is not a new finding. It is correctly documented and flagged for user awareness when `ksbod_20 > 0`. Recorded here for completeness.

**O2 -- POM `pom_doc_source_rate` is set as a side effect in `_change_with_components` but is not in `REGISTRY_DIAGNOSTICS`**

- **File:line:** `src/clearwater_modules_v3/processes/pom.py:141--146`, `pom.py:350--351`
- **Description:** The class-level comment at lines 141--146 explains the design decision: `pom_doc_source_rate` is set on `self` (for Carbon to consume via `getattr`) but is intentionally excluded from `REGISTRY_DIAGNOSTICS` because the raw dissolution rate is exposed instead under `pom_hydrolysis_rate`. This is an architectural choice, not a defect. However, if Carbon reads `self.pom_doc_source_rate` directly via `getattr`, the value is unit-converted (mg-C/L_water/d) whereas `pom_hydrolysis_rate` in the registry is in mg-D/L_sed/d, and the two are not interchangeable. Any future Carbon implementation must use the `getattr` path, not the registry path, to get the correctly unit-converted value.

  This is correctly handled in the current code. Noting it as an observation so the architecture is preserved intentionally through refactors.

**O3 -- POM `_change_with_components`: `rate_poc_settling` uses `xr.zeros_like(pom)` as the disabled-branch fallback, while the algal-settling and benthic-mortality disabled branches also use `xr.zeros_like(pom)`, but the enabled branches may return scalar `0` (from `getattr` fallback)**

- **File:line:** `src/clearwater_modules_v3/processes/pom.py:370--376`, `pom.py:382--388`
- **Description:** When `use_floating_algae` or `use_benthic_algae` is `True` but the respective process object is present and the `getattr` returns the integer fallback `0` (line 374: `getattr(..., "algal_pom_from_settling_rate", 0)`), the arithmetic `rate = ... + 0 + ...` returns a DataArray (scalar 0 broadcasts correctly in NumPy/xarray). However, if the `False` branch is taken (line 376: `rate_algal_settling = xr.zeros_like(pom)`) versus the enabled-but-uncached branch returns scalar `0`, the dtype of `rate_algal_settling` will differ across code paths (DataArray vs Python int). This is not a correctness issue under normal operations because xarray broadcasting handles it, but it is a potential source of confusion in unit tests that check the type or shape of `components`.

  Not a defect under normal v3 Model-orchestrated operation. Needs verification only if unit tests assert on container type.

---

## 3. Algorithm Parity Matrix

| v3 Term (process, line) | v1 Reference | Verdict |
|---|---|---|
| `kbod_tc = arrhenius_correction(water_temperature, kbod_20, kbod_theta)` (`cbod.py:292--294`) | `processes.py:2334--2348` | MATCH |
| `ksbod_tc = arrhenius_correction(water_temperature, ksbod_20, ksbod_theta)` (`cbod.py:295--297`) | `processes.py:2352--2366` | MATCH (theta value deviation noted in O1) |
| Oxidation: `kbod_tc * dox / (KsOxbod + dox) * cbod` when `use_DOX` (`cbod.py:302--304`) | `processes.py:2386`: `(DOX / (KsOxbod + DOX)) * kbod_tc * CBOD` | MATCH |
| Oxidation (no DOX): `kbod_tc * cbod` (`cbod.py:306`) | `processes.py:2386`: `kbod_tc * CBOD` | MATCH |
| Settling: `ksbod_tc / depth * cbod` (`cbod.py:310`) | `processes.py:2403`: `CBOD * ksbod_tc` (no depth divide) | DOCUMENTED DEVIATION -- see O1, corrections doc Section 3.5 |
| `dCBOD/dt = -oxidation_rate - settling_rate` (`cbod.py:320`) | `processes.py:2418`: `-CBOD_oxidation - CBOD_sedimentation` | MATCH |
| Sediment CBOD (SOD pathway) | v1 `SOD_tc * (DOX / (KsOxbod2 + DOX)) / depth` | NOT IN SCOPE -- handled by DOX process; correctly deferred |
| `kpom_tc = arrhenius_correction(water_temperature, kpom_20, kpom_theta)` (`pom.py:334--336`) | `processes.py:2185--2197` | MATCH (water temp, not sediment temp -- documented minor deviation per audit doc section 3) |
| Dissolution: `rate_dissolution = kpom_tc * pom` (`pom.py:340`) | `processes.py:2233`: `POM * kpom_tc` | MATCH |
| Burial: `rate_burial = self.vb * pom / self.h2` (`pom.py:355`) | `processes.py:2293`: `vb * POM / h2` | MATCH (vb unit convention: both v1/v3 in m/d; Fortran uses m/yr -- documented) |
| POC settling source: `vsoc * poc / self.h2 / self.fcom` (`pom.py:360`) | `processes.py:2252`: `vsoc * POC / h2 / fcom` | MATCH |
| Algal settling source: reads `floating_algae_process.algal_pom_from_settling_rate` (`pom.py:370--374`) | `processes.py:2216`: `vsap * Ap * rda / h2` (direct computation) | MATCH -- v3 routes via FloatingAlgae cache; algebraically equivalent |
| Benthic mortality source: reads `benthic_algae_process.balgae_pom_from_mortality_rate` (`pom.py:382--386`) | `processes.py:2275`: `Ab * kdb_tc * Fb * (1 - Fw) / h2` | MATCH -- v3 routes via BenthicAlgae cache; algebraically equivalent |
| `dPOM/dt` sign convention (`pom.py:390--395`) | `processes.py:2313`: `POM_algal_settling - POM_dissolution + POM_POC_settling + POM_benthic_algae_mortality - POM_burial` | MATCH |
| Forward Euler integration (`cbod.py:249--250`, `pom.py:286--287`) | `processes.py:2329`, `2434`: `state + dXdt * dt` | MATCH |
| `clip_negative_state` after integration (`cbod.py:253`, `pom.py:292`) | v1 uses `xr.where(state < 0, 0, state)` implicitly via the functional pipe | IMPROVED -- v3 adds diagnostic logging for clip events |

---

## 4. Parameter Default Verification

All CBOD and POM parameter defaults were verified against the three-way Fortran / v1 / v3 table in `design/clearwater_modules_v3_nsm1_audit_simple_constituents.md` (Section "Parameter defaults audit") and the corrections record in `parameter_defaults_corrections.md`.

| Parameter | v1 default | Fortran default | v3 default | Status |
|---|---|---|---|---|
| `KsOxbod` | 0.5 | 0.5 | 0.5 | MATCH |
| `kbod_20` | 0.12 | 0.12 | 0.12 | MATCH |
| `ksbod_20` | 0.0 | 0.0 | 0.0 | MATCH -- intentional dissolved-CBOD convention; corrections doc Section 2.3 |
| `kbod_theta` | 1.047 | 1.047 | 1.047 | MATCH |
| `ksbod_theta` | 1.047 | 1.024 | 1.047 | v3 follows v1; Fortran deviation dormant at `ksbod_20 = 0`; documented |
| `kpom_20` | 0.1 | 0.01 | 0.1 | v1/v3 match; Fortran 10x smaller; documented |
| `kpom_theta` | 1.047 | 1.047 | 1.047 | MATCH |
| `h2` | 0.1 | 0.1 | 0.1 | MATCH -- FIXME cleared Phase 9.F.C; Di Toro H_2 layer |
| `vsoc` (global, via `_POM_GLOBAL_DEFAULTS`) | 0.01 | 0.01 | 0.01 | MATCH |
| `fcom` (global, via `_POM_GLOBAL_DEFAULTS`) | 0.4 | 0.4 | 0.4 | MATCH |
| `vb` (global, via `_POM_GLOBAL_DEFAULTS`) | 0.01 m/d (v1 unit-convention bug) | 6.85e-6 m/d (0.0025 m/yr) | 6.85e-6 m/d | Phase 9.F.A correction applied correctly in `_POM_GLOBAL_DEFAULTS`; corrections doc Section 1.14 |

The `vsop` parameter referenced in the review mandate is a phosphorus parameter (organic-P settling velocity), not a CBOD or POM parameter; it is out of scope here.

---

## 5. Stale Comment / Marker Audit

The two stale markers required by the scope mandate were inspected.

| File | Line | Marker text | Status |
|---|---|---|---|
| `parameters/cbod.py` | 9 | `"* ``ksbod_20`` FIXME cleared (Phase 9.F.C)."` | CORRECTLY RESOLVED. The body of the comment documents the Phase 9.F.C research conclusion (zero default is intentional; literature citations present). The "FIXME cleared" phrasing correctly communicates that a prior `FIXME(phase1-audit)` was investigated and closed, with the rationale retained for traceability. No open action. |
| `parameters/pom.py` | 9 | `"* ``h2`` FIXME cleared (Phase 9.F.C)."` | CORRECTLY RESOLVED. The body documents `h2 = 0.1` as the Di Toro / QUAL2K H_2 anaerobic sediment layer thickness, replacing the prior Phase 0 "unclear physical role" note. The "FIXME cleared" phrasing is analogous to the CBOD entry. No open action. |

No stale unresolved `FIXME`, `TODO`, `BUG`, or `HACK` markers were found in any of the four files in scope.

---

## 6. Correctly Deferred Items

The following items were examined and confirmed as correctly out of scope for v3 NSM1 1.0.0.

1. **Sediment oxygen demand (SOD) within CBOD budget.** The SOD pathway (`SOD_tc * DOX / (KsOxbod2 + DOX) / depth`) is a sink term in the DOX budget, not in the CBOD budget. v3 assigns it to the DOX Process (Phase 5). The CBOD class correctly does not carry this term. Not a CBOD defect.

2. **CBOD multi-group support.** Fortran NSM1 loops over `nCBOD` groups. v3 supports a single group (`"cbod"`) with the multi-group extension path documented in the class docstring. Correctly deferred to Phase 4+ per the design specification.

3. **Full two-layer Di Toro (2001) sediment diagenesis for POM.** v3 NSM1 carries only the H_2 (lower anaerobic) layer with first-order burial/dissolution, matching v1 exactly. The aerobic H_1 layer, nutrient flux model, and SOD coupling are NSM2 sediment-diagenesis scope. Correctly deferred.

4. **Sediment temperature for POM dissolution.** Fortran applies Arrhenius to sediment temperature (`TsedC`); v1 and v3 both use water temperature. v3 correctly follows v1. The Fortran deviation is documented in the audit doc (Section 3) and is a deliberate v1/v3 architectural choice.

5. **Oxygen-weighted TDG in N2.** Out of scope for CBOD/POM review. Addressed in the N2 review.

---

## 7. xarray Refactor Checklist

| Check | CBOD result | POM result |
|---|---|---|
| Python-level cell loops in hot path | None found | None found |
| Scalar-only logic (non-broadcasting constants) | None found | None found |
| Array-truthiness `if` on DataArray (`if cbod:`, `if pom:`) | None found | None found |
| `== np.nan` (must be `.isnull()`) | None found | None found |
| Unguarded division by `depth` in thin/dry cells | `ksbod_tc / depth` at line 310 -- guarded downstream by `sanitize_rate` (line 318) which zeros inf | `self.h2 / depth` at line 351 in `pom_doc_source_rate` -- `h2 = 0.1` constant scalar, so division is `depth`-only. Unguarded by `sanitize_rate`. No `sanitize_rate` call wraps POM rates generally (only CBOD uses it). If `depth -> 0` the `pom_doc_source_rate` cache becomes inf. No wet-mask note for POM rate helpers. Practical risk is low if the orchestration-layer wet-mask gates `run` at thin cells, but the secondary defense is absent for POM. |
| `xr.where` used correctly for conditional rate gating | `use_DOX` gated by Python `if self.use_DOX` on a `bool` -- safe | `use_POC`, `use_Algae`, `use_Balgae` gated by Python `if self.use_POC`, etc. -- all `bool` scalars, safe. v1 used `xr.where(use_DOX == True, ...)` on per-column arrays; v3's scalar bool gating is equivalent for the single-group model. |

The dry-cell `pom_doc_source_rate` issue noted in the table above is a borderline Minor concern. It is not promoted to a separate finding because (a) the wet-mask orchestration layer provides the primary defense, (b) `pom_doc_source_rate` is a derived cache consumed by Carbon rather than a state variable written back to the registry, and (c) the `sanitize_rate` secondary defense is absent by design in POM (no import of `sanitize_rate` in `pom.py`). If a future refactor removes wet-mask gating or drives POM in standalone mode with thin-cell inputs, this path could silently propagate inf to Carbon. The Carbon Process should apply its own `sanitize_rate` when reading the cache.

---

## 8. Positive Notes

The following practices are worth preserving through future refactors.

The `sanitize_rate` double-defense pattern in CBOD (`cbod.py:317--318`) -- applying NaN/inf zeroing to each sub-flux independently before summing -- prevents a single bad cell in one sub-rate (for example, from a transient DOX spike) from poisoning the entire cell's net rate. This is more robust than sanitizing only the final sum.

The `getattr(..., "algal_pom_from_settling_rate", 0)` fallback pattern in POM (`pom.py:372--374`) provides graceful degradation when FloatingAlgae has not yet run or is absent. The fallback to integer `0` broadcasts correctly in xarray arithmetic, and the three-condition guard (`use_floating_algae and use_Algae and floating_algae_process is not None`) prevents any attribute access on a None process reference.

The module-level docstring of `pom.py` (lines 1--75) is exemplary: it documents the Di Toro H_2 conceptual identity of the POM state variable, the `h2` dimensional role, the `pom_doc_source_rate` unit-conversion reasoning with an explicit numerical example (`25x overcount` at depth=1 m), and the boundary between NSM1 and NSM2 scope. This level of inline physical reasoning is unusual and should be maintained.

The `vb` correction in `_POM_GLOBAL_DEFAULTS` (`pom.py:106`: `6.85e-6 m/d`) correctly applied the Phase 9.F.A fix. The inline comment records both the unit-equivalent value (0.0025 m/yr) and the match to v3 `global_vars` after the correction, providing traceability for calibration users.

The CBOD module docstring (lines 1--39) correctly documents the Phase 0 audit finding that `ksbod_20 = 0` is an intentional zero (not a placeholder), the dimensional form deviation (`ksbod_tc / depth`), and the Phase 5 DOX rate-consumer contract (`cbod_oxidation_rate` cached as `self.cbod_oxidation_rate`). This prevents the three most likely misinterpretations a reviewer or integrator would encounter.

---

## 9. Recommended Follow-up Tests or Benchmarks

1. A unit test that reads `pom_settling_rate` from the registry after a `POM.run` call and asserts that its value equals `vb * pom / h2` (not `vsoc * poc / h2 / fcom`). This test would fail today due to F1 and would serve as a regression guard after the rename fix.

2. A conservation test for the POM-DOC pathway: run a closed system with `use_POC = True`, `use_Algae = False`, `use_Balgae = False`, and nonzero initial POC. Assert that the carbon transferred from POC to POM to DOC over a simulation equals `fcom * integral(kpom_tc * pom * h2 / depth) dt` at each step, verifying the `h2 / depth` unit-conversion factor.

3. A dry-cell robustness test for `pom_doc_source_rate`: inject `depth = 1e-6 m` (near-dry) and verify that `pom_doc_source_rate` does not propagate inf to any Carbon consumer. This test would identify whether the absent `sanitize_rate` guard in POM causes downstream failure.

4. A parameter-override test for `use_DOX = False` in CBOD: construct `CBOD(parameters={"use_DOX": False})` and verify that the resulting `use_DOX` attribute is `False` (not the hardcoded `True`). This test would fail today due to F3 and serve as a regression guard after the DEFAULTS fix.

---

## 10. Open Questions

1. Is `pom_settling_rate` consumed by any existing downstream code (Carbon Process, DOX Process, any test or post-processing script) under that string name? If so, an audit of those consumers is required before the F1 rename is applied, to avoid breaking the consumer silently. A grep for `"pom_settling_rate"` across the full repository tree (excluding the four files in this review scope) will answer this.

2. The `rate()` legacy method in `pom.py` is retained for back-compat with "external tests that call `pom.rate(...)` directly." Which specific test files call `rate()` directly? Identifying them is a prerequisite for F2 remediation.

3. The `parameters` dict merge loop for CBOD (`cbod.py:149--150`) iterates over `merged` (which is `{**DEFAULTS, **user_params}`). If a user passes a key in `DEFAULTS` with a value override, it propagates correctly. But `use_DOX` is set after the merge loop (line 166: `self.use_DOX = True`) and thus cannot be overridden by the merge. Is the hardcoded `True` an intentional post-merge override, or was it simply never added to DEFAULTS? The answer determines whether the F3 fix is a one-line DEFAULTS addition or requires a deeper architecture decision.
