# v3 NSM1 Framework, Utilities, Configuration, and Temperature -- Source-Code and Science-Correctness Review

**Review date:** 2026-05-15
**Reviewer:** water-quality-model-source-code-reviewer agent
**Branch:** streaming
**Commit:** 54f2b12
**Scope reviewed line-by-line:** `model.py`, `processes/base.py`, `config/init.py`, `config/read.py`, `config/__init__.py`, `__init__.py`, `processes/riverine.py`, `utils/numerics.py`, `utils/conversions.py`, `utils/constants.py`, `utils/light.py`, `utils/sediment.py`, `utils/partitioning.py`, `utils/reaeration.py`, `utils/__init__.py`, `parameters/global_parameters.py`, `parameters/global_vars.py`, `parameters/__init__.py`, `processes/__init__.py`, `examples/nsm1_demo_setup.py`. `processes/temperature.py` spot-checked (per scope). v1 `nsm1/model.py`, v1 `nsm1/constants.py` read as parity reference.
**Languages:** Python (xarray-native `Process` framework)
**Method:** Line-by-line read of all scoped v3 files, cross-read of v1 reference, cross-check against `design/clearwater_modules_v3_review_findings.md` (2026-05-04 multi-agent review), `design/clearwater_modules_v3_nsm1_audit_utilities_params.md`, and `design/clearwater_modules_v3_tsm_audit_2026-05-05.md`. Standalone import and a numerics/conversions smoke test were executed. Benchmarks were not re-run (out of scope).

---

## 1. Summary verdict

The v3 shared scaffolding is in good shape. The 2026-05-04 review's ten CRITICAL findings (C1--C10) are all verifiably resolved in the current code, with in-line comments citing the finding ID at each fix site, and the import-blocker C8 is confirmed fixed (`PYTHONPATH=src python -c "import clearwater_modules_v3"` succeeds). The MAJOR findings are substantially resolved; the few that remain open are either correctly deferred (M4) or are low-consequence robustness items. `utils/numerics.py` (the new Q7 clip-with-log integrator safety net) and `utils/conversions.py` are correct and container-type-aware; the xarray refactor in the utilities is complete and free of the scalar-only / array-truthiness / `== np.nan` hazards the review was asked to look for.

The most significant findings are documentation-fidelity issues, not code defects. The package `README.md` and the top-level `__init__.py` docstring contain materially stale claims (Riverine/BenthicAlgae/FloatingAlgae/Nitrogen described as v2 re-exports when they are v3-native; "2 of 18 MAJOR resolved" when the authoritative triage records 17 of 18; v2 described as both "retired/removed" and "re-exported" in different files). Three authoritative design/audit documents also carry findings that the current code has already resolved (M4 Kelvin offset, `lambdam`, `vson`), which creates a real risk that a downstream reader trusts a "STILL DEFERRED" or "needs action" label that no longer reflects the code. No CRITICAL or MAJOR code defect was found in the scoped framework/utility layer. One genuinely-open framework item (Riverine cannot self-register variables; depth is a placeholder alias of wetted surface area) is correctly TODO-flagged and is a legitimate open item, not a stale comment.

**Finding counts:** CRITICAL 0, MAJOR 0, MINOR 5, OBSERVATION 9.

**Overall confidence:** High for the framework and utilities layer. The kernel-schedule, wet-mask, hotstart, and chunking semantics were traced end-to-end and match the documented v1-derived intent. Temperature was spot-checked only (per scope) and the C3/C5/C10/hotstart fixes are present; a full TSM energy-balance re-derivation was not performed and was not in scope.

---

## 2. Findings table

| ID | Severity | file:line | Category | Description | Recommended fix |
|----|----------|-----------|----------|-------------|-----------------|
| F1 | MINOR | `src/clearwater_modules_v3/README.md:17,112` | documentation | README states "remaining process classes (`Riverine`, `BenthicAlgae`, `FloatingAlgae`, `Nitrogen`) and the `Process`/`ProcessFactory` base remain re-exports from v2" and "Re-exported processes ... are bit-for-bit the v2 classes." This is false: all process modules and `processes/base.py` are v3-native in-tree (the top-level `__init__.py:11-15` says the v2 package was *removed* in the v2-retirement work). Internally inconsistent and misleading for any reader assessing what is v3-native. | Rewrite the README Status and Backward-compatibility sections to state that all processes and the `Process`/`ProcessFactory` base are v3-native in-tree; remove the "bit-for-bit v2" claim. |
| F2 | MINOR | `src/clearwater_modules_v3/README.md:12,19,30,33` | documentation | README "Status" says "Phases 0--4 complete; Phase 5 in progress" and "Phase R-3 ... In progress (2 of 18 resolved)" and lists only M6/M9 resolved. The authoritative `clearwater_modules_v3_review_findings.md:32` records 17 of 18 MAJOR resolved (only M4 deferred), and the code confirms many R-3 fixes are present (M5/M7/M8/M10/M11/M14 contracts in `model.py`; M12/M13/M15 in `config/init.py`; M17 partially done). README predates Phases 6--10 NSM1 pattern alignment entirely. | Update Status, Phase status table, and the R-3 line to match `clearwater_modules_v3_review_findings.md` (17/18 MAJOR resolved, M4 deferred); add the NSM1 pattern-alignment Phases 6--10 to the phase table. |
| F3 | MINOR | `src/clearwater_modules_v3/model.py:24-32` | documentation | Module docstring states chunking step-index comparison is "immune to floating-point drift in `current_time += time_step` arithmetic" and the kernel schedule is "exact-integer ... immune to floating-point drift" without qualification. The in-method docstring (`model.py:463-466`) correctly scopes the exactness to "the common case where `process.time_step_seconds` is an integer multiple"; the schedule modulo `delta_seconds % interval == 0` with `delta_seconds = i * time_step_seconds` is *not* exact for non-integer-second `time_step` (verified: `21 * 0.1 % 0.3 = 1.67e-16`, not 0). M18 was closed by documenting this in a test, not by making it exact. The module-level blanket claim overstates the guarantee. | Qualify the module docstring to match the precise in-method docstring: exact for integer-second `time_step`/cadence; non-integer-second `time_step` is validated/tested but relies on the cadence-multiple `ValueError` guard, not bit-exact float modulo. |
| F4 | MINOR | `design/clearwater_modules_v3_review_findings.md:203,354`; `design/clearwater_modules_v3_nsm1_audit_utilities_params.md:339` | stale-comment (design doc, not source) | M4 is labeled "STILL DEFERRED ... v3 inherits with comment 'for testing consistency with v1'", and the utilities audit Section "Documentation defects" item 3 says v3 `celsius_to_kelvin` "re-exports from v2 ... returns `T_C + 273.16`." The current code does neither: `utils/constants.py:29` `KELVIN_OFFSET = 273.15`, `utils/conversions.py:45-47` defines the function in-tree using the SI offset (verified at runtime: `celsius_to_kelvin(0.0) == 273.15`). The TSM audit (`clearwater_modules_v3_tsm_audit_2026-05-05.md:215-219`, commit `a8874e7`) resolved this one day after the multi-agent review. The "STILL DEFERRED" / "re-exports 273.16" labels are stale and will mislead a reader doing pre-LimnoTech triage. | Update `review_findings.md` M4 status to "RESOLVED 2026-05-05 (TSM audit Q5, commit a8874e7); v3 uses 273.15"; correct `audit_utilities_params.md` Section "Documentation defects" item 3. |
| F5 | MINOR | `design/clearwater_modules_v3_nsm1_audit_utilities_params.md:323-325,343-344`; `parameter_defaults_corrections.md` (Sec 1.9 referenced) | stale-comment (design doc, not source) | The utilities/params audit lists `lambdam=0.0174` as a "Likely v1 flaw NOT corrected in v3 (needs action)" and `vson_20=0.1` / `vson_theta=1.024` as "Undocumented v3 deviations." The current code has resolved all three: `parameters/global_vars.py:61` `lambdam: 0.174` (Phase 9.C), `parameters/nitrogen.py:55` `vson_20: 0.01` (Phase 9.C), and `vson_theta` removed (Phase 9.E). v1 still has `lambdam=.0174`. The audit's "Required actions before LimnoTech review" items 1--2 are satisfied in code but the doc still reads as open. | Add a resolution note to the audit doc (or a dated addendum) recording that `lambdam`, `vson_20`, and `vson_theta` were corrected in Phases 9.C/9.E; mark "Required actions" items 1--2 closed. |

---

## 3. README staleness assessment

**Verdict: the README is materially stale and should be regenerated before any LimnoTech or sponsor handoff.** It predates the NSM1 pattern-alignment Phases 6--10 and the v2-retirement refactor, and it is internally inconsistent with the package's own top-level `__init__.py` docstring.

Specific stale or contradictory claims:

1. `README.md:17` -- "The remaining process classes (`Riverine`, `BenthicAlgae`, `FloatingAlgae`, `Nitrogen`) and the `Process`/`ProcessFactory` base remain re-exports from v2 pending their own merge phases." False. `processes/base.py:1-7` is v3-native in-tree ("Class definitions in-place. Originally inherited verbatim from v2 by re-export; the v3-self-sufficient refactor moved the class bodies in-tree"). `processes/__init__.py` imports every process from in-tree v3 modules. `__init__.py:11-15` states the `clearwater_modules_v2` package was removed entirely.
2. `README.md:112` -- "Re-exported processes (`Riverine`, `BenthicAlgae`, `FloatingAlgae`, `Nitrogen`) are bit-for-bit the v2 classes, registered under the same names with the same `ProcessFactory`." False and self-contradictory with `__init__.py`'s v2-removal statement. `processes/base.py` defines a v3-owned `ProcessFactory`.
3. `README.md:12,30` -- "Phase 5 (READMEs and migration notes) is in progress." Stale relative to the much later NSM1 pattern-alignment Phases 6--10 (closeout docs `..._phase6_closeout.md` through `..._phase10b_closeout.md` exist in `design/`).
4. `README.md:19,33` -- "R-3 (MAJOR-finding cleanup) is in progress; 2 of 18 MAJOR findings (M6, M9) are resolved." Contradicts the authoritative `clearwater_modules_v3_review_findings.md:32` ("17 of 18 resolved -- only M4 deferred") and the code, which contains the M5/M7/M8/M10/M11/M14/M12/M13/M15 fixes.
5. `README.md:131` -- references the review-findings doc as "18 MAJOR" with no resolution status; acceptable as a pointer but the surrounding prose understates progress.

The "What's new in v3" section (`README.md:35-48`) is, by contrast, accurate to the code: the latent-heat fix, thin-water ramp, vectorized `mixing_ratio_air`, sediment-diffusivity 0.0432 m^2/day, dynamic sediment-temperature evolution, precomputed schedule, output-variable-scoped wet-mask, hotstart hooks, and integer-step-index chunking are all present in `temperature.py` and `model.py` as described. The YAML schema block (`README.md:90-108`) matches `config/init.py`. The staleness is confined to the Status, Phase-status, and Backward-compatibility sections.

---

## 4. 2026-05-04 review-findings cross-check (C1--C10 / M1--M18)

CRITICAL findings -- all ten verifiably resolved in the scoped code:

- **C1** (`simulation_directory` str/Path): RESOLVED. `model.py:161-164` wraps in `Path`; the explicit predicate also resolves m12.
- **C2** (`finalize_process` AttributeError): RESOLVED. `model.py:391-394` uses `getattr(process, "finalize_process", None)` + callable guard; `processes/base.py:116-117` also now defines a no-op default, so the guard is belt-and-suspenders.
- **C3** (sediment diffusivity unit/value): RESOLVED. `temperature.py:100` `sediment_diffusivity: float = 0.0432`; `temperature.py:206-217` docstring says "m^2/day"; `temperature.py:702-743` flux divides by 86400 with a comment explaining the m^2/day to m^2/s conversion.
- **C4** (`mixing_ratio_air` denom guard): RESOLVED (spot-check). `temperature.py:1445+` present with vectorized guard (consistent with README "What's new" item).
- **C5** (wet-mask masks inputs): RESOLVED. `model.py:521-562` masks only `getattr(process, "output_variables", ...)` with fallback to `variables`; `temperature.py:84` declares `output_variables`.
- **C6** (schedule timezone dependence): RESOLVED. `model.py:443-495` keys firing off `delta_seconds = i * time_step_seconds`; no `start_time.timestamp()`; cadence-multiple `ValueError` at `model.py:479-486`.
- **C7** (chunk-end membership FP-fragile): RESOLVED. `model.py:698-721` precomputes `interior_chunk_step_indices: set[int]`; `model.py:706-713` raises `ValueError` if `chunk_size` is not an integer multiple of `time_step`.
- **C8** (cannot import standalone): RESOLVED. Confirmed at runtime: `PYTHONPATH=src python -c "import clearwater_modules_v3"` succeeds without a conftest shim.
- **C9** (`_v2_init_helper` misunderstanding): RESOLVED in spirit. The v3 `config/init.py` no longer dispatches into a `clearwater_modules_v2.config.init` helper at all; `_init_processes`/`_init_model_data` are v3-native in-tree (`config/init.py:315-434`). The C9 fix evolved further (full in-tree port) than the review-findings text describes; not a defect, but the doc lags the code.
- **C10** (dynamic sediment-temperature evolution dropped): RESOLVED (spot-check). `temperature.py:103` `evolve_sediment_temperature: bool = True`; `temperature.py:498` gates the update on `use_sediment_temperature and evolve_sediment_temperature`; `temperature.py:1291-1324` `sediment_temperature_change`.

MAJOR findings within scope:

- **M5** (init_process before from_hotstart ordering): RESOLVED as a documented contract. `model.py:275-307` and `processes/base.py:31-67` document the four-step ordering; `temperature.py:395-408` `from_hotstart` overrides the fresh-start `__skip_first_time_step`. The contract is documented and exercised by Temperature, but still relies on process-author discipline (no runtime enforcement) -- acceptable per the review's "documented contract" resolution.
- **M6** (finalize only in chunked mode): RESOLVED. `model.py:663` calls `__finalize_model()` in `__process_loop_full`.
- **M7** (chunked source, unchunked run reads one step): RESOLVED. `model.py:315-318` sets `chunk_end = self.__end_time` when not chunked.
- **M8** (seed-from-hotstart first-dim fallback): RESOLVED. `model.py:586-617` recognizes only `time`/`time_step`/`datetime` dims and raises `ValueError` when a slice is requested with no time axis.
- **M9** (bare `except Exception` in wet-mask): RESOLVED. `model.py:554-558` narrowed to `except KeyError`.
- **M10** (run() callable twice): RESOLVED. `model.py:190,217-224` `__run_complete` flag raises `RuntimeError` on a second `run()`.
- **M11** (wet-mask threshold strict-> with default 0.0): RESOLVED as documented behavior. `model.py:98-122` documents the numerical hazard and the caller-side epsilon recommendation; the `>` semantic is intentional and documented, not enforced.
- **M12** (hotstart unsupported suffix): RESOLVED. `config/init.py:211-268` enumerates supported suffixes and raises a `ValueError` naming `hotstart.dataset_path` before filesystem access.
- **M13** (timestep parseability): RESOLVED. `config/init.py:277-292` eagerly `pd.to_datetime`-validates non-integer `timestep`, rejecting `bool`.
- **M14** (wet_mask variable existence): RESOLVED. `model.py:354-379` validates `wet_mask_variable in registry` at init with a helpful enumeration.
- **M15** (KeyError-to-ValueError deep path): RESOLVED. `config/init.py:189-208` `_required(d, *path)` names the offending YAML path.
- **M16** (integrator-pattern contract not enforced): RESOLVED by demotion to guideline. `processes/base.py:9-29` is now a "guideline, not enforced contract."
- **M17** (README outdated): PARTIALLY resolved. The "What's new" and YAML-schema sections were added and are accurate, but the Status/Phase-status/Backward-compatibility sections remain stale (see Section 3 and F1/F2). Recorded as still-open documentation work.
- **M18** (bare `==` mod-op on floats): RESOLVED by test + cadence-multiple validation; the residual is the module-docstring overclaim flagged as F3.
- **M1, M2, M3** (Temperature stability/NaN guards): RESOLVED (spot-check). `temperature.py:1584-1650` carries the M3 NaN-propagation and `np.errstate` RuntimeWarning suppression; M1/M2 validation present per README "What's new" and the C-series resolution notes.
- **M4** (273.16 vs 273.15): RESOLVED IN CODE (not deferred). See F4. `utils/constants.py:29` and `utils/conversions.py:45-47` use 273.15; runtime-verified. The review-findings "STILL DEFERRED" label is stale.

No CRITICAL or MAJOR finding was found to have a contradicting stale "broken/TODO" comment in the scoped source. The stale labels are all in the design/audit documents (F4, F5), not in code.

---

## 5. Framework-semantics correctness note

**Kernel compute-schedule (`model.py:443-501`).** Correct. The firing predicate is "process fires at step index `i` when `(i * time_step_seconds) % process.time_step_seconds == 0`," seeded from delta-seconds-from-`start_time`, which is timezone-independent for naive datetimes (the C6 fix). A cadence-multiple `ValueError` (`model.py:479-486`) rejects any `process.time_step_seconds` that is not an integer multiple of the model `time_step`, which is the precondition that makes the modulo exact for integer-second steps. The schedule is precomputed to `n_steps + 1` entries while the loop runs `n_steps` (`model.py:489` vs `model.py:497-501`); this one-extra-entry margin is harmless (it is never indexed past `n_steps - 1` because the loop condition is `current_time < end_time`) and is the previously-documented m14 defensive margin. The only soft spot is the module-docstring overclaim about float-time immunity for non-integer-second steps (F3); the guarantee is real for the supported integer-second case.

**Registry-level wet-mask gating (`model.py:505-562`).** Correct and matches the C5 intent. `__compute_wet_mask` returns `None` when unconfigured (legacy v2 behavior preserved), supports a `wet_mask_provider` callable override, and otherwise computes `value > threshold` (strict inequality, documented hazard for `value`-near-zero at `model.py:98-122`). `__apply_wet_mask` masks **only** the process's declared `output_variables` (with a `getattr` fallback to `variables` for unmigrated processes), so dry-cell forcings are preserved across substeps; the dtype guard (`model.py:559`) correctly skips non-floating variables; the `except KeyError` (M9) only swallows the documented "variable not yet registered" case. The mask is recomputed once per substep only when at least one process fires (`model.py:653,750`), which is the O1 profiler observation, not a defect.

**Hotstart (`model.py:255-271,326-347,566-638`; `processes/base.py:31-67`).** Correct and well-documented. Ordering invariant: load sources -> seed registry from hotstart slice -> `init_process` (fresh-start defaults) -> `from_hotstart` (override). `__seed_from_hotstart` recognizes only `time`/`time_step`/`datetime` time dims and fails loudly otherwise (M8). `to_hotstart`/`from_hotstart` are opt-in via `getattr`; `__restore_process_hotstart` shares the dataset `attrs` across processes by key-prefix convention (O2 -- a documented convention, not a defect). The remaining residual is the documented M5 hazard: a process author who adds new internal substep state must remember to handle it in both `init_process` and `from_hotstart`; this is a contract, not enforced, and is acceptable per the review resolution.

**Chunking (`model.py:665-762`).** Correct. `chunk_size` must be an integer multiple of `time_step` (`ValueError` at `model.py:706-713`). `interior_chunk_step_indices` is a `set[int]` precomputed from `chunk_size_seconds / time_step_seconds`; boundary detection is exact-integer `step_index in interior_chunk_step_indices`. The trailing partial chunk is written exactly once after the loop (`model.py:758-761`), resolving the v2 double-write. The minus-one-time-step on the next-chunk load start (`model.py:739-742`) preserves previous-step lookups across the boundary, which is the correct ghost-step handling for time-indexed reads. `__finalize_model` is called symmetrically in both loop variants (M6).

**xarray-refactor completeness (utils + model loop).** Complete in the scoped files. `utils/numerics.py` is fully container-type-aware (`xr.DataArray` / `np.ndarray` / scalar paths) and uses `rate.isnull() | np.isinf(rate)` for the DataArray branch (not `== np.nan`). `utils/conversions.py`, `utils/light.py`, `utils/sediment.py`, `utils/partitioning.py` are pure broadcasting expressions with no scalar-only logic, no cell loops, and no array-truthiness `if`. `utils/reaeration.py` correctly handles the `np.select` dim-stripping hazard via `_first_dataarray` reattachment (audit-confirmed). The model loop uses integer step indices throughout and avoids float-time arithmetic in control flow. The one unguarded-division note: `utils/reaeration.py:246` `kaw_tc / depth` and several `kah_20` formulas divide by `depth`; the v3 design routes thin/dry-cell protection through the orchestration wet-mask plus `sanitize_rate` (`utils/numerics.py:210-252`), which is the documented defense-in-depth posture, so this is by design rather than a defect at this layer.

---

## 6. Stale-comment list

No stale "broken/TODO" comments were found in the scoped **source** files. The active TODOs in `processes/riverine.py:36,78,101` are genuine open framework items, not stale (see Section 7). The `temperature.py:1580` reference to a "commented `-1` factor with a TODO" is explicitly described as **resolved** ("a TODO that was resolved per Jason Rutyna's January 2026 diff investigation") and is correctly framed as historical provenance, not an open item. `config/__init__.py:1-6` docstring says "Phase 1: re-exports ... Phase 3 will add" -- this is mildly stale phrasing (the work is done, not "will") but the module behavior is correct; folded into the documentation findings rather than listed separately.

Stale items, all in **design/audit documents** (not source), are captured as F4 and F5:

1. `clearwater_modules_v3_review_findings.md:203,354` -- M4 "STILL DEFERRED" / "v3 inherits with comment 'for testing consistency with v1'"; code uses 273.15.
2. `clearwater_modules_v3_nsm1_audit_utilities_params.md:339` -- "v3 `celsius_to_kelvin` re-exports from v2 ... returns `T_C + 273.16`"; code defines it in-tree at 273.15.
3. `clearwater_modules_v3_nsm1_audit_utilities_params.md:323-325,343-344` -- `lambdam=0.0174` "NOT corrected in v3 (needs action)"; code has `0.174`.
4. `clearwater_modules_v3_nsm1_audit_utilities_params.md:329-330` -- `vson_20=0.1` / `vson_theta=1.024` "Undocumented v3 deviations"; code has `vson_20=0.01` and `vson_theta` removed (Phase 9.E).
5. `README.md:17,112` -- processes/base "remain re-exports from v2"; code is v3-native in-tree (F1).

---

## 7. Correctly-deferred list

1. **Riverine cannot register its own variables to the registry** (`processes/riverine.py:78` "TODO: update once Riverine can register variables to the registry"). Legitimate open framework item: until `ClearwaterRiverine` exposes a registry-registration hook, the v3 `Riverine.init_process` manually mirrors `mesh.Ap/NH4/NO3/TIP/DOX` into the registry. This is an upstream-coupling dependency (clearwater_riverine API), correctly flagged, not stale. Not a defect at the v3 layer.
2. **Riverine `depth` is a placeholder alias of `wetted_surface_area`** (`processes/riverine.py:101-107` "TODO: replace this with depth calculation"). `registry.register("depth", DataArrayVariable(self.riverine_instance.mesh.wetted_surface_area.copy(deep=False)))` registers wetted surface area under the name `depth`. This is dimensionally wrong (m^2 registered as m) and would corrupt any depth-dependent kinetics, **but** it is correctly TODO-flagged and is only reached inside the `if model.has_process("FloatingAlgae")` block when Riverine drives NSM1 kinetics; the design clearly marks it as not-yet-implemented. Classify as a correctly-deferred open item rather than a finding, with the caveat that it must be resolved before any production NSM1+Riverine coupled run that consumes `depth` (worth surfacing to the riverine-coupling owner; the marker is honest about the gap).
3. **Riverine datetime-string conversion** (`processes/riverine.py:36-41` "TODO: This will be removed once Riverine is updated to use datetime objects"). The shim converts `datetime` to a `"%m-%d-%y %H:%M:%S"` string for `ClearwaterRiverine`. Legitimate upstream-API-dependency deferral, correctly flagged. Note: the `%y` (two-digit year) format is a latent Y2.1K-class ambiguity for far-future dates, but for the supported simulation horizon it is benign; mention only as a forward note, not a finding.
4. **M4 (273.16 -> 273.15) "decommission v1-parity tests" follow-up.** The code already uses 273.15; what remains deferred per the review plan is the *v1-parity test re-derivation*, which is correctly out of this framework/utils scope and is a test-suite concern.
5. **MMS energy-conservation test** (`review_findings.md:368`). Correctly deferred to v3.1 NSM1 reactive-transport work; out of scope for this framework/utilities review.
6. **NSM2 sediment diagenesis.** Any sediment-diagenesis dynamic state beyond the C10 sediment-temperature relaxation is correctly deferred to v3.x NSM2 work per project memory; nothing in the scoped files improperly claims it is implemented.

---

## 8. Positive notes (preserve through future refactors)

1. `utils/numerics.py` `clip_negative_state` / `sanitize_rate` are exemplary: container-type-aware, preserve `xr.DataArray` coords/dims/attrs, rate-limit log records with an aggregate suppressed-count stub, and degrade gracefully when `diagnostics is None`. The clip target is exactly 0 (not an epsilon), which is the correct choice for Monod `C/(C+K)` ratios. Keep this contract intact.
2. The in-line finding-ID provenance comments throughout `model.py` and `config/init.py` (each fix cites C#/M# and the review date) make regression triage straightforward and should be preserved on future edits.
3. `utils/reaeration.py` `_first_dataarray` reattachment after `np.select` is a correct and non-obvious xarray-broadcasting fix; the explanatory comment should be kept.
4. `model.py` `__seed_from_hotstart` fail-loud-on-no-time-dim and `config/init.py` `_required(d, *path)` deep-key-path errors are good scientific-software hygiene (errors name the YAML location, not a bare `KeyError`).
5. The `parameters/*.py` correction comments (e.g., `global_vars.py:42` `vb` 1460x unit bug, `:54` `q_solar` mislabeled-units note, `:60-61` `lambdas`/`lambdam`) document the Fortran/QUAL2K provenance of every non-trivial default; this is high-quality documentation-to-code fidelity at the parameter layer and contrasts favorably with the stale package README.

---

## 9. Recommended follow-up

1. Regenerate `README.md` Status / Phase-status / Backward-compatibility sections from `clearwater_modules_v3_review_findings.md` and the Phase 6--10 closeout docs (resolves F1, F2, and the still-open part of M17).
2. Add dated resolution addenda to `clearwater_modules_v3_review_findings.md` (M4) and `clearwater_modules_v3_nsm1_audit_utilities_params.md` (Kelvin offset, `lambdam`, `vson`) so downstream readers do not act on stale "deferred / needs action" labels (resolves F4, F5).
3. Tighten the `model.py` module docstring float-time wording to match the precise in-method docstring (resolves F3).
4. Before any production NSM1+Riverine coupled run, resolve `riverine.py:101` so `depth` is a true depth, not `wetted_surface_area`; add a guard or explicit `NotImplementedError` if a depth-consuming process is configured while the placeholder is active, so the dimensional error cannot silently reach kinetics.
5. Consider a lightweight runtime assertion (or a contract test) for the M5 `init_process`/`from_hotstart` parity, since the only current enforcement is author discipline and the failure mode (silent fresh-vs-resume divergence) is exactly the kind of non-reproducibility this review class is meant to catch.

---

## 10. Open questions

1. `riverine.py:101` -- is the `depth = wetted_surface_area` placeholder ever exercised by a production config today, or is NSM1+Riverine coupling still gated off pending the depth calculation? This determines whether item 4 above is a forward note or a release blocker for that coupling path.
2. Are the design/audit documents (`review_findings.md`, `audit_utilities_params.md`) considered living documents that should be patched, or point-in-time artifacts? If the latter, F4/F5 should instead be resolved by a single consolidated "current status" index that supersedes the stale labels.
