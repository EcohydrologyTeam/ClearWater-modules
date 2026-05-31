# Phase 0 Baseline Artifacts — v3 NSM1 Pattern Alignment

Originally captured 2026-05-13 against commit `186b5c4` ("Add v3 NSM1 pattern alignment specification") on branch `streaming`.

**Active baseline: `e0185de`** (re-baselined 2026-05-30 for the FloatingAlgae computed light-extinction wiring — `limit_light` now uses the optical-constituent lambda from `utils.light.L` instead of the constant `light_attenuation_coefficient`; a broad algae-light cascade through algae/DO, layered on `d530a3a`; see "Re-baseline log" below). The `d530a3a` (Pathogen `algae_floating`), `6c10f36` (terminal gold-standard through NSM1-SCI-A2), `3a8c188` (through SCI-A3), `b51df71` (CA-1+SCI-N1), `624ed7c` (CA-1 only) and `186b5c4` (pre-fix) artifacts are retained in the tree for auditability and are no longer the active reference.

These artifacts are the **gold reference** for the zero-regression contract in `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §11. Every per-Process phase commit (Phase 1 through Phase 10) of the pattern-alignment work must reproduce them bit-identically when no `REGISTRY_DIAGNOSTICS` names are pre-registered.

## Files

| File | Purpose |
|---|---|
| `baseline_coupled_trajectory_e0185de.nc` | **ACTIVE.** 4,320-substep coupled NSM1 demo trajectory at the FloatingAlgae computed light-extinction wiring (commit `e0185de`). Differs broadly from `d530a3a` (the algal light-limitation lambda is now computed from Solid/POC/Ap, cascading through algae→DO→nutrients). **Load-bearing for §11.2.** Captured + parity-verified under the conda `clearwater` test env. |
| `baseline_coupled_trajectory_d530a3a.nc` | Superseded (Pathogen canonical `algae_floating` name fix; encodes the constant light_attenuation_coefficient=1.0 algal optics). 4,320-substep coupled NSM1 demo trajectory. Retained for auditability. |
| `baseline_coupled_trajectory_6c10f36.nc` | Superseded (terminal gold-standard through NSM1-SCI-A2; encodes the Pathogen zero-algae name bug — Pathogen read the unregistered `ap`). 4,320-substep coupled NSM1 demo trajectory. Same shape/contract as `186b5c4`. Retained for auditability. |
| `baseline_coupled_trajectory_3a8c188.nc` | Superseded (through SCI-A3; encodes the SCI-A2 0.5 mortality-C mis-routing). 4,320-substep coupled NSM1 demo trajectory. Retained for auditability. |
| `baseline_coupled_trajectory_b51df71.nc` | Superseded (CA-1+SCI-N1; encodes the SCI-A3 algae total-shortwave defect). 4,320-substep coupled NSM1 demo trajectory. Retained for auditability. |
| `baseline_coupled_trajectory_624ed7c.nc` | Superseded (CA-1 only; encodes the SCI-N1 4× denitrification-alkalinity defect). Retained for auditability. |
| `baseline_coupled_trajectory_186b5c4.nc` | Superseded (pre-CA-1; encodes the raw-weight alkalinity bug). 4,320-substep coupled NSM1 demo trajectory; 20 state/forcing variables × 5 cells × 4,321 substep indices (initial condition + 4,320 substeps). Bit-identical reproducibility verified across two consecutive captures. Retained for auditability. |
| `baseline_junit_full_6c10f36.xml` | **ACTIVE (terminal).** Full repo test suite JUnit XML at the CA-1 + SCI-N1 + SCI-A3 + SCI-A2 fixes (999 passed, 2 xfailed, 0 failed, 314.59 s wall; +2 vs `3a8c188` from the new SCI-A2 regression tests). |
| `baseline_tier1_junit_6c10f36.xml` | **ACTIVE (terminal).** Tier 1 conservation tests JUnit XML at the CA-1 + SCI-N1 + SCI-A3 + SCI-A2 fixes. |
| `baseline_junit_full_3a8c188.xml` | Superseded (through SCI-A3). Full repo test suite JUnit XML (997 passed, 2 xfailed, 0 failed, 319.51 s wall). |
| `baseline_tier1_junit_3a8c188.xml` | Superseded (through SCI-A3). Tier 1 conservation tests JUnit XML (40 passed). |
| `baseline_junit_full_b51df71.xml` | Superseded (CA-1+SCI-N1). Full repo test suite JUnit XML (995 passed, 2 xfailed, 0 failed, 319.15 s wall). |
| `baseline_tier1_junit_b51df71.xml` | Superseded (CA-1+SCI-N1). Tier 1 conservation tests JUnit XML (40 passed). |
| `baseline_junit_full_624ed7c.xml` | Superseded (CA-1 only). Full repo test suite JUnit XML (993 passed, 2 xfailed, 0 failed, 327.77 s wall). |
| `baseline_tier1_junit_624ed7c.xml` | Superseded (CA-1 only). Tier 1 conservation tests JUnit XML (40 passed). |
| `baseline_junit_full_186b5c4.xml` | Superseded (pre-CA-1). Full repo test suite JUnit XML (805 passed, 2 xfailed, 0 failed, 71.36 s wall). |
| `baseline_junit_186b5c4.xml` | `tests/v3` subset JUnit XML (392 passed, 0 xfailed, 11.89 s wall). |
| `baseline_pytest_full_186b5c4.txt` | Verbose pytest text output for the full suite. |
| `baseline_pytest_v3_186b5c4.txt` | Verbose pytest text output for `tests/v3` only. |
| `baseline_pytest_summary_186b5c4.json` | Extracted summary: counts, xfailed list, slow tests >1s, wall times. **Load-bearing for §11.4.** |
| `baseline_tier1_junit_186b5c4.xml` | Tier 1 conservation tests JUnit XML (40 passed, 8.69 s wall). |
| `baseline_tier1_summary_186b5c4.json` | Tier 1 summary: per-test status and timing; asserted contract `rtol=1e-12, clip_events == {}`. **Load-bearing for §11.4.** |
| `baseline_pixi_list_186b5c4.txt` | Full conda+pypi dependency pin (329 packages) via `pixi list`. Re-capture this if dependencies change; do not silently re-baseline. |
| `capture_baseline_trajectory.py` | Script that produced `baseline_coupled_trajectory_*.nc`. Run via `pixi run --environment dev python tests/v3/nsm1/baseline/capture_baseline_trajectory.py <commit>`. |
| `check_baseline_parity.py` | Bit-identical parity check used by every subsequent phase commit. Exit code 0 = match, 1 = mismatch (prints offending variables). Run via `pixi run --environment dev python tests/v3/nsm1/baseline/check_baseline_parity.py`. |

## Reproduction

The baseline was captured in the `dev` pixi environment (only `dev` has the editable installs of `clearwater_modules`, `clearwater_data`, and `clearwater_riverine`):

```sh
pixi run --environment dev python tests/v3/nsm1/baseline/capture_baseline_trajectory.py 186b5c4
pixi run --environment dev pytest tests/ --junitxml=tests/v3/nsm1/baseline/baseline_junit_full_<commit>.xml
pixi run --environment dev pytest tests/v3/nsm1 -k tier1 --junitxml=tests/v3/nsm1/baseline/baseline_tier1_junit_<commit>.xml
```

## What "bit-identical" means here

`check_baseline_parity.py` uses `numpy.array_equal` on the raw `float64` values — `rtol=0, atol=0`, no tolerance whatsoever. A single-bit difference in a single cell at a single substep fails the check. This is the strongest invariant the test suite carries and is the only one that catches the operand-reordering / broadcast-shift class of regressions that motivated the §11.6 refactor-discipline rules.

## When to re-baseline

Only when one of the following changes:

- A dependency version (numpy, xarray, netCDF4, or anything upstream that affects floating-point evaluation order).
- A deliberate kinetics change committed outside the pattern-alignment work.
- A change to the demo's default initial conditions or default parameters that is intentional and reviewed.

Re-baselining is a separate, signed-off commit with its own short hash in the filenames. The old baseline files are not overwritten — they remain in the tree so the history of references is auditable.

## Re-baseline log

### `e0185de` — 2026-05-30 — FloatingAlgae computed light-extinction (NSM1-I port)

Trigger: a deliberate kinetics change (Fortran-port completion,
`design/clearwater_modules_v3_light_extinction.md`). Commit `e0185de` made
`FloatingAlgae.limit_light` compute the light-extinction coefficient lambda each
step from the optical constituents via `utils.light.L`
(`lambda0 + lambdas*Solid + lambdam*POC/fcom + lambda1*Ap + lambda2*Ap**2/3`, a
verified port of `fortran/NSM1/02_global/nsmi_global_params.f90:421-427`), instead
of the constant `light_attenuation_coefficient = 1.0 /m`. In the demo (Ap bloom,
POC present) the computed lambda is ~2.8 /m, so algal light limitation drops.

Scope of change vs `d530a3a`: **broad** — the algal light-limitation change
cascades through the coupled algae→DO→nutrient network; ~16 variables differ
(`algae_floating, alkalinity, ammonium, benthic_algae, cbod, dic, doc, n2,
nitrate, organic_nitrogen, organic_phosphorus, oxygen_dissolved, pathogen, poc,
pom, tip`). Captured and parity-verified under the conda `clearwater` test env
(`test_coupled_demo_parity.py`, 4 passed). All other FloatingAlgae unit tests
(v1-parity, scia3, diagnostics, phase9a1) pass unchanged. **BenthicAlgae is NOT
included** — it keeps its constant-lambda `limit_light` (own benthic form); a
parallel wiring is a follow-up. Prior artifacts (`d530a3a`, `6c10f36`,
`3a8c188`, `b51df71`, `624ed7c`, `186b5c4`) retained unmodified. References in
`test_coupled_demo_parity.py` and `check_baseline_parity.py` now point at
`e0185de`.

Scope note: as with `d530a3a`, only the load-bearing §11.2 trajectory `.nc` was
regenerated; the §11.4 JUnit/summary artifacts (no test consumes them) were not.

### `d530a3a` — 2026-05-30 — Pathogen canonical `algae_floating` name fix

Trigger: a deliberate kinetics change (bug fix). Commit `d530a3a` made the
Pathogen light-extinction / shading term read the canonical floating-algae
state `algae_floating` (falling back to the legacy `ap` only when the
canonical name is absent) instead of the unregistered name `ap`. The coupled
demo registers a floating-algae bloom (`algae_floating = 40 ug-Chla/L`), so
the pathogen die-off light term now sees algal shading instead of zero.

Scope of change vs `6c10f36`: **pathogen only** — `pathogen` differs at all
21,600 cell-substeps (4,320 × 5); all other 19 state/forcing variables are
bit-identical. Captured and parity-verified under the conda `clearwater` test
env (`pytest tests/v3/nsm1/test_coupled_demo_parity.py`, 4 passed). Prior
artifacts (`6c10f36`, `3a8c188`, `b51df71`, `624ed7c`, `186b5c4`) retained
unmodified. References in `test_coupled_demo_parity.py` and
`check_baseline_parity.py` now point at `d530a3a`.

Scope note: this re-baseline regenerated only the load-bearing §11.2
trajectory `.nc`. The §11.4 JUnit / summary / pixi-list artifacts were **not**
regenerated — no test consumes them, and the suite's test *count* has also
shifted from separately-committed work (the TSM thin-water-skip tests). A full
§11.4 re-baseline (full-suite JUnit recapture) is a separate follow-up if
desired; the `6c10f36` JUnit/summary rows below remain the last captured
§11.4 reference.

### `6c10f36` — 2026-05-16 — NSM1-SCI-A2 (terminal gold-standard baseline)

Trigger: a deliberate kinetics change (gold-standard spec Workstream C1; E1 author decision). Commit `6c10f36` fixed NSM1-SCI-A2 (MAJOR): the operative algal/benthic mortality-carbon routing fraction `f_pocp`/`f_pocb` was corrected `0.5 → 0.8` (CE-QUAL-W2 `APOM`; v1 used 0.9), routing dead algal carbon predominantly to POC rather than ~half to DOC.

Scope of change vs `3a8c188`: **broad** — the POC/DOC re-partition cascades through the carbon→DIC→DO network and onward. ~16 variables differ (`algae_floating, alkalinity, ammonium, benthic_algae, cbod, dic, doc, n2, nitrate, organic_nitrogen, organic_phosphorus, oxygen_dissolved, pathogen, poc, pom, tip`). Captured under `pixi --environment dev`; parity re-verified under the conda `clearwater` test env.

**Terminal baseline.** NSM1-CA-1, NSM1-SCI-N1, NSM1-SCI-A3 and NSM1-SCI-A2 are the trajectory-perturbing changes in the gold-standard gate set. The remaining gate items genuinely do not alter the coupled-demo trajectory (re-verified per item): **C2** CBOD settling is dormant at the shipped `ksbod_20=0`; **C3** DOX-F1 is a doc/guard with no freshwater numeric change; **C4** DOX-F2 is a warning + opt-in floor that is off by default; **C5** is documentation-only; **D** is documentation/tests. Absent a new deliberate kinetics change or dependency bump, `6c10f36` is the final gold-standard reference. References in `test_coupled_demo_parity.py` and `check_baseline_parity.py` now point at `6c10f36`; prior artifacts retained unmodified.

### `3a8c188` — 2026-05-16 — NSM1-SCI-A3

Trigger: a deliberate kinetics change (gold-standard spec Workstream B1). Commit `3a8c188` fixed NSM1-SCI-A3 (MAJOR): FloatingAlgae/BenthicAlgae now convert `solar_radiation` (total broadband shortwave) to PAR (`× Fr_PAR`, Fr_PAR=0.47) at the process boundary before light limitation, restoring the v1 convention v3 had dropped (a v1→v3 regression).

Scope of change vs `b51df71`: **broad** — unlike the alkalinity-only CA-1/SCI-N1 fixes, the algae light-limitation change cascades through the coupled network. ~16 variables differ (`algae_floating, alkalinity, ammonium, benthic_algae, cbod, dic, doc, n2, nitrate, organic_nitrogen, organic_phosphorus, oxygen_dissolved, pathogen, poc, pom, tip`). Captured under `pixi --environment dev`; parity re-verified under the conda `clearwater` test env.

Correction (logged in-record): this entry was originally written as the "terminal gold-standard baseline", asserting that **C1 carbon-routing** (along with C2–C5/D) does not perturb the trajectory. That was **wrong** — NSM1-SCI-A2 (C1, `f_pocp/f_pocb`) re-partitions mortality carbon and cascades broadly, so `3a8c188` was superseded by `6c10f36` on the same day. The over-claim is corrected here rather than rewritten away; C2–C5/D were re-verified individually as genuinely non-perturbing (see the `6c10f36` entry). References now point at `6c10f36`; `3a8c188`, `b51df71`, `624ed7c` and `186b5c4` artifacts retained unmodified.

### `b51df71` — 2026-05-16 — NSM1-SCI-N1

Trigger: a deliberate kinetics change (gold-standard spec Workstream A2). Commit `b51df71` fixed NSM1-SCI-N1 (MAJOR): the denitrification alkalinity coefficient `r_alkden` was corrected `4/14/1000` → `1/14/1000` eq/mg-N (1 eq alkalinity per mol NO₃-N reduced; CE-QUAL-W2 `water-quality.f90:3157`, Stumm & Morgan), a deliberate divergence from the upstream Fortran/v1 defect.

Scope of change vs `624ed7c`: **only `alkalinity` differs** (21,600/21,600 cell-substeps; the denitrification source term); all other 19 variables bit-identical. Captured under `pixi --environment dev`; parity re-verified under the conda `clearwater` test env.

Note: this entry was originally logged as the "terminal gold-standard baseline" on the assessment that A1/A2 were the only trajectory-perturbing gate changes. That assessment **omitted SCI-A3 (B1)**, whose Fr_PAR fix is a broad algae-driven cascade; `b51df71` was therefore superseded by `3a8c188` on the same day. Retained for auditability.

### `624ed7c` — 2026-05-16 — NSM1-CA-1 alkalinity kinetics fix

Trigger: a deliberate kinetics change committed outside the pattern-alignment work (gold-standard spec Workstream A1). Commit `624ed7c` fixed NSM1-CA-1 (CRITICAL): the alkalinity algal/benthic coupling now uses the intensive carbon ratios `rca = AWc/AWa`, `rcb = BWc/BWd` instead of the raw stoichiometric weights, correcting a 1000×/100× overstatement of the floating/benthic alkalinity flux.

Scope of change vs `186b5c4` baseline: **only `alkalinity` differs** (21,600/21,600 cell-substeps); all other 19 state/forcing variables are bit-identical. Verified bit-identical across two environments (captured under `pixi --environment dev`; parity re-verified under the conda `clearwater` test env). The `186b5c4` `.nc` and companion XML/JSON/txt artifacts are retained unmodified; references in `test_coupled_demo_parity.py` and `check_baseline_parity.py` now point at `624ed7c`.

## Pathogen warnings during capture

As of `d530a3a`, `Pathogen` prints one warning on the first capture about the optional registry variable `Solid` not being present (the demo carries no suspended solids). This is the once-only `_get_optional` warn-latch; it is part of the baseline behavior and is reproduced (silently after the first emission) on every parity run. Before `d530a3a` a second warning for `ap` was also emitted; the canonical `algae_floating` name fix removed it — the demo registers `algae_floating`, which Pathogen now reads (so the `ap` fallback is never consulted in the demo).
