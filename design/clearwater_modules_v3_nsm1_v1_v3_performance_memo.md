# v3 NSM1 — Santiam-Salem Provenance and v1-vs-v3 NSM1 Performance

**Date:** 2026-05-15
**Branch:** `streaming`
**Status:** Draft (not committed pending review)
**Scope:** (1) the verified provenance of the Santiam-Salem coupled validation run, and (2) a standalone, head-to-head v1-vs-v3 NSM1 per-step performance and memory comparison. Kinetics correctness is out of scope here; this memo is performance and provenance only.

This memo supersedes an earlier scaling memo that was withdrawn (commits `3026701` and `50a911a`, reverted by `ce4a4ec` and `ed553f3`) because its real-world anchor was sourced from a stale, unrelated run. All figures below are traced to a named source file.

---

## 1. Purpose

Two questions:

1. What are the exact parameters of the Santiam-Salem coupled run used for current validation work?
2. At the production grid sizes of interest (up to ~1,000,000 cells on the available workstation), how does v3 NSM1 per-step cost and memory footprint compare to v1 NSM1?

---

## 2. Verified Santiam-Salem coupled-run provenance

**Source:** `ClearWater-modules-phase2-ESM-streaming/case_studies/santiam_salem/output/v3_smoke_15day_wind10m_final_mumax_1_0/run_provenance.json` (git_head `a688d5b`, `run_start` 2026-05-13T00:33:27, `run_end` 2026-05-13T00:44:09). Three equivalent mu_max-sweep runs that night each completed in ≈11 minutes (`_final_sweep_20260512.log`).

| Quantity | Value |
|---|---|
| Reach | Santiam-Salem (Willamette) |
| Grid | `santiam_salem_subset_2008-09_hourly.p01.hdf` — trimmed subset prepared by the orchestrator |
| Cells | 159,634 total / 158,037 active (`nreal`) |
| Simulated period | 15 days (2008-09-01 → 2008-09-15) |
| Cadence | transport / TSM / NSM1 all hourly (`master_step_hours = 1.0`); ESM daily (`esm_step_hours = 24.0`) |
| Stack | v3 SMOKE: streaming CW-Riverine advection-diffusion + v3 NSM1 (11 Process objects on `InMemoryRegistry`) + v3 Temperature + ESM vegetation (13 constituents, 7 species seed tracers) |
| Wall clock | 641.16 s = 10.69 min (end-to-end stage, includes output writes) |

This is the configuration in active use for validation. An earlier draft of the withdrawn memo cited a different run (`Corvallis_Santiam/output/run_provenance.json`, dated 2026-04-27, 591,671 cells, daily cadence, no v3 TSM) that is a separate, diverged case study and not the run in use; that citation has been fully retracted.

---

## 3. v1-vs-v3 NSM1 performance — method

A standalone harness times **only the per-step kinetics advance**, with no transport and no orchestrator, so the comparison isolates the NSM1 implementations:

- **v1:** `clearwater_modules.nsm1.model.NutrientBudget`, constructed exactly as the orchestrator constructs it (same `global_parameters` / `global_vars` / `algae_parameters` dicts, `NSM_ACTIVE_STATES = [Ap, NH4, NO3, TIP, DOX]` with the other 11 states zeroed, `updateable_static_variables = ["q_solar", "TwaterC"]`, `track_dynamic_variables = False`, `time_dim = "days"`). Advanced via `increment_timestep(inputs)`.
- **v3:** the 11-Process NSM1 stack on `InMemoryRegistry` via `build_nsm1_demo(n_cells=N)`. Advanced via `demo.step(t)`.
- **Identical synthetic initial conditions** for both (constant per state, from the v3 demo defaults / Santiam-Salem provenance ICs: NH4 0.02, NO3 0.137, TIP 0.029, DOX 9.4, Ap 1.6, water_temp 17.35, solar 200). Constant fields are sufficient because the comparison is per-step cost, not field realism.
- **Timed region:** only the per-step advance call, after a warmup, excluding one-time construction and dataset allocation.
- **Memory:** peak resident set size (RSS) measured at 1,000,000 cells with each engine run in its own process, so the figure is not cross-contaminated by the other engine's allocations.

Harness scripts are transient (not committed); see §8.

---

## 4. Speed results (5 → 1,000,000 cells)

Per-step wall time, single run, one workstation. 2,000,000 cells was not used (exceeds available memory on the target hardware).

| Cells | v1 ms/step | v3 ms/step | v1/v3 ratio | Faster |
|---|---|---|---|---|
| 5 | 5.2 | 27.3 | 0.19 | v1 (~5×) |
| 1,000 | 5.7 | 27.5 | 0.21 | v1 (~5×) |
| 10,000 | 11.0 | 31.3 | 0.35 | v1 (~3×) |
| 100,000 | 68.8 | 70.8 | 0.97 | tied |
| 500,000 | 324.2 | 282.2 | 1.15 | v3 (~15%) |
| 1,000,000 | 630.9 | 512.2 | 1.23 | v3 (~23%) |

`v1/v3 > 1` means v1 is slower. The crossover is near ~100,000 cells: below it v1 is faster (up to ~5× at small/debug meshes); at ~100K they are tied; at and above 500K v3 is faster, and the v3 margin widens with size (15% at 500K, 23% at 1M).

Mechanism (per-step xarray dispatch is a fixed cost per operation, not per cell): at small N, v3's framework overhead dominates and v1's lower-overhead numerical path wins. As N grows, v3's fixed overhead amortizes across the array while v1's per-cell cost does not improve, so the curves cross and v3 pulls ahead in the production regime.

---

## 5. Memory results (1,000,000 cells, isolated processes)

| Engine | ms/step @ 1M | Peak RSS @ 1M |
|---|---|---|
| v3 | 512.2 | 1.69 GB |
| v1 | 630.9 | 4.46 GB |

At 1,000,000 cells v1 NSM1's peak RSS is ~2.6× v3 NSM1's. v1's footprint is driven by its `(time_steps+1, n_cells)` state-history allocation and the function-style dynamic-variable framework; v3's `InMemoryRegistry` overwrites state in place and does not accumulate per-step history. Linear extrapolation of the footprint (indicative, not measured): v1 ≈ ~9 GB at 2M, v3 ≈ ~3.4 GB at 2M.

---

## 6. Interpretation

For the production regime of interest (hundreds of thousands to ~1,000,000 cells on the available workstation), the standalone NSM1 comparison shows v3 to be both faster per step and substantially lighter on memory than v1. v1's advantage is confined to small meshes (≤ ~10,000 cells), which corresponds to the test/debug regime rather than production. On hardware where memory is the binding constraint, v3's ~2.6×-smaller footprint at 1M cells is the more decisive difference, because it sets the largest grid that fits in RAM.

---

## 7. Limitations (read before citing any number)

- **NSM1-only.** No transport, no TSM, no ESM. The coupled-stack cost (which includes CW-Riverine advection-diffusion, v3 Temperature, ESM, and per-substep mesh↔registry handling) is not measured here and is not implied by these figures.
- **Synthetic, constant initial conditions.** Kinetics are per-cell, so per-step cost should be insensitive to spatial heterogeneity, but this is an assumption, not a measured result on a real heterogeneous grid.
- **Single-run timing on one workstation.** The v1/v3 ratio moves monotonically and the crossover is unambiguous, so the direction is robust; the absolute milliseconds are indicative, not error-bounded. No multi-trial statistics.
- **v1 measured in its production configuration** (`track_dynamic_variables=False`). With dynamic-variable tracking enabled, v1 memory and time would be substantially higher; that configuration was not measured because it is not how v1 is run in the coupled pipeline.
- **Cadence.** The standalone harness steps a fixed cadence for timing; it does not model the coupled run's hourly transport/TSM/NSM1 + daily ESM scheduling.

---

## 8. Reproducibility

The harness scripts (`/tmp/cw_v1_v3_speed_sweep.py`, `/tmp/cw_1m_one.py`) are transient and not part of the repository. If a reproducible, version-controlled benchmark is wanted, a properly scoped script can be added under `tests/v3/nsm1/baseline/` as a separate, reviewed change. The verified provenance file in §2 is the authoritative source for the Santiam-Salem run parameters and is already on disk in the case-study output directory.

---

## 9. What this does and does not resolve

**Resolves:** at production NSM1 scale on the available hardware, v3 NSM1 is not slower than v1 NSM1 — it is faster and uses materially less memory. The concern that v3 represents an order-of-magnitude per-step regression is specific to small meshes and does not hold at production scale.

**Does not resolve:** the end-to-end cost of the full coupled stack (Riverine + TSM + NSM1 + ESM) at the project's target cell counts and the project's intended timestep cadence. The only recorded coupled run is the 15-day, hourly, ~158K-active-cell Santiam-Salem run in §2 (10.69 min). Coupled-stack scaling to larger grids and longer periods is a coupled-orchestrator measurement and is outside the scope of this NSM1 memo.
