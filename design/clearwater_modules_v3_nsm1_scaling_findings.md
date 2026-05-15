# v3 NSM1 Cell-Count Scaling Findings — Production Viability

**Date:** 2026-05-14
**Branch:** `streaming`
**Benchmark artifact:** `tests/v3/nsm1/baseline/scaling_benchmark.py`
**Context:** The 5-cell perf benchmark (26 ms/step) and the cProfile analysis raised an existential question: is v3 NSM1 viable at the 500K–2M-cell, 1-year, 15-min-timestep production regime the project requires? This memo answers it with measured scaling data.

---

## The question

Project requirements (from the model lead):

- **Hard requirement:** 500,000-cell models, 1–3 months of simulated time at 15-min timesteps.
- **Aspirational:** 500,000–2,000,000-cell models, 1 year at 15-min timesteps.

15-min timesteps → 8,760 steps for 3 months, 35,040 steps for 1 year.

The 5-cell test mesh (26 ms/step) extrapolates naïvely to catastrophe. It should not be extrapolated: xarray's per-operation overhead is **constant per binary op**, not per cell, so the 5-cell number is ~100% framework tax with no array work to amortise it. The 5-cell mesh is the worst-case regime for the v3 design and the least representative of production.

## Method

`build_nsm1_demo(n_cells=N)` swept N ∈ {5, 1K, 10K, 100K, 500K}, adaptive step counts (200 measured at 5 cells down to 12 at 500K), median + mean + p95 ms/step + peak RSS recorded. Same `pixi --environment dev`, same numpy-backed in-memory registry that the production path uses (the v3 `Model` runs in-memory per chunk — confirmed in the architecture spec design note and the chunked-loop source). Script: `tests/v3/nsm1/baseline/scaling_benchmark.py`.

## Measured scaling

| Cells | mean ms/step | median | p95 | peak RSS | ms/step/Mcell |
|---|---|---|---|---|---|
| 5 | 26.0 | 25.9 | 26.6 | 234 MB | (overhead-dominated) |
| 1,000 | 26.9 | 27.0 | 27.4 | 237 MB | 26,869 |
| 10,000 | 30.9 | 30.8 | 31.6 | 274 MB | 3,092 |
| 100,000 | 68.9 | 66.3 | 77.6 | 621 MB | 690 |
| 500,000 | 267.4 | 267.7 | 288.1 | 1,228 MB | 535 |

**Key observation:** going 5 → 500,000 cells (a 100,000× increase) raises step time only ~10×. The normalised `ms/step/Mcell` is still falling at 500K (535 vs 690 at 100K) — the per-cell cost has not yet reached its linear asymptote; vectorisation is still amortising the fixed overhead. The true asymptotic per-cell cost is ≤ 0.535 µs/cell/step.

## Extrapolation to production targets

Using the 500K-cell per-cell cost (conservative — the asymptote is lower):

| Scale | 3 months (8,760 steps) | 1 year (35,040 steps) |
|---|---|---|
| **500K cells** | **~0.7 hr (≈42 min)** | **~2.6 hr** |
| 1M cells | ~1.3 hr | ~5.2 hr |
| 2M cells | ~2.6 hr | ~10.4 hr |

Peak working-set memory scales ~linearly: 500K cells ≈ 1.2 GB, 2M cells ≈ ~5 GB. Feasible on a workstation. (This is the per-substep working set, independent of run length; see chunking note below.)

## Real-world anchor — Willamette Corvallis–Santiam (authoritative)

The model lead's recollection was "~150,000 cells, ~50 minutes." The authoritative record (`run_provenance.json` from the actual coupled run, git_head `5646d35`, executed 2026-04-27) corrects all three figures:

| Quantity | Recollection | **Recorded (authoritative)** |
|---|---|---|
| Reach | "Santiam-Salem" | Corvallis–Santiam–Albany (Willamette; Salem is the KSLE met station, ~30 km downstream of the mesh outlet, dropped from validation) |
| Cell count | ~150,000 | **591,671 total (586,803 active "real" cells)** |
| Wall clock | ~50 min | **82.97 min** (4,977.93 s; `run_start` → `run_end`, includes ~9.6 GB output write) |
| Simulated period | (1 year inferred) | **30 days** (2014-06-01 → 2014-06-30) |
| Cadence | (15-min inferred) | **Daily** Riverine/NSM1/ESM coupling (`n_loop = n_days − 1 = 29` outer steps) + **6-hour** TSM sub-stepping (`tsm_substeps_per_day = 4`, `tsm_dt_days = 0.25`, 120 TSM substeps) |
| Stack | (NSM1 implied) | **Full TRUE coupled stack**: CW-Riverine advection-diffusion + TSM energy balance + NSM1 kinetics + ESM vegetation (13 constituents, 7 species seed tracers) |

**This anchor does NOT validate the NSM1-only 15-min extrapolation above, and the earlier draft of this section (which accepted the ~150K / ~50-min / ~1-yr recollection) was wrong.** Three things must be kept distinct:

1. **The NSM1-only scaling benchmark** (§ above): 500K cells, 5-min substep, NSM1 kinetics only → 267 ms/step. This is what the extrapolation table is built on.
2. **The recorded coupled run**: ~592K cells, **daily** Riverine/NSM1/ESM + **6-hr** TSM, full stack, 30 days → 83 min wall (≈ 172 s per daily coupled outer step, which internally does 1 transport + 1 NSM1 + 4 TSM + 1 ESM + mesh↔registry sanitisation, plus end-to-end I/O of a 285 MB HDF in and ~9.6 GB out).
3. **The project requirement**: 500K–2M cells at **15-minute** cadence for 1–3 months / 1 year.

The recorded run is at the project's target *cell count* (~592K ≈ 500K target) but at a far *coarser cadence* (daily / 6-hr, not 15-min) and for the *full coupled stack* (not NSM1 alone). It confirms that the coupled stack at ~600K cells is tractable for a 30-day run at coarse cadence (~83 min including I/O). It does **not** confirm a 15-min-cadence coupled run at that scale — that workload is unmeasured.

### What the NSM1-only benchmark still establishes

The NSM1 kinetics component, in isolation, at 15-min cadence (the §extrapolation table) is **not** the bottleneck: ~2.6 hr for a 500K-cell 1-year NSM1-only run. The coupled-run cost is dominated by Riverine transport + TSM sub-stepping + ESM + per-substep mesh↔registry sanitisation + I/O, none of which is in scope for this NSM1 pattern-alignment work.

### The open question this surfaces

A 15-minute-cadence **coupled** run at 500K–2M cells is the actual project workload and **has never been benchmarked**. Extrapolating the one recorded daily-cadence coupled point to 15-min cadence is not defensible without measurement: going from daily to 15-min Riverine/NSM1 cadence is 96× more transport+kinetics outer steps. Whether that is 96× more wall clock depends entirely on how much of the 172 s/daily-step is per-step transport+kinetics vs fixed per-day overhead (mask construction, sanitisation, I/O) — which is not decomposed in the provenance. **This is the measurement the project actually needs, and it is a coupled-orchestrator question, not an NSM1 pattern-alignment question.**

## Verdict against requirements

| Requirement | Result | Status |
|---|---|---|
| 500K cells, 1 month | ~14 min | **Met comfortably** |
| 500K cells, 3 months | ~42 min | **Met comfortably** |
| 500K cells, 1 year | ~2.6 hr | **Met** (coffee-break run) |
| 1M cells, 1 year | ~5.2 hr | Met (long but feasible) |
| 2M cells, 1 year | ~10.4 hr | Feasible as an overnight job |

**v3 NSM1 kinetics, in isolation, is fit for purpose at the project's production scale.** The 5-cell panic was a measurement artifact. The xarray-DataArray design's implicit bet — that vectorisation amortises the per-operation framework overhead at production scale — pays off precisely in the 500K–2M-cell regime, for NSM1 kinetics. It looks catastrophic only at the 5-cell test mesh, the worst case for xarray.

**Caveat (load-bearing): this verdict covers NSM1 kinetics only.** The table above is NSM1-only at 15-min cadence. The actual project workload is the *full coupled stack* (Riverine + TSM + NSM1 + ESM) at 15-min cadence. The only recorded coupled run (§ below) is at *daily/6-hr* cadence and has not been re-measured at 15-min. The end-to-end viability of a 15-min-cadence coupled run at 500K–2M cells is **unmeasured and remains an open question** — it is a coupled-orchestrator measurement, outside the NSM1 pattern-alignment scope.

## On the chunked / dask path

The v3 `Model.__process_loop_chunked` iterates the **time axis**: it loads forcing inputs and writes state outputs per time-window. This is a **memory-management feature for long runs** — a 35,040-step output timeseries at 500K cells is written to disk in windows rather than held in RAM. It does **not** reduce per-substep compute cost; each substep does the same per-cell xarray work the benchmark measured. Consequences:

- The benchmark's ms/step numbers are production-representative regardless of chunking.
- Chunking makes arbitrarily long runs memory-feasible without changing wall-clock.
- The production registry is numpy-backed (not dask); the benchmark path matches production. A future dask-backed chunked path could parallelise across chunks for additional speed, but that is not the current design and is not required to meet the stated requirements.

## On the pattern-alignment +15%

At 500K cells, 15% of 267 ms ≈ 40 ms/step. Over a 1-year run that adds 35,040 × 0.04 s ≈ 23 minutes to a 2.6-hour run. Noticeable but not decisive. The profiling memo's Option A (drop per-sub-flux `sanitize_rate` in DOX/CBOD; sanitise once at the net rate) recovers ~10% (≈ 15 min off the 1-year 500K run) at low-medium risk under the §11 contract — available if a calibration campaign wants it, not required to meet requirements.

## Recommendation

1. **v3 NSM1 meets the project's hard requirement and the 500K-cell aspirational target with comfortable margin.** No architectural rework is required to be production-viable.
2. **The xarray→numpy/numba hot-path rewrite (profiling memo Option C, ~30–40% recoverable) is NOT needed to meet requirements.** It would turn the 2M-cell 1-year run from ~10.4 hr to ~6-7 hr — a nice-to-have for the largest aspirational case, not a survival item. High risk, 1-2 week scope; defer unless 2M-cell annual runs become routine and the overnight turnaround is binding.
3. **Pattern-alignment Option A** (~10%, low-medium risk) is the only optimisation worth considering near-term, and only if a calibration campaign reports the wall-clock as binding.
4. **Validate against the real Willamette run parameters** when convenient: confirm the exact cell count and timestep count of the ~50-min Santiam-Salem run to tighten the extrapolation constant. The interpolation is consistent but a measured production data point is stronger than a synthetic-demo extrapolation.

## What this resolves

- The "NSM is an order of magnitude slower / I am screwed" concern: **the order-of-magnitude gap is v1→v3 (numba vs xarray), present before pattern-alignment, and it does not matter at production scale because the per-operation overhead amortises away.** v3 at 500K cells is ~0.5 µs/cell/step — within a small factor of raw numpy.
- The pattern-alignment +15%: real, small at scale, recoverable, not decisive.
- The design question (NSM1 only): the v3 (xarray) design is appropriate for the production regime *for NSM1 kinetics*. It is the wrong design only if you primarily run small meshes — not the project's use case.

## What this does NOT resolve (escalated)

- **The full-coupled-stack 15-min-cadence cost at 500K–2M cells is unmeasured.** The only recorded coupled run (Willamette Corvallis–Santiam, `run_provenance.json`) is 591,671 cells / 30 days / **daily Riverine-NSM1-ESM + 6-hr TSM** / 82.97 min — *not* the 15-min cadence the project requires. The model lead's "~150K cells / ~50 min / ~1 yr" recollection was inaccurate on every figure (real: ~592K / ~83 min / 30 days / coarse cadence).
- The next decisive measurement is **a 15-min-cadence coupled run benchmark at ~500K cells** — owned by the coupled-orchestrator workstream (`ClearWater-modules-phase2-ESM-streaming/.../08_run_coupled*.py`), not by NSM1 pattern-alignment. Until that exists, statements about whether the project's 1-year, 15-min, 500K–2M-cell target is feasible are extrapolation, not measurement.
