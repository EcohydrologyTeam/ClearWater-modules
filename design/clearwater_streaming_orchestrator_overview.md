# ClearWater Streaming Stack: Orchestrator Overview

**Date:** 2026-05-06
**Audience:** Todd (lead designer) and collaborators reorienting to the streaming stack
**Scope:** Cross-cutting summary of how a model is run today across `ClearWater-Riverine-streaming`, `ClearWater-modules-streaming` (TSM v3, NSM v3), `ClearWater-modules-phase2-ESM-streaming` (ESM + coupled exemplar), and the v2 dependencies that NSM v3 and SSM still rely on.

## Two orchestrators at different layers

There are **two distinct orchestrators** in the streaming stack. Confirm which one is meant when collaborators say "the orchestrator."

### 1. Module-side orchestrator: `clearwater_modules_v3.model.Model`

`ClearWater-modules-streaming/src/clearwater_modules_v3/model.py:68` — the framework piece that TSM v3 and NSM v3 plug into for module-only (no Riverine) runs.

What it owns:

- **Time loop** (`run()` at `model.py:208`) with `__process_loop_full()` and `__process_loop_chunked()` variants.
- **Per-process firing schedule** precomputed at init and indexed by integer step number (`model.py:349`). Drift-immune; replaces per-step modulo arithmetic.
- **Wet-mask gating** applied automatically after each process writes its outputs. Configured via `wet_mask_variable` and `wet_mask_threshold` (`model.py:83-122`); dry cells receive NaN.
- **Hotstart** via `to_hotstart()` / `from_hotstart()` Process hooks. Processes that don't implement them are silently skipped (`model.py:268-270`), so v2 processes stay compatible.
- **Chunked I/O** via `ChunkedZarrDataStore`; `chunk_size` controls window size, decoupled from the kinetic timestep. Chunk loop at `model.py:660-754`.
- **Shared `VariableRegistry`** as the blackboard — no more per-module `.dataset` / `.mesh` objects to thread together.

User-facing change vs. the Willowbend-era `init_from_file(modules.yml)` entry: the YAML entry-point is still supported, but the `Process` tuple, registry, and data-source dict can also be passed directly. Adding a new module is "add a Process to the tuple" rather than "write a custom coupling function."

### 2. Coupled-stack orchestrator: the `08_run_coupled.py` exemplar

`ClearWater-modules-phase2-ESM-streaming/case_studies/santiam_salem/scripts/08_run_coupled.py:1-1242`. **Not yet a framework class** — a hand-written script that directly instantiates and drives:

- `ClearwaterRiverine` (transport solver) — `:528`
- `NutrientBudget` (NSM1 kinetics) — `:549`
- `EcohydrologyModel` (ESM vegetation) — `:716`

Loop structure (`:816-1060`): outer per-day loop, inner per-substep loop that calls `transport.update()` → extracts state from the mesh → `reaction.increment_timestep()` → sanitizes and writes back → ESM runs once per day at the day boundary.

**Flagship feedback** is the seed-tracer pathway (`:997-1060`): ESM writes mature-cell density into a CW-Riverine source term, transport advects it, and ESM reads downstream concentration to drive germination.

A reusable `CoupledOrchestrator` / `ESMOrchestrator` class is on the roadmap but not built. Today the coupled stack lives as an exemplar script that case studies copy and adapt.

## What "streaming" means

From the CW-Riverine README (lines 5–14): **memory-bounded execution for multi-month, 100k+ cell runs on workstation hardware**. Three concrete mechanisms:

1. **Zarr time-slab flushing** — `streaming_interval` controls cadence; flushed timesteps are released from the in-memory mesh (`transport.py:759-761`). `finalize()` consolidates the chunks (`transport.py:1012-1102`).
2. **On-disk mesh cache** — HEC-RAS HDF geometry is parsed once, cached, then reused on subsequent runs.
3. **Checkpoint / resume** — preserves time index, per-constituent mass-flux arrays, boundary state, and streaming bookkeeping.

On the modules side, the analogous concept is the v3 chunked process loop: load a window of data sources, run the schedule for that chunk, write outputs, release, advance.

"Streaming" is **not** about reactive/event-driven kinetics or per-cell streaming computation — it's about not materializing the full simulation domain in memory at once.

## Where TSM v3, NSM v3, and SSM fit

- **TSM v3 / NSM v3** — `src/clearwater_modules_v3/processes/`. Each subclasses `Process` (re-exported from v2 via `processes/base.py:67`, which is the v2 dependency in NSM v3) and implements `run(time, registry)`. The v3 `Model` dispatches them via the firing schedule. Two endorsed update patterns in `processes/base.py:9-28`:
  - **(a)** Compute rate of change and apply Forward Euler: `state_new = state_old + rate * dt_seconds`.
  - **(b)** Compute per-substep delta directly and add to state (e.g., Temperature's depth-ramp guard).
- **SSM** — Sediment Simulation Module. The Python port lives in `clearwater_modules_v2`. Its integration into the coupled run is one of the open questions on the follow-up list (see § Open questions below). In the Santiam-Salem coupled exemplar today, sediment forcing is **read from a precomputed NetCDF** at `08_run_coupled.py:739` and passed to ESM via the `HydraulicData` namedtuple — meaning v2 SSM Python is not invoked inside the per-substep loop on this branch, even though the v2 SSM source still exists.

## The user-workflow arc, condensed

| Generation | User instantiates | Time loop | State exchange |
|---|---|---|---|
| Original (`ClearWater-modules` + `ClearWater-riverine`) | Two objects, separately | User writes the `for` loop and a custom coupling function | Manual xarray indexing + name-mapping dicts |
| Willowbend publications (intermediate) | `init_from_file(modules.yml)` returns one coupled `model` | `model.run()` — hidden | Internal registry |
| Streaming (Willamette + Santiam-Salem) | v3 `Model(processes=(...), registry=..., data_sources=...)` for modules; hand-wired script for the full coupled stack | `model.run()` for v3 modules; explicit per-day script for the coupled stack | `VariableRegistry` blackboard for v3; mesh ↔ registry sanitization in the coupled script |

## Follow-up findings (resolved 2026-05-06)

### 1. SSM v2 wiring — dormant in the coupled stack

**Verdict:** v2 SSM Python exists but is not invoked by any coupled driver on streaming. The previous-section claim that sediment is "read from a precomputed NetCDF" is correct, and the precomputation does not involve any SSM Python code.

- v2 SSM source lives at `ClearWater-modules-streaming/src/clearwater_modules_v2/processes/sediment/` — main class `SSM(Process)` at `ssm.py:183` with `run()` at `ssm.py:697`, registered via `@ProcessFactory.register("sediment")` at `ssm.py:292`, plus 13 supporting modules (armoring, bed, bedload, classes, consolidation, contracts, coupling, deposition, erosion, settling, shear, io/*).
- Zero `from clearwater_modules_v2.processes.sediment` or `SSM(` imports across `ClearWater-modules-phase2-ESM-streaming/`, `ClearWater-Riverine-streaming/`, `clearwater_modules_v3/`, or `LargeProjects/ClearWater-riverine-case-study-Willamette/`.
- In the coupled exemplars, sediment fields are generated by `05_synthesize_sediment.py` using **synthetic physics** (velocity-driven scour, deposition proportional to sediment supply) and written to `synthetic_sediment.nc`, then read at `08_run_coupled.py:739`. Comment in stage 5 explicitly notes the USGS HEC-RAS plan was run with `Run Sediment=0`, so all sediment information is synthesized.
- `design/ssm_improvement_plan.md:11-14` flags four CRITICAL defects preventing production use:
  - **C1** (`bed.py:630-712`) — `ssm_bed_elevation` allocated but never written → ESM reads zeros for scour mortality.
  - **C2** (`ssm.py:759`) — `ssm_bed_d50_surface` computed but never persisted → ESM reads zeros for biostabilization.
  - **C3** (`bed.py:327`) — per-layer `tau_crit` set at t=0 only and never propagated → erosion threshold collapses.
  - **C4** (`ssm.py:1024-1028`) — sediment source injection into Riverine has no consumer; loop is open.
- `design/ssm_improvement_plan.md:64`: "SSM is **not** production-ready... two of ESM's five expected inputs are zeros, the riverine source-injection loop is open."
- Consistent with prior memory: CSM/MSM/NSM2 prototypes were retired 2026-05-04 (commit 8deca2f); SSM is **not** in the retirement list because it was never activated in production in the first place. v3 SSM port is deferred until the four CRITICAL defects are fixed.

### 2. Willamette case study driver state — split across two repos, drivers exist

The "Willamette case study" lives in **two complementary locations**:

- **Data + design + observations:** `/Users/todd/LargeProjects/ClearWater-riverine-case-study-Willamette/` — HEC-RAS model files (`Corvallis_Santiam.{g01,p01,u01-u03}.hdf`), boundary conditions, observed data, design docs (`master-case-study-plan.md`, `reach-selection.md`, `site-assessment.md`, `Case_study_demonstration_needs_and_repositories.md`). The only Python file here is `scripts/make_study_area_map.py`. **No coupled drivers.**
- **Drivers + outputs:** `ClearWater-modules-phase2-ESM-streaming/case_studies/corvallis_santiam_albany/` — full pipeline (stages 01–13) for the same Corvallis-Santiam reach of the Willamette. This **is** the Willamette driver set; it was forked from `santiam_salem` and has diverged.

Three coupled drivers exist in `corvallis_santiam_albany/scripts/`:

| Driver | Lines | Cadence | Stack |
|---|---|---|---|
| `08_run_coupled.py` | 1647 | daily | CW-Riverine + NSM1 (v1/v2) + ESM; baseline that ships today |
| `08_run_coupled_15min.py` | 2484 | 15-min substeps inside a daily outer loop | Same as baseline + chunked per-day NetCDF writes from `_orchestrator_helpers.py` |
| `08_run_coupled_v3_smoke.py` | 1191 | multi-rate (master = HDF cadence; TSM subcycles; ESM macro-steps) | **`clearwater_modules_v3` + `ClearWater-Riverine-streaming`** — first driver against the new stack |

`_orchestrator_helpers.py` (551 lines) already factors out the chunked-I/O abstractions: `find_latest_checkpoint`, `allocate_day_buffer`, `flush_day_chunk`, `DailyAggregator`, `concatenate_per_day_chunks`. So abstraction has already started — incrementally, not big-bang.

`design/refactor_strategy.md` documents the approved (2026-04-30) workstream split: A = CW-Riverine `transport.mesh` streaming, B = ESM in-memory streaming, C = NSM1 optimization (LimnoTech handoff). The 15-min and v3-smoke drivers are direct outputs of workstreams A and B as they progress.

The diverged-copy problem: `corvallis_santiam_albany/scripts/08_run_coupled.py` differs from `santiam_salem/scripts/08_run_coupled.py` (confirmed via `diff -q`). Changes do not auto-propagate.

### 3. Promote `08_run_coupled.py` to a `CoupledOrchestrator` class? — recommend incremental, not big-bang

**Recommendation: defer the class extraction; keep widening `_orchestrator_helpers.py` first.**

Evidence for promotion now:
- Three drivers, 5,322 lines combined, with substantial overlap.
- Two case studies (`santiam_salem`, `corvallis_santiam_albany`) hold diverged copies of `08_run_coupled.py` — manual propagation is already a tax.
- `_orchestrator_helpers.py` proves the appetite for shared abstraction exists.

Evidence to defer:
- The v3 smoke test (`08_run_coupled_v3_smoke.py`) introduces a multi-rate scheduling contract (master + TSM subcycling + ESM macro-stepping) that is **not yet validated at scale** — only a 4-day smoke run.
- ESM streaming (Workstream B) is in flight. The orchestrator API will need to change to handle ESM's new per-day on-disk slice.
- The seed-tracer / ESM-feedback contract is still evolving.
- Three drivers is a thin sample for choosing the abstraction shape.

**Suggested phasing:**

1. **Now:** continue widening `_orchestrator_helpers.py` to cover multi-rate scheduling. Make `08_run_coupled_15min.py` and `08_run_coupled_v3_smoke.py` both call into the same scheduling helpers.
2. **After the v3 smoke test passes a 30-day Albany run AND ESM streaming (Workstream B) lands:** extract a `CoupledOrchestrator` class into `phase2-ESM-streaming/src/esm/orchestration/`. Two configuration variants to start: daily-cadence (Santiam-Salem) and multi-rate (v3). Keep the case-study scripts as thin wrappers that build the configuration object.
3. **After a third independent case study lands using the class:** consider promoting `CoupledOrchestrator` to `clearwater_modules_v3.coupling` so it becomes part of the framework rather than the ESM repo.

Decision signal to watch: how often a change to one of the three drivers requires the same edit in another. That ratio rising past ~50% is the right moment to extract — the drivers will be telling you the contract has stabilized.

## Key file references

| Concern | File |
|---|---|
| v3 Model orchestrator | `ClearWater-modules-streaming/src/clearwater_modules_v3/model.py:68` |
| v3 Process base | `ClearWater-modules-streaming/src/clearwater_modules_v3/processes/base.py` |
| Transport solver | `ClearWater-Riverine-streaming/src/clearwater_riverine/transport.py:54` (class), `:618` (`update()`), `:1012-1102` (`finalize()`) |
| Coupled exemplar | `ClearWater-modules-phase2-ESM-streaming/case_studies/santiam_salem/scripts/08_run_coupled.py:1-1242` |
| Seed-tracer feedback | `08_run_coupled.py:997-1060` |
| Sediment forcing handoff | `08_run_coupled.py:739` |
