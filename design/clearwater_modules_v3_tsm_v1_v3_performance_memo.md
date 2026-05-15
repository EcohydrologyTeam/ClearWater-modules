# v3 TSM — v1-vs-v3 TSM Performance

**Date:** 2026-05-15
**Branch:** `streaming`
**Status:** Draft (not committed pending review)
**Scope:** Standalone, head-to-head v1-vs-v3 **TSM** (Temperature Simulation Module) per-step energy-balance timing and peak memory. Kinetics/numerical-correctness is out of scope; performance only. Companion to `clearwater_modules_v3_nsm1_v1_v3_performance_memo.md` (the NSM1 equivalent).

**Read this alongside the NSM1 memo — the conclusions are not the same.** The NSM1 result (v3 faster than v1 at production scale) does **not** transfer to TSM. They are reported in separate memos precisely so the findings are not conflated.

---

## 1. Purpose

At the production grid sizes of interest (up to ~1,000,000 cells on the available workstation), how does v3 TSM per-step cost and memory footprint compare to v1 TSM?

---

## 2. Method

A standalone harness times **only the per-step energy-balance advance**, no transport, no orchestrator, isolating the TSM implementations:

- **v1:** `clearwater_modules.tsm.model.EnergyBudget`, constructed exactly as the dedicated v1 TSM driver constructs it (`ClearWater-modules-phase2-ESM-streaming/case_studies/corvallis_santiam_albany/scripts/12_run_tsm_hourly.py`): state = `{water_temp_c, volume, surface_area}`; the six meteo forcings (`air_temp_c, q_solar, wind_speed, eair_mb, pressure_mb, cloudiness`) declared `updateable_static_variables`; `use_sed_temp=False`; `track_dynamic_variables=False`; `time_dim="hours"`. Advanced via `increment_timestep(inputs)`.
- **v3:** `clearwater_modules_v3.processes.temperature.Temperature` on an in-memory registry seeded with the 11 input variables; `use_sediment_temperature=False` to match v1's `use_sed_temp=False`. Advanced via `run(time, registry)`. v3 `Temperature` skips its first `run` call by design (`__skip_first_time_step`); the warmup absorbs the skipped step so it is outside the timed region.
- **Identical synthetic constant IC / forcing** for both (water_temp 17.35 °C, air_temp 20 °C, q_solar 400 W/m², wind 3 m/s, eair 1.0 mb, pressure 1013 mb, cloudiness 0.1, volume 1.5, surface_area 1.0, hourly dt). Constant fields are sufficient because the comparison is per-step cost, not field realism.
- **Timed region:** only the per-step advance call, after warmup, excluding one-time construction.
- **Memory:** peak RSS at 1,000,000 cells, each engine in its own process (so the figure is not cross-contaminated).

Harness: `tests/v3/nsm1/baseline/v1_v3_tsm_benchmark.py` (modes `sweep` and `isolated`, default ceiling 1,000,000 cells; 2,000,000 excluded — exceeds the target workstation's memory for v1).

> **Note on the v1 construction reference.** `12_run_tsm_hourly.py` is the only dedicated *standalone* v1 TSM driver in the repository and is used here solely as the authoritative recipe for how the v1 `EnergyBudget` object is constructed (state set, updateable statics, `use_sed_temp`, `track_dynamic_variables`, `time_dim`). It is **not** a validation-provenance source. Santiam-Salem runs TSM only in coupled mode (`08_run_coupled*.py`); there is no standalone Santiam-Salem TSM script, so none could be cited. The construction recipe is reach-independent — it specifies object configuration, not run inputs. The validation anchor for this work remains the Santiam-Salem coupled provenance (this memo §7; NSM1 memo §2).

---

## 3. Speed results (5 → 1,000,000 cells)

Single run, one workstation.

| Cells | v1 ms/step | v3 ms/step | v1/v3 ratio | Faster |
|---|---|---|---|---|
| 5 | 0.7 | 8.6 | 0.08 | v1 (~12×) |
| 1,000 | 0.7 | 9.2 | 0.08 | v1 (~12×) |
| 10,000 | 1.3 | 9.7 | 0.13 | v1 (~8×) |
| 100,000 | 7.0 | 17.1 | 0.41 | v1 (~2.4×) |
| 500,000 | 36.1 | 53.8 | 0.67 | v1 (~1.5×) |
| 1,000,000 | 73.9 | 112.6 | 0.66 | v1 (~1.5×) |

`v1/v3 > 1` means v1 slower. **v1 TSM is faster than v3 TSM at every cell count tested, including 1,000,000.** The gap narrows monotonically with scale (~12× at small meshes to ~1.5× at 500K–1M) but **no crossover was observed within the 1,000,000-cell ceiling**.

---

## 4. Memory results (1,000,000 cells, isolated processes)

| Engine | ms/step @ 1M | Peak RSS @ 1M |
|---|---|---|
| v3 | 100.4 | 1.04 GB |
| v1 | 69.7 | 2.25 GB |

At 1,000,000 cells v1 TSM is ~1.4× faster per step but uses ~2.2× the memory of v3 TSM. v1's footprint is driven by its `(time_steps+1, n_cells)` state-history allocation and the function-style framework; v3's registry overwrites in place.

(The §3 and §4 v3 timings differ — 112.6 vs 100.4 ms — because §3 runs both engines in one process and §4 is isolated; this is expected measurement variance, not a discrepancy. The direction and ratio are stable across both.)

---

## 5. Interpretation, and the contrast with NSM1

For NSM1 (separate memo) v3 overtakes v1 near ~100,000 cells and is faster at production scale. **TSM behaves differently: v1 stays faster through 1,000,000 cells.**

The mechanism is the per-operation nature of the xarray overhead. v3's fixed framework cost is paid per array operation and amortizes as cell count grows. The 11-Process NSM1 stack performs on the order of hundreds of array operations per step, so it carries a large fixed overhead that amortizes substantially at scale — enough for v3 to overtake v1. TSM is a single state variable with a comparatively light energy-balance computation and far fewer array operations per step, so there is much less fixed overhead to amortize; v1's lower-overhead numerical path retains its advantage across the tested range. The monotonic narrowing of the v1/v3 gap (12× → 1.5×) is consistent with a crossover at some larger cell count, but no crossover was measured and none is claimed.

Memory direction is the same for both modules (v3 leaner), with a smaller magnitude for TSM (~2.2×) than for NSM1 (~2.6×).

---

## 6. Limitations (read before citing any number)

- **TSM-only.** No transport, no NSM1, no ESM. The coupled-stack cost is not measured here and is not implied.
- **Synthetic, constant IC / forcing.** Per-step cost should be insensitive to spatial heterogeneity (the energy balance is per-cell), but this is an assumption, not a measured result on a real heterogeneous grid.
- **Single-run timing on one workstation.** The v1/v3 ratio moves monotonically, so the direction is robust; absolute milliseconds are indicative, not error-bounded. In-process sweep RSS is unreliable (see §4 note); memory figures use isolated processes.
- **v1 measured in its production configuration** (`track_dynamic_variables=False`, `use_sed_temp=False`); v3 set `use_sediment_temperature=False` to match. Enabling sediment-temperature evolution or dynamic-variable tracking would change both engines and was not measured.
- **No crossover claim.** v1 is faster at all tested sizes; whether v3 would overtake beyond 1,000,000 cells is not measured and is not asserted.

---

## 7. What this does and does not resolve

**Resolves:** at the production TSM scale on the available hardware, v1 TSM is faster per step than v3 TSM (~1.5× at 500K–1M), while v3 TSM uses ~2.2× less memory. The NSM1 finding (v3 faster at scale) is module-specific and does not generalize to TSM.

**Does not resolve:** the end-to-end coupled-stack cost; whether a speed crossover exists beyond 1,000,000 cells; and whether the v3 TSM per-step cost matters in practice given that, in the recorded Santiam-Salem coupled run, water temperature is computed by v3 TSM at the same hourly cadence as NSM1 and the whole 15-day coupled run completes in ~11 minutes (see the Santiam-Salem provenance in the NSM1 memo §2). Per-module micro-benchmarks bound the components; they do not substitute for a coupled-stack measurement at the project's target cadence and scale.
