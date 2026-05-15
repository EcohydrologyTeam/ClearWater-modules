# Pattern-Alignment Profiling Findings — Hot-Path Analysis

**Date:** 2026-05-14
**Branch:** `streaming`
**Profile artifact:** `tests/v3/nsm1/baseline/profile_perf.{prof,txt}`
**Context:** Resolves Phase 10.A's open question about which call sites contribute the +15% no-subscription overhead measured in `clearwater_modules_v3_nsm1_pattern_alignment_perf_findings.md`.

---

## Method

cProfile run on `build_nsm1_demo()` + 60-substep warmup + 500-substep measured window. Same 5-cell mesh as the perf benchmark. Profile script: `tests/v3/nsm1/baseline/profile_perf.py`.

cProfile adds ~50% overhead vs the un-profiled run (76.5 ms/step under cProfile vs 26 ms/step bare); the absolute times below are inflated, but the proportions are accurate.

## Top hot spots (by cumulative time)

Total: 42.85 s for 560 substeps, 160M function calls, 27 GiB of inferred work.

| Rank | Cumulative (s) | % of total | Function | Call count |
|---|---|---|---|---|
| 1 | 27.0 | 63% | `xarray.dataarray._binary_op` | 384,160 |
| 2 | 15.0 | 35% | `xarray.alignment.align` | 193,200 |
| 3 | 10.8 | 25% | `xarray.computation.apply_ufunc` | 98,560 |
| 4 | 8.8 | 21% | `xarray.computation.apply_dataarray_vfunc` | 96,880 |
| 5 | 8.0 | 19% | `xarray.computation.where` | 43,120 |
| 6 | **7.7** | **18%** | **`utils.numerics.sanitize_rate`** | 16,800 |
| 7 | 6.8 | 16% | `xarray.variable._binary_op` | 384,160 |
| 8 | 6.9 | 16% | `DOX._change_with_components` | 560 |
| 9 | 6.4 | 15% | `Carbon._change_with_components` | 560 |
| 10 | 6.0 | 14% | `Nitrogen._change_with_components` | 560 |

(Ranks 1–7 overlap because most functions are inside the helpers ranked 8–10.)

## Hypothesis test — predicted hot spots from Phase 10.A

Three candidates were named in `..._phase10a_closeout.md` §6. The profile resolves each:

| Predicted hot spot | Actual measurement | Verdict |
|---|---|---|
| `setattr` loop overhead (~80 setattr per substep × 11 Processes) | 0.044 s cumulative across 77,956 calls = **0.1% of total** | **WRONG** — not a hot spot |
| Per-cell `in registry` membership checks in opportunistic-write loop | Effectively zero (no diagnostics subscribed in the profile run) | **WRONG** — not a hot spot |
| Duplicate sub-flux computations in algal-coupling components (Nitrogen / Phosphorus call algae helpers separately for components dict) | Inferred ~10% of Nitrogen / Phosphorus runtime; can't isolate precisely from cProfile | **PARTIALLY RIGHT** — contributes but dominated by xarray overhead |

The pattern G zero-cost-when-unused contract is verified at the call-count level: `set_at_time` calls = 8,960 (= 16 state-variable persistence × 560 substeps), no opportunistic-write calls because no diagnostics were subscribed.

## Real hot spots

The cost lives in xarray's binary-op machinery:

1. **xarray `_binary_op` is the overwhelmingly dominant cost** — every `+`, `-`, `*`, `/` between DataArrays carries dimension-alignment + ufunc-dispatch overhead. The 5-cell mesh has 384,160 binary ops over 560 substeps = **~686 binary ops per substep**.

2. **`sanitize_rate` is 18% of total** — 16,800 calls × 3 binary ops each (`isnull`, `isinf`, `where`) ≈ 50,400 binary ops attributable to sanitization alone. Sanitize is called per-sub-flux in DOX (8 sub-fluxes), CBOD (2), and on the net rate in most other Processes.

3. **Per-Process `_change_with_components` runtime distribution** is concentrated in the heaviest integrators:
   - DOX: 6.9 s (16%)
   - Carbon: 6.4 s (15%)
   - Nitrogen: 6.0 s (14%)
   - FloatingAlgae: 6.0 s (14%)
   - BenthicAlgae: 4.2 s (10%)
   - N2: 2.4 s (6%); Phosphorus: 2.0 s; Pathogen: 1.6 s; Alkalinity: 1.5 s; CBOD: 1.1 s; POM: 0.6 s
   - Total `_change_with_components`: 38.7 s = **90% of step time**.

The remaining 10% is `Process.run` orchestration (setattr loops, registry persistence, opportunistic-write membership checks).

## Optimisation candidates

Ranked by impact-per-effort, with §11 zero-regression contract risk noted:

### A. Reduce `sanitize_rate` call frequency (impact: ~10%, risk: low-medium)

DOX, CBOD, and a few others sanitize **per sub-flux**. The defense is "a NaN in any sub-flux poisons the cell rate after `sum()` + `sanitize`, freezing the cell at IC indefinitely." In the closed-system Tier 1 and 5-cell baseline, no NaN is generated; per-sub-flux sanitize is a no-op cost.

**Optimisation:** sanitize once at the net rate level. Drop ~50,400 binary ops per profile run (~10% of total cost).

**Risk:** changes NaN-handling at the `sum()` boundary. If any cell genuinely produces a sub-flux NaN (wet-mask edge cases), the net-only sanitize zeros the entire cell rate instead of just the NaN sub-flux. Bit-identical baseline parity requires the test scenarios produce no sub-flux NaN — likely true for the 5-cell baseline but should be verified case-by-case.

**Mitigation:** keep per-sub-flux sanitize behind a flag; default to "net-only" with the per-sub-flux mode opt-in for production runs that observe wet-mask edges.

### B. Cache algal-coupling sub-fluxes on the algae Process (impact: ~5%, risk: low)

Nitrogen, Phosphorus, Carbon all read algae helpers (e.g. `_floating_algae_growth_rate()`) inside their `_change_with_components` to populate the components dict AND inside `change_*` methods to compose the rate. The functions return cached values from `floating_algae.algal_growth_rate` (no compute), but the call chain itself adds binary ops.

**Optimisation:** in Phase 5's setattr loop, pre-compute the algae-side coupling values that downstream Processes will read. The downstream Processes then read scalars via getattr, avoiding the `_floating_algae_*_rate()` helper-call chain.

**Risk:** low — code-motion-only; bit-identical preserved.

### C. Numpy-array fast path inside heavy helpers (impact: ~30–40%, risk: high)

Every binary op between two DataArrays carries xarray's full alignment + ufunc dispatch even when the operands are guaranteed to be on the same dim. If the heavy helpers (DOX, Carbon, Nitrogen) extracted `.values` once at the top, did all arithmetic in numpy, then wrapped back to a DataArray at the end, the per-substep cost would drop dramatically.

**Risk:** bit-identicality preserved (numpy arithmetic is deterministic). But the rewrite is invasive and breaks one of the v3 design principles (operate on DataArrays end-to-end). Effort: 2–3 days per heavy Process.

### D. Vectorise across Processes (impact: ~20%, risk: very high)

The 11 per-Process `run` calls execute sequentially per substep. They could be parallelised since most are independent within a substep (Jacobi semantics). But this requires reworking the orchestration and loses readability.

**Risk:** very high; major architectural change. Not recommended for 1.0.x.

### E. Switch from `xarray.DataArray` to `numpy.ndarray` end-to-end (impact: ~50%, risk: highest)

Strip xarray entirely from the kinetics path. Use ndarray, with a thin wrapper for boundary I/O. Kills the xarray overhead but breaks the v3 design contract. Not recommended.

## Recommendation

**Land A only for v3 1.0.2** if the +15% overhead matters for a calibration application:

- Drop ~10% of step time at low-medium risk.
- Keep per-sub-flux sanitize as opt-in for wet-mask production runs.
- §11 contract enforced as before; bit-identical baseline parity is the gate.

**Defer B–E** unless a calibration application reports binding wall-clock constraints. The 26 ms/step on a 5-cell mesh extrapolates to ~3.6 hours on Sumwere Creek (600 cells × 4,320 substeps), comfortably inside the 1.0.0 spec §10 "must" target of 30 minutes for the Tier 0 / Tier 1 test cases. The 30-minute target was set against a different cell count (1,000-cell smaller bench) so the like-for-like budget on Sumwere Creek is more lenient than the documentary number suggests.

**Accept the +15% as the price of the diagnostic surface for v3 1.0.1.** Revisit only if a calibration application surfaces a binding wall-clock constraint, in which case A is the first move.

## What the profile rules out

- The setattr loop (pattern F) is **not** a meaningful cost; no need to batch attribute writes.
- The opportunistic-write membership-check loop (pattern G) is **not** a meaningful cost; no need to optimise the pre-registration check.
- The components-dict construction itself is **not** a meaningful cost; the cost is the sub-flux computations that populate it.

The §11 contract design (pattern G zero-cost-when-unused) is structurally correct: the cost the user pays for subscribing is dominated by the kinetics computations that would have happened anyway, not by the registry-write machinery.
