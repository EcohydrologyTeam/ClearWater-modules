# v3 TSM thin-water heat-flux skip

**Date:** 2026-05-30
**Status:** Specification + implementation note. Documents the opt-in
`q_net_depth_skip_threshold` parameter added to
`clearwater_modules_v3.processes.temperature.Temperature` on 2026-05-27.
The feature is **disabled by default** (`q_net_depth_skip_threshold = 0.0`)
and preserves byte-identity with prior runs when left at the default.

## Summary

The v3 TSM already carries two thin-water stability regularizations
ported from v1 (implemented in `_temperature_change_with_factors`; see
`design/clearwater_modules_v3_tsm_design_specification.md` and the tests
in `tests/v3/test_tsm_stability_ramp_v3.py` — note those tests cite a
`design/tsm_stability_thin_water.md` memo that is not currently present
in the repo):

1. **Depth ramp** on the net flux:
   `ramp = min(1, depth / q_net_depth_ramp_ref)` (default
   `q_net_depth_ramp_ref = 0.3 m`).
2. **Rate cap** on the per-substep delta T:
   `|dT| ≤ dTdt_max_per_hour · dt_hours` (default `5.0 K/hr`).

This memo adds a third, harder regularization — a **hard depth skip** —
to address a specific artifact observed in the Corvallis–Santiam
September 2008 application: **newly-wet cells** (cells that transition
from dry to wet within a substep, leaving a film of water only a few
centimetres deep) produce spurious temperature spikes from the heat-flux
kinetics before the cell either deepens or dries again.

The new `q_net_depth_skip_threshold` parameter lets an application skip
the heat-flux kinetics **entirely** for cells thinner than a chosen
threshold, so those cells evolve by transport only (their temperature
tracks neighbours through CE-QUAL-W2-Riverine / CWR's transport LHS) and
the kinetics resume once the cell deepens past the threshold.

## Motivation — the newly-wet-cell artifact

The per-substep water-temperature kinetics have the base form

```
dT = q_net · SA · dt / (V · ρ · cp)  =  q_net · dt / (depth · ρ · cp)
```

The `1 / depth` factor makes the kinetics stiff as `depth → 0`: the same
surface flux is absorbed by a vanishing heat capacity, so a thin film
heats or cools far faster than physically meaningful. At `depth < 5 cm`
the surface-area-to-volume ratio is large enough that even moderate solar
(≈ 200 W/m²) drives multi-K/hr temperature deltas **before the dT/dt cap
fires**, and such a cell typically dries up within a substep or two, so
its kinetics contribution is dominated by numerical artifact rather than
real physics.

The depth ramp damps but does **not** eliminate this — at `depth = 0.05 m`
with `q_net_depth_ramp_ref = 0.3 m` the ramp factor is `0.05 / 0.3 ≈ 0.167`,
so 1/6 of an already-stiff flux still drives the spike. For applications
with frequent wetting/drying along a fluctuating shoreline (the
Corvallis–Santiam case), the residual artifact is large enough to bias
shallow-margin temperatures.

## Design

### Parameter

```python
Temperature(..., q_net_depth_skip_threshold: float = 0.0)
```

- **Units:** metres (m).
- **Default `0.0`:** the skip is disabled. The branch
  `if self.q_net_depth_skip_threshold > 0.0` is never entered, so output
  is **byte-identical** to a run without the feature. This is the
  back-compatibility contract — existing applications and the existing
  test baselines are unaffected.
- **Recommended non-zero value `0.05` m (5 cm):** sits below the linear
  ramp's natural shutoff and below any depth at which the kinetics term
  is physically meaningful, per the motivation above. Applications that
  exhibit newly-wet-cell temperature artifacts should set this; others
  should leave it at the default.

### Mechanism

The skip is applied at the end of
`Temperature._temperature_change_with_factors`, after the depth ramp and
the rate cap have produced `delta_clipped` and the per-cell `ramp` factor:

```python
if self.q_net_depth_skip_threshold > 0.0:
    thin_water_skip = depth < self.q_net_depth_skip_threshold
    delta_clipped = xr.where(thin_water_skip, 0.0, delta_clipped)
    if not isinstance(ramp, xr.DataArray):       # promote the disabled-ramp scalar
        ramp = xr.ones_like(depth) * ramp
    ramp = xr.where(thin_water_skip, 0.0, ramp)
```

Two things are zeroed for a skipped (thin) cell:

1. **The water-side kinetics delta** (`delta_clipped → 0`): the cell's
   water temperature does not change from the heat-flux kinetics this
   substep.
2. **The `ramp` factor** (`ramp → 0`): this is what propagates the skip
   to the sediment side. `Temperature.run` scales the sediment-side delta
   by `ramp · clip_ratio` (audit finding F2, 2026-05-05 — the per-cell
   water↔sediment energy pair-cancellation invariant). With `ramp = 0`,
   the sediment-side delta is also zeroed for the skipped cell.

The comparison `depth < threshold` is **strict**, so a cell exactly at
the threshold is **not** skipped. Dry cells (`depth = 0` from the
`surface_area > 0` guard) are flagged as thin and zeroed; this is
harmless and consistent with `Model.__apply_wet_mask`, which overwrites
dry cells with NaN at the orchestration layer regardless.

### Energy conservation

Because **both** the water-side and sediment-side kinetics deltas are
zeroed for a skipped cell, the per-cell water↔sediment heat exchange is
`0 + 0 = 0` — trivially conservative. The skip does not introduce the
one-sided energy sink that audit finding F2 warned against (which arises
when only one of the two reservoirs is scaled). A skipped cell is frozen
on the **kinetics** side only; it still evolves by **transport** (CWR's
advective–dispersive LHS moves heat in and out of the cell), so the cell
is not artificially isolated — it relaxes toward its wet neighbours.

### Diagnostics are NOT zeroed

The per-component flux dictionary (`components`: `q_sensible`, `q_latent`,
`q_longwave_up`, `q_longwave_down`, `q_solar`, `q_sediment`, `q_net`) is
**not** zeroed by the skip. These remain available as a registry-recorded
diagnostic of what the kinetics *would* have applied absent the skip,
which is useful for auditing how much heat the skip suppressed in a given
run. Only the applied `delta` and the `ramp` factor are zeroed.

### Validation

The constructor validates the parameter with the same finite-and-
non-negative predicate used for `q_net_depth_ramp_ref`:

```python
if not (np.isfinite(q_net_depth_skip_threshold)
        and q_net_depth_skip_threshold >= 0.0):
    raise ValueError(...)
```

`0.0` (the disable value) passes; negatives, `+inf`, and `NaN` are
rejected at construction time so a misconfiguration surfaces immediately
rather than silently mis-skipping cells.

## Relationship to the depth ramp and rate cap

The three thin-water regularizations form an escalating stack, applied in
order within `_temperature_change_with_factors`:

| Tier | Parameter | Default | Effect on a thin cell |
| --- | --- | --- | --- |
| 1. Depth ramp | `q_net_depth_ramp_ref` | `0.3 m` (on) | Linearly damps `q_net` by `depth / ref` |
| 2. Rate cap | `dTdt_max_per_hour` | `5.0 K/hr` (on) | Clips `|dT|` to `cap · dt_hours` |
| 3. Hard skip | `q_net_depth_skip_threshold` | `0.0` (**off**) | Zeroes kinetics + ramp below threshold |

The skip is a **generalisation of the linear ramp** to a hard cutoff:
where the ramp asymptotically reduces the flux but never removes it, the
skip removes it entirely below the threshold. The two compose — with the
default ramp on and the skip set to `0.05 m`, a cell at `0.10 m` is ramped
(`0.10 / 0.3 ≈ 0.33`) but not skipped, while a cell at `0.02 m` is skipped
outright. The skip is shipped **off by default** because, unlike the ramp
and cap (which are conservative regularizations that barely touch deep
cells), a hard skip is an application-specific modelling choice.

## Configuration and hotstart

- **Config-driven runs:** `Temperature.from_config` is
  `Temperature(**config)`, so `q_net_depth_skip_threshold` is settable
  directly from a config dict / YAML with no additional plumbing.
- **Hotstart:** the parameter is a constructor-time setting, not
  per-substep state, so it is not part of `to_hotstart` / `from_hotstart`
  (which only persist `__skip_first_time_step`). A hotstarted run must be
  reconstructed with the same constructor argument.

## Testing

Coverage lives in `tests/v3/test_tsm_thin_water_skip_v3.py`:

- **Default off / byte-identity:** with `q_net_depth_skip_threshold = 0.0`,
  even an extremely thin cell is **not** zeroed — output equals the
  ramp+cap result.
- **Skip fires below threshold:** thin cell → `delta == 0` and
  returned `ramp == 0`.
- **Strict boundary:** a cell exactly at the threshold is **not** skipped.
- **Above threshold unaffected:** a cell above the threshold is
  bit-identical with the skip on vs off.
- **Multi-cell mixed depths:** per-cell behaviour on a heterogeneous mesh
  (some cells skipped, others ramped/unguarded).
- **Ramp promotion when ramp disabled:** with `q_net_depth_ramp_ref = 0`
  the helper promotes the scalar `ramp = 1.0` to an array and zeroes it on
  skipped cells.
- **Validation:** constructor rejects negative / `inf` / `NaN`; accepts
  `0.0` and positive finite values.
- **Sediment-side propagation (energy pair-cancellation):** end-to-end
  through `Temperature.run`, a skipped thin cell freezes **both** the
  water and sediment temperatures (both deltas exactly `0`), while a deep
  cell is unaffected.
- **Diagnostics preserved:** `q_net` written to the registry on a skipped
  cell is the unzeroed flux.

## Provenance

- Implementation: `src/clearwater_modules_v3/processes/temperature.py`
  (constructor parameter, validation, and the skip block at the end of
  `_temperature_change_with_factors`), added 2026-05-27.
- Motivation: Corvallis–Santiam September 2008 newly-wet-cell temperature
  investigation.
- Energy-conservation contract (the `ramp · clip_ratio` sediment-side
  scaling the skip relies on): audit finding F2, 2026-05-05
  (`design/clearwater_modules_v3_tsm_audit_2026-05-05.md`).
- Related thin-water regularizations: the depth ramp and rate cap, v1
  port (`design/clearwater_modules_v3_tsm_design_specification.md`,
  `tests/v3/test_tsm_stability_ramp_v3.py`).
