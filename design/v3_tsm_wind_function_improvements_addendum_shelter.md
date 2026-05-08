# v3 TSM wind-function improvements — addendum: wind sheltering and wind input height

**Date:** 2026-05-08
**Status:** Addendum to `v3_tsm_wind_function_improvements.md` in this
repo. Items 1, 2, 3, and 5 of the original hand-off note (wind_c
default 3.0 → 1.0, wind_c validator, 2 m wind reference height
docstring, YAML config exposure) have been implemented. **Items 4
(`wind_input_height`) and 6 (wind shelter) are not yet implemented.**
This addendum revisits item 6 with the authoritative CE-QUAL-W2
source evidence and elevates it from "optional v3+ feature" to a
recommended near-term enhancement, then proposes a refined design.

## Why revisit wind shelter now

The Sep 2008 Santiam-Salem case-study sensitivity sweep produced:

| Configuration | Salem-cell T bias | RMSE |
|---|---|---|
| `c = 3` (old default), no fixes | -2.68 °C | 2.68 |
| `c = 3` + `--wind-height-factor 0.78` + `all_edges` | -1.65 °C | 1.67 |
| **`c = 2` + same wind-height + `all_edges`** | **-0.28 °C** | **0.35** |
| `c = 1` + same wind-height + `all_edges` (queued) | TBD | TBD |

The `c = 2` run essentially calibrates the case study, but with two
inheritances the principled formulation should not have:

1. **`c = 2` is not the literature default** — `c = 1` is, per
   the analysis in `v3_tsm_wind_function_improvements.md`. The new
   v3 default is correctly 1.0.

2. **The `--wind-height-factor 0.78` covers only the over-water
   log-law height correction** — it ignores both:
   * the airfield-to-water roughness change at the land/water
     boundary (about a further factor of 0.85–0.90), and
   * canopy/topographic sheltering of the channel from open-fetch
     wind (CE-QUAL-W2's `WSC(I)`).

A combined effective wind-reduction factor for typical
ASOS-over-airfield → 2 m over a wide vegetated river runs roughly:

| Component | Factor | Cumulative |
|---|---|---|
| Height: 10 m → 2 m, log-law over water (`z₀ = 0.001 m`) | 0.825 | 0.825 |
| × airfield-to-water roughness change | ~0.875 | ~0.722 |
| × open river WSC | ~0.85 | ~0.614 |
| × canopied channel WSC | ~0.65 | ~0.469 |

The case-study runner's `--wind-height-factor 0.78` captures only the
first row. The implicit absorption of the missing factors by `c = 2`
(rather than `c = 1`) is what made the `c = 2` configuration land
numerically calibrated; with `c = 1` the principled answer requires
explicit shelter handling.

## CE-QUAL-W2 wind shelter — full reference

Repo: `/Users/todd/GitHub/CE-QUAL-W2-ERDC/CE-QUAL-W2-ERDC-dev/src/W2_v2026.02/`

### Variable

`WSC(I)`, dimension `(IMX,)`, where `IMX` is the total mesh segment
count. Per-segment scalar; **time-varying**.

### Where it's read

`time-varying-data.f90`:
* Lines 22, 34, 98: `WSCNX` (next-time-step buffer), `NXWSC` (next-
  time-stamp), allocation.
* Lines 185–196: open the wind-shelter file `WSCFN` (logical unit
  `WSH`) and read the first two records — typical W2 time-varying
  file pattern with header consumption.
* Lines 1512–1517: in the time-step advance loop, when
  `JDAY ≥ NXWSC`, copy `WSCNX` into the active `WSC` array and read
  the next record.

So `WSC(I)` is a **per-segment, file-driven, time-varying** array.
Common file content patterns: monthly values for seasonal vegetation,
or a single time-invariant value per segment if shelter is treated as
geometry-only.

### Where it's applied

`w2_4_unix.f90:480, 487` — primary application:

```fortran
WIND2(I) = WIND(JW) * WSC(I) &
         * DLOG(2.0D0 / Z0(MetRegWB(JW))) / DLOG(WINDH(MetRegWB(JW)) / Z0(MetRegWB(JW)))
```

Composition order:
1. `WIND(JW)` — raw wind from input file (any height, user-specified
   `WINDH(JW)`).
2. `× WSC(I)` — per-segment shelter, scales the wind magnitude.
3. `× log_law` — height correction from `WINDH(JW)` to 2 m.

The order matters mathematically only if log_law depends on wind
magnitude (it doesn't, so the operations commute), but the W2 form
puts shelter first because it's a wind-magnitude reduction and log_law
is purely a height transform.

`w2_4_unix.f90:661, 663` — also applies `WSC(I)` to a separate `WIND10(I)`
computed at 10 m for use in models that want 10 m wind (e.g., the
diagenesis sub-model, lines 1802 of
`Diagenesis Sediment Flux Model 05.f90`).

### Typical values (W2 manual)

* `1.0` — no sheltering (default if no `WSCFN` file is supplied; some
  W2 setups omit the shelter file and rely on this).
* `0.85–1.0` — open lake or wide river without significant canopy
  influence.
* `0.5–0.85` — narrow channel, partial riparian canopy, mild
  topographic shelter.
* `0.3–0.5` — heavily shaded backwater, narrow steep canyon.

Calibration is typically done against a measured ratio between the
station anemometer and an in-channel observation. Where in-channel
wind isn't available, a digital surface model (canopy heights from
LiDAR, channel width) feeds engineering judgement — there's no
closed-form derivation analogous to the log-law height correction.

## Current state of v3 TSM (after the recent default/validator
implementation)

`src/clearwater_modules_v3/processes/temperature.py`:

* `wind_c` default = `1.0` (linear, consensus default). ✓
* `wind_c` validator: warns outside `{1.0, 2.0}`, rejects outside
  `(0.0, 3.0]`. ✓
* Docstring documents 2 m above water as the wind reference height
  convention; says the application is responsible for height
  conversion. ✓

What's still missing relative to the W2 reference and the original
hand-off note:

* **No mechanism for wind sheltering** in the TSM module itself or
  in the registry forcing schema. The `wind_speed` registry variable
  is read as-is (line 257 in `temperature.py:run`).
* **No `wind_input_height` parameter** for in-module log-law height
  correction — applications must pre-correct.
* **No `wind_shelter_coefficient`** as a standard forcing variable
  the user can register per cell.

## Recommendations for v3 TSM

Three options ordered by intrusiveness. Items A and B together would
match W2's first-class-feature posture for shelter; item C is the
fully-equivalent file-driven variant.

### Option A — Scalar shelter parameter on the constructor (lightest)

Add `wind_shelter: float = 1.0` to `Temperature.__init__`. In
`wind_function`, multiply `wind_speed` by `self.wind_shelter` before
the exponentiation:

```python
def wind_function(self, wind_speed, richardson_function):
    sheltered_wind = wind_speed * self.wind_shelter
    return richardson_function * (
        self.wind_a / 1e6
        + (self.wind_b / 1e6) * sheltered_wind ** self.wind_c
    )
```

* **Pro:** trivial change; no new registry variable; no breaking
  change (`1.0` is a no-op).
* **Con:** scalar — same shelter applied to every cell. Loses W2's
  per-segment granularity. Reasonable for single-reach simulations
  where shelter is roughly homogeneous within the simulation domain.

### Option B — Per-cell `wind_shelter_coefficient` registry forcing (recommended)

Register a new optional forcing variable
`wind_shelter_coefficient` (default 1.0 if not registered). In
`wind_function`, look it up alongside `wind_speed`:

```python
shelter = registry.get_at_time("wind_shelter_coefficient", time)
sheltered_wind = wind_speed * shelter
```

* **Pro:** matches W2's per-segment posture (per cell here). The
  application can supply a 1D or 2D field of shelter values from a
  canopy-height raster or hand-tuned per-segment table.
* **Con:** a new registry variable; downstream graphs/tests need to
  tolerate the optional presence. The pattern is already established
  in v3 (e.g., several optional registry forcings in NSM1 modules).
* **Default:** 1.0 (no shelter) — same back-compat posture as item A.

### Option C — Time-varying file-driven shelter (W2 parity)

If the application needs seasonal vegetation effects (leaf-on vs.
leaf-off), the registry value can simply be set at each TSM substep
to a different number — that's already supported by the
`registry.set_at_time` API. No code change in TSM is required beyond
option B; the time-variation lives in the application layer's
forcing code, not in TSM. So **option B subsumes option C** in the
v3 architecture; mention both for completeness but recommend B as
the implementation path.

### `wind_input_height` parameter (independent of shelter)

The original item 4 in `v3_tsm_wind_function_improvements.md`
proposed surfacing a `wind_input_height` parameter and applying the
log-law internally. That hasn't been implemented; the current state
expects the application to pre-correct.

**Two reasons to add `wind_input_height` after all:**

1. **Symmetry with W2.** W2 takes a raw `WIND(JW)` at user-specified
   `WINDH(JW)` and corrects internally. v3's "responsibility lives
   at the application layer" posture is documented but not
   conventional in the broader heat-balance modeling community.

2. **Couples cleanly with shelter.** If options A or B above are
   adopted, the wind function's input transform pipeline becomes
   `raw_wind → × shelter → log_law(input_height → 2 m)` —
   a clean composition that mirrors W2's `WIND2(I)` line.

### Concrete proposed signatures

```python
class Temperature(Process):
    def __init__(
        self,
        wind_a: float = 0.3,
        wind_b: float = 1.5,
        wind_c: float = 1.0,
        # NEW — height handling moved into the module:
        wind_input_height: float = 2.0,   # m above water (or land)
        wind_input_z0: float = 0.001,     # m, water surface roughness
        # NEW — scalar shelter (option A); set to 1.0 as no-op default:
        wind_shelter: float = 1.0,
        ...
    ) -> None:
        ...
```

In `wind_function`, the input transform becomes:

```python
def wind_function(self, wind_speed, richardson_function):
    # Per-cell shelter (from registry, optional). Falls back to the
    # constructor scalar `wind_shelter` if no registry variable is set
    # (or to 1.0 if both are absent).
    shelter_per_cell = self._maybe_get_shelter_from_registry(...)
    shelter = shelter_per_cell if shelter_per_cell is not None else self.wind_shelter

    # Height correction from input height to 2 m (log-law over water).
    # No-op if wind_input_height == 2.0.
    if self.wind_input_height != 2.0:
        height_factor = (
            np.log(2.0 / self.wind_input_z0)
            / np.log(self.wind_input_height / self.wind_input_z0)
        )
    else:
        height_factor = 1.0

    effective_wind = wind_speed * shelter * height_factor
    return richardson_function * (
        self.wind_a / 1e6
        + (self.wind_b / 1e6) * effective_wind ** self.wind_c
    )
```

This composition mirrors W2's `WIND2(I)` line exactly:

| W2 (`w2_4_unix.f90:480`) | v3 TSM (proposed) |
|---|---|
| `WIND(JW)` | `wind_speed` (registry) |
| `× WSC(I)` | `× shelter` (per-cell registry or scalar) |
| `× DLOG(2/Z0)/DLOG(WINDH/Z0)` | `× height_factor` (constructor params) |
| `= WIND2(I)` | `= effective_wind` |

The `effective_wind` then enters the Edinger formula at `c=1` (the
new v3 default).

### Test plan

1. Unit test: scalar `wind_shelter < 1.0` produces a smaller `f(W)`
   than the no-shelter case; mass-balance and Edinger structure
   unchanged.
2. Unit test: `wind_input_height = 10.0` over `wind_input_z0 = 0.001`
   produces effective wind = `0.825 × wind_speed`, matching the
   log-law factor used in CE-QUAL-W2.
3. Unit test: registry `wind_shelter_coefficient` (option B) overrides
   constructor scalar when both are set.
4. Regression test: pre-existing tests with `wind_shelter` and
   `wind_input_height` left at defaults reproduce the prior numerical
   output to machine precision.

## Recommendation for the Santiam-Salem case study

Independent of the v3 TSM enhancements above, the case-study runner
(`08_run_coupled_v3_smoke.py` in the case-study repo) should grow a
`--wind-shelter-factor` flag analogous to the existing
`--wind-height-factor`, applied at met-data load time. That allows
the calibration sweep to use literature-default `c = 1` with a
combined `(height × shelter)` reduction factor that absorbs the
missing roughness/canopy components, rather than the empirical
`c = 2` workaround the case study currently relies on. Once the v3
TSM features above are in, `--wind-shelter-factor` becomes redundant
(and the case-study flag can be retired in favor of the in-module
parameter).

## Cross-references

* `v3_tsm_wind_function_improvements.md` (this repo, original
  hand-off; items 1/2/3/5 implemented, items 4/6 deferred).
* `edinger_wind_exponent_audit.md` (this repo, CE-QUAL-W2 source
  check that established `c=1` as the consensus default).
* `clearwater_modules_v3_tsm_audit_2026-05-05.md` (this repo, the
  earlier audit; Q4/Q6 documented the v1 inheritance of `c=3`).
* `/Users/todd/GitHub/ecohydrology/ClearWater-modules-phase2-ESM-streaming/design/`
  (case-study repo): `tsm_heat_balance_audit.md` and
  `bc_inflow_continuity_findings.md` for the investigation chain
  that surfaced the wind-height and wind-exponent issues.
* `/Users/todd/GitHub/CE-QUAL-W2-ERDC/CE-QUAL-W2-ERDC-dev/src/W2_v2026.02/`:
  `heat-exchange.f90`, `w2_4_unix.f90`, `time-varying-data.f90`
  for the W2 wind shelter reference implementation.
