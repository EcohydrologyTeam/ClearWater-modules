# DOX wind-shelter consumption (closing the reaeration shelter gap)

**Date:** 2026-05-24
**Status:** Design note for a small, surgical fix in `clearwater_modules_v3.processes.dox` (and the supporting `utils/reaeration.kaw_20`). Parallels the existing `wind_shelter_coefficient` consumption in `processes.temperature`.
**Scope:** Make the wind-driven oxygen-reaeration path consume the same per-cell `wind_shelter_coefficient` registry forcing that the TSM heat-budget wind function already consumes, so that sheltering reduces both evaporation and gas transfer consistently.

**Status update — IMPLEMENTED 2026-05-30 (DOX):** `utils.reaeration.kaw_20` gained a `wind_shelter` keyword (default `1.0`), applied to the raw wind before the height rescale (§2/§3.2). `DOX.run` reads the optional `wind_shelter_coefficient` registry forcing (default `1.0`) and threads it through `_change_with_components` → `kaw_20` (§3.1). No re-baseline: the demo registers no `wind_shelter_coefficient`, so the default `1.0` leaves the coupled trajectory bit-identical. Tests: `tests/v3/test_reaeration_wind_shelter.py` (12). **N2 and Carbon (§3.3) are NOT yet wired** — they call the now-shelter-aware `kaw_20` but do not pass `wind_shelter`; extending them is the tracked follow-up for full pipeline consistency.

Companion CWR document:
[`wind_sheltering_design.md`](../../ClearWater-riverine/design/wind_sheltering_design.md)
specifies how ClearWater-riverine constructs and publishes the per-cell `wind_shelter_coefficient` field that this APL change consumes. The CWR side can land independently; this APL change is what closes the consistency gap.

---

## 1. Background — current state

The v3 TSM wind function reads `wind_shelter_coefficient` from the registry and applies it inside `_compute_effective_wind` together with the log-law height correction (`clearwater_modules_v3_tsm_wind_function_specification.md` §3–§4):

```
processes/temperature.py:425-436        # cache wind_shelter_coefficient per substep
processes/temperature.py:1394-1459      # _compute_effective_wind + wind_function
```

The composition is:

```
effective_wind = raw_wind * shelter * height_factor
```

mirroring `w2_4_unix.f90:480` exactly.

The DOX reaeration path **does not consume the shelter forcing**. `processes/dox.py:744-810` reads `wind_speed` from the registry and syncs `wind_input_height` from the Temperature process (`dox.py:478-492`), then calls `utils/reaeration.kaw_20(wind_speed, wind_reaeration_option, wind_input_height)` (`utils/reaeration.py:131-205`). `kaw_20` internally rescales to 10 m via `(10 / wind_input_height) ** 0.143` and evaluates the empirical formulas — but no shelter is applied at any step.

**Symptom.** With everything wired as specified in the CWR side, a sheltered cell sees reduced wind for latent + sensible heat but full wind for gas transfer. Cold-season DO and downstream-of-canopy reaches will compute reaeration consistent with open-water gauge wind even though the heat budget correctly registers the canopy. This is a quiet inconsistency: no error, no warning, just systematically high reaeration in sheltered cells.

This is the only outstanding wind-shelter gap in the v3 pipeline. Closing it makes shelter behave the same way W2 does (W2 applies `WSC(I)` once into `WIND2(I)` and that single sheltered wind drives both heat exchange and reaeration).

---

## 2. Decision

Apply `wind_shelter_coefficient` inside the reaeration path with the **same composition order** as Temperature: `effective = raw × shelter × height_factor`. The shelter must be applied **before** the `(10 / wind_input_height) ** 0.143` height rescale, since the W2 reference (`w2_4_unix.f90:480`) applies the shelter to the raw wind in the same step that rescales height.

**No new public parameter.** The hook is the same registry variable that Temperature already consumes (`wind_shelter_coefficient`). Applications that don't register it see the existing behavior. Applications that do register it (e.g., ClearWater-riverine per the companion CWR design) get sheltered reaeration automatically.

---

## 3. Implementation sketch

### 3.1 `processes/dox.py` — cache shelter per substep

Parallel to `temperature.py:425-436`. Inside `DOX.run()`, after the `wind_speed` registry read (currently at `dox.py:744-758`), add:

```python
# Phase X: read wind_shelter_coefficient from the registry to mirror
# Temperature's _compute_effective_wind composition. Without this,
# sheltered cells get reduced wind for heat exchange but full wind
# for reaeration -- silent inconsistency.
if "wind_shelter_coefficient" in registry:
    wind_shelter = registry.get_at_time("wind_shelter_coefficient", time)
else:
    wind_shelter = 1.0  # back-compat: no shelter registered
```

Pass `wind_shelter` through the `kaw_20` call site (currently `dox.py:810`, `dox.py:930`).

### 3.2 `utils/reaeration.py:kaw_20` — accept and apply shelter

Add a new keyword argument `wind_shelter` (default `1.0` for back-compat), and multiply before the height rescale:

```python
def kaw_20(
    kaw_20_user,
    wind_speed,
    wind_reaeration_option,
    wind_input_height: float = 2.0,
    wind_shelter: ArrayLike | float = 1.0,
):
    ...
    # Composition order matches Temperature._compute_effective_wind
    # (raw * shelter * height_factor) and W2 w2_4_unix.f90:480.
    # Shelter is applied BEFORE the height rescale; if the user
    # is at 10 m, the rescale is a no-op and only shelter matters.
    sheltered_wind = wind_speed * wind_shelter
    if abs(wind_input_height - 10.0) < 1e-12:
        Uw10 = sheltered_wind
    else:
        Uw10 = sheltered_wind * (10.0 / wind_input_height) ** 0.143
    ...
```

The `wind_shelter` argument accepts a scalar `1.0` (back-compat default) or an `xr.DataArray` aligned with `wind_speed`. The xarray broadcast machinery handles per-cell shelter automatically.

### 3.3 N2 and Carbon (if applicable)

N2 (`processes/n2.py`) and Carbon (`processes/carbon.py`) also have a wind-driven reaeration path (see `clearwater_modules_v3_nsm1_audit_c_dox.md` if it exists for cross-references). Audit those files for the same `kaw_20` (or analogous) call sites and apply the identical pattern. Out of scope for the minimum fix, but should be tracked.

---

## 4. Validation

### 4.1 Unit (`tests/v3/test_reaeration_wind_shelter.py`, new)

- `kaw_20(wind_speed=U, wind_shelter=1.0)` produces the same output as the pre-change `kaw_20(wind_speed=U)` — back-compat.
- `kaw_20(wind_speed=U, wind_shelter=0.5)` produces the same output as `kaw_20(wind_speed=0.5 * U, wind_shelter=1.0)` — verifies the shelter is mathematically equivalent to pre-multiplying.
- Per-cell shelter: `kaw_20` with `wind_shelter` as an `xr.DataArray` of shape `(nface,)` broadcasts correctly over a `wind_speed` array of shape `(nface,)`, producing per-cell reaeration velocities.
- For `wind_reaeration_option == 3` (the piecewise formula with the 3.5 m/s threshold), shelter must be applied before the threshold check. Test that a sheltered wind crossing the threshold lands in the correct branch: e.g., `wind_speed = 4.0`, `wind_shelter = 0.5` → effective 2.0 → low-wind branch.

### 4.2 Integration (`tests/v3/test_dox_wind_shelter_integration.py`, new)

- Two-cell DOX run with `wind_shelter_coefficient = [1.0, 0.3]` registered. Reaeration velocity in cell 1 should equal the unsheltered cell 0 velocity evaluated at `0.3 * U`.
- Coupled smoke: run Temperature + DOX together with the same shelter forcing; assert that both `_cached_shelter` (Temperature) and the `kaw_20` call (DOX) see the same per-cell array.

### 4.3 Regression

- The full `tests/v3` suite must remain green. Where existing tests bind `kaw_20` via positional args, no behavior changes (shelter defaults to `1.0`).

---

## 5. Backward compatibility

- `kaw_20` gains a keyword-only argument with a default value of `1.0`; existing callers that don't pass it get unchanged numerical output.
- `DOX.run()` gains an `if "wind_shelter_coefficient" in registry` branch; applications that never register the variable get unchanged numerical output.
- The TSM Temperature path is unchanged — it already consumes the variable.
- No changes to parameter defaults, registry semantics, or the `wind_reaeration_option` enum.

---

## 6. Out of scope

- Adding `wind_shelter` as a *constructor parameter* on `DOX` (paralleling `Temperature(wind_shelter=...)`). The current Temperature constructor parameter exists for the standalone-APL use case where there is no registry; CWR coupling always provides per-cell shelter via the registry, so the DOX-side scalar would only be exercised in standalone DOX runs. Defer pending demand.
- Refactoring `_compute_effective_wind` into a shared utility. Two call sites with identical three-line composition is below the abstraction threshold.
- N2 and Carbon wind-reaeration paths (see §3.3). Track as follow-ups; identical pattern.

---

## 7. References

- `clearwater_modules_v3_tsm_wind_function_specification.md` (2026-05-08) — the spec the Temperature side implements; §4 defines the shelter parameter, §3 the log-law height correction, and the W2 reference at `w2_4_unix.f90:480` for composition order.
- `clearwater_modules_v3/processes/temperature.py:425-436, 1394-1459` — the pattern to mirror.
- `clearwater_modules_v3/processes/dox.py:744-810, 861-945` and `utils/reaeration.py:131-205` — the call sites this fix touches.
- `ClearWater-riverine/design/wind_sheltering_design.md` (2026-05-24) — the consumer-side design that drives the need for this fix.
- CE-QUAL-W2 source `w2_4_unix.f90:480` — `WIND2(I) = WIND(JW) * WSC(I) * log(2/Z0)/log(WINDH/Z0)`.
