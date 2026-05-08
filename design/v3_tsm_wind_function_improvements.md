# v3 TSM wind-function improvements — design hand-off

**Date:** 2026-05-08
**Audience:** v3 TSM maintainer working in
`/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/`
(separate Claude session).
**Source motivation:** Sep 2008 Santiam-Salem case-study investigation
in this repo isolated a multi-degree cool bias in TSM water
temperature to two compounding causes: (a) silent passthrough of
10 m ASOS wind into a 2 m-calibrated formula, and (b) a non-standard
cubic wind exponent that has no primary literature source.
Authoritative CE-QUAL-W2 source confirms only `c=1` and `c=2` are
supported there; `c=3` is explicitly flagged as undetermined.
Broader investigation of legacy Fortran TSM, CE-QUAL-W2, and QUAL2K
establishes that **`c=1` (linear in wind) is the consensus default**
across these references (Fortran TSM uses `1.0 / 1.0 / 1.0`
placeholders, QUAL2K is linear, W2's Ryan-Harleman alternative is
linear, W2 Edinger supports both `CFW=1` and `CFW=2`). See
`design/edinger_wind_exponent_audit.md` and
`design/tsm_heat_balance_audit.md` in this repo for full evidence.

This note proposes targeted changes to v3 TSM module defaults and
APIs to fix both issues at the modules layer rather than via
case-study workarounds.

## Items to change

### 1. Change `wind_c` default from `3.0` to `1.0`

**File:** `src/clearwater_modules_v3/processes/temperature.py:88–102`
(`Temperature.__init__`).

**Current:**

```python
wind_a: float = 0.3,
wind_b: float = 1.5,
wind_c: float = 3.0,
```

**Proposed:**

```python
wind_a: float = 0.3,
wind_b: float = 1.5,
wind_c: float = 1.0,   # consensus default across legacy Fortran TSM,
                       # QUAL2K, and W2 Ryan-Harleman; W2 Edinger
                       # also supports CFW=1
```

`c=1.0` is the consensus default across the reference family:

* **Legacy Fortran TSM:** placeholder defaults `1.0 / 1.0 / 1.0`
  (the audit-quoted "placeholder values intended to be overridden"
  use `c=1`).
* **QUAL2K:** wind function is linear in wind, `f(W) = vw_a + vw_b · W`.
* **CE-QUAL-W2 Ryan-Harleman alternative** (`heat-exchange.f90:141`):
  linear in wind, `FW = 3.59·DTV^(1/3) + 4.26·W`.
* **CE-QUAL-W2 Edinger** (`heat-exchange.f90:143`): supports both
  `CFW=1` (linear) and `CFW=2` (quadratic). `CFW=1` is the more
  basic / standard case.
* **`c=3.0`** (v3 / v1 inheritance) is not supported by CE-QUAL-W2
  (`heat-exchange.f90:78` comment "CFW not determined for other
  values of CFW") and has no primary literature citation.

Keep `wind_a` and `wind_b` for now to avoid compounding the change;
recommend a follow-up calibration study to revisit those once the
exponent is back in line with the literature.

The 2026-05-05 v3 TSM audit memo (Q4/Q6) should be updated to record
that `c=3.0` was the v1 inheritance with no primary source, and
`c=1.0` is the corrected default per the reference family.

Test impact: tests that pin behavior at `c=3` will need either an
explicit `wind_c=3.0` override (to preserve the prior numerical
result) or a re-baselined expected value. Note that the change is
quantitatively large at typical wind speeds (1–3 m/s) — at W=2 m/s,
`b·Wᶜ` drops from `1.5·8 = 12` (c=3) to `1.5·2 = 3` (c=1), a 4×
reduction in the wind-driven part of `f(W)`.

### 2. Add a validator on `wind_c`

In the `Temperature.__init__` validation block (after the existing
`q_net_depth_ramp_ref` and `dTdt_max_per_hour` checks):

```python
if wind_c not in (1.0, 2.0):
    warnings.warn(
        f"wind_c = {wind_c} is outside the values supported by the "
        f"reference family. CE-QUAL-W2 Edinger explicitly supports "
        f"CFW=1.0 (linear, the consensus default across legacy "
        f"Fortran TSM, QUAL2K, and W2 Ryan-Harleman) and CFW=2.0 "
        f"(quadratic, also valid). W2 heat-exchange.f90:78 flags "
        f"other values as 'CFW not determined'. Coefficient `b` is "
        f"unit-coupled to `c`; using a non-standard exponent without "
        f"re-calibrating `b` will produce unphysical heat fluxes.",
        UserWarning,
    )
```

Reject `c <= 0` and `c > 3` outright (those are physically
indefensible regardless of `b`):

```python
if not (0.0 < wind_c <= 3.0):
    raise ValueError(
        f"wind_c must be in (0.0, 3.0]; got {wind_c!r}"
    )
```

(Allow `c=3.0` as the upper bound to preserve back-compat for any
existing run that explicitly opted in. Future versions may tighten
to `(0.0, 2.0]`.)

### 3. Document the wind reference-height assumption

The `wind_function` and `Temperature.__init__` docstrings currently
describe the form but not the height assumption. Add an explicit
note:

> The wind speed `wind_speed` registered into the v3 registry is
> assumed to be at **2 m above the water surface**, matching the
> Edinger 1974 / CE-QUAL-W2 convention. If the application's
> meteorological source reports wind at a different height (e.g.,
> 10 m for ASOS / METAR / GridMET / NLDAS), the application is
> responsible for converting before registering. Standard log-law
> correction over open water (`z₀ ≈ 0.001 m`):
> `U_2m / U_10m ≈ 0.78`. See `design/edinger_wind_exponent_audit.md`
> in the case-study repo for the CE-QUAL-W2 reference and a
> derivation of the log-law factor.

This is a soft change (documentation only) but clarifies the contract
that the case-study layer must satisfy.

### 4. Optional: surface a `wind_height` parameter and apply log-law internally

Heavier-touch alternative to (3): add a `wind_input_height: float`
parameter (default `2.0`, meaning "registry value is already at 2 m")
to `Temperature.__init__`. When non-2.0, apply the log-law correction
internally before the wind function. Pseudocode:

```python
def wind_function(self, wind_speed, richardson_function):
    if self.wind_input_height != 2.0:
        # log-law from input height to 2 m, water z0 ~ 0.001 m
        wind_speed = wind_speed * (
            np.log(2.0 / self.surface_z0)
            / np.log(self.wind_input_height / self.surface_z0)
        )
    return richardson_function * (
        self.wind_a / 1e6 + (self.wind_b / 1e6) * wind_speed**self.wind_c
    )
```

with `surface_z0` also exposed as a parameter (default `0.001 m`).

**Trade-off:** centralises the height conversion in one place (no
per-application case-study code needed), but moves a responsibility
into the module that has historically lived at the application
layer. The user can decide which posture to take. CE-QUAL-W2 keeps
the height correction out of the wind function (it's done at
`WIND2(I)` precomputation in `w2_4_unix.f90`) and v3 could follow
that pattern by leaving (3) alone.

If you choose (4), the case-study runner's `--wind-height-factor`
flag becomes redundant and should be removed in a follow-up.

### 5. Expose AFW / BFW / CFW (or `wind_a / wind_b / wind_c`) overrides via the `init_from_file` YAML path

The `Temperature.__init__` already accepts these as constructor
kwargs but the YAML config path may not. Confirm and (if missing)
add the keys to `from_config` so users can override per-waterbody
without code edits, mirroring CE-QUAL-W2's per-waterbody convention
(`AFW(JW)`, `BFW(JW)`, `CFW(JW)`).

### 6. Optional: accept a per-cell wind-shelter coefficient

CE-QUAL-W2 multiplies `WIND(JW)` by a per-segment `WSC(I)` shelter
coefficient (`w2_4_unix.f90:480`) to account for vegetation,
topography, riparian canopy reducing wind speed reaching the water
surface. This is meaningful in narrow channels with riparian
forest cover (Salem reach is relevant — Minto Brown is well-shaded).

Implementation: register `wind_shelter_coefficient` as an optional
forcing variable (default 1.0 if not registered). In
`wind_function`, multiply `wind_speed` by `wind_shelter_coefficient`
before the exponentiation. Document that the application is
responsible for setting per-cell shelter values.

This is a v3+ feature, not a defect fix. Defer unless a calibration
study shows it's needed.

## Test plan

1. Update unit tests for `wind_function` to use the new default
   (`c=2.0`) and add a regression test that pins the prior `c=3.0`
   result behind an explicit override.
2. Add a unit test for the validator: assert `c=3.0` warns,
   `c=2.0` does not.
3. Re-run any v2-parity tests; any failures will identify places
   where the old `c=3` was implicitly relied on.

## Out of scope

* Recalibrating `wind_a` and `wind_b`. They're inherited the same
  way `wind_c` was; a focused calibration study against a small set
  of water-temperature observations would be the principled way to
  resolve them. For now, hold them at the v3 defaults and only fix
  the exponent.
* Replacing the Edinger family with Ryan-Harleman as a default. Both
  are in CE-QUAL-W2; the choice is application-dependent. The
  `RH_EVAP(JW)` flag in W2 is the W2-side switch. v3 TSM doesn't
  currently expose this option.

## Trace back to the case-study finding

The case-study evidence is in
`/Users/todd/GitHub/ecohydrology/ClearWater-modules-phase2-ESM-streaming/`
under:

* `design/bc_inflow_continuity_findings.md` — isolates the cool bias
  to TSM.
* `design/tsm_heat_balance_audit.md` — identifies wind-height
  passthrough as the leading candidate (closes ~38% of bias).
* `design/edinger_wind_exponent_audit.md` — CE-QUAL-W2 source check;
  c=3 is not standard; predicts c=2 closes the rest.
* `case_studies/santiam_salem/output/v3_smoke_optC_calibrated_15day_obsT_wind2m_alledges/`
  — bias measurement at -1.65 °C with `c=3` after wind-height
  correction.

Once the v3 TSM defaults are updated, the case-study runner's
`--wind-c` override becomes optional rather than required, and the
Stage F validation should converge on the obs IQR for water
temperature without any case-study-side compensation.
