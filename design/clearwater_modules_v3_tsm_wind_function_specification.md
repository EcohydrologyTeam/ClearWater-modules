# v3 TSM wind-function specification

**Date:** 2026-05-08
**Status:** Specification document. Replaces and supersedes the
prior wind-function design memos in this folder (see Provenance
at the end). Describes the target state for the v3 TSM Edinger-form
wind function in `clearwater_modules_v3.processes.temperature`.

## Summary

The v3 TSM wind function has the form
`f(W) = (a + b · Wᶜ) / 1e6 × Ri(stability)` and is used in the
Edinger family latent and sensible heat flux computations. This
specification:

1. Adopts `wind_c = 2.0` as the v3 default, matching the explicit
   default in the CE-QUAL-W2 user manual and the consensus across
   the Edinger family of parameterisations.
2. Documents 2 m above the water surface as the wind-speed reference
   height required by the Edinger 1974 calibration.
3. Surfaces an optional `wind_input_height` parameter that lets v3
   apply a log-law height correction internally when the application
   passes wind from a different reference height (e.g. 10 m ASOS).
4. Adds an optional `wind_shelter` parameter (scalar, with a future
   per-cell registry-forcing extension) to represent canopy and
   topographic wind reduction, mirroring CE-QUAL-W2's `WSC(I)`.
5. Validates `wind_c` with a warn-on-non-standard / reject-on-out-of-
   range pattern.
6. Holds `wind_a` and `wind_b` at the v3 inheritance values
   (`0.3 / 1.5`) for now and tags a follow-up calibration study as
   future work. Note the v3 normalisation `/1e6` is part of the
   coefficient convention — the v3 numerical values are not directly
   interchangeable with CE-QUAL-W2's SI calibration values
   (`9.2 / 0.46`).

## Authoritative reference

### CE-QUAL-W2

**Manual** (`cequalw2-docs-dev/docs/part3/02_input_and_output_files_data_description.md:512–514`):

| Coefficient | Value | Units |
|---|---|---|
| `AFW` | **9.2** | W m⁻² mmHg⁻¹ |
| `BFW` | **0.46** | W m⁻² mmHg⁻¹ (m s⁻¹)⁻ᶜᶠʷ |
| `CFW` | **2.0** | dimensionless |

> *"The generalized wind function is: `f(Wz) = AFW + BFW · Wz^CFW`.
> The default formulation is from Edinger et al. (1974). For
> thermally loaded systems such as cooling lakes, use the
> Ryan-Harleman formulation with `[RHEVAP]` set to `ON`."*

**Source** (`heat-exchange.f90:143`):

```fortran
FW = AFW(JW) + BFW(JW) * WIND2(I)**CFW(JW)
```

with British-to-SI unit conversion factors `BCONV` defined only
for `CFW = 1.0` and `CFW = 2.0` (`heat-exchange.f90:76–78`); the
trailing comment `! SW Issues: CFW not determined for other values
of CFW` documents that values outside `{1, 2}` are not part of the
designed range.

**Wind reference height** (`w2_4_unix.f90:480`):

```fortran
WIND2(I) = WIND(JW) * WSC(I) &
         * DLOG(2.0D0 / Z0(JW)) / DLOG(WINDH(JW) / Z0(JW))
```

`WIND(JW)` is the raw wind at user-specified `WINDH(JW)`; `WSC(I)`
is the per-segment time-varying wind-shelter coefficient; the log-
law converts to 2 m above water with surface roughness `Z0(JW)`.

**Empirical uniformity across W2 example cases.** All seven example
case studies in `CE-QUAL-W2/examples/` use the same defaults:

| Case | Type | AFW | BFW | CFW | WINDH (m) | WSC |
|---|---|---|---|---|---|---|
| Spokane River | River | 9.2 | 0.46 | 2.0 | 2 | per-segment file |
| Columbia Slough | Estuary | 9.2 | 0.46 | 2.0 | 10 | 0.80 |
| Bonneville Dam | Tailrace | 9.2 | 0.46 | 2.0 | 2 | — |
| Detroit Reservoir | Reservoir | 9.2 | 0.46 | 2.0 | 6 | 1.00 |
| DeGray Reservoir | Reservoir | 9.2 | 0.46 | 2.0 | 10 | 0.90 |
| Long Lake | Reservoir | 9.2 | 0.46 | 2.0 | 2 | — |
| Berlin Milton | River/lake | 9.2 | 0.46 | 2.0 | 10 | — |

No river-vs-reservoir differentiation in `AFW / BFW / CFW`. River-
specific tuning is done via `WSC(I)` (wind shelter coefficient),
not via different wind-function coefficients.

### QUAL2K

**Default formulation:** Brady, Graves, Geyer (1969), same Edinger
family.

```
f(Uw) = 19.0 + 0.95 · Uw²    [m d⁻¹ mmHg⁻¹]
```

with `Uw` at **7 m** above water. QUAL2K internally adjusts wind
speeds taken at other heights using the TVA 1972 exponential law.

**Coefficient convention differs from W2:** different magnitude
(`19.0 / 0.95` vs `9.2 / 0.46`) reflecting different unit systems
(`m d⁻¹ mmHg⁻¹` vs `W m⁻² mmHg⁻¹`) and different calibration
datasets (BGG cooling-pond vs Edinger 1974). **Exponent is
identical: `c = 2.0`.**

### Consensus

Both W2 (Edinger default) and QUAL2K (BGG default) use `c = 2.0`
quadratic. The Edinger-family wind function is universally
quadratic in the modern reservoir/river WQ modelling literature.
`c = 1.0` linear is reserved for some derivative parameterisations
(W2's Ryan-Harleman, intended for "thermally loaded systems such
as cooling lakes") and is not the natural-waterbody default.

## Current v3 state

`src/clearwater_modules_v3/processes/temperature.py:88-103`:

```python
def __init__(
    self,
    wind_a: float = 0.3,
    wind_b: float = 1.5,
    wind_c: float = 1.0,    # PROPOSED CHANGE: 1.0 → 2.0
    sediment_density: ArrayLike = 1600.0,
    sediment_specific_heat: float = 1673.0,
    air_diffusivity_ratio: float = 1.0,
    sediment_diffusivity: float = 0.0432,
    time_step: timedelta = timedelta(minutes=5),
    use_sediment_temperature: bool = True,
    evolve_sediment_temperature: bool = True,
    q_net_depth_ramp_ref: float = 0.3,
    dTdt_max_per_hour: float = 5.0,
) -> None:
```

The wind function:

```python
def wind_function(self, wind_speed, richardson_function):
    return richardson_function * (
        (self.wind_a / 1_000_000)
        + (self.wind_b / 1_000_000) * (wind_speed**self.wind_c)
    )
```

A validator checks `wind_c ∈ (0.0, 3.0]` (raises `ValueError`
outside) and warns when `wind_c not in {1.0, 2.0}`.

## Specification

### 1. Default `wind_c = 2.0`

**File:** `src/clearwater_modules_v3/processes/temperature.py`,
`Temperature.__init__`.

```python
wind_c: float = 2.0,    # CE-QUAL-W2 default; consensus across the
                        # Edinger family (W2 Edinger, QUAL2K
                        # Brady-Graves-Geyer)
```

### 2. Wind-function-coefficient docstring

Replace the `wind_a, wind_b, wind_c` block in
`Temperature.__init__`'s docstring with:

```
wind_a, wind_b, wind_c: Wind-function parameters in the
    Edinger, Brady & Geyer (1974) form
    ``f(W) = (a + b * W^c) / 1e6`` (multiplied internally
    by the Richardson stability function in
    :py:meth:`wind_function`). Defaults are
    ``a = 0.3, b = 1.5, c = 2.0``.

    The exponent ``c = 2.0`` (quadratic in wind) matches the
    explicit default in the CE-QUAL-W2 user manual
    (``AFW / BFW / CFW = 9.2 / 0.46 / 2.0``) and the
    QUAL2K Brady-Graves-Geyer default
    (``19.0 / 0.95 / 2.0``, different unit system, same
    exponent). All seven example case studies shipped with
    CE-QUAL-W2 — including Spokane River, Columbia Slough
    Estuary, Detroit, DeGray, Long Lake, Bonneville Dam,
    and Berlin Milton — use ``CFW = 2.0`` regardless of
    waterbody type. River-specific tuning in W2 is done
    via the per-segment wind shelter coefficient
    (``WSC(I)``), not via different wind-function
    exponents.

    The v3 magnitude coefficients ``a = 0.3, b = 1.5`` are
    inherited from v1's ``clearwater_modules.tsm.constants``
    and use the v3 ``/1e6`` normalisation. They are **not**
    directly interchangeable with W2's SI values (``9.2 /
    0.46``) — the units bake the normalisation into the
    coefficients. A focused calibration study to revisit
    ``a`` and ``b`` against observation is tracked as future
    work.

    ``wind_c`` is validated: values not in ``{1.0, 2.0}``
    emit a ``UserWarning``, and values outside ``(0.0, 3.0]``
    raise ``ValueError``. ``c = 3.0`` is allowed at the
    upper bound for back-compat with explicit opt-ins from
    runs that have already been calibrated against the
    prior v3 default.

    Pass any subset to override per-instance; YAML configs
    may also override via ``wind_a / wind_b / wind_c`` keys
    at ``init_from_file`` time.

    References:
    * Edinger, J.E., D.K. Brady, and J.C. Geyer (1974),
      *Heat exchange and transport in the environment*,
      Report 14, Cooling Water Discharge Research Project
      (RP-49), Electric Power Research Institute, Palo
      Alto, CA, 125 pp.
    * Brady, D.K., W.L. Graves, and J.C. Geyer (1969),
      *Surface heat exchange at power plant cooling
      lakes*, Cooling Water Discharge Research Project
      Report 5, Edison Electric Institute, New York.
    * CE-QUAL-W2 User Manual, AFW/BFW/CFW default entries.
    * Chapra, S.C. (2008), *QUAL2K User Manual*,
      §4.1.4 wind function f(Uw).
```

### 3. Wind reference height

Add a `wind_input_height` constructor parameter (and `surface_z0`
for the log-law correction):

```python
def __init__(
    self,
    ...
    wind_input_height: float = 2.0,    # m above water; matches
                                       # Edinger 1974 / W2 default
                                       # convention
    surface_z0: float = 0.001,         # m, water surface roughness
                                       # for log-law correction
    ...
) -> None:
```

Apply inside `wind_function`:

```python
def wind_function(self, wind_speed, richardson_function):
    # If the application registers wind at a non-2-m height, apply
    # the log-law correction internally. No-op when the input is
    # already at 2 m.
    if self.wind_input_height != 2.0:
        height_factor = (
            np.log(2.0 / self.surface_z0)
            / np.log(self.wind_input_height / self.surface_z0)
        )
        wind_speed = wind_speed * height_factor

    return richardson_function * (
        self.wind_a / 1e6
        + (self.wind_b / 1e6) * wind_speed ** self.wind_c
    )
```

Docstring entry:

```
wind_input_height: Height (m) above the water surface at
    which the application's ``wind_speed`` registry
    variable was measured. Default 2.0 (no correction),
    matching the Edinger 1974 / CE-QUAL-W2 convention. Set
    to 10.0 when registering raw ASOS / METAR / GridMET /
    NLDAS wind without external pre-correction. The
    log-law conversion uses ``surface_z0`` (default
    0.001 m, typical for open water) per ``U_2 / U_z = ln(2 / z0)
    / ln(z / z0)``. For ASOS over a flat airfield (i.e.,
    not over the water surface itself) the user should
    additionally apply an airfield-to-water roughness
    correction at the application layer; the in-module
    log-law assumes the wind was measured above the same
    surface as the receiving water.

surface_z0: Roughness length (m) of the water surface, used
    by the in-module log-law height correction when
    ``wind_input_height != 2.0``. Default 0.001 m. Larger
    values (e.g., 0.003 m) are sometimes used for
    smoother-water regimes; sensitivity is small.
```

### 4. Wind sheltering

Add an optional scalar `wind_shelter` parameter to the constructor:

```python
def __init__(
    self,
    ...
    wind_shelter: float = 1.0,    # 1.0 = no shelter; values < 1
                                  # represent canopy / topographic
                                  # reduction
    ...
) -> None:
```

Apply inside `wind_function`, before the height correction:

```python
def wind_function(self, wind_speed, richardson_function):
    # 1. Per-cell wind shelter from the registry, if registered;
    #    falls back to the constructor scalar.
    shelter_per_cell = self._maybe_get_shelter_from_registry(...)
    shelter = (
        shelter_per_cell
        if shelter_per_cell is not None
        else self.wind_shelter
    )

    # 2. Height correction (no-op when wind_input_height == 2.0).
    if self.wind_input_height != 2.0:
        height_factor = (
            np.log(2.0 / self.surface_z0)
            / np.log(self.wind_input_height / self.surface_z0)
        )
    else:
        height_factor = 1.0

    # 3. Apply both factors before the Edinger formula.
    effective_wind = wind_speed * shelter * height_factor

    return richardson_function * (
        self.wind_a / 1e6
        + (self.wind_b / 1e6) * effective_wind ** self.wind_c
    )
```

The composition `raw_wind × shelter × height_factor` mirrors the
CE-QUAL-W2 `WIND2(I) = WIND × WSC × log_law` line in
`w2_4_unix.f90:480` exactly.

Optional follow-on: register `wind_shelter_coefficient` as an
optional registry forcing variable; the helper
`_maybe_get_shelter_from_registry` checks for its presence and
returns `None` when absent (so the constructor scalar applies).
This matches W2's per-segment posture for `WSC(I)` while remaining
back-compatible with applications that don't supply per-cell
shelter values.

Docstring entry:

```
wind_shelter: Scalar wind shelter coefficient applied to
    ``wind_speed`` before the Edinger formula and before
    any height correction. Default 1.0 (no shelter). Values
    < 1 represent canopy or topographic wind reduction;
    typical W2 values are 1.0 (open lake), 0.85–0.90 (open
    river / reservoir), 0.5–0.85 (narrow channel with
    riparian canopy), 0.3–0.5 (heavily shaded backwater).
    Mirrors CE-QUAL-W2's ``WSC(I)`` per-segment shelter
    coefficient (see ``w2_4_unix.f90:480``). For per-cell
    shelter, register ``wind_shelter_coefficient`` as an
    optional forcing variable; the registry value overrides
    this constructor scalar when present.
```

### 5. Validator update

Replace the warning message text to reflect the corrected
consensus:

```python
if wind_c not in (1.0, 2.0):
    warnings.warn(
        f"wind_c = {wind_c} is outside the values supported by "
        f"the Edinger family of wind-function parameterisations. "
        f"CE-QUAL-W2 explicitly defaults to CFW = 2.0 and supports "
        f"CFW = 1.0; the v3 default is c = 2.0. QUAL2K's "
        f"Brady-Graves-Geyer default is also c = 2.0. Other "
        f"values are flagged as 'CFW not determined' "
        f"(W2 heat-exchange.f90:78). Coefficient `b` is "
        f"unit-coupled to `c`; using a non-standard exponent "
        f"without re-calibrating `b` will produce unphysical "
        f"heat fluxes.",
        UserWarning,
        stacklevel=2,
    )
```

The reject-on-out-of-range guard `if not (0.0 < wind_c <= 3.0)`
stays as is — it preserves explicit opt-ins to the prior `c = 3`
default for back-compat while rejecting unphysical values.

### 6. YAML config exposure

Confirm `wind_a / wind_b / wind_c / wind_input_height /
surface_z0 / wind_shelter` are accepted by the
`Temperature.from_config` path (already exists for `wind_a / wind_b
/ wind_c`; extend for the three new parameters).

## Test plan

1. **Default-change regression.** Verify that `Temperature()` with
   no kwargs constructs with `wind_c = 2.0`. Update any test that
   pinned the prior `c = 3.0` numerical result behind an explicit
   `Temperature(wind_c=3.0)` override; the test name should record
   that `c = 3.0` is a back-compat case, not the recommended path.

2. **Validator behavior.**
   - `Temperature(wind_c=1.0)` and `Temperature(wind_c=2.0)`
     construct without warning.
   - `Temperature(wind_c=1.5)` constructs and emits a single
     `UserWarning` whose message references the W2 / QUAL2K
     consensus.
   - `Temperature(wind_c=3.0)` constructs and emits the warning
     (it's outside `{1.0, 2.0}` but inside the legal upper
     bound).
   - `Temperature(wind_c=0.0)`, `Temperature(wind_c=-1.0)`,
     and `Temperature(wind_c=3.5)` raise `ValueError`.

3. **Wind-input-height correction.**
   - `Temperature(wind_input_height=2.0)` produces `f(W)` identical
     to the no-correction case (no-op).
   - `Temperature(wind_input_height=10.0, surface_z0=0.001)`
     produces effective wind = `0.825 × wind_speed`, matching
     the log-law factor used in CE-QUAL-W2.
   - Sensitivity test: at `W = 2 m/s`, verify that
     `wind_input_height = 10.0` reduces the wind term `b × Wᶜ` by
     the expected factor (`0.825² ≈ 0.68` for `c = 2`).

4. **Wind shelter.**
   - `Temperature(wind_shelter=1.0)` produces `f(W)` identical
     to the no-shelter case.
   - `Temperature(wind_shelter=0.5)` reduces the wind term `b × Wᶜ`
     by a factor of `0.5² = 0.25` for `c = 2`.
   - When a `wind_shelter_coefficient` registry forcing is set,
     it overrides the constructor scalar. When absent, the scalar
     applies. (Once the optional registry-forcing path is
     implemented.)

5. **Composition.** Combined `wind_input_height = 10`,
   `wind_shelter = 0.65`, `c = 2` produces the same `f(W)` as a
   single multiplicative factor `0.825 × 0.65 = 0.536` applied to
   `wind_speed` and the result squared. Verify analytically.

6. **Mass / energy conservation regression.** Existing F2
   sediment-water energy-conservation tests, MMS-style flux
   reconstruction tests, and any existing regression tests should
   continue to pass under the new default. Differences in absolute
   numerical output are expected and are the intentional product
   of the default change.

## Out of scope

- **Recalibrating `wind_a` and `wind_b`.** They keep their v3
  inheritance values (`0.3 / 1.5`) for now. A focused calibration
  study against a small set of water-temperature observations is
  the principled way to revisit them. The v3 magnitude convention
  (with the `/1e6` normalisation) is not directly comparable to
  CE-QUAL-W2's SI convention, so a literal substitution of
  W2's `9.2 / 0.46` would not be correct without unit-conversion
  work that's out of scope here.

- **Replacing Edinger as the default.** CE-QUAL-W2 reserves
  Ryan-Harleman for "thermally loaded systems such as cooling
  lakes." The Edinger family is the natural-waterbody default
  there and remains the v3 default.

- **File-driven time-varying shelter.** v3's `registry.set_at_time`
  API already supports time-varying values; an application can
  drive a per-cell shelter value at any cadence by calling that
  method each TSM substep. A separate file-driven mechanism
  (analog of W2's `WSCFN` file) is not needed at the module
  level.

- **Per-segment shelter file format.** Out of scope for this
  spec; the per-cell registry forcing is the mechanism.
  Applications that want a file-driven workflow can write a
  small loader at the application layer.

## Provenance

This specification supersedes two prior design memos in this
folder:

- `edinger_wind_exponent_audit.md` (2026-05-08).
- `v3_tsm_wind_function_improvements.md` (2026-05-08).

Both were written based on a reading of CE-QUAL-W2 source that
weighted the Fortran placeholder defaults (`1.0 / 1.0 / 1.0`) and
the W2 Ryan-Harleman alternative (linear in wind, cited but
intended for "thermally loaded systems") more heavily than the
explicit W2 manual default `CFW = 2.0`. The corrected reading —
the W2 user-manual entry, the uniform `9.2 / 0.46 / 2.0` across
all seven W2 example cases, and the QUAL2K Brady-Graves-Geyer
default also at `c = 2.0` — is in §"Authoritative reference"
above.

The case-study investigation that surfaced the underlying
wind-function questions is documented in
`/Users/todd/GitHub/ecohydrology/ClearWater-modules-phase2-ESM-streaming/design/`:

- `bc_inflow_continuity_findings.md` — boundary-condition
  mass-flux mechanism in CW-Riverine; explains why early
  attempts to drive the case study with observation-based
  BCs produced no measurable change at the validation point
  until the `continuity_correction='all_edges'` mode was
  enabled.
- `tsm_heat_balance_audit.md` — case-study audit that
  identified the wind-input height as the leading
  contributor to the cool-T bias once BCs were verified to
  be flowing.
- `sep_2008_temperature_drift_findings.md`,
  `sep_2008_observed_bc_plan.md` — earlier investigation
  documents that frame the broader case-study workflow.

A 3-point sensitivity sweep at the case-study level
(`wind_c ∈ {3.0, 2.0, 1.0}` with the wind-height correction
factor 0.78 applied at the case-study layer for the still-pending
v3 `wind_input_height` parameter) measured Salem-cell
water-temperature residuals at -1.65, -0.28, and +0.70 °C
respectively; `c = 2.0` lands within the year-to-year
climatological variability for the validation site, consistent
with the literature default.

## File index

Recommended location for the changes described in this spec:

- `src/clearwater_modules_v3/processes/temperature.py`
  — `Temperature.__init__` signature, validator block, and
  `wind_function` body.
- `src/clearwater_modules_v3/processes/temperature.py`
  — docstring revisions.
- `tests/processes/test_temperature.py` (or equivalent)
  — new validator tests, regression-baseline updates,
  composition tests.
- `src/clearwater_modules_v3/processes/temperature.py`
  — `from_config` path: confirm `wind_input_height`,
  `surface_z0`, `wind_shelter` are accepted from YAML.

This spec stands alone; any prior wind-function memos in this
folder can be removed after the changes here are implemented and
the corresponding tests pass.
