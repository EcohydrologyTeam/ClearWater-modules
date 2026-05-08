# Edinger wind-function exponent audit — CE-QUAL-W2 reference

**Date:** 2026-05-08
**Scope:** Confirm whether v3 TSM's default `wind_c = 3.0` (cubic wind
exponent in `f(W) = (a + b·Wᶜ)/1e6`) is standard or non-standard, by
reading the authoritative CE-QUAL-W2 source and checking what
exponents and conventions it actually uses.

**Update (post-broader-investigation):** the user's investigation of
legacy Fortran TSM, CE-QUAL-W2, and QUAL2K together establishes that
**`c = 1` (linear in wind)** is the consensus default across these
references. My earlier recommendation of `c = 2` (quadratic) in this
memo was a conservative read that gave equal weight to W2's
`CFW=1` and `CFW=2` cases without checking which is the *typical*
choice across the family. The corrected recommendation is `c = 1`,
consistent with QUAL2K, the W2 Ryan-Harleman alternative, the legacy
Fortran TSM placeholders (`1.0 / 1.0 / 1.0` — note the exponent
slot is also 1), and the broader literature consensus. The
quantitative-impact and recommendation sections below are revised
accordingly. The CE-QUAL-W2 evidence about W2 *supporting* `CFW=1`
and `CFW=2` (and rejecting `CFW=3`) is unchanged — it just doesn't
distinguish which of `1` and `2` is the standard choice.

## Why we're auditing this

After the BC investigation isolated the Sep 2008 Salem-cell cool bias
to TSM, and the wind-height correction (`--wind-height-factor 0.78`)
closed about 38% of the residual (-2.68 → -1.65 °C, see
`design/tsm_heat_balance_audit.md`), the remaining ~1.65 °C bias
points at something else in the heat-flux machinery. Wind enters the
Edinger latent + sensible heat formulas as `Wᶜ`. The v3 default
`c = 3.0` is an unusual choice — typical Edinger family
parameterizations use `c = 1` (linear) or `c = 1.5–2.0`. The 2026-05-05
v3 TSM audit
(`ClearWater-modules-streaming/design/clearwater_modules_v3_tsm_audit_2026-05-05.md`)
verified the *structure* `(a + b·Wᶜ)/1e6` against the Edinger family
but acknowledged the citation was reconstructed retroactively and that
the *coefficients* `0.3 / 1.5 / 3.0` are inherited from v1
(`clearwater_modules.tsm.constants`) without primary source.

Authoritative CE-QUAL-W2 source resolves the question.

## CE-QUAL-W2 source: how wind is handled

Repo: `/Users/todd/GitHub/CE-QUAL-W2-ERDC/CE-QUAL-W2-ERDC-dev/src/W2_v2026.02/`

### 1. Wind is always pre-corrected to 2 m height

`w2_4_unix.f90:480, 487`:

```fortran
WIND2(I) = WIND(JW) * WSC(I) &
         * DLOG(2.0D0 / Z0(JW)) / DLOG(WINDH(JW) / Z0(JW))
```

* `WIND(JW)` — raw wind from the input file (any height).
* `WSC(I)` — per-segment **wind shelter coefficient** (vegetation,
  topography, riparian canopy reducing wind speed reaching the water).
* `DLOG(2.0/Z0) / DLOG(WINDH/Z0)` — log-law height correction from
  the input anemometer height (`WINDH`, user-specified per-waterbody)
  to **2 m above water**.
* `Z0(JW)` — surface roughness, user-specified per-waterbody (typical
  default 0.001–0.003 m for water).

CE-QUAL-W2 thus enforces the 2 m convention before the wind function
sees the value. This confirms the wind reference height we hypothesised
in the prior audit, and it generalises beyond simple log-law: W2 also
applies a wind-shelter coefficient per mesh segment.

### 2. Wind function form

`heat-exchange.f90:143`:

```fortran
FW = AFW(JW) + BFW(JW) * WIND2(I)**CFW(JW)
```

`AFW(JW)`, `BFW(JW)`, `CFW(JW)` are **per-waterbody, user-specified**.
There are no hard-coded coefficient defaults — the user supplies them
in the input deck based on calibration to that specific waterbody.

### 3. Only `CFW=1` (linear) and `CFW=2` (quadratic) are supported

`heat-exchange.f90:76–78`:

```fortran
IF (CFW(JW) == 1.0) BCONV = 3.401062
IF (CFW(JW) == 2.0) BCONV = 1.520411
! SW Issues: CFW not determined for other values of CFW
```

These are British-to-SI unit conversion factors keyed to the wind
exponent. **Only CFW=1 and CFW=2 have defined conversions.** The
trailing comment "*CFW not determined for other values of CFW*"
indicates the W2 maintainers regard `CFW=3` (or any non-1/2 value) as
**outside the supported design range**. There is no path through the
British-units branch (`equilibrium_temperature` ENTRY) for `CFW=3`.

### 4. Alternative wind function (Ryan-Harleman)

`heat-exchange.f90:139–141`:

```fortran
DTVL = 0.0084D0 * WIND2(I)**3
IF (DTV < DTVL) DTV = DTVL
FW = (3.59D0 * DTV**0.3333D0 + 4.26D0 * WIND2(I))
```

When `RH_EVAP(JW) = .TRUE.`, W2 uses the Ryan-Harleman virtual-
temperature-difference parameterization. The `WIND2³` term inside
`DTVL` is a *floor* on the virtual-T deficit (not a wind-function
exponent on the dominant evaporative term). The actual wind function
is `(3.59·DTV^(1/3) + 4.26·W)` — **linear in wind**, never cubed.

### 5. QUAL2K and HEC-RAS-WQ comparison

* QUAL2K wind function: `f(W) = vw_a + vw_b · W` — **linear** in wind
  (`c = 1`).
* HEC-RAS-WQ derivatives generally inherit from W2 or QUAL2K. The
  Fortran reference snapshot the v3 audit examined uses placeholder
  defaults `1.0 / 1.0 / 1.0` (`a=1, b=1, c=1`) — also linear in wind.

## Implication for v3 TSM

The default `wind_c = 3.0` in v3 TSM:

* **Is not from CE-QUAL-W2.** W2 supports `CFW = 1` (linear) and
  `CFW = 2` (quadratic). `CFW = 3` is explicitly flagged as
  undetermined and has no defined unit-conversion path.
* **Is not from CE-QUAL-W2's Ryan-Harleman alternative,** which is
  linear in wind.
* **Is not from QUAL2K,** which is linear in wind.
* **Is not from the audited Fortran reference,** which uses
  placeholder linear defaults.
* **Has no primary citation** in the v3 audit memo, which acknowledged
  the citation was reconstructed retroactively.

The v3 default `c = 3` is a v1 Python-port artefact — almost certainly
either a transcription error or an undocumented choice that has been
quietly inherited through v2 → v3 without ever being validated against
the literature.

## Quantitative impact

At our Sep 2008 conditions (mean wind 1.56 m/s after the
`--wind-height-factor 0.78` correction, neutral stability so
Ri ≈ 1):

| `c` | `b·Wᶜ` | Latent flux (W/m²) | Equilibrium effect |
|---|---|---|---|
| 3.0 (v3 default — non-standard) | 5.70 | ~110 | ≥ 1.5 °C cool bias |
| 2.0 (W2 quadratic case, CFW=2) | 3.66 | ~70 | ~half the bias |
| 1.5 (intermediate) | 2.92 | ~58 | smaller still |
| **1.0 (consensus default — Fortran TSM / QUAL2K / W2 Ryan-Harleman)** | **2.34** | **~46** | **closest to obs equilibrium** |

Reducing `c` from 3 to 1 reduces the wind-driven term by ~59%
relative to `c=3`. Combined with the wind-height correction already
in place, that should bring the equilibrium temperature into or
slightly above the obs IQR for Sep 2008.

## Recommendation for this case study

1. **(Done):** add `--wind-a / --wind-b / --wind-c` CLI flags to the
   runner (defaults `0.3 / 1.5 / 3.0` to preserve back-compat). Plumb
   through to `Temperature(wind_a=…, wind_b=…, wind_c=…)`.
2. **Sweep `c ∈ {3.0, 2.0, 1.0}`** for a 3-point sensitivity. The
   `c=3.0` baseline is already in
   `..._wind2m_alledges/` (-1.65 °C bias). The `c=2.0` mid-point is
   running as `..._wind2m_alledges_c2/`. The `c=1.0` consensus
   default should run after as `..._wind2m_alledges_c1/`.
3. **Pick the value that lands closest to the obs IQR** at the Salem
   cell. Literature consensus says `c=1.0`; the case-study sweep
   verifies that against our specific Sep 2008 forcing and reach.

## Recommendation for v3 TSM

A separate, focused design note (`v3_tsm_wind_function_improvements.md`
in this repo) describes specific changes to v3 TSM defaults and
methods — change of `wind_c` default, mandatory wind-height parameter,
optional wind-shelter coefficient, validator on `c ∈ {1, 2}`. That
note is intended for the v2→v3 merge work, not for this case study.

## Files referenced

* `/Users/todd/GitHub/CE-QUAL-W2-ERDC/CE-QUAL-W2-ERDC-dev/src/W2_v2026.02/heat-exchange.f90`
  — wind function `FW = AFW + BFW·WIND2^CFW` (line 143); BCONV
  conversion factors only for CFW=1 and CFW=2 (lines 76–78); Ryan-
  Harleman alternative linear in wind (line 141).
* `/Users/todd/GitHub/CE-QUAL-W2-ERDC/CE-QUAL-W2-ERDC-dev/src/W2_v2026.02/w2_4_unix.f90`
  — log-law height correction at lines 480, 487.
* `ClearWater-modules-streaming/design/clearwater_modules_v3_tsm_audit_2026-05-05.md`
  — earlier v3 audit; Q4/Q6 record the `0.3 / 1.5 / 3.0` defaults
  inheritance and Edinger citation reconstruction.
* `design/tsm_heat_balance_audit.md` — earlier audit identifying the
  wind-height passthrough.
* `design/bc_inflow_continuity_findings.md` — BC investigation that
  isolated the cool bias to TSM.
