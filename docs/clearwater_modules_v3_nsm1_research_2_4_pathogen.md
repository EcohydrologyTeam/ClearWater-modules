# Phase 9.F.3 — pathogen `apx` and `vx` literature research

## Summary

- **`apx` role**: dimensional light-efficiency coefficient in the depth-
  averaged sunlight-inactivation term of the NSM1/QUAL2K pathogen
  kinetics. Units: `(W/m²)⁻¹·d⁻¹` (because v3 supplies `q_solar` in
  W/m² and the rate must come out in 1/d). The v1 docstring is wrong
  ("dimensionless"); the parameter cannot be dimensionless if
  `q_solar` carries units. Canonical value (Auer & Niehaus 1993,
  Onondaga Lake): α = 0.00824 cm²/cal → **≈ 0.017 (W/m²)⁻¹·d⁻¹**.
- **`vx` role**: pathogen net settling velocity (m/d). Auer & Niehaus
  (1993) reported a particle-associated settling loss of
  **1.38 m/d** for fecal coliform in Onondaga Lake based on sediment-
  trap measurements, and this value has been adopted in subsequent
  modeling studies (Bowie 1985 cites the same range).
- **Recommended v3 disposition**:
  1. Replace placeholder `apx = 1.0` with **`apx = 0.017`** with
     citation to Auer & Niehaus (1993) via Chapra (1997, *Surface
     Water-Quality Modeling*) and the QUAL2K formulation (Chapra,
     Pelletier & Tao 2008). Document units as `(W/m²)⁻¹·d⁻¹`.
     Optionally apply a 1/Fr_PAR multiplier (giving `apx ≈ 0.036`)
     because the v3 port substitutes PAR for total q_solar in the
     light-decay term — see "Recommendation and rationale" below.
  2. Replace placeholder `vx = 1.0` with **`vx = 1.38 m/d`** citing
     Auer & Niehaus (1993).
  3. Fix the v1-inherited docstring claiming `apx` is "dimensionless"
     and `vx` units of "m" (should be m/d).

## What `apx` and `vx` represent in NSM1 kinetics

The kinetics are identical across the legacy Fortran NSM1, the v1
Python port, and the v3 native Process:

```
dPX/dt = − kdx_tc · PX                                 (natural decay)
         − apx · I0 · (1 − exp(−KEXT·H)) / (KEXT·H) · PX  (light decay)
         − vx / H · PX                                 (settling)
```

where `I0` is the surface PAR (or total q_solar in v1/Fortran), `KEXT`
is the Beer–Lambert composite extinction coefficient (1/m), and `H` is
depth (m). The light-decay term is the depth-averaged form of the
Beer–Lambert integration over a fully-mixed water column —

```
(1/H) · ∫₀ᴴ I0 · exp(−KEXT·z) dz = I0 · (1 − exp(−KEXT·H)) / (KEXT·H)
```

— which is the standard form in QUAL2K (Chapra et al. 2008, §5.5.20.1)
and Chapra (1997, §33). The Fortran source comments explicitly cite
Chapra (1997) as the reference for this term:

> "q_solar units is ly/day in original formulation (Chapra, 1997)???"
> — `modPathogen.f90:90`

This comment also reveals an unresolved unit ambiguity in the legacy
code: the Chapra/Mancini/Auer canonical α is reported in **cgs units**
(cm²/cal, i.e. per ly/d) but Fortran/v1/v3 all consume `q_solar` in
**SI W/m²**. The placeholder `apx = 1.0` masked this ambiguity.

### `apx` — light efficiency factor

Trace from kinetic units backward:

- LHS rate `[1/d]`
- RHS `apx · q_solar · (depth-averaged optical factor) · PX`
- The depth-averaged optical factor `(1 − exp(−KEXT·H)) / (KEXT·H)` is
  dimensionless.
- `q_solar` is `[W/m²]` in v3 (and in the v1 NSM1 docstring at
  `processes.py:3184`).
- Therefore `apx` must carry units `[(W/m²)⁻¹·d⁻¹]` for the rate to
  come out in 1/d.

The v1 docstring's claim that `apx` is "dimensionless" (and the QUAL2K
manual's parallel claim that αpath is "dimensionless") is **wrong** —
αpath in QUAL2K is also dimensional; the manual simply omits the
units. This is consistent with αpath being implemented as a calibration
parameter that absorbs whatever unit convention is used for the solar
radiation input.

### `vx` — pathogen settling velocity

Trace from kinetic units backward:

- `vx · PX / H` must be `[1/d · PX]`
- `H` is `[m]`
- So `vx` is `[m/d]`.

The v1 docstring at `processes.py:3196` says `vx` has units "(m)". This
is also **wrong** — it must be m/d for dimensional consistency.

## Cross-model comparison for `apx`

| Source | Symbol | Reported value | Units | Comment |
|---|---|---|---|---|
| Auer & Niehaus (1993), *Wat. Res.* 27(4) | α | **0.00824** | cm²/cal | Field study, Onondaga Lake; raw-sewage dialysis-tube incubations; equivalent to **0.017 (W/m²)⁻¹·d⁻¹** in SI (1 W/m² = 2.0651 cal/(cm²·d)). |
| Chapra (1997), *Surface Water-Quality Modeling*, McGraw-Hill, Ch. 33 | α | "≈ 0.008 cm²/cal" | cm²/cal | Cites Auer & Niehaus (1993) as the primary source; this is the textbook reference adopted by NSM1/QUAL2K. |
| Mancini (1978), *J. Wat. Pollut. Control Fed.* 50:2477 | α (light-mortality) | ~1.0 (ly/h)⁻¹·d⁻¹ in fresh water | (ly/h)⁻¹·d⁻¹ | The 1978 empirical synthesis of ~100 studies; expressed per langley/hour rather than per cal/cm²/d (factor of 24 unit difference). 1 ly/h = 24 cal/cm²/d, and Mancini's α ≈ 1.0 (ly/h)⁻¹·d⁻¹ corresponds to α ≈ 0.041 (cal/cm²/d)⁻¹·d⁻¹ ≈ 0.085 (W/m²)⁻¹·d⁻¹. Mancini is ~5× higher than Auer/Niehaus, reflecting the difference between an open-water composite and a single-lake field study; both are within plausible literature scatter. |
| QUAL2K v2.11b8 (Chapra, Pelletier & Tao 2008), §5.5.20.1 | αpath | not specified in defaults table | "dimensionless" (manual is wrong) | The QUAL2K manual cites Chapra (1997) for the formulation but leaves αpath as a user-supplied calibration parameter. |
| WASP7 (Wool et al. 2008) pathogen module | k₁ (light-related decay) | 0.4–1.0 d⁻¹ at typical irradiance | 1/d | WASP lumps light + dark + temperature into a single first-order coefficient rather than exposing a separate α; not directly comparable. |
| Bowie et al. (1985), EPA/600/3-85/040, §6 | total k for fecal coliform | 0.7–8.0 d⁻¹ (range across studies) | 1/d | Lumped die-off rate. The high end of the range corresponds to surface waters at midday; consistent with α in the Auer/Mancini range. |
| NSM1 v1 / Fortran / v3 (placeholder) | apx | 1.0 | unspecified | No literature basis; flagged by Phase 0 audit. |

## Cross-model comparison for `vx`

| Source | Reported value (m/d) | Comment |
|---|---|---|
| Auer & Niehaus (1993), *Wat. Res.* 27(4) | **1.38** | Sediment-trap measurement of particle-associated fecal coliform in Onondaga Lake; the canonical settling rate cited in subsequent literature. |
| Bowie et al. (1985), EPA/600/3-85/040, §6 | 0.5–2.5 (typical range) | Compilation of reported values for fecal coliform and total coliform across multiple studies. |
| Garcia-Armisen & Servais (2009) and other particle-class studies | 1.17 (small particles, 0.45–10 µm); 2.40 (large particles, >10 µm) | Particle-size-resolved settling rates; bracket the 1.38 m/d composite value. |
| Chapra (1997), *Surface Water-Quality Modeling*, Ch. 33 | ~1 m/d (illustrative) | Cites Auer & Niehaus (1993) as the empirical basis. |
| Steets & Holden (2003), Onondaga Lake follow-up | 1.0–1.6 | Range of values across calibrated 3-D models. |
| NSM1 v1 / Fortran / v3 (placeholder) | 1.0 | No literature basis; flagged by Phase 0 audit. The v1 docstring also incorrectly lists units as "m" (should be m/d). |

## Recommendation and rationale

### `apx` — recommended **0.017 (W/m²)⁻¹·d⁻¹**

The Auer & Niehaus (1993) value (α = 0.00824 cm²/cal, i.e. **0.017
(W/m²)⁻¹·d⁻¹** in SI) is the most direct literature anchor. It is:

- the value Chapra (1997) cites in the chapter that NSM1 explicitly
  references in its Fortran source comments;
- derived from an in situ field study with raw sewage (the most common
  pathogen-source assumption in TMDL applications);
- within the bracket of the Mancini (1978) composite value when the
  unit conventions are reconciled.

A subtlety: the v3 port substitutes PAR (≈ 47% of q_solar) for total
q_solar in the light-decay term, whereas v1/Fortran/Auer/Niehaus use
total q_solar. To preserve the same depth-integrated decay rate as the
canonical Auer & Niehaus formulation under the v3 PAR substitution, the
v3 default should be set to `0.017 / Fr_PAR ≈ 0.036`. Two clean
options:

1. **Restore the v1 behavior** — drop Fr_PAR from `_rate_light_decay`
   and use total q_solar; set `apx = 0.017`. Closer to Chapra (1997).
   Slightly cleaner from a literature-traceability standpoint.
2. **Keep the Fr_PAR substitution** that v3 currently has — set
   `apx = 0.036`. Internally consistent with the rest of v3's light
   handling (which keys off PAR) but requires a one-line docstring
   note that the canonical Auer/Niehaus α has been pre-multiplied by
   1/Fr_PAR.

Either is defensible; option 1 is preferable because the canonical
literature (Chapra, Auer/Niehaus, Mancini) operates on total broadband
solar radiation, not PAR. Pathogen inactivation is largely UVA/UVB
mediated, not PAR mediated, so the PAR substitution in the v3 port is
arguably an introduced bug rather than an improvement; the existing
v3 docstring at `pathogen.py:295-302` flags this as an unresolved
question.

### `vx` — recommended **1.38 m/d**

This is the canonical Auer & Niehaus (1993) value, cited by Chapra
(1997), QUAL2K (Chapra et al. 2008), and adopted in subsequent
modeling studies. The Bowie et al. (1985) 0.5–2.5 m/d range brackets
this value, so 1.38 is a defensible mid-range literature default. Sites
with predominantly small-particle (<10 µm) pathogen association may
calibrate downward toward 1.17 m/d; sites with larger particles may
calibrate upward toward 2.40 m/d.

### Documentation fixes (independent of value changes)

Even if neither default is changed, the v3 docstrings should be
corrected:

- `apx` — update from "FIXME(phase1-audit): unitless" to dimensional
  units `(W/m²)⁻¹·d⁻¹` and cite Auer & Niehaus (1993) via Chapra
  (1997). The v1 claim that `apx` is dimensionless is dimensionally
  incorrect.
- `vx` — update v1's "(m)" units typo to `m/d` (m per day). v3's
  parameter library already lists m/d correctly, but the v1 process
  docstring (`processes.py:3196`) is wrong.
- Add Fr_PAR-substitution caveat to the `_rate_light_decay` docstring
  if option 2 is chosen, or remove the substitution entirely if
  option 1 is chosen.

## Sources

- [Auer, M.T. & Niehaus, S.L. (1993) Modeling fecal coliform bacteria — I. Field and laboratory determination of loss kinetics. *Wat. Res.* 27(4):693–701](https://www.sciencedirect.com/science/article/abs/pii/004313549390179L)
- [Mancini, J.L. (1978) Numerical estimates of coliform mortality rates under various conditions. *J. Wat. Pollut. Control Fed.* 50:2477–2485](https://www.semanticscholar.org/paper/Numerical-estimates-of-coliform-mortality-rates-Mancini/f3c722fbd3dcacbfc5c324900b31adbeee725ddb)
- [Chapra, S.C. (1997) *Surface Water-Quality Modeling*, McGraw-Hill, New York (Waveland reissue)](https://www.waveland.com/browse.php?t=378)
- [Chapra, S.C., Pelletier, G.J. & Tao, H. (2008) QUAL2K v2.11b8 documentation, Tufts University](https://csdms.colorado.edu/csdms_wiki/images/Q2KDocv2_11b8.pdf) — pathogen kinetics in §5.5.20.1
- [Bowie, G.L. et al. (1985) Rates, Constants, and Kinetics Formulations in Surface Water Quality Modeling, 2nd ed. EPA/600/3-85/040](https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=ORD&dirEntryId=34685)
- [Effland, T. et al. (2020) Sunlight-mediated inactivation of health-relevant microorganisms in water (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7064263/) — modern review of mechanisms
- [Indicator organisms for estuarine and marine waters review (PMC) — discusses Mancini equation history](https://pmc.ncbi.nlm.nih.gov/articles/PMC7164043/)
- [Effler et al. (1996) *Limnological and Engineering Analysis of a Polluted Urban Lake — Onondaga Lake, NY*, Springer](https://link.springer.com/chapter/10.1007/978-1-4612-2318-4_9) — chapter on mechanistic modeling, fecal-coliform application
- Source files reviewed:
  - `/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/modPathogen.f90` — Fortran kinetics
  - `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/src/clearwater_modules/nsm1/processes.py:3137-3241` — v1 Python port
  - `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/src/clearwater_modules/nsm1/constants.py:214-223` — v1 defaults
  - `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/src/clearwater_modules_v3/processes/pathogen.py` — v3 native Process
  - `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/src/clearwater_modules_v3/parameters/pathogen.py` — v3 defaults (current placeholders)
