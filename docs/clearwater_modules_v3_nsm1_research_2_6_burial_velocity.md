# Phase 9.F.5 — `vb` (burial velocity) literature research

## Summary

- **Recommended v3 disposition: change the default from `0.01 m/d` to `6.85e-6 m/d`
  (i.e., `2.5e-3 m/yr` ≡ `0.25 cm/yr`), and re-document the units as `m/d`.**
  This matches the WASP7/WASP8 sediment-diagenesis default (verbatim, to 3
  significant figures), the legacy Fortran NSM1 default (after its built-in
  `/365` conversion), and the broad consensus of the Di Toro (2001) sediment-
  flux literature for low-to-moderate-deposition aquatic systems.
- **Authoritative cited values** (all converge on the same physical rate):
  - WASP7/WASP8 diagenesis manual (EPA), Layer-2 burial velocity default:
    `6.85e-6 m/day` (≈ `2.5e-3 m/yr` = `0.25 cm/yr`).
  - Fortran NSM1 (`modGlobalParam.f90:138`): `vb = 0.0025 m/yr`, with the
    `m/yr → m/d` conversion `vb / 365.0` performed at the use site
    (`modPOM.f90:114`); effective per-day rate `≈ 6.85e-6 m/d`.
  - Di Toro (2001) *Sediment Flux Modeling*: the two-layer aerobic/anaerobic
    sediment formulation that QUAL2K, WASP, and CE-QUAL-W2 all derive from
    treats `w₂` ("burial velocity") as a slow geological-scale process; the
    canonical Chesapeake Bay calibration (Di Toro & Fitzpatrick 1993) and
    subsequent applications (Brady, Testa, Di Toro 2013) all use values in
    the `0.1–0.5 cm/yr` range.
  - General sediment-accumulation literature (Baskaran 2015 review and the
    Springer 2023 freshwater-reservoir review): typical lake/reservoir
    accumulation rates `0.5–13 mm/yr` (`0.05–1.3 cm/yr`); fluvial systems
    `1.6–13 mm/yr` (`0.16–1.3 cm/yr`); both ranges bracket the
    Di Toro / WASP / Fortran-NSM1 default of `0.25 cm/yr`.
- **Headline finding:** v3's `vb = 0.01 m/d` (= `3.65 m/yr` = `365 cm/yr`)
  is **~1460× larger than the canonical value** and is well outside the
  range reported in any of the authoritative sources surveyed. The error
  was inherited from v1, which in turn was introduced by a partial unit
  conversion: v1's `processes.py:2293` removed the Fortran `/365` factor
  ("note removed 365 from FORTRAN") **without** also converting the
  default value from `m/yr` to `m/d`, so the parameter was relabeled as
  `m/d` while keeping its `m/yr` magnitude. The Fortran value of
  `0.0025` is, in this case, the correct canonical reference.

## Cross-model comparison

| Source | `vb` (m/d) | `vb` (m/yr) | `vb` (cm/yr) | Comment |
|---|---|---|---|---|
| **Di Toro (2001)** *Sediment Flux Modeling* | `~6.85e-6` | `~2.5e-3` | `~0.25` | Foundational two-layer sediment model; canonical Chesapeake Bay calibration; described as a slow geological process. |
| **QUAL2K** v2.11 (Chapra, Pelletier, Tao 2008) | (not tabulated as a default) | — | — | Section 5.6 implements the Di Toro framework verbatim; symbol table p. 88 documents `w₂ = burial velocity, m/d` but the user manual leaves the numerical default to the user (parameter file `Sediment` sheet inherits Di Toro values). |
| **WASP7 / WASP8** sediment-diagenesis module (EPA) | **`6.85e-6`** | `2.50e-3` | `0.25` | Documented default in the WASP diagenesis input table: `Burial velocity for layer 2 to inactive sediments (m/day) = 6.85E-06`. Active sediment thickness `H₂ = 0.1 m` (10 cm), matching NSM1. |
| **CE-QUAL-W2** sediment-diagenesis module (Berger, Wells et al.) | (per-application calibration) | `~2e-3 to 5e-3` | `~0.2 to 0.5` | Adopts Di Toro framework; per-reservoir calibration but stays in the Di Toro range. |
| **Stream/river accumulation literature** (Britannica; Baskaran 2015) | `~3e-6 to 4e-5` | `~1e-3 to 1.3e-2` | `~0.1 to 1.3` | Net long-term sediment accumulation in fluvial systems. |
| **Lake/reservoir accumulation literature** (Springer review 2023; Mendes 2023) | `~1.4e-6 to 4e-5` | `~5e-4 to 1.3e-2` | `~0.05 to 1.3` | Higher in agriculture-impacted reservoirs; lower in oligotrophic lakes. |
| **v1 NSM1** (`constants.py:325`) | `0.01` | `3.65` | `365` | Default labeled `m/d`; consumer `processes.py:2293`: `vb * POM / h2` with comment "removed 365 from FORTRAN". |
| **Fortran NSM1** (`modGlobalParam.f90:138`) | `~6.85e-6` (effective) | `0.0025` | `0.25` | Declared `m/a` (m/yr) at line 39; used as `vb / 365.0 * POM2 / h2(r)` at `modPOM.f90:114` so the per-day rate is `0.0025 / 365 ≈ 6.85e-6 m/d`. |
| **v3 NSM1** (`parameters/global_vars.py:42`) | `0.01` | `3.65` | `365` | Inherited verbatim from v1, including the units error. Already flagged with `FIXME(phase1-audit)`. |

> **Unit verification across the table.** The five conversions used:
> `m/d × 365 = m/yr`; `m/yr × 100 = cm/yr`; `m/d × 365 × 100 = cm/yr`. The
> WASP default `6.85e-6 m/d` × 365 = `2.5e-3 m/yr` = `0.25 cm/yr`. The
> Fortran default `0.0025 m/yr` ÷ 365 = `6.849e-6 m/d` (matches WASP to
> 3 sig fig — these are the same number expressed in different units).

## Trace through v3 POM consumer

The burial term in `src/clearwater_modules_v3/processes/pom.py:265` is

```python
rate_burial = self.vb * pom / self.h2
```

with the inline default `_POM_GLOBAL_DEFAULTS["vb"] = 0.01` (line 78,
labeled `m/d`) and `h2 = 0.1` (m, from `parameters/pom.py:11`).

### Dimensional check

| Symbol | Units | Numerical | Notes |
|---|---|---|---|
| `vb` | `m/d` | `0.01` (current v3) | Should be `~6.85e-6 m/d` per WASP/Fortran. |
| `pom` | `mg/L` | state variable | Bed-sediment POM concentration (Fortran `POM2`). |
| `h2` | `m` | `0.1` | Active anaerobic sediment layer thickness (Di Toro `H₂`). |
| `rate_burial` | `[m/d] · [mg/L] / [m] = mg/L/d` | — | Correctly resolved to a per-day rate. |

The dimensional formula is correct. The only defect is the **magnitude of
`vb`**: at the current default, the bed-POM compartment loses ~10% of its
mass per day to burial alone (`rate_burial / pom = 0.01 / 0.1 = 0.1 d⁻¹`),
which is geologically implausible. Replacing `vb` with the WASP default
`6.85e-6 m/d` yields `6.85e-5 d⁻¹` ≈ a `(1/6.85e-5) ≈ 14,600 d ≈ 40 yr`
e-folding burial timescale, consistent with the geological-scale timescales
described by Di Toro (2001) and observed in lake/reservoir cores.

### Comparison to other terms

The other POM sinks in v3 `pom.py`:

- `POM_dissolution = kpom_tc * POM` with `kpom = 0.1 d⁻¹` (typical) — this
  is the dominant loss term and is correctly first-order on the order
  of days.
- `vb / h2 = 0.01 / 0.1 = 0.1 d⁻¹` (current v3) accidentally **equals
  the dissolution rate**, which is the smell-test failure: burial should
  be many orders of magnitude slower than microbial mineralization in
  the active sediment layer, not equal to it.
- `vb / h2 = 6.85e-6 / 0.1 = 6.85e-5 d⁻¹` (WASP/Fortran value): about
  **`1500×` slower than dissolution**, which is physically reasonable
  (mineralization recycles labile POM back to the water column on a
  timescale of weeks; net burial is the small residual that escapes to
  the geological record).

## Recommendation and rationale

**Change v3's default `vb` from `0.01 m/d` to `6.85e-6 m/d`** in
`src/clearwater_modules_v3/parameters/global_vars.py:42` and the inline
fallback in `src/clearwater_modules_v3/processes/pom.py:78`. Replace the
existing `FIXME(phase1-audit)` comment with a citation pointing to WASP7
and Fortran NSM1.

**Rationale:**

1. **Three independent canonical sources agree:** WASP (`6.85e-6 m/d`),
   Fortran NSM1 (`0.0025 m/yr` ≡ `6.85e-6 m/d`), and the Di Toro (2001)
   sediment-flux literature converge on the same value. v3's `0.01 m/d`
   is an outlier by ~3 orders of magnitude.
2. **The error has a documented unit-conversion provenance.** v1's
   `processes.py:2293` author removed the Fortran `/365` conversion
   factor and relabeled the units from `m/yr` to `m/d`, but did not
   adjust the default value's magnitude. The Fortran source is
   internally consistent (`m/yr` declared, `/365` applied at use); v1
   silently broke that consistency. v3 inherited the broken default
   verbatim, and the Phase 0 audit's `FIXME` comment correctly flagged
   the magnitude as suspect even without identifying the root cause.
3. **Dimensional smell-test fails at v3's current value.** A burial
   velocity numerically equal to the active-sediment-layer dissolution
   rate (`vb/h2 ≈ kpom`) implies a "burial" process that competes with
   microbial mineralization — physically wrong. Burial should be a small
   residual on a much longer timescale (decades to centuries) than
   diagenesis (days to weeks).
4. **The Phase 0 audit explicitly anticipated this.** The Phase 0
   `FIXME(phase1-audit)` comment in `global_vars.py:42` reads
   "magnitude not validated"; this research validates that the magnitude
   is wrong by ~1500× and identifies the canonical replacement.

The fix is a single-line value change (no formula or unit change), no
test-suite changes (since current tests do not exercise long-time-horizon
mass-balance closure on POM2), and no API impact. The change should be
landed alongside an updated parameter-corrections document entry citing
this research and noting that the Fortran value was, in this case,
correct.

> **Note on Fortran-as-truth:** the [Phase 0 audit framing
> doc](./project_legacy_fortran_validation_status.md) and the user's
> guidance both emphasize that the legacy Fortran is *not* independently
> validated. This finding is therefore not "match Fortran to match
> Fortran"; it is "match Fortran *because* the Fortran value happens to
> agree with the WASP7 default and the Di Toro 2001 literature, all
> three of which agree with each other." If the Fortran value disagreed
> with the literature it would be a different question; here it does
> not.

## Sources

- [WASP Sediment Diagenesis Routines: Model Theory and User's Guide (EPA, 2017)](https://www.epa.gov/sites/default/files/2018-05/documents/wasp8_sod_module_v1.pdf) — local extracted copy of this PDF in the conversation tool-results cache shows the verbatim parameter table entry `Burial velocity for layer 2 to inactive sediments (m/day) (0.00000685) Y 6.85E-06` (p. 9 of the input-file specification, Appendix A).
- [QUAL2K v2.11 Documentation and User's Manual (Chapra, Pelletier, Tao 2008)](https://csdms.colorado.edu/csdms_wiki/images/Q2KDocv2_11b8.pdf) — Section 5.6 "Sediment Water Flux Model" pp. 65-75; symbol table p. 88 documents `w₂` units as `m/d`.
- [Di Toro, D.M., 2001. *Sediment Flux Modeling*. Wiley-Interscience.](https://www.wiley.com/en-us/Sediment+Flux+Modeling-p-9780471135357) (book reference; Chesapeake Bay calibration discussed in Chapter 14)
- [Di Toro, D.M., Fitzpatrick, J.J., 1993. *Chesapeake Bay Sediment Flux Model*. USACE Contract Report EL-93-2.](https://www.chesapeakebay.net/channel_files/35712/cb_sediment_flux_model_1993.pdf) — original calibration of the two-layer model that WASP and QUAL2K derive from.
- [Brady, D.C., Testa, J.M., Di Toro, D.M., et al. 2013. Sediment flux modeling: Calibration and application for coastal systems. *Estuarine, Coastal and Shelf Science* 117:107–124.](https://www.sciencedirect.com/science/article/abs/pii/S0272771412004374) — recent recalibration of the Di Toro model for multiple coastal systems.
- [Mendonça, R., et al. 2023. *A review of sedimentation rates in freshwater reservoirs* (Aquatic Sciences).](https://link.springer.com/article/10.1007/s00027-023-00960-0) — meta-analysis of sediment accumulation rates in reservoirs.
- [Baskaran, M. 2015. Sediment accumulation rates and sediment dynamics using five different methods.](https://huw.wayne.edu/research/baskaran-2015-sediment-accumulation-rates.pdf) — survey of measurement methods and typical values.
- [Bowie, G.L., et al. 1985. *Rates, Constants, and Kinetics Formulations in Surface Water Quality Modeling*, EPA/600/3-85/040.](https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=ORD&dirEntryId=34685) — historical EPA reference for kinetic parameter ranges.
- [v1 NSM1 `clearwater_modules/nsm1/constants.py`](../src/clearwater_modules/nsm1/constants.py) (line 325, default `vb = 0.01`).
- [v1 NSM1 `clearwater_modules/nsm1/processes.py`](../src/clearwater_modules/nsm1/processes.py) (line 2293, `vb * POM / h2 #note removed 365 from FORTRAN`).
- [Fortran NSM1 `modGlobalParam.f90`](file:///Users/todd/Downloads/NSM_comparison/NSM1/Source%20Files/modGlobalParam.f90) (line 39 declares units as `m/a`; line 138 sets `vb = 0.0025`; line 201 applies `vb = vb / 365.0` on input-file override).
- [Fortran NSM1 `modPOM.f90`](file:///Users/todd/Downloads/NSM_comparison/NSM1/Source%20Files/modPOM.f90) (line 114, `POM2_Burial = vb(r) / 365.0 * POM2 / h2(r)`).
- [v3 NSM1 `parameters/global_vars.py`](../src/clearwater_modules_v3/parameters/global_vars.py) (line 42, current default `vb = 0.01 m/d` flagged with `FIXME`).
- [v3 NSM1 `processes/pom.py`](../src/clearwater_modules_v3/processes/pom.py) (line 265, `rate_burial = self.vb * pom / self.h2`).
- [Phase 9.F.4 sibling research note on `h2`](./clearwater_modules_v3_nsm1_research_2_5_pom_h2.md) — pre-flagged the v1/Fortran `vb` units mismatch as a related finding.
