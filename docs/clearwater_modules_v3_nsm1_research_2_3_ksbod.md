# Phase 9.F.2 — `ksbod_20` (CBOD settling) literature research

## Summary

- **Recommended v3 disposition**: **Keep `ksbod_20 = 0.0` as the default**, but
  fix three documentation/units bugs (see "Recommendation and rationale"
  below). The zero default is consistent with the modern peer-reviewed
  convention (QUAL2K/QUAL2Kw treat CBOD as fully **dissolved**, with
  particulate organic matter tracked separately as *detritus*); it is
  also consistent with how legacy Fortran NSM1 ships defaults and with
  EPA TMDL guidance practice for treated/dissolved CBOD ("`Ks` is
  assumed to be zero for the secondary effluent").
- **Authoritative cited values**:
  - **Bowie et al. (1985), Brown & Barnwell (1987) / QUAL2E**, EPA
    TMDL Technical Guidance (Book II, Eq. 2-4 + Sec. A.4): `Ks` ("BOD
    settling rate") in **/day**, defined as `Kr − Kd` (total CBOD
    removal rate minus deoxygenation rate). Typical ranges *implied*
    by the Kd-vs-depth correlations are ≤ a few tenths /day; specific
    Ks values are not tabulated authoritatively because Ks is
    fundamentally site-, source-, and treatment-specific.
  - **Chapra, Pelletier & Tao (2008), QUAL2K v2.11**: **No CBOD
    settling parameter exists in QUAL2K**. Both fast CBOD (`cf`) and
    slow CBOD (`cs`) are explicitly defined as the *dissolved*
    ultimate CBOD ("CBODFNU = ultimate dissolved CBOD"). Particulate
    organic matter is carried separately as detritus (`mo`) with
    settling velocity `vdt` [m/d].
  - **Yamuna River, Delhi (Parmar & Keshari, citing Kazmi & Agrawal
    2005)**: in a heavily polluted urban stretch where CBOD removal is
    dominated by settling of organic matter, calibrated K3 = **0.9
    /day** for all 16 reaches. This is the high end of the realistic
    range.
- **Headline cross-model comparison**: Modern QUAL2K and WASP CBOD
  modules contain **no CBOD settling rate at all**; QUAL2E (the older
  Brown & Barnwell 1987 framework) and its derivatives treat CBOD as a
  bulk pool that may include particulates and provide a `K3` settling
  rate in /day with a default of 0 and case-specific calibrated values
  spanning 0 to ≈ 1 /day. NSM1 is a QUAL2E-lineage formulation; its
  zero default is correct *for the dissolved-CBOD convention*, but the
  parameter is exposed for users who treat CBOD as a bulk pool.

## What `ksbod_20` represents in NSM1 kinetics

The NSM1 CBOD mass balance (legacy Fortran `modCBOD.f90` lines 109–116,
v1 `clearwater_modules/nsm1/processes.py`, and v3
`clearwater_modules_v3/processes/cbod.py`) is

```
dCBOD/dt = − kbod_tc · DOX/(DOX + KsOxbod) · CBOD     (oxidation)
         − ksbod_tc · CBOD                            (settling/sedimentation)
```

where `ksbod_tc = arrhenius_correction(ksbod_20, ksbod_theta, T)`. This
is the QUAL2E `K3` first-order settling rate, in **/day** (1/d). It
represents the apparent first-order loss of CBOD from the water column
due to settling of particulate carbonaceous material (i.e., it bundles
the *particulate fraction × particulate settling velocity / depth* into
a single lumped 1/d coefficient). It is therefore meaningful only when
CBOD is treated as a *bulk* pool that includes both dissolved and
particulate organic carbon.

### Units inconsistency among the three NSM1 implementations

| Implementation | Declared units of `ksbod_20` | Equation form |
|---|---|---|
| Legacy Fortran (`modCBOD.f90`) | "(1/day) Range {-0.36-0.36}" (line 13) | `CBOD_Sediment = ksbod_tc · CBOD` (line 114), no depth division |
| v1 Python (`clearwater_modules/nsm1/constants.py` line 12) | comment says `m/d at 20 C` | matches Fortran (no depth division) per Phase 0 audit |
| v3 Python (`clearwater_modules_v3/parameters/cbod.py` line 12; `processes/cbod.py` line 240) | comment says `m/d at 20 C` | code divides by depth: `settling_rate = ksbod_tc / depth · cbod` |

This is a real bug: **v3 and v1 disagree with the Fortran equation form**
on whether `ksbod_20` is a 1/d rate constant (Fortran convention) or an
m/d settling velocity (the v3 implementation form). With the default
`ksbod_20 = 0.0`, the bug is silent — both forms produce zero — but
any non-zero user value will produce values that differ from
Fortran/QUAL2E by a factor of `1/depth`. This needs a separate fix
ticket regardless of the disposition of the default value.

## Cross-model comparison

| Source | Parameter | Default value | Range | Units | Comment |
|---|---|---|---|---|---|
| **QUAL2E** (Brown & Barnwell 1987, EPA/600/3-87/007) | `K3` | 0.0 | -0.36 to +0.36 (per Fortran NSM1 comment, which appears to inherit this from QUAL2E) | /day | First-order CBOD settling. CBOD is a bulk pool. Negative values represent benthic resuspension. |
| **QUAL2K v2.11** (Chapra, Pelletier & Tao 2008) | (none) | — | — | — | Slow CBOD (`cs`) and fast CBOD (`cf`) are explicitly *dissolved* (CBODFNU). No settling term. Particulate organic carbon is carried separately as detritus (`mo`) with settling velocity `vdt` [m/d]. |
| **QUAL2Kw** (Pelletier, Chapra & Tao 2006; Bagmati R. case study, Kannel et al. 2007) | (none) | — | — | — | Same as QUAL2K — calibrated parameter set lists slow/fast CBOD oxidation and slow CBOD hydrolysis only; no CBOD settling. |
| **WASP7 EUTRO / Eutrophication** (Wool et al., EPA) | (none in EUTRO CBOD module) | — | — | — | WASP CBOD is treated as dissolved. POM (PON, POP, POC) settles via solids settling velocities (typically order 0.1–1 m/d) in the diagenesis routines, but CBOD itself does not have an independent settling parameter. |
| **CE-QUAL-W2** (Cole & Wells, ERDC) | (none for dissolved CBOD groups) | — | — | — | CBOD is a dissolved state variable; particulate OM (LPOM, RPOM) has its own settling velocity (POMS, RPOMS) typically 0.1–0.35 m/d. |
| **EPA TMDL Technical Guidance Book II** (Eq. 2-4; App. A.4) | `Ks` | 0 (assumed for secondary effluent in worked example B-3) | site-specific | /day | Defined as `Kr − Kd`. Sample text: "*The CBOD removal rate by settling (Ks) is assumed to be zero for the secondary effluent, and the CBOD oxidation rate (Kd) equals the total removal rate (Kr).*" |
| **Bowie et al. (1985)** EPA/600/3-85/040 | `Ks` | (no canonical default) | (not tabulated separately from `Kr/Kd`) | /day | Bowie's Table for CBOD kinetics (TMDL guidance Table A-20) summarizes `Kd` (deoxygenation) values 0.01 – 5.6 /day across 22 rivers. Ks is implicit in `Kr − Kd` only where field data resolved both. Bowie does provide phytoplankton settling velocities (Table A-29: 0.0–30 m/d) and POM settling, but no consolidated `Ks` table. |
| **Yamuna River QUAL2E case** (Parmar & Keshari, citing Kazmi & Agrawal 2005) | `K3` | — | 0.9 (uniform across 16 reaches) | /day | Heavily polluted urban stretch where particulate-laden CBOD settling dominated removal. Cited as the upper end of the practical range. |
| **Hudson, Patuxent, Wilsons Cr., etc.** (TMDL guidance Table A-21) | `Kd` | varies | 0.10 – 0.61 (pre-improvement); 0.15 – 0.35 (post-improvement) | /day | Treatment-plant-discharge-dominated streams. The implied `Ks = Kr − Kd` is small (single tenths /day) where data resolve both, and cannot be reliably extracted from this table alone. |
| **v1 NSM1** (`clearwater_modules/nsm1/constants.py` `DEFAULT_CBOD`) | `ksbod_20` | **0.0** | — | comment says "m/d", code matches Fortran (1/d) | Inherited from Fortran NSM1. |
| **Fortran NSM1** (`modCBOD.f90` ll. 13, 36) | `ksbod%rc20` | **0.0**; theta=1.024 | "{-0.36-0.36}" | /day | Inherits QUAL2E `K3` default and theta. Equation `CBOD_Sediment = ksbod_tc · CBOD` is dimensionally consistent with /day. |
| **v3 NSM1** (`clearwater_modules_v3/parameters/cbod.py`) | `ksbod_20` | **0.0**; theta=1.047 | — | comment says "m/d at 20 C" | Inherited from v1; `theta` differs from Fortran (1.047 vs 1.024), and code form `ksbod_tc/depth · CBOD` *implements* the parameter as if it were m/d, which conflicts with the Fortran convention. |

## Recommendation and rationale

### Disposition: keep `ksbod_20 = 0.0`, but fix three documentation/code bugs.

**Why keep the zero default:**

1. **Modern peer-EPA models (QUAL2K, QUAL2Kw, WASP, CE-QUAL-W2) do
   not have a CBOD settling parameter at all.** They treat CBOD as
   dissolved-only and route particulate organic carbon through a
   separate detritus / POM state variable with an explicit settling
   velocity. Setting NSM1's `ksbod_20` to zero by default is the
   modern-convention-equivalent behavior.
2. **The legacy QUAL2E `K3` default (which NSM1 inherits) is also 0**
   in the original Brown & Barnwell (1987) documentation. Brown &
   Barnwell describe `K3` as a user-supplied calibration parameter for
   cases where the modeler wants to lump particulate CBOD into the
   bulk CBOD pool — explicitly noting that 0 is the correct default
   when CBOD is dissolved or after secondary treatment.
3. **EPA TMDL guidance (Book II, sample calc B-3) explicitly assumes
   `Ks = 0` for treated effluent**, with the caveat that for
   particulate-laden, untreated discharges the value can be calibrated
   from the difference `Kr − Kd`. NSM1's approach (provide the knob,
   default it to zero) is therefore consistent with current EPA
   guidance.
4. **Sediment coupling supersedes lumped settling.** v3's roadmap
   includes proper PSDM/HEC-RAS sediment exchange and (eventually) a
   v3-native diagenesis module (NSM2 path). When particulate
   organic matter is properly tracked, an additional lumped CBOD
   settling term would double-count.
5. **Conservatism for unvalidated Fortran.** The Phase 0 audit
   correctly flagged `ksbod_20 = 0` as suspect, but the literature
   review confirms zero is the *intentional, defensible* default. The
   FIXME comment should be cleared with citation, not the value
   changed.

### Required documentation/code fixes (separate ticket — out of scope here):

1. **Resolve the units inconsistency.** Decide whether `ksbod_20` is
   a /day rate constant (Fortran/QUAL2E convention, recommended for
   v3 to match the canonical literature) or an m/d settling velocity
   (current v3 code, but no published reference uses this form for
   CBOD). The v3 code in `processes/cbod.py` divides by depth, which
   is *not* what Fortran NSM1 does. Recommended fix: change v3 to
   match Fortran (`settling_rate = ksbod_tc · cbod`, no depth
   division) and update the parameter docstring to "1/d at 20 C".
2. **Align `ksbod_theta` with Fortran.** v3 uses 1.047 while Fortran
   uses 1.024. The Bowie (1985) and QUAL2E recommendation for
   sedimentation/settling-type coefficients is 1.024; for biological
   reactions it is 1.047. Change v3 to 1.024 to match Fortran and
   QUAL2E.
3. **Replace the FIXME comment with a citation block** noting that 0
   is the correct default per QUAL2K (no settling), Brown & Barnwell
   1987 (default 0), and EPA TMDL guidance Book II Sec. B-3 (`Ks = 0`
   for secondary effluent), and pointing users at QUAL2E case studies
   (e.g., Yamuna River K3 = 0.9 /day) for guidance when they need to
   lump particulate CBOD.

## Sources

- [Chapra, S.C., Pelletier, G.J., & Tao, H. (2008). *QUAL2K: A Modeling Framework for Simulating River and Stream Water Quality, Version 2.11: Documentation and Users Manual.* Tufts University / Washington Dept. of Ecology / EPA.](https://csdms.colorado.edu/csdms_wiki/images/Q2KDocv2_11b8.pdf) — Sections 5.5.9–5.5.10 define slow/fast CBOD as dissolved with no settling term; Section 5.5.8 defines detritus settling separately.
- [Brown, L.C. & Barnwell Jr., T.O. (1987). *The Enhanced Stream Water Quality Models QUAL2E and QUAL2E-UNCAS: Documentation and User Manual.* USEPA, Athens GA. EPA/600/3-87/007.](https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=NERL&dirEntryId=41777) — Original `K3` definition; cited by NSM1 Fortran as the source of the {-0.36, +0.36} /day range.
- [Bowie, G.L., et al. (1985). *Rates, Constants, and Kinetics Formulations in Surface Water Quality Modeling*, 2nd ed. EPA/600/3-85/040.](https://cfpub.epa.gov/si/si_public_record_report.cfm?Lab=ORD&dirEntryId=34685) — The canonical EPA reference for water-quality kinetics. Provides Kd values across 22+ rivers (TMDL Guidance Table A-20 reproduces this).
- [USEPA (1997). *Technical Guidance Manual for Performing Wasteload Allocations Book II: Streams and Rivers, Part 1: BOD/DO and Nutrients/Eutrophication.*](https://19january2021snapshot.epa.gov/sites/static/files/2019-12/documents/technical-guidance-tmdl-book2.pdf) — Eq. 2-4 (`dL/dt = -(Kd+Ks)L`), Sec. A.4 ("CARBONACEOUS DEOXYGENATION RATE"), Tables A-20, A-21, and Sample Calculation B-3 ("`Ks` is assumed to be zero for secondary effluent").
- [Kannel, P.R., et al. (2007). Application of automated QUAL2Kw for water quality modeling and management in the Bagmati River, Nepal. *Ecological Modelling* 202(3-4): 503–517.](https://apps.ecology.wa.gov/publications/documents/0703035.pdf) — QUAL2Kw calibrated parameter set; lists slow/fast CBOD oxidation rates and slow CBOD hydrolysis; **no CBOD settling parameter**.
- [Parmar, D.L. & Keshari, A.K. (n.d.). *Calibration and Validation of QUAL2E model on the Delhi stretch of the Yamuna River.*](https://swat.tamu.edu/media/56955/h1-3-parmar.pdf) — Uses QUAL2E formulation with `K3` (BOD settling) = 0.9 /day uniform across 16 reaches; cites Kazmi & Agrawal (2005) as the source.
- [Wool, T.A., Ambrose, R.B., Martin, J.L., et al. *Water Quality Analysis Simulation Program (WASP) v8 User Documentation*, US EPA.](https://www.epa.gov/hydrowq/wasp-model-documentation) — WASP EUTRO carries CBOD as a dissolved state variable with no independent settling term; particulate organic matter (PON/POP/POC) settles via solids velocities.
- [Cole, T.M. & Wells, S.A. *CE-QUAL-W2: A Two-Dimensional, Laterally Averaged Hydrodynamic and Water Quality Model.* US Army Corps of Engineers / Portland State University.](https://www.ce.pdx.edu/w2/) — Treats CBOD groups as dissolved; LPOM/RPOM particulate OM has its own settling velocities POMS/RPOMS.
- [Chapra, S.C. (1997). *Surface Water-Quality Modeling.* McGraw-Hill / Waveland Press.](https://www.waveland.com/browse.php?t=378) — Standard textbook reference; bottle-test rates for sewage-derived CBOD = 0.05–0.3 /day (cited in QUAL2K manual p. 38). Confirms CBOD framework treats CBOD as dissolved.
- Legacy Fortran NSM1 source: `/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/modCBOD.f90` (lines 13, 36, 114) — `ksbod%rc20 = 0.0`, `theta = 1.024`, units stated as `(1/day)`, range `{-0.36-0.36}`.
- v1 Python NSM1 source: `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/src/clearwater_modules/nsm1/constants.py` line 172, `DEFAULT_CBOD` at line 178 (`ksbod_20 = 0.0`, `ksbod_theta = 1.047`).
- v3 NSM1 sources: `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/src/clearwater_modules_v3/parameters/cbod.py` and `.../processes/cbod.py`.
