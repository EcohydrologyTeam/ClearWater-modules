# v3 NSM1 Science-Correctness Validation: CBOD, POM/Organic Matter, Phosphorus Partitioning, Pathogen

**Review date:** 2026-05-15
**Reviewer:** Claude Code (water-quality model source-code review agent, claude-opus-4-7)
**Repository:** ClearWater-modules-streaming, branch `streaming`, commit `54f2b12`
**Scope axis:** Science correctness against *validated* references. Parity (v3↔v1) was
established by prior agents (`review_cbod_pom.md`, `review_phosphorus_pathogen.md`); this
review tests whether the *shared* v1/v3 kinetic forms and defaults are themselves defensible,
with emphasis on terms where v1 and v3 may be *jointly* wrong.

**Authoritative references consulted:**

- CE-QUAL-W2 v2026.02 (LOCAL, AUTHORITATIVE): `CE-QUAL-W2-ERDC-dev/src/W2_v2026.02/water-quality.f90`.
  Entries read: `KINETIC_RATES` (283), `BACTERIA` (1183), `PHOSPHORUS` (1409), `LABILE_POM` (1699),
  `REFRACTORY_POM` (1740), `BIOCHEMICAL_O2_DEMAND` (1813), `DISSOLVED_OXYGEN` (1859),
  CBODD coefficient (440), PARTP PO4 sorption (447--464), OM decay rates (342--381),
  PO4 budget incl. `PO4NS` sorbed-settling (1448, 1472).
- QUAL2K v2.11b8 (Chapra, Pelletier & Tao 2008); QUAL2E (Brown & Barnwell 1987, EPA/600/3-87/007);
  WASP7 EUTRO; Bowie et al. 1985 (EPA/600/3-85/040); Chapra 1997 (*Surface Water-Quality Modeling*);
  Auer & Niehaus 1993 (*Wat. Res.* 27(4):693--701); Mancini 1978; Thomann & Mueller 1987.

---

## 1. Verdict

The shared v1/v3 NSM1 kinetic forms for CBOD oxidation, POM hydrolysis/settling, organic-P
hydrolysis/settling, and pathogen die-off are scientifically defensible first-order
formulations consistent with QUAL2E, QUAL2K, WASP, and Bowie et al. (1985). None of the
five adjudicated terms exhibits a *form* error that is jointly wrong in v1 and v3 under the
shipped default parameters. The single most important science finding is that the **v3
phosphorus sorption isotherm is correct and the v1 NSM1 runtime stub (`fdp ≡ 1.0`) is
wrong**: v3 restores the linear-equilibrium isotherm that CE-QUAL-W2 (`PARTP`, water-quality.f90:449)
and QUAL2K both implement. This is the headline V3-CORRECT-V1-WRONG ruling and the prior
parity reviews correctly identified it (parity F1/F2) but, being parity-scoped, did not
adjudicate it against the validated W2/QUAL2K isotherm. This review supplies that adjudication
(Section 3.2).

Findings: zero terms are **V1&V3-BOTH-WRONG in form**. One term is **V3-CORRECT-V1-WRONG**
(phosphorus `fdp` sorption). Three items are documented v3 default corrections that this
review confirms as literature-defensible improvements over v1 placeholders (`apx`, `vx`,
`vsop`). Two items are **STRUCTURAL-OK** (NSM1 lumped CBOD vs W2 split LDOM/RDOM/LPOM/RPOM;
NSM1 single-layer POM2 burial vs W2 settling cascade). Two carry **science-finding caveats
the parity reviews could not raise because they only compared v3 to v1**: (a) the v3 CBOD
DO-attenuation uses a Monod half-saturation form whereas W2 uses a binary oxic/anoxic switch
and QUAL2E/Bowie use a *first-order* (no DO limitation) decay -- the v3/v1 Monod form is the
*more* defensible modern choice but its `KsOxbod = 0.5 mg/L` default and the consequence at
low DO deserve a note (Section 3.1); (b) the v3 `ksbod_theta = 1.047` for CBOD *settling*
contradicts the canonical Bowie/QUAL2E settling-Arrhenius value of 1.024 that CE-QUAL-W2
does not even apply (W2 has no settling-rate Arrhenius). This is latent at `ksbod_20 = 0`
but is a genuine *defaults-vs-literature* defect inherited from v1, not just a v1↔Fortran
parity note (Section 3.1, Section 5).

No defect blocks v3 1.0.0 at the validated default-parameter regime, because every divergent
term is gated to zero by a shipped default (`ksbod_20 = 0`, `kdpo4 = 0`, `rpo4_20 = 0`).
The science risk is entirely in the *calibration* regime: a user who activates any of these
terms inherits a form or coefficient that, while parity-faithful to v1, is not the value a
QUAL2K/W2/Bowie-literate modeler would expect.

---

## 2. Cross-Model Validation Matrix

Verdict legend: CONSISTENT = v3 form/default agrees with the validated reference family.
V3-DIVERGENT = v3 differs from validated references (regardless of v1). V1&V3-BOTH-WRONG =
the shared v1/v3 form or default is wrong against validated references. V3-CORRECT-V1-WRONG
= v3 fixes a genuine v1 defect, confirmed against a validated reference. STRUCTURAL-OK =
v3/v1 model structure differs from W2 but the comparable kinetic form and coefficients are
defensible.

### CBOD

| Term / default | v3 | v1 | CE-QUAL-W2 (file:line) | QUAL2K / QUAL2E | WASP/Bowie | Literature value | VERDICT |
|---|---|---|---|---|---|---|---|
| CBOD oxidation form | `kbod_tc·DOX/(KsOxbod+DOX)·cbod` (cbod.py:302-304) | same (processes.py:2386) | `KBOD·TBOD^(T-20)·DO3·CBOD`, `DO3` = binary oxic switch (water-quality.f90:440, 1819; DO3 at 298) | QUAL2E/QUAL2K: first-order `k·CBOD`, *optional* DO half-sat in QUAL2K-CBODs | WASP EUTRO: optional `k·CBOD·DO/(K+DO)` | Monod DO limitation is the modern WASP/QUAL2K-fast option | CONSISTENT (Monod form is defensible; differs in mechanism from W2 binary switch -- see 3.1) |
| `kbod_20` | 0.12 /d | 0.12 | `KBOD` read from control (typ. 0.02--0.30 /d) | QUAL2E `K_1` 0.02--3.4 /d, typ 0.1--0.5 | Bowie Table: 0.05--0.5 /d carbonaceous | 0.1--0.3 /d typical | CONSISTENT |
| `kbod_theta` | 1.047 | 1.047 | `TBOD` (typ. 1.047) | QUAL2E `K_1` θ = 1.047 | Bowie: 1.047 (BOD oxidation) | 1.047 | CONSISTENT |
| `KsOxbod` (DO half-sat) | 0.5 mg-O2/L | 0.5 | n/a (W2 uses binary `DO3`, not half-sat) | QUAL2K CBOD-fast `Ksocf` ~0.5--2 mg/L | WASP half-sat ~0.5 mg/L | 0.5--2.0 mg-O2/L | CONSISTENT (low end of literature range) |
| CBOD settling form | `ksbod_tc/depth·cbod` (cbod.py:310; velocity m/d) | `ksbod_tc·cbod` (rate 1/d; processes.py:2403) | `CBODS·(CBOD_{k-1}-CBOD_k)·BI/BH2` (settling *velocity*, water-quality.f90:1817) | QUAL2K: CBOD is dissolved-only, no settling | Bowie: CBOD settling rare; particulate carried as detritus | velocity convention if used | V3-DIVERGENT vs v1 form; v3's velocity form is closer to W2's velocity convention -- but see 3.1 |
| `ksbod_20` | 0.0 | 0.0 | `CBODS` default 0 | QUAL2E `K_3` default 0 | 0 (modern dissolved-CBOD convention) | 0 | CONSISTENT (intentional zero) |
| `ksbod_theta` | 1.047 | 1.047 | W2 applies **no Arrhenius to CBODS settling** | QUAL2E settling `K_3` θ = **1.024** | Bowie: settling-coeff θ = 1.024 | 1.024 (settling), 1.047 (decay) | **V1&V3-BOTH-WRONG** (latent at ksbod_20=0; see 3.1, Section 5) |
| CBOD→O2 stoich (`roc`) | reciprocal of OM C:O ratio, applied in DOX process | same | `RBOD(JCB)` per-group, default ~1.0--1.2 (water-quality.f90:1868) | QUAL2E `roc` 2.69 gO2/gC if CBOD in C; 1.0 if CBOD in O2-equiv | Bowie: CBOD usually in O2-equiv → roc=1 | depends on CBOD basis | CONSISTENT (NSM1 CBOD is O2-equivalent; roc handled in DOX) -- see 3.5 |
| CBOD→DIC stoich | via OM C:O in Carbon process | same | `BODC(JCB)·CBODD·CBOD` (water-quality.f90:1942) | QUAL2K: CBOD-fast → DIC via `rcco` | WASP: CBOD → DIC | 0.375 gC/gO2 (32/12 reciprocal) | CONSISTENT (deferred to Carbon/DIC; see 3.5) |

### POM / Organic Matter

| Term / default | v3 | v1 | CE-QUAL-W2 (file:line) | QUAL2K / QUAL2E | WASP/Bowie | Literature value | VERDICT |
|---|---|---|---|---|---|---|---|
| POM hydrolysis/dissolution | `kpom_tc·pom` (pom.py:340) | same (processes.py:2233) | `LPOMD = OMTRM·LPOMDK·LPOM·DO3` (water-quality.f90:354) | QUAL2K detritus dissolution first-order | Bowie: POM hydrolysis first-order | first-order | STRUCTURAL-OK (NSM1 lumps; W2 splits LPOM/RPOM; same first-order kinetic) |
| `kpom_20` | 0.1 /d | 0.1 | `LPOMDK` typ. 0.08--0.12 /d (labile) | QUAL2K detritus dissolution 0.01--0.25 /d | Bowie POM: 0.03--0.25 /d | 0.05--0.12 /d labile | CONSISTENT |
| `kpom_theta` | 1.047 | 1.047 | W2 `OMTRM` (Arrhenius-like temp multiplier) | QUAL2K θ = 1.047 | Bowie: 1.047 | 1.047 | CONSISTENT |
| POM settling/burial | `vb·pom/h2` burial (pom.py:355) | same (processes.py:2293) | `POMS·(LPOM_{k-1}-LPOM_k)·BI/BH2` water-column settling (water-quality.f90:1709) | QUAL2K detritus settling velocity | Bowie POM settling 0.1--2.5 m/d | settling 0.1--2 m/d | STRUCTURAL-OK (NSM1 POM is Di Toro layer-2 bed POM with burial; W2 POM is water-column with settling -- comparable velocity convention) |
| `vsoc` (POC→POM settling) | 0.01 m/d | 0.01 | W2 `POMS` water-column | QUAL2K particulate-C settling | Bowie: 0.1--2.5 m/d range | 0.01--2 m/d | CONSISTENT (low end; defensible for fine POC) |
| `vb` (burial velocity) | 6.85e-6 m/d | 6.85e-6 (corrected) | n/a (W2 has no bed-POM burial in this layer model) | Di Toro 2001 burial w2 ~0.25--0.5 cm/yr | Di Toro: 0.25 cm/yr | 6.85e-6 m/d = 0.0025 m/yr | CONSISTENT (Di Toro burial; v1 unit bug corrected in v3) |
| `h2` (sediment layer) | 0.1 m | 0.1 | W2 `BH2` is geometric, not Di Toro H2 | Di Toro/QUAL2K H2 = 0.1 m | Di Toro 2001 §: 0.1 m | 0.1 m | CONSISTENT |
| `fcom` (OM C fraction) | 0.4 | 0.4 | W2 `ORGC(JW)` ~0.45 | QUAL2K `rcd` ~0.40 gC/gD | Redfield C:D ~0.40--0.45 | 0.40--0.45 | CONSISTENT |

### Phosphorus partitioning and kinetics

| Term / default | v3 | v1 (as run) | CE-QUAL-W2 (file:line) | QUAL2K / QUAL2E | WASP/Bowie | Literature value | VERDICT |
|---|---|---|---|---|---|---|---|
| **Dissolved-P fraction `fdp`** | `1/(1+kdpo4·Solid·1e-6)` (partitioning.py:50) | `xr.where(use_TIP,1,0)` ≡ **1.0, no isotherm** (nsm1/processes.py:290) | `FPSS = PARTP·TISS/(PARTP·TISS+PARTP·Fe·DO1+1)`; dissolved = `1/(1+PARTP·TISS+PARTP·Fe·DO1)` (water-quality.f90:449) | QUAL2K: inorganic-P sorbed fraction via linear `Kd·m` isotherm | Bowie: linear equilibrium `Kd` partition | linear `1/(1+Kd·SS)` | **V3-CORRECT-V1-WRONG** (see 3.2) |
| `kdpo4` default | 0.0 L/kg | 0.0 (ignored by stub) | `PARTP` control default 0 | QUAL2K default 0; if used ~10^2--10^5 L/kg | Bowie: PO4 Kd 10^2--10^4 L/kg | 100--10000 L/kg | CONSISTENT (zero default; isotherm correct when activated) |
| TIP settling | `vs/depth·(1-fdp)·tip` (phosphorus.py:401) | same | `PO4NS = (SSSI·PO4_{k-1}-SSSO·PO4_k)·BI/BH2` (water-quality.f90:1448) | QUAL2K: sorbed-P settles with solids | WASP: `(1-fdp)` settles | sorbed fraction settles | CONSISTENT (form correct; v1 fdp stub makes it inert -- v3 restores) |
| `vs` (TIP settling vel.) | 0.1 m/d | 999 sentinel (v1 bug) | tied to `SSS`/`FeSetVel` in W2 | QUAL2K sorbed-P settling 0.1--2 m/d | Bowie: 0.05--2 m/d | 0.1--1 m/d | V3-CORRECT-V1-WRONG (v1 999 sentinel is unusable; v3 0.1 defensible) |
| OrgP hydrolysis | `kop_tc·orgp` (phosphorus.py:394) | same | `LPOMPD = ORGP·LPOMD` etc. (water-quality.f90:379) | QUAL2K OrgP→DIP first-order | Bowie: 0.01--0.7 /d | 0.05--0.3 /d | CONSISTENT |
| `kop_20` | 0.1 /d | 0.1 | `LPOMPDK`/`ORGP·LPOMDK` | QUAL2K OrgP hydrolysis 0.02--0.4 /d | Bowie: 0.1--0.7 /d | 0.1--0.3 /d | CONSISTENT |
| OrgP settling | `vsop/depth·orgp` (phosphorus.py:408) | same | OrgP carried in POM settling (W2 has no separate OrgP state) | QUAL2K OrgP settles with detritus | Bowie: 0.1--2 m/d | 0.1--1 m/d | CONSISTENT (no Arrhenius -- correct, matches W2/QUAL2K) |
| `vsop` | 0.1 m/d | 999 sentinel | n/a separate state | QUAL2K 0.1--2 m/d | Bowie 0.1--1 m/d | 0.1--1 m/d | V3-CORRECT-V1-WRONG (999 sentinel → 0.1; documented Phase 9.E) |
| Sediment-P release `rpo4_20` | 0.0 g-P/m²/d | 0.0 | `PO4SR = PO4R·SODD·DO2` (water-quality.f90:1447) | QUAL2K `kdip` sed P flux | Bowie: 0--50 mg/m²/d | site-specific | CONSISTENT (zero default = NSM2 deferral; W2 ties to SOD) |

### Pathogen

| Term / default | v3 | v1 | CE-QUAL-W2 (file:line) | QUAL2K / Mancini | WASP/Bowie | Literature value | VERDICT |
|---|---|---|---|---|---|---|---|
| Natural die-off | `kdx_tc·px`, θ^(T-20) (pathogen.py:356-359) | same | `DK1_BACT = BACT1DK·BACTQ10^(T-20)·BACT·DO3` (water-quality.f90:1215) | Chick's law first-order; Mancini 1978 temp/salinity | Bowie: Chick die-off | first-order Chick | CONSISTENT (W2 adds DO3 oxic gate; NSM1/QUAL2K omit -- defensible, see 3.4) |
| `kdx_20` | 0.8 /d | 0.8 | `BACT1DK` typ. 0.5--2 /d | Mancini freshwater base ~0.8 /d (T,S=0) | Bowie fecal coli 0.5--2.0 /d | 0.5--2.0 /d | CONSISTENT |
| `kdx_theta` | 1.07 | 1.07 | `BACTQ10` (Q10-style, typ ~1.07) | Mancini: 1.07 | Bowie: 1.07 | 1.07 | CONSISTENT |
| Light decay form | `apx·q_solar·(1-e^{-kd})/kd·px` (pathogen.py:398-424) | same | `PHOTO_BACT = BACTLDK·LIGHT·BACT`, `LIGHT` = depth-avg Beer-Lambert (water-quality.f90:1196,1220) | Chapra/Auer&Niehaus depth-averaged irradiance | Bowie/Mancini: light proportional to I0 | depth-avg I0 Beer-Lambert | CONSISTENT (identical depth-averaging structure to W2) |
| `apx` (light efficiency) | 0.017 (W/m²)⁻¹d⁻¹ | 1.0 placeholder | `BACTLDK` calibrated, no canonical default | Auer & Niehaus 1993: α=0.00824 cm²/cal ≈ 0.017 (W/m²)⁻¹d⁻¹ | Chapra 1997 Ch.33: canonical | 0.017 (W/m²)⁻¹d⁻¹ | V3-CORRECT-V1-WRONG (v1 placeholder 1.0; v3 canonical, see 3.4) |
| `vx` (settling vel.) | 1.38 m/d | 1.0 placeholder | `BACTS` calibrated | Auer & Niehaus 1993 sediment-trap 1.38 m/d; Chapra Ch.33 | Bowie: 0.5--2.5 m/d range | 1.38 m/d | V3-CORRECT-V1-WRONG (v1 placeholder 1.0; v3 canonical, see 3.4) |
| Settling form | `vx/depth·px` (pathogen.py:433) | same | `SET_BACT = BACTS·(BACT_{k-1}-BACT_k)·BI/BH2` (water-quality.f90:1200) | velocity/depth first-order | Bowie: settling velocity | velocity/depth | CONSISTENT |

---

## 3. The Five Adjudications

### 3.1 CBOD oxidation: rate, θ, DO half-saturation; and the CBOD-settling Arrhenius defect

The v3/v1 CBOD oxidation form `kbod_tc · DOX/(KsOxbod + DOX) · CBOD` is a Monod DO-limited
first-order decay. Validated-reference comparison:

- **CE-QUAL-W2** (`CBODD = KBOD·TBOD^(T-20)·DO3·CBOD`, water-quality.f90:440 and the
  `CBODSS = -CBODD·CBOD` sink at 1819) applies a *binary* oxic switch `DO3 = (1+sign(1,O2-1e-10))·0.5`
  (water-quality.f90:298): full first-order decay when O2 > ~0, zero when anoxic. W2 does
  *not* use a Monod half-saturation.
- **QUAL2E** (Brown & Barnwell 1987) uses pure first-order `K_1·CBOD` with *no* DO
  limitation at all. **Bowie et al. (1985)** Table compilations of `K_1` (0.05--0.5 /d,
  θ = 1.047) are first-order coefficients.
- **QUAL2K v2.11b8** and **WASP7 EUTRO** introduce the optional Monod attenuation
  `DO/(K_socf + DO)` that v3/v1 use. The half-saturation in QUAL2K-fast (`Ksocf`) is
  typically 0.5--2 mg O2/L.

Ruling: the v3/v1 Monod form is the *most modern and most defensible* of the three families
and is CONSISTENT with QUAL2K/WASP. `kbod_20 = 0.12 /d` and `kbod_theta = 1.047` are squarely
within Bowie/QUAL2E ranges. `KsOxbod = 0.5 mg/L` sits at the low end of the QUAL2K range; it
is defensible but a calibration user should know that at DO = 0.5 mg/L the oxidation rate is
halved, and that this Monod attenuation is *more aggressive at moderate DO* than W2's binary
switch (which keeps full rate down to ~0 mg/L). This is a behavioral difference from W2, not
an error. No joint v1/v3 form defect.

The genuine *defaults-vs-literature* defect that the parity reviews could only flag as a
v1↔Fortran discrepancy is the **CBOD settling Arrhenius coefficient `ksbod_theta = 1.047`**.
The canonical settling-coefficient temperature correction in Bowie et al. (1985) and QUAL2E
is **θ = 1.024**, distinct from the 1.047 used for *oxidation*. CE-QUAL-W2 applies **no
Arrhenius at all** to `CBODS` settling (water-quality.f90:1817 is a bare velocity term). So
v3/v1 are jointly wrong against the literature here: they apply an oxidation-θ to a settling
process. This is **V1&V3-BOTH-WRONG**, fully latent at the shipped `ksbod_20 = 0` default,
and already partly noted in `parameters/cbod.py:50` as a "follow-up." This review elevates
it from a parity footnote to a confirmed science-correctness defect against Bowie/QUAL2E:
if `ksbod_20 > 0` is calibrated, the settling rate will carry a 2x-too-steep temperature
sensitivity (1.047 vs 1.024 over a 30 °C span is roughly a 1.9x vs 1.4x multiplier). The
v3 velocity-form (`ksbod_tc/depth`) is, separately, *closer* to W2's velocity convention
than v1's rate-form, so on the units axis v3 is the better choice; but the θ value is the
deeper science error and it is shared by both versions.

### 3.2 Phosphorus sorption `fdp`: the headline ruling (V3-CORRECT-V1-WRONG)

This is the single most consequential science adjudication in scope. The question posed:
is v3's `fdp = 1/(1 + kdpo4·Solid·1e-6)` the *correct* equilibrium sorption isotherm, is
the `1e-6` unit factor and default `kdpo4` dimensionally/numerically correct, and is the
v1 stub a v1 bug that v3 fixed?

**The correct isotherm.** Linear-equilibrium (Kd) partitioning of orthophosphate between
dissolved and solid-sorbed phases gives a dissolved fraction

  f_d = 1 / (1 + Kd · m)

where Kd is the solid-water partition coefficient and m is the particulate concentration,
with Kd·m dimensionless. This is the textbook form (Chapra 1997; Thomann & Mueller 1987;
Bowie et al. 1985 §3 on adsorption) and is exactly what QUAL2K uses for its inorganic-P
sorbed fraction.

**The validated reference (CE-QUAL-W2).** W2 implements precisely this isotherm, generalized
to two sorbents (suspended solids and ferric iron). At water-quality.f90:449:

  FPSS = PARTP·TISS / (PARTP·TISS + PARTP·Fe·DO1 + 1.0)

so the W2 *dissolved* fraction is `1/(1 + PARTP·TISS + PARTP·Fe·DO1)`, and dropping the Fe
term (NSM1 has no Fe sorption) this is `1/(1 + PARTP·TISS)` -- structurally identical to
v3's `1/(1 + kdpo4·Solid·1e-6)`. W2's `PARTP·TISS` is dimensionless because `PARTP` is in
m³/g and `TISS` is g/m³. The W2 sorbed-P then settles via `PO4NS` (water-quality.f90:1448),
exactly analogous to v3's `vs/depth·(1-fdp)·tip` (phosphorus.py:401). The structural and
dimensional correspondence is exact. **v3's isotherm is the validated form.**

**Dimensional/numerical check of the `1e-6` factor.** v3 takes `kdpo4` in L/kg and `Solid`
in mg/L. The product `kdpo4·Solid` has units (L/kg)(mg/L) = mg/kg. To render it the
dimensionless mass ratio required by the isotherm, multiply by `1e-6 kg/mg`. So
`kdpo4·Solid·1e-6` is dimensionless and correct. This is the same magnitude as W2's
`PARTP[m³/g]·TISS[g/m³]` once unit systems are reconciled (L/kg = 1e-3 m³/kg = 1e-6 m³/g,
mg/L = g/m³, so `kdpo4[L/kg]·Solid[mg/L]·1e-6 = kdpo4[m³/g equiv]·Solid[g/m³]`). The factor
is dimensionally and numerically correct, and matches the Fortran NSM1 `modGlobalParam.f90:228`
`/1.0E6` direction. The v1 NSM1 *shared-library* copy with `/0.000001` (i.e. ×1e6) is
inverted and would drive `fdp → 0` at trivially small `kdpo4·Solid`; v3 corrected the
direction.

**Is the v1 stub a v1 bug?** Yes. The v1 NSM1 model *as wired* (`nsm1/processes.py:290`,
referenced by `dynamic_variables.py:99`) computes `fdp = xr.where(use_TIP, 1, 0)` -- a
constant 1.0 that ignores `kdpo4` and `Solid` entirely. With `fdp ≡ 1.0`, the particulate
fraction `(1 - fdp) = 0`, so TIP sorbed-settling is identically zero *for all parameter
values*, including any nonzero `kdpo4` a user might calibrate. This is a degenerate stub,
not a model: it silently disables a process the model claims to support. Against the
validated W2 `PARTP` isotherm and QUAL2K, this is a **v1 defect**. v3 replaces it with the
correct linear-equilibrium isotherm with the correct unit factor.

**Ruling: V3-CORRECT-V1-WRONG, confirmed against CE-QUAL-W2 (`PARTP`, water-quality.f90:449)
and QUAL2K.** The v3 `fdp` form, the `1e-6` factor, and the `kdpo4 = 0` default are all
correct. The only residual issue is documentation (the prior parity review's F1/F2 already
capture the stale audit text and the need to declare the deliberate departure from v1
runtime behavior). At the shipped `kdpo4 = 0` default both reduce to `fdp = 1.0` so v1↔v3
benchmark parity holds; the divergence is latent and v3 is the scientifically correct side
of it. This review additionally records that the v1 stub is a *true science defect against
a validated reference*, not merely a v3-vs-v1 parity gap -- the parity reviews could not
make that ruling because their reference was v1 itself.

### 3.3 POM settling and hydrolysis rates and θ

v3/v1 POM dissolution `kpom_tc·POM` with `kpom_20 = 0.1 /d`, `kpom_theta = 1.047` is a
first-order hydrolysis consistent with CE-QUAL-W2's `LPOMD = OMTRM·LPOMDK·LPOM·DO3`
(water-quality.f90:354) and with QUAL2K detritus dissolution and Bowie et al. (1985) POM
hydrolysis ranges (0.03--0.25 /d labile; θ = 1.047). The structural difference -- NSM1
carries a single lumped POM as Di Toro layer-2 *bed* sediment with a *burial* sink
(`vb·POM/h2`), whereas W2 carries water-column LPOM/RPOM with a *settling* velocity
(`POMS·ΔLPOM·BI/BH2`, water-quality.f90:1709) -- is a documented STRUCTURAL-OK difference,
not an error. The comparable kinetic forms (first-order temperature-corrected hydrolysis;
velocity-based vertical transport) and the coefficients (`kpom_20 = 0.1`, `vsoc = 0.01 m/d`,
`fcom = 0.4`, `h2 = 0.1 m` Di Toro convention, `vb = 6.85e-6 m/d` Di Toro burial) are all
within validated ranges. The v1 `vb` unit bug (0.01 m/d, off by ~1460x from the Di Toro
0.0025 m/yr) was corrected in v3 to 6.85e-6 m/d; this is a V3-CORRECT-V1-WRONG on the burial
coefficient, already documented and regression-pinned. No joint form defect; defaults
defensible.

### 3.4 Pathogen die-off: form and `apx`/`vx`/`kdx` defaults

The v3/v1 pathogen budget is `dPX/dt = -(kdx_tc·PX + apx·I0·[depth-avg Beer-Lambert]·PX
+ vx/depth·PX)`: Chick first-order natural die-off with Arrhenius θ^(T-20), light-induced
inactivation proportional to depth-averaged irradiance, and settling. This is the canonical
Mancini (1978) / Chapra (1997, Ch. 33) / QUAL2K structure. CE-QUAL-W2 `BACTERIA`
(water-quality.f90:1183--1224) implements the same three terms: `DK1_BACT = BACT1DK·BACTQ10^(T-20)·BACT·DO3`
(first-order with Q10 temperature and an oxic `DO3` gate), `PHOTO_BACT = BACTLDK·LIGHT·BACT`
with `LIGHT` the depth-averaged Beer-Lambert irradiance `LAM1·(1-exp(-γH))/(γH)`
(water-quality.f90:1196), and `SET_BACT` settling. The v3 depth-averaging
`(1 - e^{-kd})/kd` (pathogen.py:421) is *algebraically identical* to W2's
`(1-exp(-γH))/(γH)`. The only W2 addition is the `DO3` oxic gate on natural die-off, which
NSM1/QUAL2K/Mancini omit; omitting it is the standard pathogen-modeling choice (die-off is
not redox-limited in the Mancini formulation), so this is CONSISTENT, not a defect.

Defaults. `kdx_20 = 0.8 /d`, `kdx_theta = 1.07` are the Mancini (1978) freshwater base
coefficients and within Bowie's 0.5--2.0 /d fecal-coliform range. The parity round's claim
that v3's `apx` and `vx` are *improvements over v1 placeholders* is **confirmed against the
canonical literature**: v1 ships `apx = 1.0` and `vx = 1.0` as dimensionless placeholders.
v3 corrects `apx = 0.017 (W/m²)⁻¹d⁻¹` (Auer & Niehaus 1993, α = 0.00824 cm²/cal in cgs,
the canonical value reproduced in Chapra 1997 Ch. 33 and QUAL2K v2.11b8 §5.5.20.1) and
`vx = 1.38 m/d` (Auer & Niehaus 1993 Onondaga Lake sediment-trap measurement, within
Bowie's 0.5--2.5 m/d range). These are **V3-CORRECT-V1-WRONG**: the v1 placeholders are
not physical (a dimensionless `apx = 1.0` against a W/m² irradiance gives a nonsensical
rate magnitude; `vx = 1.0` is an unsourced placeholder). v3's values are the canonical
literature calibration and are correctly tied to the broadband `q_solar` after the
Phase 9.F.B revert of the PAR substitution -- which this review confirms is the correct
choice, since Auer & Niehaus's α is defined against broadband irradiance and pathogen
inactivation is UV-mediated, not PAR-mediated. Ruling: pathogen form CONSISTENT with W2
and Mancini; `apx`/`vx` corrections are confirmed literature-grounded improvements.

### 3.5 CBOD→DIC and CBOD→O2 stoichiometry consistency

NSM1 carries CBOD in oxygen-equivalent units (mg-O2/L), the WASP/QUAL2E convention. The
CBOD oxidation sink therefore consumes DO 1:1 with the CBOD decayed (the oxygen-demand
*is* the state), so the effective CBOD→O2 stoichiometry is roc = 1.0 in O2-equivalent
basis -- consistent with Bowie et al. (1985) and QUAL2E when CBOD is expressed as
ultimate oxygen demand. CE-QUAL-W2 carries a per-group `RBOD(JCB)` multiplier on the DO
sink (`DOBOD = RBOD·CBODD·CBOD`, water-quality.f90:1868) to allow CBOD expressed on a
mass-substrate basis; with CBOD already in O2-equivalent, `RBOD = 1` reproduces the NSM1
behavior. The CBOD→DIC path in W2 is `TICBOD = BODC·CBODD·CBOD` (water-quality.f90:1942),
the carbon released per unit O2-demand oxidized. v3 NSM1 routes the CBOD→DIC carbon return
through the Carbon/DIC process via the OM C:O stoichiometry (out of scope for the CBOD
class itself, correctly deferred -- prior parity review §6). The two stoichiometric
constants are reciprocals of the OM C:O mass ratio (≈0.375 gC/gO2 from 12/32) and are
applied consistently in the respective consumer processes. No stoichiometric inconsistency
or double-counting was found in the CBOD class; the cross-process contract is CONSISTENT
with W2's `RBOD`/`BODC` decomposition.

---

## 4. Science Findings the Parity Reviews Necessarily Missed

The two prior reviews are parity-scoped: their reference is v1, so by construction they
cannot find a defect that v1 and v3 *share*. The following are science-correctness items
that fall in that blind spot. None is a v3 implementation defect; all are inherited shared
v1/v3 issues, surfaced here against validated references.

1. **SF1 (Major, latent): CBOD settling Arrhenius θ is the wrong coefficient.**
   `ksbod_theta = 1.047` (parameters/cbod.py:50) applies an oxidation-process temperature
   coefficient to a *settling* process. Bowie et al. (1985) and QUAL2E specify θ = 1.024
   for settling coefficients; CE-QUAL-W2 applies no Arrhenius to `CBODS` at all. Shared by
   v1 and v3. Latent at the shipped `ksbod_20 = 0`. Consequence: any calibration with
   `ksbod_20 > 0` carries a settling temperature sensitivity that is too steep by the
   1.047/1.024 ratio (compounding to ~1.4x error in the rate multiplier over a 30 °C span).
   Recommendation: set `ksbod_theta = 1.024` (the Bowie/QUAL2E settling value) or document
   the deviation as a deliberate, calibrated choice with a literature basis. This is the
   highest-value finding in this review because it is a *form-coefficient error against a
   validated reference that both versions share*, exactly the class the task prioritized.

2. **SF2 (Observation): CBOD Monod DO attenuation is more aggressive than CE-QUAL-W2's
   binary switch at moderate DO.** Not a defect -- the v3/v1 Monod form is the modern
   QUAL2K/WASP choice and is the more defensible of the available formulations -- but a
   calibration user porting a W2 application should know that, at DO in the 0.5--3 mg/L
   range, NSM1 will predict slower CBOD oxidation than W2 (which keeps full first-order
   rate until DO ≈ 0). With `KsOxbod = 0.5 mg/L` the half-rate point is 0.5 mg/L; raising
   it toward the QUAL2K upper bound (2 mg/L) would widen the divergence. Recommendation:
   document the mechanistic difference in the CBOD module docstring so cross-model
   calibration is not silently miscompared.

3. **SF3 (Observation): NSM1 lumped CBOD/OrgP vs W2 split labile/refractory DOM+POM is
   STRUCTURAL-OK but the single `kpom_20 = 0.1 /d` implicitly assumes a labile-dominated
   pool.** W2 distinguishes `LPOMDK` (labile, ~0.08--0.12 /d) from `RPOMDK` (refractory,
   ~0.001--0.01 /d). NSM1's single `kpom_20 = 0.1 /d` is the labile value; applying it to
   a pool that is in reality refractory-dominated would over-predict hydrolysis. This is
   an inherent limitation of the NSM1 lumped structure, not a code defect, and is the
   correct design per the NSM1 specification. Recorded so a calibration user does not treat
   `kpom_20 = 0.1` as universal.

---

## 5. Defaults vs Bowie / Literature Reconciliation Table

| Parameter | v3 default | Bowie et al. 1985 range | QUAL2E/QUAL2K | Validated W2 | Disposition |
|---|---|---|---|---|---|
| `kbod_20` | 0.12 /d | 0.05--0.5 /d (K_1 carbonaceous) | QUAL2E K_1 0.02--3.4, typ 0.1--0.5 | KBOD control input | Within range. OK. |
| `kbod_theta` | 1.047 | 1.047 (BOD oxidation) | 1.047 | TBOD ~1.047 | Canonical. OK. |
| `KsOxbod` | 0.5 mg-O2/L | DO half-sat 0.5--2 mg/L (WASP) | QUAL2K Ksocf 0.5--2 | n/a (binary) | Low end; defensible. SF2. |
| `ksbod_20` | 0.0 /d (or m/d) | 0 (modern dissolved-CBOD) | QUAL2E K_3 default 0 | CBODS default 0 | Intentional zero. OK. |
| `ksbod_theta` | **1.047** | **1.024 (settling)** | **QUAL2E 1.024 (settling)** | **no Arrhenius on CBODS** | **WRONG coefficient. SF1.** |
| `kpom_20` | 0.1 /d | 0.03--0.25 /d (POM labile) | QUAL2K detritus 0.01--0.25 | LPOMDK ~0.08--0.12 | Labile value; within range. SF3. |
| `kpom_theta` | 1.047 | 1.047 | 1.047 | OMTRM temp fn | Canonical. OK. |
| `vsoc` (POC→POM) | 0.01 m/d | 0.1--2.5 m/d (POM settling) | QUAL2K detritus settle | POMS control | Low end; defensible for fine POC. OK. |
| `vb` (burial) | 6.85e-6 m/d | Di Toro 0.25 cm/yr | Di Toro/QUAL2K w2 | n/a | = 0.0025 m/yr Di Toro. OK (v1 bug fixed). |
| `kdpo4` | 0.0 L/kg | PO4 Kd 10²--10⁴ L/kg if active | QUAL2K Kd | PARTP control 0 | Zero default; isotherm correct when set. OK. |
| `vs` (TIP settle) | 0.1 m/d | 0.05--2 m/d | QUAL2K sorbed-P 0.1--2 | tied to SSS | Within range (v1 999 sentinel fixed). OK. |
| `kop_20` (OrgP) | 0.1 /d | 0.1--0.7 /d | QUAL2K 0.02--0.4 | ORGP·LPOMDK | Within range. OK. |
| `vsop` (OrgP settle) | 0.1 m/d | 0.1--1 m/d | QUAL2K 0.1--2 | n/a sep state | Within range (v1 999 sentinel fixed). OK. |
| `rpo4_20` (sed P) | 0.0 g-P/m²/d | 0--50 mg/m²/d | QUAL2K kdip | PO4R·SODD | Zero default = NSM2 deferral. OK. |
| `kdx_20` (pathogen) | 0.8 /d | 0.5--2.0 /d fecal coli | Mancini base ~0.8 | BACT1DK | Mancini freshwater base. OK. |
| `kdx_theta` | 1.07 | 1.07 | Mancini 1.07 | BACTQ10 ~1.07 | Canonical. OK. |
| `apx` | 0.017 (W/m²)⁻¹d⁻¹ | (via Auer&Niehaus α) | Chapra Ch.33 / QUAL2K | BACTLDK calibrated | Auer&Niehaus 1993 canonical (v1 placeholder 1.0 fixed). OK. |
| `vx` | 1.38 m/d | 0.5--2.5 m/d | Chapra Ch.33 | BACTS calibrated | Auer&Niehaus 1993 (v1 placeholder 1.0 fixed). OK. |

Single defaults-vs-literature defect: `ksbod_theta` (SF1). All other defaults are within
Bowie/QUAL2E/QUAL2K validated ranges or are intentional gated zeros.

---

## 6. Summary of Verdicts by Adjudication

1. CBOD oxidation form and `kbod`/θ/`KsOxbod`: **CONSISTENT** with QUAL2K/WASP (Monod is
   the modern defensible form; more aggressive than W2 binary switch at moderate DO -- SF2).
2. Phosphorus `fdp` sorption: **V3-CORRECT-V1-WRONG**, confirmed against CE-QUAL-W2 `PARTP`
   (water-quality.f90:449) and QUAL2K linear-equilibrium isotherm. v3 form, `1e-6` unit
   factor, and `kdpo4 = 0` default are all dimensionally and numerically correct. The v1
   NSM1 runtime stub (`fdp ≡ 1.0`) is a genuine science defect against a validated
   reference, not merely a parity gap.
3. POM settling/hydrolysis/θ: **STRUCTURAL-OK / CONSISTENT**; v1 `vb` unit bug is
   V3-CORRECT-V1-WRONG and already fixed.
4. Pathogen die-off form and `apx`/`vx`/`kdx`: form **CONSISTENT** with W2 `BACTERIA` and
   Mancini 1978; `apx`/`vx` corrections **V3-CORRECT-V1-WRONG** confirmed against Auer &
   Niehaus 1993 / Chapra 1997 canonical values. The parity round's improvement claim is
   confirmed.
5. CBOD→DIC and CBOD→O2 stoichiometry: **CONSISTENT** with W2 `RBOD`/`BODC` decomposition;
   no double-counting; cross-process contract sound.

Net: one shared v1/v3 form-coefficient defect against validated references (`ksbod_theta`,
SF1, Major-latent), one confirmed V3-CORRECT-V1-WRONG science fix (phosphorus `fdp`),
multiple confirmed V3-CORRECT-V1-WRONG default corrections (`apx`, `vx`, `vsop`, `vs`, `vb`),
zero V1&V3-BOTH-WRONG defects that are active at shipped defaults.
