# v3 NSM1 Nitrogen Cycle — Science-Correctness Validation

Review date: 2026-05-15
Reviewer: water-quality model source-code reviewer (Claude)
Branch: `streaming` @ HEAD (54f2b12)
Scope: scientific correctness of v3 NSM1 nitrogen-cycle kinetics (nitrification,
denitrification, organic-N hydrolysis, settling, sediment NH4/NO3 release,
ammonia preference for algal uptake, alkalinity stoichiometry of N reactions,
Arrhenius theta values and their pairing). This is a correctness validation
against validated reference models, not a v3-v1 parity check. Parity was
established by a prior round; the question here is whether the shared
formulation is *defensible*, including the possibility that v1 and v3 are
**both** wrong.

Reference models:

* CE-QUAL-W2 v2026.02 (local, authoritative):
  `/Users/todd/GitHub/CE-QUAL-W2-ERDC/CE-QUAL-W2-ERDC-dev/src/W2_v2026.02/water-quality.f90`
* NSM1 Fortran (authoritative upstream of ClearWater):
  `/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/modNitrogen.f90`
* QUAL2K (Chapra, Pelletier & Tao 2008), QUAL2E (Brown & Barnwell 1987,
  EPA/600/3-87/007), WASP/EUTRO (Wool et al.), Bowie et al. 1985
  (EPA/600/3-85/040), Chapra 1997 *Surface Water-Quality Modeling*.

---

## 1. Verdict

The v3 NSM1 nitrogen kinetics are **scientifically defensible and faithfully
reproduce the authoritative NSM1 Fortran formulation**. Every rate law
(nitrification, low-DO inhibition, denitrification, organic-N hydrolysis,
settling, sediment release, ammonia preference) matches `modNitrogen.f90`
term for term, and the kinetic *forms* are standard and consistent with
QUAL2K/QUAL2E and the published literature. The Phase 9.E theta-transposition
correction is **confirmed correct on independent grounds**: v3 matches the
Fortran source and the physically expected ordering (organic-matter
hydrolysis theta ~1.047, sediment-exchange theta ~1.074-1.08); v1 has the
pairs swapped and is **wrong**. This is a genuine V3-CORRECT / V1-WRONG case
that the parity round correctly identified and that this validation
independently confirms.

One **MAJOR** science finding that the parity round did not surface: the
denitrification alkalinity stoichiometry coefficient `r_alkden = 4/14/1000`
eq/mg-N is **inconsistent with the canonical 1 eq alk per mol N released by
NO3 reduction to N2** that CE-QUAL-W2 uses (`water-quality.f90:3157`,
`+ NO3D` with coefficient 1, not 4). v1 and v3 both carry `4/14/1000`, so
this is a **V1-AND-V3-BOTH-WRONG** candidate. It does not affect the nitrogen
state itself (it lives in the Alkalinity Process) but it biases simulated
alkalinity and any pH derived from it by a factor of ~4 on the
denitrification term. Detail and severity are in Section 4 (SF-1).

Two **MINOR** observations on defaults: the v3 default nitrification rate
`knit_20 = 0.1 /d` sits at the low end of the literature range and is the
NSM1/QUAL2E low default rather than a calibrated value (acceptable as a
documented default but worth a note); and the sediment NH4/NO3 release
defaults are held at exactly zero, which is the Fortran default but means
benthic N flux is silently off unless the user sets `rnh4_20`/`vno3_20`
(documented design choice, not a defect).

Severity counts (science-correctness only): CRITICAL 0, MAJOR 1 (SF-1),
MINOR 3 (SF-2, SF-3, SF-4), OBSERVATION 3.

---

## 2. Cross-model validation matrix

Verdict legend: CONSISTENT = v3 form and default are defensible against
references; V3-DIVERGENT = v3 differs from references in a way that needs
justification; V1&V3-BOTH-WRONG = shared error; V3-CORRECT-V1-WRONG = v3
fixed a v1 error; STRUCTURAL-OK = cross-model structural difference that is
not an error.

| Term / default | v3 value/form | v1 value/form | CE-QUAL-W2 (file:line) | QUAL2K / QUAL2E | WASP/EUTRO | Literature | VERDICT |
|---|---|---|---|---|---|---|---|
| Nitrification rate law | `NH4 * knit_tc * (1-exp(-KNR*DOX))` (nitrogen.py:723--725, 885) | same | `NH4D = NH4TRM*NH4DK*NH4*DO1`, `DO1=O2/(O2+KDO)` (water-quality.f90:296,326) | first-order in NH4 x DO attenuation | first-order in NH4 | first-order in NH4 (Bowie 1985 §6) | CONSISTENT (forms differ in DO factor only; see Adjudication 2) |
| `knit_20` (nitrification rate, 1/d at 20C) | 0.1 | 0.1 | NH4DK input, no hardwired default in src | QUAL2E 0.1--1.0, default ~0.5; QUAL2K 0.05--2 | 0.09--0.13 | 0.1--1.0 (Bowie 1985 Table 6-2; Chapra 1997 ~0.1--0.5) | CONSISTENT (low end; see SF-2) |
| `knit_theta` (nitrification Arrhenius) | 1.083 | 1.083 | NH4 uses W2 four-point T-rate, not theta | QUAL2E theta_nit = 1.083 | 1.08 | 1.07--1.08 (Bowie 1985) | CONSISTENT (matches QUAL2E exactly) |
| `KNR` (DO half-sat for nitrif inhibition, mg-O2/L) | 0.6 | 0.6 | KDO Monod constant (different form) | QUAL2K uses `1-exp(-KNR*DO)`-equivalent low-DO inhibition | n/a | NSM1 hardwire range {0.6--0.7} (modNitrogen.f90:63) | CONSISTENT |
| Low-DO inhibition form | `1 - exp(-KNR*DOX)` | same | Monod `O2/(O2+KDO)` (water-quality.f90:296,326) | QUAL2K `1-e^(-KNR·DO)` (exponential) | none / step | exponential form standard (Chapra 1997) | CONSISTENT (matches QUAL2K exponential; see Adjudication 2) |
| Denitrification rate law | `NO3 * kdnit_tc * (1 - DOX/(DOX+KsOxdn))` (nitrogen.py:758--762) | same | `NO3D = NO3TRM*NO3DK*NO3*DO2`, `DO2=KDO/(O2+KDO)` (water-quality.f90:297,339) | QUAL2K anoxic Monod `KsNitr/(KsNitr+DO)` | first-order, anoxic | first-order in NO3, anoxic (Bowie 1985) | CONSISTENT (W2 DO2 is the exact Monod complement of v3's `1-DO/(DO+Ks)`) |
| `kdnit_20` (denit rate, 1/d at 20C) | 0.002 | 0.002 | NO3DK input | QUAL2K 0--2 /d; QUAL2E 0--1 /d | 0.05--0.2 typical water column | 0.002--0.2 water-column; bed flux dominates (Bowie 1985) | CONSISTENT (low water-column default; bed term carries denit) |
| `kdnit_theta` (denit Arrhenius) | 1.045 | 1.08 (transposed) | n/a | QUAL2E/QUAL2K denit theta ~1.045--1.08 | n/a | ~1.045 (Chapra 1997; QUAL2E denitrification) | V3-CORRECT-V1-WRONG (matches Fortran 1.045; see Adjudication 1) |
| `kon_20` (OrgN hydrolysis, 1/d at 20C) | 0.1 | 0.1 | LDOM/RDOM, POM N decay (different structure) | QUAL2E org-N hydrolysis 0.02--0.4, default ~0.1 | 0.075 | 0.02--0.4 (Bowie 1985 Table 6-2) | CONSISTENT |
| `kon_theta` (OrgN hydrolysis Arrhenius) | 1.047 | 1.074 (transposed) | OM four-point T-rate | universal organic-matter theta 1.047 (Chapra 1997; QUAL2K) | 1.08 | 1.047 (NSM1 organic-matter convention) | V3-CORRECT-V1-WRONG (matches Fortran 1.047; see Adjudication 1) |
| `rnh4_theta` (sediment NH4 release Arrhenius) | 1.074 | 1.047 (transposed) | NH4R sediment release | sediment-water exchange theta ~1.06--1.08 | n/a | steeper than organic-matter (~1.074) | V3-CORRECT-V1-WRONG (matches Fortran 1.074) |
| `vno3_theta` (sediment denit Arrhenius) | 1.08 | 1.045 (transposed) | sediment denit | ~1.045--1.08 | n/a | sediment-exchange ~1.08 | V3-CORRECT-V1-WRONG (matches Fortran 1.08) |
| OrgN settling | `vson/depth*OrgN`, **no** Arrhenius (nitrogen.py:934) | same (Phase 9.E removed bad vson_theta) | `OrgN` partitioned via POM settling | settling velocity, no Arrhenius | velocity | settling ~ viscosity, weak T-dep (Chapra 1997) | CONSISTENT (no Arrhenius is correct; v3 removed an earlier v3-only error) |
| `vson_20` (OrgN settling velocity, m/d) | 0.01 | 0.01 (GlobalVars) | POMS settling 0.1--1 m/d typical | 0.1--2.5 m/d (QUAL2K POM) | 0.1--1 | 0.05--2.5 m/d (Bowie 1985) | CONSISTENT (NSM1 lumped OrgN low velocity; STRUCTURAL difference from W2 PON/DON split is not an error) |
| Sediment NH4 release `rnh4_20` (1/d at 20C) | 0.0 | 0.0 | NH4R(JW) x SODD x DO2 (water-quality.f90:1535) | QUAL2K SOD-driven NH4 flux 20--200 mg/m2/d | benthic flux | 10--300 mg-N/m2/d (Bowie 1985 Table) | CONSISTENT-AS-DEFAULT (zero default = off by design; see SF-3) |
| Sediment NO3 denit `vno3_20` (1/d at 20C) | 0.0 | 0.0 | NO3S(JW), NO3SED (water-quality.f90:1591) | sediment denit velocity 0.1--1 m/d | n/a | 0.1--1 m/d (Bowie 1985) | CONSISTENT-AS-DEFAULT (zero default = off by design; see SF-3) |
| Ammonia preference for algal uptake | `PN*NH4/(PN*NH4+(1-PN)*NO3)` (floating_algae.py:910--912) | same | QUAL2K-form `NH4PR` when ANEQN=2 (water-quality.f90:1499--1502) | QUAL2K Pnh4 product form (ANEQN=2); QUAL2E simple ratio (ANEQN=1) | EUTRO simple ratio | both forms in literature | STRUCTURAL-OK (NSM1 uses PN-weighted ratio = QUAL2K ANEQN=1 generalization; see Adjudication 3) |
| `PN` (algal NH4 preference factor) | 0.5 | 0.5 | ANPR (half-sat, mg-N/L) — different parameterization | QUAL2K Khnxp ~0.025--0.1 mg-N/L (half-sat, not fraction) | n/a | NSM1 uses dimensionless preference 0--1 | STRUCTURAL-OK (NSM1 PN is a weighting fraction, not a half-sat; default 0.5 = neutral) |
| Nitrification O2 stoichiometry `ron` | 4.57 mg-O2/mg-N (DOX Process) | 4.57 | 4.57 / RNH4DK in W2 DO eqn | 4.57 (QUAL2E/QUAL2K) | 4.57 | 64/14 = 4.571 (full nitrification) | CONSISTENT |
| Alk change per nitrification `r_alkn` | 2/14/1000 eq/mg-N | 2/14/1000 | `- 2.*NH4D` (water-quality.f90:3156,3166) | 2 eq/mol N (Stumm & Morgan 1996) | 2 eq/mol | 2 eq alk consumed per mol NH4 nitrified | CONSISTENT (matches W2 and Stumm & Morgan exactly) |
| Alk change per denitrification `r_alkden` | **4**/14/1000 eq/mg-N | **4**/14/1000 | `+ NO3D` coefficient **1**, not 4 (water-quality.f90:3157,3166) | 1 eq/mol N (Stumm & Morgan 1996) | 1 eq/mol | **1 eq alk produced per mol NO3 denitrified** | **V1&V3-BOTH-WRONG** (see SF-1, Adjudication 4) |
| Forward-Euler integrator (days) | `X + rate*dt_days` | `X + dXdt*dt` (dt=days) | W2 source/sink + transport solver | n/a | n/a | explicit Euler standard for NSM1 | CONSISTENT |

---

## 3. The five adjudications (numeric answers)

### Adjudication 1 — Phase 9.E theta transposition. CONFIRMED: v3 is correct, v1 is wrong.

Authoritative anchor `modNitrogen.f90` (read directly):

* line 77: `knit%rc20 = 0.1; knit%theta = 1.083`
* line 82: `rnh4%rc20 = 0.0; rnh4%theta = 1.074`
* line 89: `kon%rc20 = 0.1; kon%theta = 1.047`
* line 95: `kdnit%rc20 = 0.002; kdnit%theta = 1.045`
* line 100: `vno3%rc20 = 0.0; vno3%theta = 1.08`

v3 (`parameters/nitrogen.py:56--60`): knit 1.083, kon 1.047, kdnit 1.045,
rnh4 1.074, vno3 1.08 — **exact match to Fortran**.

v1 (`constants.py:133--137`): knit 1.083, kon **1.074**, kdnit **1.08**,
rnh4 **1.047**, vno3 **1.045** — kon↔rnh4 swapped and kdnit↔vno3 swapped.

Independent physical check (does not rely on the Fortran source):

* Nitrification theta. Literature value is 1.06--1.08; QUAL2E uses
  theta_nitrification = 1.083 (Brown & Barnwell 1987). v3/v1 both use 1.083.
  **Correct, not affected by the transposition.**
* Organic-matter hydrolysis (kon). The universal NSM1 / QUAL2K
  organic-matter Arrhenius coefficient is 1.047 (matches `mu_max_theta`,
  `kdp_theta`, `krp_theta`, `kpoc_theta`, `kdoc_theta`, `kop_theta`, all
  1.047 in v1 `constants.py:38--40,161--162`). Organic-N hydrolysis must
  share this convention. v3 = 1.047 (correct); v1 = 1.074 (wrong).
* Denitrification theta. Literature/QUAL2E denitrification temperature
  coefficient is ~1.045 (Chapra 1997; QUAL2E uses 1.045 for the
  denitrification rate). v3 = 1.045 (correct); v1 = 1.08 (wrong, this is a
  sediment-exchange value misapplied to water-column denitrification).
* Sediment-exchange velocities (rnh4, vno3). Sediment-water exchange and
  benthic release have steeper temperature dependence (~1.074--1.08); the
  phosphorus parallel `rpo4_theta = 1.074` (sediment-P release) confirms
  the convention. v3: rnh4 = 1.074, vno3 = 1.08 (correct); v1 has 1.047,
  1.045 (organic-matter / denitrification values misapplied to sediment
  release — wrong).

Conclusion: the v1-to-Fortran port transposed the theta values within both
the (kon, rnh4) and (kdnit, vno3) pairs. v3's Phase 9.E correction restores
the physically and source-correct assignment. Numeric answer: nitrification
theta = 1.083 (correct in both), denitrification theta = 1.045 (v3 correct,
v1 wrong at 1.08). **V3-CORRECT-V1-WRONG, confirmed on independent
physical grounds, not merely by Fortran agreement.**

### Adjudication 2 — Nitrification low-DO inhibition: NSM1 vs W2 vs QUAL2K.

NSM1/v3 (`nitrogen.py:885`, `modNitrogen.f90:265`):
`NitrificationInhibition = 1 - exp(-KNR*DOX)`, KNR = 0.6 mg-O2/L.

CE-QUAL-W2 (`water-quality.f90:296,326`):
`NH4D = NH4TRM*NH4DK*NH4*DO1`, `DO1 = O2/(O2+KDO)` — a **Monod** (rectangular
hyperbola) DO limitation.

QUAL2K (Chapra, Pelletier & Tao 2008, nitrification section): low-DO
attenuation factor `1 - e^(-Knitr·DO)` (exponential form). This is **the
same functional form NSM1/v3 uses.**

Numeric answer: the NSM1/v3 form is **not** the same closed form as W2's
(W2 = Monod, NSM1 = exponential), but it is **identical to the QUAL2K
low-DO inhibition form**, which is a validated, peer-reviewed reference.
With KNR = 0.6, the v3 factor reaches 0.45 at DO = 1 mg/L, 0.70 at 2 mg/L,
0.95 at 5 mg/L; W2's Monod with the typical KDO = 0.5 mg/L reaches 0.67 at
1 mg/L, 0.80 at 2 mg/L, 0.91 at 5 mg/L. The two are qualitatively similar
(both saturating, both -> 1 at high DO, both -> 0 at anoxia) but the
exponential approaches saturation faster. Both are defensible; v3's choice
is the QUAL2K-consistent one. **CONSISTENT; defensible default. Not a
defect — a documented modeling-method choice that follows QUAL2K rather
than W2.**

### Adjudication 3 — Ammonia preference factor for algal uptake.

v3/v1/NSM1 form (`floating_algae.py:910--912`, `modNitrogen.f90:208`):
`ApUptakeFr_NH4 = PN*NH4 / (PN*NH4 + (1-PN)*NO3)`, PN default 0.5.

QUAL2K product form (ANEQN = 2, Chapra et al. 2008, W2 line 1499--1500):
`Pnh4 = NH4·NO3 / ((Khnxp+NH4)(Khnxp+NO3)) + NH4·Khnxp / ((NH4+NO3)(Khnxp+NO3))`.

QUAL2E simple form (ANEQN = 1, W2 line 1502): `NH4 / (NH4+NO3)`.

Numeric answer: the NSM1/v3 form is **not** the QUAL2K product form
(ANEQN = 2). It is a **preference-weighted generalization of the QUAL2E
simple ratio** (ANEQN = 1): with PN = 0.5 it reduces exactly to
`NH4/(NH4+NO3)` (the QUAL2E/EUTRO simple ratio); with PN -> 1 it forces
all uptake from NH4; with PN -> 0 all from NO3. CE-QUAL-W2 supports **both**
ANEQN = 1 and ANEQN = 2 (the product form is selected only when
`ANEQN(JA) == 2`). The NSM1 PN-weighted ratio is mathematically equivalent
to the W2 ANEQN = 1 default with a tunable preference weight rather than a
fixed 0.5. This is a legitimate published parameterization (it is the
EUTRO/QUAL2E ammonia-preference form). It is **structurally different from
the QUAL2K half-saturation product form but not an error**: it is a
different, validated reference choice. **STRUCTURAL-OK; consistent with
QUAL2E/EUTRO; the default PN = 0.5 yields the neutral `NH4/(NH4+NO3)`
ratio.** Note: NSM1's `PN` is a dimensionless 0--1 weighting, whereas
QUAL2K's `Khnxp` is a half-saturation concentration (mg-N/L). They are not
interchangeable parameters; the v3 default 0.5 is correct for the NSM1
form.

### Adjudication 4 — Nitrification O2 stoichiometry and alkalinity change.

Nitrification O2 demand `ron`: v3/v1 use 4.57 mg-O2/mg-N (in the DOX
Process; outside the nitrogen-state scope but verified consistent). Full
nitrification NH4+ + 2 O2 -> NO3- + H2O + 2 H+ requires 64 g O2 / 14 g N
= 4.571. **CONSISTENT** with QUAL2E, QUAL2K, W2, and stoichiometry.

Alkalinity change, nitrification (`r_alkn`): v3 = `2/14/1000` eq/mg-N
(`parameters/alkalinity.py:16`). W2 (`water-quality.f90:3156,3166`):
nitrification term is `- 2.*NH4D` with the comment "Nitrification of
ammonium results in an alkalinity decrease: 2 eq. alk per 1 mole ammonium".
Stumm & Morgan (1996) Table 4.5: 2 eq alk consumed per mol NH4 nitrified
(the 2 H+ produced). **v3 r_alkn = 2/14/1000 is CORRECT and matches W2 and
Stumm & Morgan exactly.**

Alkalinity change, denitrification (`r_alkden`): v3 = `4/14/1000` eq/mg-N
(`parameters/alkalinity.py:17`). W2 (`water-quality.f90:3157,3166`):
denitrification term is `+ NO3D` with the explicit comment
"Denitrification of nitrate (to nitrogen gas) results in an alkalinity
increase: **1 eq. alk per 1 mole nitrate**". The canonical denitrification
half-reaction (NO3- + 6 H+ + 5 e- -> 1/2 N2 + 3 H2O, coupled to organic-C
oxidation) produces **1 equivalent of alkalinity per mole of N reduced**
(Stumm & Morgan 1996; Chapra 1997). The correct coefficient is therefore
`1/14/1000`, not `4/14/1000`. **Numeric answer: v3 (and v1) overstate the
denitrification alkalinity production by a factor of 4. r_alkn is correct
(2/14/1000); r_alkden is WRONG and should be 1/14/1000.** See SF-1.

### Adjudication 5 — Sediment / benthic NH4 and NO3 release.

v3 defaults: `rnh4_20 = 0.0`, `vno3_20 = 0.0` (`parameters/nitrogen.py:53--54`).
NSM1 Fortran defaults: `rnh4%rc20 = 0.0`, `vno3%rc20 = 0.0`
(`modNitrogen.f90:82,100`). v3 matches Fortran exactly.

CE-QUAL-W2 sediment NH4 release (`water-quality.f90:1535`):
`NH4SR = NH4R(JW) * SODD(K,I) * DO2(K,I)` — release is SOD-driven and
DO-gated, with `NH4R` a site-specific input (no hardwired default in src;
W2 control-file driven). W2 sediment NO3 loss (`water-quality.f90:1591`):
`NO3SED = NO3*NO3S(JW)*NO3TRM*BI/BH2` — a settling/uptake velocity `NO3S`,
again site-specific input.

Typical literature benthic fluxes (Bowie et al. 1985, Table; QUAL2K SOD
section): sediment NH4 release 10--300 mg-N/m2/d; sediment denitrification
velocity 0.1--1 m/d.

Numeric answer: the v3 zero defaults are **consistent with the NSM1
Fortran defaults and are a deliberate "off-unless-configured" design**.
This is defensible because (a) it matches the authoritative source, and
(b) W2 itself has no universal default for these (they are always
site-specific). The risk is operational, not a formulation error: a user
running v3 with defaults gets **zero benthic N flux**, which for many
systems (eutrophic lakes, low-gradient rivers with organic sediments)
materially understates the NH4 budget. This is the same behavior as NSM1
Fortran, so it is not a v3 regression, but it is a defaults-vs-reality gap
worth flagging (SF-3). **CONSISTENT-AS-DEFAULT, with a usability caveat.**

---

## 4. Science findings the parity round missed

The parity review (`review_nitrogen_n2.md`) reported CRITICAL 0, MAJOR 0
and concentrated on comment hygiene. Because parity was its frame, it did
not adjudicate whether the *shared* v1/v3 formulation is physically
correct. The following are science-correctness findings on the shared
formulation.

### SF-1 (MAJOR) — Denitrification alkalinity coefficient is 4x too large.

Location: `src/clearwater_modules_v3/parameters/alkalinity.py:17`
(`'r_alkden': 4.0 / 14.0 / 1000.0`), consumed at
`src/clearwater_modules_v3/processes/alkalinity.py:325`
(`return self.r_alkden * denit_flux * EQ_TO_MG_CACO3`).

Issue: denitrification of NO3 to N2 produces 1 equivalent of alkalinity
per mole of N (Stumm & Morgan 1996 Table 4.5; Chapra 1997). CE-QUAL-W2
encodes exactly this: `water-quality.f90:3157` comment "1 eq. alk per 1
mole nitrate" and `:3166` uses `+ NO3D` with coefficient 1 (compared with
`- 2.*NH4D` for nitrification, coefficient 2). v3 (and v1) use a
coefficient of 4. The 2:1 alkalinity ratio between nitrification and
denitrification that W2 enforces (2 eq consumed per mol nitrified, 1 eq
produced per mol denitrified) is broken in v3/v1, which has a 2:4 ratio.

Consequence: simulated alkalinity gains from denitrification are
overstated by a factor of 4. In systems where denitrification is a
significant alkalinity term (anoxic hypolimnia, organic-rich sediments,
wetland-influenced reaches), this biases the alkalinity state and any
pH/CO2 speciation derived from it. The nitrogen state variables (NH4,
NO3, OrgN, N2) are **not** affected — `r_alkden` lives only in the
Alkalinity Process — so this is not a nitrogen-mass-balance defect, but
it is a charge-balance / alkalinity-correctness defect driven by a
nitrogen reaction.

Classification: MAJOR. It is a definite formulation error against the
authoritative reference (W2) and standard geochemistry, it affects a
core derived quantity (alkalinity, hence pH), and it is V1-AND-V3-BOTH-
WRONG (the parity round could not have caught it because v3 faithfully
reproduces the v1 error). It is not CRITICAL only because it does not
corrupt the nitrogen state or violate nitrogen mass balance, and
alkalinity is a secondary derived field; a site without significant
denitrification is unaffected.

Origin confirmed: the error is **upstream in NSM1 Fortran itself**.
`/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/modAlkalinity.f90:54`
sets `ralkden = 4.0 / 14.0 / 1000.0` and `:103` applies it as
`Alk_Denit = ralkden * NO3_Denit * 50000.0` with
`dAlkdt = ... + Alk_Denit ...` (`:109`). NSM1 Fortran's `ralkn` at line 53
is `2.0/14.0/1000.0` (correct, matches W2). So NSM1 Fortran has the same
2:4 nitrification:denitrification ratio defect; v1 and v3 faithfully
reproduce it. This is therefore an **upstream NSM1 Fortran error
propagated unchanged through v1 to v3**, not a v1 port mistake.

Recommended fix: change `parameters/alkalinity.py:17` to
`'r_alkden': 1.0 / 14.0 / 1000.0` and update the corresponding docstring
in `processes/alkalinity.py:49`. Because this is an upstream NSM1 Fortran
error, the v3 change is a **deliberate, reference-anchored divergence
from NSM1 Fortran**: document it as such in the fix commit and the
corrections doc, citing CE-QUAL-W2 `water-quality.f90:3157` (explicit
"1 eq. alk per 1 mole nitrate") and Stumm & Morgan (1996) Table 4.5.
Report the defect upstream to the NSM1 Fortran maintainers
(`modAlkalinity.f90:54`). Add a regression/benchmark test: closed-system
NO3 -> N2 with no other alkalinity terms, assert dAlk/dt = 1 eq per mol N
(50.044 g CaCO3 per 14 g N per unit denit flux), matching W2.

### SF-2 (MINOR) — Nitrification rate default `knit_20 = 0.1 /d` is at the low end of the literature range and is not flagged as a screening default.

Location: `src/clearwater_modules_v3/parameters/nitrogen.py:50`.

Issue: `knit_20 = 0.1 /d` matches the NSM1 Fortran default
(`modNitrogen.f90:77`) and is within the published range (Bowie et al.
1985 Table 6-2: 0.1--1.0 /d; QUAL2E typical 0.5 /d; Chapra 1997
0.1--0.5 /d for rivers, lower in lakes). It is the *lowest* literature
value. For warm, well-oxygenated, biofilm-rich streams, measured
nitrification rates are often 0.3--1.0 /d. A user accepting the default
will systematically under-predict NH4 -> NO3 conversion.

Consequence: not a formulation error and not a v3 regression (it matches
the authoritative Fortran default), but a screening-default-vs-typical-
value gap. NH4 will tend to persist longer and NO3 production lag relative
to a calibrated system.

Classification: MINOR. The form is correct and the value is defensible as
a documented default; the issue is purely that the default is at the
extreme low end and the comment does not say so.

Recommended fix: keep the value (changing it would break Fortran parity)
but expand the inline comment at `parameters/nitrogen.py:50` to state the
literature range (0.1--1.0 /d, Bowie et al. 1985) and that 0.1 is the
low-end NSM1 screening default requiring site calibration.

### SF-3 (MINOR) — Zero sediment-release defaults silently disable benthic N flux.

Location: `src/clearwater_modules_v3/parameters/nitrogen.py:53--54`
(`rnh4_20 = 0.0`, `vno3_20 = 0.0`).

Issue: matches NSM1 Fortran (`modNitrogen.f90:82,100`) so not a
regression, but for many eutrophic or organic-sediment systems the
benthic NH4 flux (typically 10--300 mg-N/m2/d, Bowie et al. 1985) is a
first-order term in the NH4 budget. With the default, `ammonium_from_bed`
and `nitrate_bed_denitrification` are identically zero
(`nitrogen.py:769--785` and the `rnh4_tc/depth` path), so a defaulted run
omits a process that W2 includes by default whenever the user provides
`NH4R`/`NO3S` (W2 lines 1535, 1591). The zero default plus the
`use_SedFlux=True -> NotImplementedError` guard means there is no
turnkey benthic-flux path in v3 1.0.0.

Consequence: usability/completeness gap, not a formulation error. A user
must know to set `rnh4_20`/`vno3_20` to get any benthic N exchange.

Classification: MINOR (documented design choice; matches authoritative
default). Worth surfacing because the parity round classified the same
fact as a non-issue purely on parity grounds without assessing the
science consequence.

Recommended fix: no code change required; add a user-facing note in the
parameter docstring and the model documentation that benthic N flux is
off by default and must be enabled via `rnh4_20`/`vno3_20`, with the
typical literature range. Consider a startup INFO-level diagnostic when
both are zero and `use_OrgN`/algae are active.

### SF-4 (MINOR) — Nitrification/denitrification DO factors differ from the authoritative W2 form without an in-code reference.

Location: `nitrogen.py:885` (`1 - exp(-KNR*DOX)`) and
`nitrogen.py:758--762` (`1 - DOX/(DOX+KsOxdn)`).

Issue: these forms are correct and match NSM1 Fortran and QUAL2K (see
Adjudications 2), but CE-QUAL-W2 uses Monod `O2/(O2+KDO)` for
nitrification and its complement for denitrification. The v3 code cites
the v1/Fortran provenance but not the QUAL2K reference that makes the
exponential form defensible against the W2 alternative. A future reviewer
comparing v3 to a W2 application could mistake the form difference for an
error.

Consequence: documentation-to-reference traceability gap; no numerical
defect.

Classification: MINOR.

Recommended fix: add a one-line citation in the `nitrification_inhibition`
and `nitrate_denitrification` docstrings noting the exponential / Monod
low-DO forms follow QUAL2K (Chapra, Pelletier & Tao 2008) and NSM1
Fortran, and differ deliberately from CE-QUAL-W2's Monod DO factor.

### Observations (not defects)

* O-SCI-1: The denitrification NaN fallback differs between v3 (`-> 0`,
  `nitrogen.py:765--767`) and NSM1 Fortran (`-> kdnit_tc*NO3`,
  `modNitrogen.f90:266`). Reachable only at the unphysical point
  `DOX = -KsOxdn`. No practical consequence; recorded for completeness.
* O-SCI-2: `nitrification_inhibition` gates on `use_nitrate`
  (`nitrogen.py:878`) where NSM1 Fortran gates on `use_DOX`
  (`modNitrogen.f90:264`). Under all v3 1.0.0 supported configurations
  DOX is present, so the behavior is identical. Needs verification only
  if a `use_DOX=False` configuration becomes supported.
* O-SCI-3: NSM1 lumps organic N into a single `OrgN` pool while CE-QUAL-W2
  splits labile/refractory DON and PON (water-quality.f90:2472, 2505,
  2518). This is a STRUCTURAL-OK cross-model difference, not an error;
  the v3 lumped hydrolysis (`kon_tc*OrgN`) and settling (`vson/depth*OrgN`)
  kinetic forms are the standard QUAL2E single-pool organic-N
  formulation and are individually defensible.

---

## 5. Defaults vs literature table

| Parameter | v3 default | Literature / reference range | Source | Assessment |
|---|---|---|---|---|
| `knit_20` nitrification rate | 0.1 /d | 0.1--1.0 /d (typ. 0.5) | Bowie 1985 Table 6-2; QUAL2E; Chapra 1997 | In range, low end (SF-2) |
| `knit_theta` | 1.083 | 1.06--1.08 | QUAL2E (=1.083); Bowie 1985 | Exact match QUAL2E |
| `KNR` nitrif DO half-sat | 0.6 mg-O2/L | 0.6--2.0 mg-O2/L | NSM1 hardwire {0.6--0.7}; QUAL2K KNR | Defensible (low end) |
| `kdnit_20` denit rate | 0.002 /d | 0.002--0.2 /d water column | Bowie 1985; QUAL2K | In range; water-column denit is minor, bed term intended to dominate |
| `kdnit_theta` | 1.045 | 1.045--1.08 | QUAL2E; Chapra 1997 | Correct (v3); v1's 1.08 wrong |
| `kon_20` OrgN hydrolysis | 0.1 /d | 0.02--0.4 /d | Bowie 1985 Table 6-2; QUAL2E | In range, central |
| `kon_theta` | 1.047 | 1.047 (universal OM) | Chapra 1997; QUAL2K | Correct (v3); v1's 1.074 wrong |
| `rnh4_theta` | 1.074 | 1.06--1.08 (sediment exchange) | sediment-water exchange convention | Correct (v3); v1's 1.047 wrong |
| `vno3_theta` | 1.08 | ~1.045--1.08 | sediment-exchange convention | Correct (v3); v1's 1.045 wrong |
| `vson_20` OrgN settling | 0.01 m/d | 0.05--2.5 m/d (PON) | Bowie 1985; QUAL2K POM | Low; NSM1 lumped-OrgN convention (STRUCTURAL-OK) |
| `rnh4_20` sediment NH4 | 0.0 /d | 10--300 mg-N/m2/d when active | Bowie 1985 | Off by design = NSM1 default (SF-3) |
| `vno3_20` sediment denit | 0.0 /d | 0.1--1 m/d when active | Bowie 1985 | Off by design = NSM1 default (SF-3) |
| `PN` algal NH4 preference | 0.5 | 0--1 (0.5 = neutral) | QUAL2E/EUTRO preference form | Correct neutral default |
| `KsOxdn` denit DO half-sat | 0.1 mg-O2/L | 0.1--0.5 mg-O2/L | QUAL2K anoxic Monod | In range (low end, sharp anoxic switch) |
| `ron` nitrif O2 demand | 4.57 mg-O2/mg-N | 4.57 (64/14) | stoichiometry; QUAL2E/W2 | Exact |
| `r_alkn` nitrif alk | 2/14/1000 eq/mg-N | 2 eq/mol N | W2 line 3156; Stumm & Morgan 1996 | Exact, correct |
| `r_alkden` denit alk | **4**/14/1000 eq/mg-N | **1** eq/mol N | W2 line 3157; Stumm & Morgan 1996 | **WRONG, 4x high (SF-1)** |

---

## 6. Recommended follow-up benchmarks

1. Closed-system denitrification alkalinity test (validates SF-1 fix):
   single cell, NO3 only, denitrification on, all other alkalinity terms
   off; assert dAlk = `(50.044/14.00674) * 1 * (NO3 consumed)` mg-CaCO3,
   i.e. 1 eq per mol N, matching W2 `water-quality.f90:3166`.
2. Nitrification alkalinity regression: same setup with NH4 only,
   nitrification on; assert dAlk = `- (50.044/14.00674) * 2 * (NH4
   nitrified)` (already correct; lock it in so an SF-1 fix does not
   regress r_alkn).
3. Low-DO inhibition curve test: tabulate `nitrification_inhibition(DO)`
   at DO = 0, 0.5, 1, 2, 5, 10 mg/L and assert against the analytic
   `1 - exp(-0.6*DO)`; document the QUAL2K reference (Adjudication 2).
4. Ammonia-preference reduction test: assert that with PN = 0.5 the v3
   `_ap_uptake_fr_nh4` equals the QUAL2E simple ratio `NH4/(NH4+NO3)`
   exactly across a NH4/NO3 grid (validates Adjudication 3).
5. Theta-pairing regression: assert
   `(kon_theta, kdnit_theta, rnh4_theta, vno3_theta) ==
   (1.047, 1.045, 1.074, 1.08)` to lock in the Phase 9.E correction
   against future re-transposition.

## 7. Open questions

1. SF-1 origin: RESOLVED. NSM1 Fortran `modAlkalinity.f90:54` itself uses
   `ralkden = 4.0/14.0/1000.0`. The defect is upstream in NSM1 Fortran
   and propagated unchanged to v1 and v3. The v3 fix is therefore a
   deliberate, reference-anchored divergence from upstream NSM1 and
   should be reported upstream. Remaining question is the product
   decision on whether to break NSM1-Fortran parity to correct the
   science (recommended: yes, with documentation).
2. Is a `use_DOX=False` NSM1 configuration in scope for v3 1.x? If so the
   `nitrification_inhibition` gating predicate (O-SCI-2) must change from
   `use_nitrate` to `use_DOX` to match NSM1 Fortran.
3. Should v3 ship a turnkey constant-benthic-flux path (non-zero
   `rnh4_20`/`vno3_20` example or a documented enable recipe), given that
   `use_SedFlux=True` raises `NotImplementedError`? Product decision
   (SF-3).
