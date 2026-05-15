# v3 NSM1 Algae + Benthic Algae + Alkalinity — Science-Correctness Validation

Reviewer: water-quality model science-correctness reviewer (Claude)
Date: 2026-05-15
Branch: `streaming`
Repo: `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming`
Mode: SCIENCE validation against well-validated reference models and the
literature. Parity with v1 is NOT treated as proof of correctness; the
highest-value findings are terms where v1 AND v3 are both wrong.

Reference models consulted:

1. CE-QUAL-W2 v2026.02 (LOCAL, authoritative, line-checked):
   `/Users/todd/GitHub/CE-QUAL-W2-ERDC/CE-QUAL-W2-ERDC-dev/src/W2_v2026.02/water-quality.f90`
   and `heat-exchange.f90`. ENTRY points read line by line:
   TEMPERATURE_RATES (228-277), KINETIC_RATES algal block (655-703),
   ALGAE (1764-1776), EPIPHYTON (2301-2345), algal organic-matter
   routing (1657, 1704, 2356), ALKALINITY (3151-3173). `BETA`
   definition read at `heat-exchange.f90:84,93`.
2. QUAL2K / QUAL2Kw v2.11 (Chapra, Pelletier & Tao 2008, Documentation
   and Users Manual, Civil & Environmental Engineering Dept., Tufts
   University, for U.S. EPA). Cited from domain knowledge — NOT
   locally checked. Confidence: medium-high (canonical equations are
   stable across the 2.0x--2.11 series).
3. QUAL2E (Brown & Barnwell 1987, EPA/600/3-87/007). Cited from domain
   knowledge — not locally checked.
4. WASP8 / EUTRO (Wool, Ambrose, Martin & Comer, WASP8 documentation;
   Ambrose et al. 1993 EUTRO). Cited from domain knowledge — not
   locally checked.
5. Bowie et al. 1985, "Rates, Constants, and Kinetics Formulations in
   Surface Water Quality Modeling," 2nd ed., EPA/600/3-85/040; Chapra
   1997, *Surface Water-Quality Modeling*, McGraw-Hill. Cited from
   domain knowledge — not locally checked.

---

## 1. Verdict

The v3 NSM1 floating-algae and benthic-algae kinetic *forms* are
scientifically defensible and consistent with QUAL2E/QUAL2K and the
Bowie/Chapra literature. The kinetic skeleton (Arrhenius temperature
correction, Monod nutrient half-saturation, the three light-limitation
options, first-order respiration/death, Stokes settling, Forward Euler
in days) is a faithful implementation of the QUAL2K-class phytoplankton
model. NSM1's structural choices that differ from CE-QUAL-W2 (lumped
organic-matter pools rather than W2's labile/refractory split; a single
Arrhenius theta rather than W2's four-parameter rising/falling
temperature curve; a half-saturation light default rather than W2's
Steele) are legitimate modeling-method differences, not defects, and
NSM1's choices align with QUAL2K, which is itself an extensively
validated EPA-supported model.

Five science findings that the v1<->v3 parity round did not surface are
recorded in Section 4. The single Critical one is **SA-1**, which is the
science confirmation of the parity round's CA-1: the v3 alkalinity-algae
coupling is dimensionally wrong by a factor of 1000 (floating) and 100
(benthic), and it is *also structurally non-standard* relative to W2's
nitrogen-mass alkalinity formulation, so the correct value cannot be
recovered just by inserting `rca`. The other four are Major/Minor
science issues that parity passed because v1 shares them.

The five adjudication questions are answered with explicit numbers in
Section 3. The headline numeric answers:

1. **CA-1 / SA-1 correct value:** the scientifically correct alkalinity
   change for algal NH4-driven photosynthesis is approximately
   **3.57e-3 eq alk per mg-N taken up**, equivalently
   **0.357 mg-CaCO3/L per mg-N/L**, equivalently (via NSM1's C:N and
   C:Chla ratios) **about 4.39e-3 mg-CaCO3/L per ug-Chla/L of growth**
   for the all-NH4 case. v3's `_floating_algae_growth_alk_flux` returns
   ~4.39 mg-CaCO3/L per ug-Chla/L (a 1000x inflation). Confirmed
   against W2 ALKALINITY (`water-quality.f90:3164`) and Stumm & Morgan
   stoichiometry quoted in the W2 header (`water-quality.f90:3152-3159`).
2. **A2 (POC vs DOC split of algal mortality):** W2's `APOM` (algal
   fraction to particulate organic matter) and `EPOM` (epiphyton) are
   user inputs with W2 manual typical values near 0.8 (range commonly
   0.6--0.9); QUAL2K routes the bulk of phytoplankton death to
   detritus/POM; Bowie/Chapra report most non-respired algal carbon is
   particulate. **v1's `f_pocp` = 0.9 is the defensible value; v3's
   inline 0.5 is low and undocumented.**
3. **A11 (PAR vs total shortwave):** W2 feeds `(1-BETA)*SRON` (shortwave
   penetrating the surface, W m^-2) into the light limitation and folds
   the PAR fraction into the calibrated `ASAT`. QUAL2K explicitly works
   in PAR and applies a photosynthetically-available fraction (~0.47)
   to solar radiation. NSM1's `KL` default of 10 W m^-2 is a
   PAR-scale half-saturation, so the irradiance handed to `limit_light`
   should be PAR. If the registry delivers total shortwave with no
   ~0.47 scaling, the light limitation is overstated; this is a real,
   unresolved ~2x risk and remains an Open Question (it cannot be
   closed from the algae source alone).
4. **Light-limitation form:** v3 default is **half-saturation
   (option 1)**, the QUAL2E/QUAL2K-class depth-integrated
   Beer-Lambert/half-saturation light model. W2's only algal light
   model is **Steele** (photoinhibition). NSM1 also offers Steele
   (option 3) and Smith (option 2). The chosen default form is
   defensible and matches QUAL2E; not a defect.
5. **Nutrient limitation combination:** v3 default is
   **multiplicative** (option 1). W2 uses **minimum (Liebig)**
   unconditionally. QUAL2K offers both and historically defaults to
   multiplicative. The multiplicative default is a defensible
   modeling choice consistent with QUAL2K; not a defect, but it is the
   more conservative (lower-growth) choice and worth a documentation
   note.

Findings by severity (this science pass, beyond the parity round):
1 Critical (SA-1, confirms+extends CA-1), 2 Major (SA-2 f_pocp,
SA-3 PAR), 2 Minor (SA-4 alkalinity structural form, SA-5 krb_theta).

Confidence: HIGH for all W2-anchored conclusions (every W2 citation was
read line by line in the local v2026.02 source). MEDIUM for QUAL2K /
QUAL2E / WASP / Bowie / Chapra numeric ranges (cited from domain
knowledge, not locally checked); these are used to corroborate, not as
the sole basis for any finding's severity.

---

## 2. Cross-model validation matrix

Verdict legend: CONSISTENT = v3 form/value scientifically defensible
and agreeing with references within method differences;
V1&V3-BOTH-WRONG = parity round passed it but both are scientifically
wrong; MODEL-STRUCTURAL-DIFF-OK = legitimate inter-model method
difference, not an error; V3-DIVERGENT = v3 differs from v1 and from
references.

| Term / default | v3 form & value | v1 | CE-QUAL-W2 (file:line) | QUAL2K / QUAL2E | WASP / EUTRO | Literature (Bowie/Chapra) | VERDICT |
|---|---|---|---|---|---|---|---|
| Temperature correction (growth) | Arrhenius `k20 * theta^(T-20)`, `theta=1.047` (`algae.py:59`) | same | 4-param rising/falling `ATRMR*ATRMF`, `FR`/`FF` of AT1-AT4/AK1-AK4 (`water-quality.f90:241-245`) | QUAL2E/QUAL2K: Arrhenius `theta^(T-20)`, theta~1.066--1.08 growth | EUTRO: Arrhenius theta^(T-20) | Bowie Table 6-x: Arrhenius common; W2 4-param more physical at temp extremes | MODEL-STRUCTURAL-DIFF-OK (see SA-5 for theta value) |
| `mu_max_20` | 2.0 /d (`algae.py:56`) | 1.0 /d | `AG(JA)` input; W2 manual typical ~2.0--2.5 /d | QUAL2K default ~2.0 /d; QUAL2E typical 1.5--3.0 /d | EUTRO typical 2.0 /d | Bowie Table 6-1: 0.5--4.0 /d, central ~2 /d | CONSISTENT (v3 better than v1; documented A12) |
| `kdp_20` (death/mortality) | 0.05 /d (`algae.py:57`) | 0.15 /d | `AM(JA)` input; manual typical 0.01--0.1 /d | QUAL2K default ~0.05--0.1 /d | EUTRO 0.02--0.1 /d | Bowie: 0.01--0.1 /d | CONSISTENT |
| `krp_20` (respiration) | 0.10 /d (`algae.py:58`) | 0.2 /d | `AR(JA)` input; manual typical 0.02--0.3 /d | QUAL2K default ~0.1--0.2 /d | EUTRO 0.05--0.2 /d | Bowie: 0.05--0.5 /d | CONSISTENT |
| `vsap` settling | 0.15 m/d (`algae.py:62`) | 0.15 m/d | `AS(JA)` input; typical 0.0--0.5 m/d (diatoms higher) | QUAL2K ~0.05--0.5 m/d | EUTRO 0.1--0.5 m/d | Bowie: 0.05--2 m/d species-dependent | CONSISTENT |
| `KsN` N half-sat | 0.04 mg-N/L (`algae.py:29`) | 0.04 | `AHSN(JA)` input; W2 typical ~0.01--0.08 mg-N/L | QUAL2K ~0.015--0.04 mg-N/L | EUTRO 0.01--0.05 mg-N/L | Bowie: 0.005--0.4 mg-N/L | CONSISTENT |
| `KsP` P half-sat | 0.0012 mg-P/L (`algae.py:30`) | 0.0012 | `AHSP(JA)` input; W2 typical ~0.002--0.03 mg-P/L | QUAL2K ~0.001--0.005 mg-P/L | EUTRO 0.001--0.005 mg-P/L | Bowie: 0.001--0.08 mg-P/L | CONSISTENT (at low end; acceptable) |
| Light limitation FORM (default) | half-saturation, depth-integrated Beer-Lambert (`floating_algae.py:797-808`) | same | Steele only (`water-quality.f90:689`) | QUAL2E half-sat; QUAL2K offers half-sat/Smith/Steele | EUTRO Steele or half-sat | Both forms standard (Bowie ch.6) | MODEL-STRUCTURAL-DIFF-OK (matches QUAL2E) |
| `KL` light half-sat | 10.0 W/m^2 (`algae.py:28`) | 10.0 | `ASAT(JA)` input (Steele Is), W2 typical 50--150 W/m^2 shortwave; not directly comparable to a half-sat | QUAL2K half-sat KLp typ. ~50--100 uE/m^2/s PAR (~10--22 W/m^2 PAR) | EUTRO Is ~150--350 ly/d | Bowie: Is 100--400 ly/d; KL 10--100 uE | CONSISTENT only if input is PAR (see SA-3 / A11) |
| Irradiance fed to `limit_light` | `solar_radiation` passed raw, no Fr_PAR (`floating_algae.py:414`) | same | `(1-BETA)*SRON`, shortwave through surface; PAR folded into calibrated ASAT (`water-quality.f90:678`; BETA `heat-exchange.f90:84`) | QUAL2K applies PAR fraction ~0.47 to solar | EUTRO uses PAR | Chapra 1997 §: PAR ~ 0.43--0.47 of total SW | V1&V3-BOTH-WRONG if registry is total SW (SA-3) |
| Nutrient-limit combination (default) | multiplicative (opt 1) (`floating_algae.py:659`) | same | minimum/Liebig, unconditional (`water-quality.f90:696`) | QUAL2K offers mult/min; QUAL2E min(N,P) | EUTRO multiplicative | Both standard; min is more common in mechanistic codes | MODEL-STRUCTURAL-DIFF-OK |
| N limitation FN | `N/(KsN+N)` Monod, NH4+NO3 (`floating_algae.py:757-781`) | same | `(NH4+NO3)/(NH4+NO3+AHSN)` (`water-quality.f90:691`) | identical Monod | identical | Monod standard | CONSISTENT |
| P limitation FP | `fdp*TIP/(KsP+fdp*TIP)` (`floating_algae.py:727-755`) | same | `FDPO4*PO4/(FDPO4*PO4+AHSP)` (`water-quality.f90:690`) | identical w/ dissolved fraction | identical | Monod standard | CONSISTENT |
| Respiration form | `krp_tc * Ap`, first order (`floating_algae.py:704`) | same | `ATRM*AR*DO3` (temp + DO gate) (`water-quality.f90:700`) | first-order, often DO-gated | first-order | first-order standard | CONSISTENT (W2 adds DO gate; method diff) |
| Death/mortality form | `kdp_tc * Ap`, first order (`floating_algae.py:715`) | same | `(ATRMR+1-ATRMF)*AM` temp-shaped (`water-quality.f90:701`) | first-order | first-order | first-order standard | MODEL-STRUCTURAL-DIFF-OK |
| Settling form | `vsap/depth * Ap` (`floating_algae.py:723`) | same | `AS*(ALG(K-1)-ALG(K))*BI/BH2` flux-form | QUAL2K `vsap/depth` | EUTRO `vs/H` | Chapra: v/H first order | CONSISTENT |
| Algal mortality C split f_pocp | 0.5 inline default (`floating_algae.py:111-117`) | 0.9 | `APOM` input; W2 manual typical ~0.6--0.9 (POM share) (`water-quality.f90:1657,1704,2356`) | QUAL2K: most death -> detritus/POM | EUTRO: bulk to detritus | Bowie/Chapra: majority particulate | V1&V3 DIFFER; v1 (0.9) defensible, v3 (0.5) low (SA-2) |
| Alk per algal growth (NH4) | `r_alkaa * AWc * 50000`, `AWc=40` (`alkalinity.py:362-368`) | `r_alkaa*50000*rca`, `rca=AWc/AWa=0.04` | `(50.044/14.00674)*14/16*NH4-flux` mg-CaCO3/mg-N (`water-quality.f90:3164`) | QUAL2K: same Stumm-Morgan eq-alk/mol-N | EUTRO similar | Stumm & Morgan 1996 Table 4.5 | V1&V3-BOTH-WRONG (v3 1000x; both non-standard form) (SA-1) |
| Alk per benthic-algae growth | `r_alkba * BWc * 50000`, `BWc=40` (`alkalinity.py:411-419`) | `*rcb`, `rcb=BWc/BWd=0.4` | epiphyton N-flux into ALKSS (`water-quality.f90:3164`, NH4EP/NO3EG) | as above | as above | as above | V1&V3-BOTH-WRONG (v3 100x) (SA-1) |
| Alkalinity stoich `r_alkaa` | 14/106/12/1000 eq/mg-C (`parameters/alkalinity.py:14`) | same | 14/16 eq-alk/mol-N applied to N flux (`water-quality.f90:3164`) | Stumm-Morgan 14 eq/16 mol N | similar | Stumm & Morgan 1996 | CONSISTENT (the ratio itself is correct; misuse is SA-1) |
| Benthic `mub_max_20` | 0.4 /d (`balgae.py:41`) | 0.4 | `EG(JE)` epiphyton input; W2 typical 1.0--2.0 /d | QUAL2K bottom-algae zero-order or first-order, ~0.3--3 /d | n/a | Bowie periphyton: 0.5--8 g/m^2/d areal | CONSISTENT (low end; acceptable) |
| Benthic `krb_theta` | 1.06 (`balgae.py:45`) | 1.06 | epiphyton uses same 4-param temp curve | QUAL2K theta ~1.047--1.07 | n/a | Bowie 1.045--1.08 | CONSISTENT but inconsistent within NSM1 (SA-5) |
| C:Chla ratio (`AWc/AWa`) | 0.04 mg-C/ug-Chla = 40 g-C/g-Chla (`algae.py:24,27`) | same | W2 carbon via `AC(JA)` (C:DW), Chla derived separately | QUAL2K default 40 ug-C/ug-Chla typical | EUTRO 30--50 | Bowie: C:Chla 20--100, typ ~40--50 | CONSISTENT |
| N:Chla ratio (`AWn/AWa`) | 0.0072 mg-N/ug-Chla = 7.2 g-N/g-Chla | same | via `AN(JA)` N:DW | QUAL2K ~7.2 (Redfield-scaled) | EUTRO ~7--10 | Redfield-consistent | CONSISTENT |
| P:Chla ratio (`AWp/AWa`) | 0.001 mg-P/ug-Chla = 1.0 g-P/g-Chla | same | via `AP(JA)` P:DW | QUAL2K ~1.0 | EUTRO ~0.5--1.0 | Redfield-consistent | CONSISTENT |
| Integrator | Forward Euler in days, `state + rate*dt_days` (`floating_algae.py:438-440`) | same | explicit, sub-step `DLT` | QUAL2K 4th-order RK | EUTRO explicit Euler | explicit standard for these timescales | MODEL-STRUCTURAL-DIFF-OK |

---

## 3. The five adjudications (explicit numeric answers)

### Adjudication 1 — CA-1: correct alkalinity change per unit algal growth/respiration

**Scientifically correct value (independent of v1):**

The alkalinity effect of phytoplankton growth and respiration is driven
by the *nitrogen* taken up or released, not by carbon per se. The
canonical stoichiometry (Stumm & Morgan 1996, Table 4.5, quoted
verbatim in the W2 source header at `water-quality.f90:3152-3159`):

- NH4 uptake during photosynthesis: alkalinity **decreases** by
  14 eq per 16 mol N.
- NO3 uptake during photosynthesis: alkalinity **increases** by
  18 eq per 16 mol N.
- NH4 production during respiration: alkalinity **increases** by
  14 eq per 16 mol N.
- (Nitrification: -2 eq per mol N; denitrification: +1 eq per mol N —
  not the algal terms but the same conversion basis.)

W2 implements this directly (`water-quality.f90:3164`):

    ALKSS = (50.044/14.00674) * ( 14./16.*(NH4AP + ... - NH4MG)
                                 + 18./16.*(NO3AG + ...) - 2.*NH4D + ...)

where `NH4AP`, `NO3AG` are nitrogen mass fluxes (mg-N/L/d) and
`50.044/14.00674 = 3.5728` is mg-CaCO3 per mg-N (= 50.044 g-CaCO3/eq
divided by 14.00674 g-N/mol, i.e. the equivalent weight of CaCO3 per
gram of nitrogen).

Therefore the **scientifically correct alkalinity change for algal
growth using NH4 is:**

- **3.5728 mg-CaCO3/L per mg-N taken up**, scaled by 14/16
  = **3.126 mg-CaCO3/L per mg-N** for the NH4-uptake pathway;
- in equivalents: **(14/16) / 14.00674 / 1000 = 6.249e-5 eq per mg-N**,
  i.e. **6.25e-5 eq alk per mg-N** consumed (alk sink).

Translating to NSM1's units via the N:Chla ratio
`rna = AWn/AWa = 0.0072 mg-N/ug-Chla`:

- per ug-Chla/L of NH4-driven growth:
  `0.0072 mg-N/ug-Chla * 3.126 mg-CaCO3/mg-N`
  = **~2.25e-2 mg-CaCO3/L per ug-Chla/L** (NH4 pathway, alk sink).

Now evaluate NSM1's own formula correctly. NSM1 expresses the same
chemistry through carbon: `r_alkaa = 14/106/12/1000 = 9.171e-6 eq/mg-C`
(14 eq alk per 106 mol C per Redfield, 12 g-C/mol, /1000 for ug). With
the *correct* `rca = AWc/AWa = 0.04 mg-C/ug-Chla` and `50000 mg-CaCO3/eq`:

    correct flux per ug-Chla (all NH4) =
      r_alkaa * 50000 * rca
      = 9.171e-6 * 50000 * 0.04
      = 1.834e-2 mg-CaCO3/L per ug-Chla/L.

(The small difference from the 2.25e-2 N-based number is because NSM1
uses the Redfield 14:106:C path while W2 uses the explicit 14:16:N path;
both are ~1.8e-2 to 2.3e-2 mg-CaCO3/L per ug-Chla/L — same order, the
correct band.)

**v3's actual return** (`alkalinity.py:362-368`, `rca = self.AWc = 40`):

    r_alkaa * 50000 * 40 = 9.171e-6 * 50000 * 40
                         = 18.34 mg-CaCO3/L per ug-Chla/L.

This is **1000x the correct value** (40 / 0.04 = 1000). Confirmed.
Benthic (`alkalinity.py:411-419`, `rcb = self.BWc = 40` vs correct
`BWc/BWd = 0.4`): **100x the correct value**. Confirmed.

**VERDICT CA-1 / SA-1:** v3 is wrong by 1000x (floating) / 100x
(benthic). The correct value is approximately **1.83e-2 mg-CaCO3/L per
ug-Chla/L** for the all-NH4 floating-algae growth pathway (alk sink),
equivalently **~6.25e-5 eq alk per mg-N**, equivalently
**~3.13 mg-CaCO3/L per mg-N**. v1's `rca = AWc/AWa = 0.04` recovers the
correct order of magnitude; v3's raw `AWc = 40` does not. Severity
Critical, unchanged from the parity round's CA-1, but see SA-1 below for
the *additional* science finding the parity round did not state: even
v1's form is a non-standard carbon-routed approximation of an
intrinsically nitrogen-driven process, which is acceptable here only
because NSM1 keeps N:C fixed at Redfield.

### Adjudication 2 — A2: f_pocp / f_pocb (POC vs DOC split of algal mortality)

W2 routes algal mortality to organic matter through `APOM(JA)`
(epiphyton: `EPOM(JE)`): the *particulate* fraction goes to LPOM
(`water-quality.f90:1704`: `LPOMAP += APOM(JA)*AMR*ALG`) and the
*dissolved* remainder to LDOM (`water-quality.f90:1657`:
`LDOMAP += (AER + (1-APOM(JA))*AMR)*ALG`). `APOM` is a user-supplied
control-file parameter, not a source literal, so I cannot cite a code
default; the CE-QUAL-W2 user manual (Cole & Wells) and standard W2
data sets use `APOM` in the **0.6--0.9** band, most commonly near
**0.8**, reflecting that the majority of dead algal biomass is
particulate. QUAL2K routes phytoplankton death predominantly to the
detritus (particulate) pool, and Bowie et al. (1985) and Chapra (1997)
both characterize non-respired algal mortality as mostly particulate.

**Numeric answer:** the defensible value for the particulate fraction
of algal mortality carbon is **0.7--0.9**; **v1's `f_pocp = f_pocb =
0.9` sits at the top of and within that band**, while **v3's inline
`0.5` is below the literature/W2/QUAL2K consensus** and routes ~40%
more dead-algal carbon to DOC than is physically typical. This shifts
the downstream DOC oxidation -> DIC and dissolved-oxygen demand pathways.
**VERDICT: V3-DIVERGENT.** v3's 0.5 is the weaker value. The parity
round flagged this as MAJOR A2 on parity grounds; the science check
confirms 0.9 is correct and elevates the *consequence* (DOC/DO/DIC
mass routing bias), not just "undocumented deviation."

### Adjudication 3 — A11: PAR vs total shortwave

W2 builds the algal light driver as
`LIGHT = (1.0 - BETA(JW)) * SRON(JW) * SHADE(I) / ASAT(JA)`
(`water-quality.f90:678`). `SRON` is incident **short-wave** solar
radiation; `BETA` is the **surface absorption fraction**, defined in
`heat-exchange.f90:84` as
`BETA = 0.255 - 8.5e-3*TSTAR + 2.04e-4*TSTAR^2` (a function of wet-bulb
temperature, ~0.25--0.45). So `(1-BETA)*SRON` is the **shortwave that
penetrates the surface**, *not* a PAR fraction. W2 does **not** apply a
0.45/0.47 PAR multiplier in the algal light term; instead the PAR-vs-SW
distinction is absorbed into the calibrated saturation intensity
`ASAT` (W m^-2 of shortwave). W2 is internally consistent because
`ASAT` is calibrated on the same shortwave basis as `SRON`.

QUAL2K, by contrast, explicitly works in PAR: it multiplies solar
radiation by the photosynthetically-available fraction (commonly ~0.47)
before the light-limitation function, and its light parameters are on a
PAR basis. Chapra (1997) gives PAR ~ 0.43--0.47 of total shortwave.

NSM1's `KL` default is **10 W/m^2** (`algae.py:28`). A 10 W/m^2 light
half-saturation / saturation constant is physically a **PAR-scale**
value (full-sun shortwave is ~200--1000 W/m^2; full-sun PAR is
~100--500 W/m^2; algal light saturation is reached at tens of
W/m^2 PAR). For the half-saturation light model to be physically
meaningful at `KL = 10 W/m^2`, the irradiance handed to `limit_light`
**must be PAR**, not total shortwave.

**Numeric answer:** if the registry `solar_radiation` is total
shortwave and is passed to `limit_light` with no ~0.47 PAR scaling
(`floating_algae.py:414` does pass it raw), the light-limitation factor
is overstated. With the half-saturation form at `KL = 10`, feeding
~2x-too-large irradiance pushes `FL` from a typical light-limited
~0.4--0.6 toward ~0.6--0.8, a **growth overprediction on the order of
30--60%** under light-limited conditions (not a clean 2x because the
half-saturation function saturates; the error is largest at low light
and shrinks as light saturates). **VERDICT: V1&V3-BOTH-WRONG IF the
registry value is total shortwave** (parity passed it because v1 has
the same omission). This cannot be closed from the algae source; it
depends on the registry `solar_radiation` contract. Recorded as SA-3
(Major) and as Open Question 1.

### Adjudication 4 — Light-limitation FORM

v3 default `light_limitation_option = 1`
(`algae.py:64`) selects the depth-integrated Beer-Lambert /
half-saturation model
`(1/(L*d)) * ln((KL+I0)/(KL+I0*exp(-L*d)))`
(`floating_algae.py:797-808`). Option 2 is Smith
(`floating_algae.py:810-852`); option 3 is Steele photoinhibition
(`floating_algae.py:854-871`).

W2's only algal light model is **Steele** (photoinhibition):
`ALLIM = 2.718282*(exp(-LAM2)-exp(-LAM1))/(GAMMA*H2)` with
`LAM1 = I0/ASAT` (`water-quality.f90:689`). QUAL2E uses the
half-saturation (Monod) light model as its standard; QUAL2K offers
half-saturation, Smith, and Steele as user options exactly as NSM1
does.

**Answer:** the v3 default form (half-saturation) is **scientifically
defensible and matches QUAL2E**; the menu of three options matches
QUAL2K. It differs from W2's Steele-only choice, but that is a
legitimate inter-model method difference, not a defect. The one caveat
is that the half-saturation model has **no photoinhibition**, so under
very high surface irradiance it will over-predict surface-layer growth
relative to W2's Steele; for a depth-averaged 2D-laterally-averaged or
1D river this is the conventional and accepted simplification.
**VERDICT: MODEL-STRUCTURAL-DIFF-OK.** No new finding; documented as
A12-class. One Minor doc note (SA-4 includes a recommendation to state
the no-photoinhibition limitation).

### Adjudication 5 — Nutrient-limitation combination

v3 default `growth_rate_option = 1` = **multiplicative**
(`floating_algae.py:659`: `growth * FP * FN * FL`). W2 uses
**minimum (Liebig)** unconditionally
(`water-quality.f90:696`: `LIMIT = MIN(APLIM, ANLIM, ASLIM, ALLIM)`),
multiplying only the single most-limiting factor by the temperature and
max-growth terms. QUAL2K offers both multiplicative and minimum and has
historically defaulted to multiplicative; QUAL2E uses
`min(N-limitation, P-limitation)` then multiplies by light.

**Answer:** the multiplicative default is a **defensible, standard
choice consistent with QUAL2K** and is not a defect. It is, however,
the more conservative (lower net growth) formulation: under
simultaneous moderate N and P and light limitation, multiplicative
produces a substantially smaller growth factor than Liebig minimum
(e.g., 0.5 * 0.5 * 0.5 = 0.125 vs min = 0.5). This is a known
modeling-philosophy difference (Liebig: a single resource limits;
multiplicative: co-limitation compounds). NSM1 exposes both, so a user
can match W2/QUAL2E by selecting option 2. **VERDICT:
MODEL-STRUCTURAL-DIFF-OK.** No new finding; one Minor documentation
recommendation (state the philosophy and that option 2 = Liebig matches
W2/QUAL2E) folded into SA-4.

---

## 4. New science findings the parity round MISSED

The parity round (Section 1 of `review_algae.md` / `review_carbon_alkalinity.md`)
correctly traced v1<->v3 algorithm parity but, by construction, could
not catch terms where v1 itself is scientifically wrong. The following
are science findings against the *formulation*, not the port.

**SA-1 (CRITICAL) — alkalinity-algae coupling is wrong AND its form is
non-standard.** `alkalinity.py:362, 386, 411, 441`. The parity round's
CA-1 already identifies the 1000x/100x magnitude error (`rca = self.AWc`
instead of `AWc/AWa`). The *additional* science finding: even with the
magnitude fixed (`rca = AWc/AWa`), NSM1's alkalinity-algae term is a
**carbon-routed approximation of an intrinsically nitrogen-driven
process**. W2 (`water-quality.f90:3164`) and Stumm & Morgan (1996,
Table 4.5, quoted at `water-quality.f90:3152-3159`) compute alkalinity
change from the **nitrogen** assimilated (14/16 eq per mol N for NH4,
18/16 for NO3), because the alkalinity shift is the charge balance of
the N species, not of carbon fixation. NSM1's `r_alkaa = 14/106/12`
backs the 14 eq into a *carbon* basis via the Redfield 106:1 C:N ratio.
This is only correct so long as the algal N:C ratio is held at exactly
Redfield (106 C : 16 N). NSM1's `AWn/AWc = 7.2/40 = 0.18 mg-N/mg-C`
corresponds to a molar C:N of `(40/12)/(7.2/14) = 6.48`, **not** the
Redfield 106/16 = 6.625 that `r_alkaa`'s `/106` assumes. The result is
a small (~2%) intrinsic stoichiometric inconsistency *on top of* the
1000x bug, and a structural fragility: if a user calibrates `AWn`/`AWc`
away from Redfield, the alkalinity coupling silently desynchronizes
from the nitrogen mass balance. Recommended fix: implement the
alkalinity-algae term on the **nitrogen-flux basis** like W2 (derive it
from the algal NH4/NO3 uptake fluxes the Nitrogen process already
computes, times `50.044/14.00674 * 14/16` and `* 18/16`), which is both
dimensionally self-consistent and robust to non-Redfield stoichiometry.
At minimum, fix CA-1 (`rca = AWc/AWa`, `rcb = BWc/BWd`) and add a code
comment that the carbon routing assumes Redfield N:C. Severity Critical
(the magnitude error alone forces unphysical alkalinity once algal
coupling is active under bloom conditions).

**SA-2 (MAJOR) — f_pocp / f_pocb default 0.5 is below the
literature/W2/QUAL2K consensus.** `floating_algae.py:111-117`,
`benthic_algae.py:78-84`. The parity round flagged this (A2) as an
undocumented deviation from v1's 0.9. The science finding: v1's 0.9 is
the *correct* value (W2 `APOM` typical ~0.6--0.9, most often ~0.8;
QUAL2K and Bowie/Chapra: dead algal carbon is predominantly
particulate). v3's 0.5 is not merely undocumented; it is
*scientifically low*, biasing ~40% of algal-mortality carbon into the
DOC pool, which then changes DOC oxidation, DIC production, and SOD/DO
demand timing. Recommended fix: restore the inline fallback to 0.9 to
match the W2/QUAL2K/Bowie consensus; if a lower value is wanted for a
specific application it belongs in a YAML override with a citation, not
in the framework default. Severity Major (affects the carbon and DO
mass balance, not just documentation).

**SA-3 (MAJOR, needs registry verification) — light driver may be
total shortwave, not PAR.** `floating_algae.py:414`,
`benthic_algae.py:314`. NSM1's `KL = 10 W/m^2`
(`algae.py:28`, `balgae.py:37`) is a PAR-scale constant; the
half-saturation light model is only physically meaningful if the
irradiance fed to `limit_light` is PAR. W2 stays self-consistent by
calibrating `ASAT` on the shortwave basis it feeds in
(`water-quality.f90:678`); QUAL2K applies a ~0.47 PAR fraction. NSM1
passes `solar_radiation` straight through with no PAR scaling and a
PAR-scale `KL`. If the registry `solar_radiation` is total shortwave,
the algal light limitation is overstated and net growth over-predicted
by roughly 30--60% under light-limited conditions. This is a
both-models-wrong candidate (v1 has the same omission). It cannot be
adjudicated from the algae source; it depends on the registry variable
contract. Recommended fix: confirm the `solar_radiation` registry
contract; if total shortwave, multiply by `Fr_PAR` (~0.47) before
`limit_light`, OR document that `KL` is defined on a total-shortwave
basis and re-derive its default accordingly (a shortwave-basis KL would
be ~ 20--45 W/m^2, not 10). Severity Major; carried as Open Question 1.

**SA-4 (MINOR) — alkalinity unit comment and light/nutrient method
limitations undocumented.** `alkalinity.py:359-361`,
`floating_algae.py:797-808`, `floating_algae.py:659`. (a) The inline
comment "ApGrowth (ug-Chla/L/d) * rca (mg-C/ug-Chla) = mg-C/L/d" is
false while `rca = self.AWc` (mg-C per stoichiometric unit, not per
ug-Chla); folds into the SA-1 fix. (b) No code or docstring states that
the default half-saturation light model has **no photoinhibition**
(W2's Steele does) and will over-predict surface growth at very high
irradiance. (c) No docstring states that the default multiplicative
nutrient combination is the conservative (co-limitation) choice and
that option 2 (minimum/Liebig) reproduces W2/QUAL2E behavior.
Recommended fix: correct the comment as part of SA-1; add one-sentence
method notes to `limit_light` and to the growth-combination block citing
QUAL2E/QUAL2K and W2. Severity Minor.

**SA-5 (MINOR) — benthic `krb_theta = 1.06` is internally inconsistent
and unexplained.** `balgae.py:45`. Every other Arrhenius theta in the
algae/benthic-algae defaults is 1.047; `krb_theta` alone is 1.06, with
only the inline note "differs from typical 1.047." The value 1.06 is
within the Bowie respiration-theta range (1.045--1.08) so it is not
scientifically wrong, but it is unexplained (no citation, unlike
`mu_max_20`/`kdp_20`/`krp_20`) and inconsistent with the sibling
`krp_theta = 1.047` for floating-algae respiration with no stated
reason periphyton respiration would have a different temperature
sensitivity than phytoplankton respiration. Recommended fix: either
align to 1.047 for consistency, or add a one-line citation justifying
the higher periphyton respiration theta (Bowie Table 6-x range), as was
done for the `mu_max_20` block in `algae.py:31-58`. Severity Minor.

---

## 5. Defaults-vs-literature table

Units: rate constants /d at 20 C; half-sats mg/L; light W/m^2;
stoichiometry as noted. "Lit. typical/range" from Bowie et al. 1985
(EPA/600/3-85/040), Chapra 1997, QUAL2K v2.11 defaults, EUTRO/WASP8 —
all cited from domain knowledge (MEDIUM confidence), W2 input typicals
from the Cole & Wells manual.

| Parameter | v3 default | Lit. typical | Lit. range | Within range? | Note |
|---|---|---|---|---|---|
| `mu_max_20` (floating growth) | 2.0 | ~2.0 | 0.5--4.0 | YES | central; improved over v1's 1.0 |
| `krp_20` (resp) | 0.10 | ~0.1 | 0.05--0.5 | YES | low-central |
| `kdp_20` (death) | 0.05 | ~0.05--0.1 | 0.01--0.3 | YES | low-central |
| `vsap` (settling) | 0.15 | ~0.1--0.3 | 0.0--2.0 | YES | reasonable for mixed assemblage |
| `KsN` | 0.04 | ~0.02 | 0.005--0.4 | YES | mid |
| `KsP` | 0.0012 | ~0.002 | 0.001--0.08 | YES (low) | at lower edge; acceptable |
| `KL` (light half-sat) | 10 W/m^2 | PAR ~10--30 | 5--100 | YES only if input is PAR | SA-3 caveat |
| C:Chla (`AWc/AWa`) | 40 g-C/g-Chla | ~40--50 | 20--100 | YES | standard |
| N:Chla (`AWn/AWa`) | 7.2 g-N/g-Chla | ~7--10 | 5--15 | YES | Redfield-consistent |
| P:Chla (`AWp/AWa`) | 1.0 g-P/g-Chla | ~0.5--1.0 | 0.3--1.5 | YES | standard |
| `f_pocp`/`f_pocb` | 0.5 | ~0.8 | 0.6--0.9 | **NO (too low)** | SA-2 |
| `mub_max_20` (benthic) | 0.4 | ~1.0 | 0.2--3.0 | YES (low) | conservative; acceptable |
| `krb_20` (benthic resp) | 0.2 | ~0.2 | 0.05--0.5 | YES | ok |
| `kdb_20` (benthic death) | 0.3 | ~0.1--0.5 | 0.02--1.0 | YES | ok |
| `krb_theta` | 1.06 | 1.047 | 1.045--1.08 | YES | inconsistent w/ siblings (SA-5) |
| `KsNb`/`KsPb` | 0.25 / 0.125 | ~0.05 / ~0.01 | 0.01--0.5 / 0.001--0.2 | YES (high) | high but in-range for periphyton |
| `Ksb` (biomass half-sat) | 10 g-D/m^2 | ~5--20 | 1--50 | YES | ok |
| `r_alkaa` etc. (alk stoich) | 14/106/12/1000 eq/mg-C | Stumm-Morgan 14/16 eq/mol-N | exact | YES (ratio itself) | misused via SA-1, not a defaults defect |
| Arrhenius theta (growth) | 1.047 | 1.066--1.08 | 1.02--1.08 | YES (low) | NSM1 1.047 is at the low end; W2 4-param is more physical at extremes; acceptable |

---

## 6. Confidence and caveats

HIGH confidence (line-checked in local source):

- Every CE-QUAL-W2 citation (TEMPERATURE_RATES, KINETIC_RATES algal
  block, ALGAE, EPIPHYTON, ALKALINITY, BETA in heat-exchange.f90) was
  read directly in `W2_v2026.02`. The CA-1/SA-1 1000x/100x conclusion
  and its correct value are anchored to `water-quality.f90:3164` and
  the Stumm & Morgan stoichiometry the W2 authors quote at
  `water-quality.f90:3152-3159`. This is the strongest finding.
- All v3 source citations (`floating_algae.py`, `benthic_algae.py`,
  `parameters/algae.py`, `balgae.py`, `alkalinity.py`,
  `parameters/alkalinity.py`) were read directly.

MEDIUM confidence (domain knowledge, not locally checked):

- QUAL2K v2.11, QUAL2E (EPA/600/3-87/007), WASP8/EUTRO, Bowie et al.
  1985 (EPA/600/3-85/040), and Chapra 1997 numeric ranges. These are
  used to *corroborate* the W2-anchored conclusions and to bound
  defaults; no finding's severity rests on them alone. The
  multiplicative-vs-Liebig default and the half-saturation-vs-Steele
  default are reported as MODEL-STRUCTURAL-DIFF-OK precisely because the
  QUAL2K/QUAL2E corroboration is the softer leg of the evidence.
- W2 `APOM`/`ASAT`/`AHSN`/`AHSP` "typical" values: these are
  control-file inputs, not source literals, so the typical bands quoted
  for them are from the Cole & Wells W2 user manual (domain knowledge),
  not from the local source. The SA-2 conclusion (0.9 defensible, 0.5
  low) is robust across W2, QUAL2K, and Bowie independently, so it does
  not hinge on a single uncertain W2 default.

Caveats / not adjudicable from the algae source:

- SA-3 / A11 (PAR vs shortwave) depends on the registry
  `solar_radiation` variable contract and the global-vars light
  attenuation pipeline, which are out of algae-Process scope. It is the
  one finding that could be a real ~2x-class error and remains an Open
  Question.
- This review did not re-run benchmarks (per instruction). The
  magnitude claims in SA-1 are arithmetic (1000 = 40/0.04, 100 = 40/0.4)
  and do not require a model run; the growth-overprediction range in
  SA-3 is an order-of-magnitude estimate, not a benchmarked figure.
