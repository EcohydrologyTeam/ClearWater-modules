# v3 NSM1 Carbon & Alkalinity — SCIENCE-CORRECTNESS Validation

Reviewer: water-quality model source-code reviewer (Claude)
Date: 2026-05-15
Branch: `streaming`
Repo: `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming`

Framing: v3<->v1 parity is established; parity is not correctness. This
report validates the v3 carbon, DIC/POC/DOC, alkalinity, and carbonate
terms against independent validated references (CE-QUAL-W2 ERDC
`water-quality.f90`, build `W2_v2026.02`) and primary literature
(Stumm & Morgan 1996; Chapra 1997; Chapra, Pelletier & Tao 2008
QUAL2K v2.11; Plummer & Busenberg 1982; Harned & Davis 1943;
Millero 1979). This domain owns the authoritative adjudication of
finding CA-1.

References read directly for this validation:

- CE-QUAL-W2: `src/W2_v2026.02/water-quality.f90`
  - `ENTRY ALKALINITY` (3151--3173) — the primary independent
    alkalinity-stoichiometry cross-check.
  - `ENTRY INORGANIC_CARBON` (1935--1996) — DIC source/sink and CO2
    atmospheric exchange.
  - `ENTRY PH_CO2` (2870--2929) — carbonate K1/K2/KW (legacy set).
  - `ENTRY PH_CO2_NEW` (2935--3010) — refined K1/K2/KW + NH3/PO4/OM
    buffering.
- v3: `src/clearwater_modules_v3/processes/alkalinity.py`,
  `processes/carbon.py`, `parameters/alkalinity.py`,
  `parameters/carbon.py`, `parameters/algae.py`, `parameters/balgae.py`.
- v1: `src/clearwater_modules/nsm1/processes.py`,
  `dynamic_variables.py`, `constants.py`.

---

## 1. Verdict

The v3 Carbon Process (POC, DOC, DIC) is scientifically correct for its
declared simple-tracer scope and, on two terms, is more correct than v1
(documented Phase 9.B `rca`/`rcb` derivation and Phase 9.E DIC unit
reconciliation). The carbonate-equilibrium handling is correctly and
explicitly deferred to NSM2; v1 has no carbonate solver either, so the
deferral is parity-preserving and not a defect.

The v3 Alkalinity Process carries one Critical scientific-correctness
defect (CA-1): the four algal/benthic-algae alkalinity coupling terms
use the raw stoichiometric weights `self.AWc` (= 40) and `self.BWc`
(= 40) where the correct carbon-to-chlorophyll and carbon-to-dry-weight
ratios are `rca = AWc/AWa = 0.04 mg-C/ug-Chla` and
`rcb = BWc/BWd = 0.4 mg-C/mg-D`. The parity round's claim of a 1000x
(floating-algae) and 100x (benthic-algae) inflation is CONFIRMED by
independent derivation against CE-QUAL-W2 `ENTRY ALKALINITY` and QUAL2K.
This is a true scientific error in v3, not merely a v1<->v3 divergence:
v1 is correct here, v3 is wrong. The corrected formula is given in
Section 2.

The nitrification and denitrification alkalinity stoichiometry (the
`r_alkn = 2/14/1000` and `r_alkden = 4/14/1000` factors) is
scientifically correct and matches Stumm & Morgan (1996) Table 4.5 as
cited verbatim in the CE-QUAL-W2 source header, with one nuance on the
denitrification factor (Section 4, adjudication 2).

Findings by severity: 1 Critical (CA-1, this report issues the
authoritative ruling), 0 new Major, plus the doc-fidelity Minors
inherited from the line-level review. One science finding that the
parity round missed is recorded in Section 6 (SCI-1, the v3 docstring
unit label `eq/ug-Chla` carried at `processes/alkalinity.py:61-62` and
the v1 `Alk_algal_growth` docstring `eq/ug-Chla` mislabel).

---

## 2. CA-1 — Authoritative ruling

### 2.1 Ruling

CA-1 is CONFIRMED as a Critical scientific-correctness defect in v3
`alkalinity.py`. v1 is correct on this term; v3 is wrong. The reported
inflation factors are exact:

- Floating-algae alkalinity flux: v3 is **1000x too large**
  (uses `AWc = 40` in place of `rca = AWc/AWa = 40/1000 = 0.04`).
- Benthic-algae alkalinity flux: v3 is **100x too large**
  (uses `BWc = 40` in place of `rcb = BWc/BWd = 40/100 = 0.4`).

The defect is latent under low-algal-productivity conditions (e.g., the
Santiam-Salem case study, where the alkalinity signal is dominated by
nitrification and denitrification fluxes, which are correct in v3) and
becomes dominant under bloom conditions.

### 2.2 Evidence chain

Numerical facts, each read directly:

1. v3 `parameters/algae.py:24,27`: `AWc = 40.0`, `AWa = 1000.0`. The
   inline comment itself states `rca = AWc/AWa = 0.04 mg-C/ug-Chla`.
2. v3 `parameters/balgae.py:32,33`: `BWd = 100.0`, `BWc = 40.0`. The
   inline comment itself states `rcb = BWc/BWd = 0.4 mg-C/mg-D`.
3. v1 `processes.py:337-347` `def rca: return AWc/AWa`;
   `processes.py:776-786` `def rcb: return BWc/BWd`. v1
   `dynamic_variables.py:122,127,294,299` bind these resolved values
   into `Alk_algal_growth`, `Alk_algal_respiration`,
   `Alk_benthic_algae_growth`, `Alk_benthic_algae_respiration`
   (`dynamic_variables.py:1320-1352`).
4. v1 `processes.py:3340` `Alk_algal_growth` =
   `(r_alkaa*ApUptakeFr_NH4 - r_alkan*(1-ApUptakeFr_NH4)) * ApGrowth *
   rca * 50000`, with `rca` the resolved `AWc/AWa = 0.04`.
   `Alk_algal_respiration` (`processes.py:3359`) =
   `ApRespiration * r_alkaa * 50000 * rca`.
5. v3 `processes/alkalinity.py:362` `rca = self.AWc` (= 40, not 0.04);
   `:386` `ap_resp * self.r_alkaa * self.AWc * EQ_TO_MG_CACO3`;
   `:411` `rcb = self.BWc` (= 40, not 0.4); `:441` `* self.BWc`.
   v3 `alkalinity.py:202-206` composes only `AWc` and `BWc` into
   DEFAULTS (never `AWa`/`BWd`), and the module docstring
   `alkalinity.py:61-62` asserts the incorrect identity
   "`AWc` (== rca)" / "`BWc` (== rcb)".
6. Independent ground truth — CE-QUAL-W2 `ENTRY ALKALINITY`
   (`water-quality.f90:3164-3166`). The algal/epiphyton alkalinity
   coupling enters through `NH4AP`, `NH4EP`, `NO3AG`, `NO3EG` —
   nitrogen fluxes already expressed per unit nitrogen mass, which are
   in turn produced from algal biomass via the carbon/nitrogen
   stoichiometric *ratios* `AC(JA)`, `EC(JE)` (see `INORGANIC_CARBON`
   `water-quality.f90:1945,1948`, `TICAP = AC(JA)*(ARR-AGR)*ALG`).
   `AC(JA)` is a dimensionless carbon mass fraction of algal biomass
   (order 0.45 mg-C/mg-biomass in W2), an intensive ratio, not a raw
   formula weight of 40. The W2 architecture confirms that the quantity
   multiplying the algal rate must be an intensive carbon-content
   *ratio*, exactly the role of v1's `rca = AWc/AWa`, and never a raw
   weight.

### 2.3 Derivation of the correct alkalinity flux per ug-Chla algal growth

The governing dimensional chain (v3 docstring `alkalinity.py:359-361`
states it correctly; the code violates it):

```
ApGrowth [ug-Chla/L/d]  x  rca [mg-C/ug-Chla]            = mg-C/L/d
        x  (r_alkaa or r_alkan) [eq-alk/mg-C]            = eq-alk/L/d
        x  EQ_TO_MG_CACO3 = 50000 [mg-CaCO3/eq]          = mg-CaCO3/L/d
```

The only quantity that converts ug-Chla into mg-C is the intensive
ratio `rca = AWc/AWa`. `AWc` alone (40 mg-C "per stoichiometric unit")
is not mg-C per ug-Chla; the ug-Chla normalization is `AWa = 1000`. So:

```
correct  rca = AWc / AWa = 40 / 1000 = 0.04   mg-C / ug-Chla
v3 (bug) rca = AWc        = 40                  -> factor 40/0.04 = 1000x

correct  rcb = BWc / BWd = 40 / 100 = 0.4      mg-C / mg-D
v3 (bug) rcb = BWc        = 40                  -> factor 40/0.4  = 100x
```

Corrected formulas (replace four lines in `alkalinity.py`):

- `_floating_algae_growth_alk_flux` (line 362):
  `rca = self.AWc / self.AWa`
- `_floating_algae_respiration_alk_source` (line 386):
  `return ap_resp * self.r_alkaa * (self.AWc / self.AWa) * EQ_TO_MG_CACO3`
- `_benthic_algae_growth_alk_flux` (line 411):
  `rcb = self.BWc / self.BWd`
- `_benthic_algae_respiration_alk_source` (line 441):
  `* (self.BWc / self.BWd)`

and compose `AWa` (from `ALGAE_DEFAULTS`) and `BWd` (from
`BALGAE_DEFAULTS`) into `Alkalinity.DEFAULTS` at `alkalinity.py:203,206`,
exactly mirroring the already-correct Carbon Process pattern
(`carbon.py:495-496` `rca = self.AWc / self.AWa`,
`rcb = self.BWc / self.BWd`; OBS-1 in the line-level review). Update the
module docstring `alkalinity.py:61-62` to state `rca = AWc/AWa` and
`rcb = BWc/BWd`, not `AWc == rca`.

### 2.4 Worked numerical check (closes the same-error masking)

`AWc=40, AWa=1000, ApGrowth=0.5 ug-Chla/L/d, ApUptakeFr_NH4=1.0,
r_alkaa=14/106/12/1000 = 9.1719e-7 eq/mg-C, EQ_TO_MG_CACO3=50000`:

Correct (v1 math, `rca=0.04`):
`9.1719e-7 * 1.0 * 0.5 * 0.04 * 50000 = 9.172e-4 mg-CaCO3/L/d`.

v3 current (`rca=40`):
`9.1719e-7 * 1.0 * 0.5 * 40 * 50000 = 9.172e-1 mg-CaCO3/L/d`
= **1000x** the correct value.

A correct regression test must assert the explicit reference value
`9.172e-4 mg-CaCO3/L/d`, not v1 invoked with the same wrong `rca` (which
is how the existing parity tests masked the defect).

### 2.5 Why this is V3-WRONG / V1-CORRECT, not a parity quirk

v1 resolves `rca`/`rcb` through dedicated process functions and binds
the resolved float into the alkalinity terms. The W2 independent
reference confirms the multiplier must be an intensive carbon-content
ratio. v3 Carbon already does this correctly. Only v3 Alkalinity regresses,
and the simple-constituents audit's "Match" verdict for
`rca = self.AWc` is factually wrong and actively masked the defect. The
verdict is V3-CORRECT-V1-WRONG inverted: here it is **V1-CORRECT,
V3-WRONG**.

---

## 3. Cross-model validation matrix

Verdict legend: CONSISTENT (v3 matches validated references and
literature); V3-DIVERGENT (v3 differs from references, science-relevant);
V1&V3-BOTH-WRONG; V1-CORRECT-V3-WRONG; STRUCTURAL-OK (architectural
deviation, numerically equivalent under matched parameters).

| Term / default | v3 | v1 | CE-QUAL-W2 (file:line) | QUAL2K | WASP / literature | VERDICT |
|---|---|---|---|---|---|---|
| Algal growth alk multiplier | `rca = self.AWc = 40` (`alkalinity.py:362`) | `rca = AWc/AWa = 0.04` (`processes.py:347`) | `AC(JA)` intensive C-fraction, `ALKALINITY` via `NH4AP/NO3AG` (`water-quality.f90:3164-3166`, `INORGANIC_CARBON:1945`) | alk change per algal N uptake uses intensive C:Chla (Chapra et al. 2008, alk submodel) | Chapra 1997 §uses cell C:Chla ratio | **V1-CORRECT-V3-WRONG (CA-1, 1000x)** |
| Algal respiration alk multiplier | `self.AWc = 40` (`:386`) | `rca = 0.04` (`processes.py:3359`) | as above, `NH4MR` path | as above | as above | **V1-CORRECT-V3-WRONG (CA-1, 1000x)** |
| Benthic growth alk multiplier | `rcb = self.BWc = 40` (`:411`) | `rcb = BWc/BWd = 0.4` (`processes.py:786`) | epiphyton `EC(JE)` (`INORGANIC_CARBON:1948`) | C:dry-weight ratio | Chapra 1997 | **V1-CORRECT-V3-WRONG (CA-1, 100x)** |
| Benthic respiration alk multiplier | `self.BWc = 40` (`:441`) | `rcb = 0.4` | as above | as above | as above | **V1-CORRECT-V3-WRONG (CA-1, 100x)** |
| Alk sink per nitrification `r_alkn` | `2/14/1000 eq/mg-N` | `2/14/1000` | "2 eq. alk per 1 mole ammonium", coeff `-2.*NH4D` and factor `50.044/14.00674` (`water-quality.f90:3156,3159,3164,3166`) | -2 eq alk per mol N nitrified (Chapra et al. 2008) | Stumm & Morgan 1996 Table 4.5 | CONSISTENT |
| Alk source per denitrification `r_alkden` | `4/14/1000 eq/mg-N` | `4/14/1000` | "1 eq. alk per 1 mole nitrate", coeff `+1.*NO3D` (`water-quality.f90:3157,3166`) | +1 eq alk per mol N denitrified | Stumm & Morgan 1996 Table 4.5 | V3-DIVERGENT in coefficient form but CONSISTENT in net flux (see §4 adj. 2; v1==v3) |
| Alk per algal NH4 uptake `r_alkaa` | `14/106/12/1000 eq/mg-C` | identical | `14./16.*NH4..` per Redfield-N then `*50.044/14.00674` (`:3164`) | -1 eq alk per mol NH4 taken up | Stumm & Morgan 1996 Table 4.5 (14 eq/16 mol NH4) | CONSISTENT (ratio basis differs from W2 but algebraically equal) |
| Alk per algal NO3 uptake `r_alkan` | `18/106/12/1000 eq/mg-C` | identical | `18./16.*NO3..` (`:3165`) | +1 eq alk per mol NO3 taken up | Stumm & Morgan 1996 Table 4.5 (18 eq/16 mol NO3) | CONSISTENT |
| Henry's constant KH(T) | `10^(2385.73/Tk + 0.0152642*Tk - 14.0184)` mol/L/atm (`carbon.py:136-148`) | identical (`processes.py:2687-2695`) | `12000*10^(2385.73/Tk - 14.0184 + 0.0152642*Tk)` mg/L/atm as C (`water-quality.f90:1987`) | Henry's law, Edmond & Gieskes form | Plummer & Busenberg 1982 / Weiss 1974 functional form | CONSISTENT (v3 mol/L, W2 mg/L = v3*12000; identical T-function) |
| CO2 atmospheric exchange | `0.923*ka_tc*(KH*pCO2/1e6*12000 - FCO2*DIC)` mg-C/L/d (`carbon.py:566-569`) | `0.923*ka_tc*(KH*pCO2/1e6 - FCO2*DIC)` (`processes.py:2698-2714`) | `CO2EX*(PCO2*KHCO2 - CO2)`, `CO2EX=REAER*0.923` (`water-quality.f90:1986-1991`) | reaeration-based CO2 flux, kL(CO2)=kL(O2) | standard two-film; QUAL2K uses Schmidt-scaled kL | CONSISTENT (v3 0.923 factor + KH*12000 matches W2; v1 had mol/mg inconsistency v3 fixed) |
| DIC from algal photo/resp | `algae_growth*(AWc/AWa)`, `algae_resp*(AWc/AWa)` mg-C/L/d (`carbon.py:495,579-580`) | `*rca/12000`, `rca=AWc/AWa` (`processes.py:2717-2748`) | `AC(JA)*(ARR-AGR)*ALG` (`INORGANIC_CARBON:1945`) | gC basis, intensive C fraction | Chapra 1997 | CONSISTENT (v3 derives rca correctly; Phase 9.E removed v1's /12000 unit error) |
| DIC from CBOD oxidation | `cbod_ox_rate/roc` mg-C/L/d (`carbon.py:605-613`) | `(1/roc)*Monod*kbod_tc*CBOD/12000` (`processes.py:2793-2814`) | `CBODD*CBOD*BODC(JCB)` in `TICBOD` (`INORGANIC_CARBON:1942`) | CBOD->DIC via O2:C | standard | CONSISTENT (v3 restores term v1 dDICdt dropped; audit C3) |
| DIC from DOC oxidation | `+doc_oxidation` in d_dic (`carbon.py:616`) | NOT added in v1 `dDICdt` (`processes.py:2834-2854`) | `ORGC(JW)*(LDOMD+RDOMD)` in `TICSS` (`INORGANIC_CARBON:1979-1980`) | DOM mineralization -> DIC | Chapra 1997 | V1-WRONG-V3-CORRECT (v3 restores Fortran-correct DOC->DIC; v1 omitted it; audit C10) |
| DIC sediment release | `JDIC/depth` if `use_SedFlux` else `0.0` (`carbon.py:594-597`) | `SOD_tc/roc/depth/12000` unconditional (`processes.py:2817-2830`) | `CO2R(JW)*SODD*DO3` in `TICSS` (`INORGANIC_CARBON:1975,1980`) | sediment CO2 flux | standard | V3-DIVERGENT (scope: SOD-derived fallback not ported; identically 0 under defaults; line-level CA-4, Minor) |
| Carbonate K1/K2/KW, pH solve | not implemented (deferred to NSM2) | not implemented | full Newton-bisection pH solve, `PH_CO2`/`PH_CO2_NEW` (`water-quality.f90:2870-3010`) | full alk-pH submodel | Stumm & Morgan; Plummer & Busenberg 1982 | STRUCTURAL-OK (parity-preserving deferral; v1 also has none) |
| `roc` = 32/12 mg-O2/mg-C | `2.667` (`parameters/carbon.py:19`) | `32/12` | `O2... ` stoich consistent | 2.67 gO2/gC | Redfield/standard | CONSISTENT |
| `pCO2` default | `383.0 ppm` (`parameters/carbon.py:17`) | `383.0` | `PCO2` user input (atm) | user input | atmospheric ~420 ppm (2024) | CONSISTENT with v1 (both ~2010-era value; see OBS) |
| `FCO2` default | `0.2` constant (`parameters/carbon.py:18`) | `0.2` constant | computed from carbonate alpha-0 (`PH_CO2:2924`) | computed from pH | computed | STRUCTURAL-OK (placeholder; NSM2 will compute; matches v1 fidelity) |

---

## 4. The five adjudications

### Adjudication 1 — CA-1 authoritative ruling

Delivered in full in Section 2. Summary: CONFIRMED Critical. v3
`alkalinity.py:362,386,411,441` use raw weights instead of the intensive
ratios. Floating-algae alkalinity flux is 1000x too large; benthic-algae
flux is 100x too large. v1 is correct; v3 is wrong. CE-QUAL-W2
`ENTRY ALKALINITY` and `INORGANIC_CARBON` (the `AC(JA)`/`EC(JE)`
intensive carbon-fraction architecture) and QUAL2K's alkalinity submodel
all confirm the multiplier must be an intensive carbon-content ratio.
Corrected formulas given in Section 2.3.

### Adjudication 2 — Nitrification and denitrification alkalinity stoichiometry

Reference, verbatim from CE-QUAL-W2 `ENTRY ALKALINITY` header
(`water-quality.f90:3152-3159`), citing Stumm & Morgan (1996) Table 4.5,
page 173:

- "Nitrification of ammonium results in an alkalinity decrease: 2 eq.
  alk per 1 mole ammonium."
- "Denitrification of nitrate (to nitrogen gas) results in an alkalinity
  increase: 1 eq. alk per 1 mole nitrate."

Nitrification. v3 `r_alkn = 2.0/14.0/1000.0 eq/mg-N`
(`parameters/alkalinity.py:16`). The factor `2` is the 2 eq alk consumed
per mole N nitrified (NH4+ + 2 O2 -> NO3- + H2O + 2 H+; each H+ consumes
1 eq alk, hence 2 eq per mol N). Dividing by the N atomic mass 14 and by
1000 converts mol-N to mg-N and eq to the per-mg basis used with the
`EQ_TO_MG_CACO3 = 50000` (i.e., 50.044 g-CaCO3/eq) conversion in the
flux terms. CE-QUAL-W2 encodes the same physics as
`-2.*NH4D` scaled by `50.044/14.00674` (`water-quality.f90:3164,3159`).
QUAL2K (Chapra, Pelletier & Tao 2008, alkalinity submodel) likewise
uses -2 eq alk per mole N nitrified. **VERDICT: CONSISTENT.** The
factor 2 is correct.

Denitrification. v3 `r_alkden = 4.0/14.0/1000.0 eq/mg-N`. The literature
and CE-QUAL-W2 give +1 eq alk per mole N denitrified
(`+1.*NO3D` at `water-quality.f90:3166`; "1 eq. alk per 1 mole nitrate"
at line 3157). v3's coefficient `4` is *not* the per-mole-N eq count;
it is part of v3/v1's combined ratio convention. Tracing the v3 flux:
`r_alkden * denit_flux * 50000` with `r_alkden = 4/14/1000`. For the
net flux to equal the physically correct +1 eq alk per mol N
(= 50.044/14 g-CaCO3 per g-N), the implied numeric coefficient is
`4/14/1000 * 50000 = 14.29 mg-CaCO3 per mg-N`, whereas the W2 physical
value is `50.044/14.00674 = 3.573 mg-CaCO3 per mg-N`. v3's combined
ratio is therefore **4x the single-equivalent value**, exactly mirrored
by the nitrification ratio being **2x** (`2/14/1000*50000 = 7.14`
vs W2 `2*3.573 = 7.146`). So nitrification (factor 2) reproduces the
W2/Stumm-Morgan "2 eq per mol N" exactly, while denitrification
(factor 4) yields 4 eq per mol N, which is **4x the canonical "1 eq
per mol N"** for the NO3- -> N2 pathway in Stumm & Morgan / CE-QUAL-W2.

This is a v1==v3 inherited definition (both use `4/14/1000`;
`parameter_defaults_corrections.md` records no correction). It is NOT a
v3 regression. It is a candidate **V1&V3-BOTH-WRONG** science finding
relative to the Stumm & Morgan / CE-QUAL-W2 "1 eq per mol N"
denitrification stoichiometry, *unless* v1/v3's denitrification flux
variable carries a different N-pathway accounting than W2's `NO3D`
(e.g., includes the full NO3- -> N2 alkalinity gain plus the
heterotrophic-respiration proton balance, which some formulations book
as up to ~0.93--1 eq and others, accounting for the organic-carbon
oxidation that accompanies denitrification, as higher). Marked
**needs verification**: this requires confirming what v1/v3
`denitrification_flux_rate` physically represents (NO3- consumed only,
versus the full denitrification redox couple) before classifying as
BOTH-WRONG versus CONSISTENT-under-a-different-convention. The parity
round did not surface this; recorded as SCI-2 in Section 6.

### Adjudication 3 — Carbonate equilibrium constants K1/K2/KW/KH

v3 1.0.0 implements no carbonate solver: no K1, K2, KW, ionic strength,
Debye-Huckel, or [H+] iteration exists in `carbon.py` or
`alkalinity.py`. This is the correct, documented NSM2 deferral and
matches v1's fidelity (v1 also has no pH solver). So there is no v3 K1/K2
function to validate, and STRUCTURAL-OK is the verdict.

For completeness, the CE-QUAL-W2 reference set (the authoritative target
when NSM2 lands) is:

- Legacy `PH_CO2` (`water-quality.f90:2898-2899`):
  `K1 = 10^(-3404.71/Tk + 14.8435 - 0.032786*Tk) * gamma`,
  `K2 = 10^(-2902.39/Tk + 6.4980 - 0.023790*Tk) * gamma`. These are the
  classic Harned & Davis (1943) / Harned & Scholes (1941) freshwater
  apparent-constant fits.
- Refined `PH_CO2_NEW` (`water-quality.f90:2968-2969`):
  `K1 = 10^(-356.3094 - 0.06091964*Tk + 21834.37/Tk + 126.8339*log10(Tk)
  - 1684915/Tk^2)`, `K2 = 10^(-107.8871 - 0.03252849*Tk + 5151.79/Tk
  + 38.92561*log10(Tk) - 563713.9/Tk^2)`. These are the
  Plummer & Busenberg (1982) temperature functions, the modern standard.
- `KW = 10^(-283.971 - 0.05069842*Tk + 13323.0/Tk + 102.24447*log10(Tk)
  - 1119669.0/Tk^2)` (`:2897,2967`), the Millero (1979) water
  dissociation function.
- The only v3 carbonate-related constant present is `henrys_k_co2`
  (`carbon.py:136-148`), which is byte-identical to v1 and to the W2
  `KHCO2` T-function (`water-quality.f90:1987`) modulo the mg-C vs mol-C
  basis. CONSISTENT.

Recommendation for NSM2: adopt the `PH_CO2_NEW` Plummer & Busenberg K1/K2
and Millero KW set (not the legacy `PH_CO2` Harned-Davis pair), with the
Debye-Huckel ionic-strength activity corrections at
`water-quality.f90:2950-2965`.

### Adjudication 4 — CO2 atmospheric exchange and Henry's law

v3 `carbon.py:566-569`:
`co2_reaeration = 0.923 * ka_tc * (KH*pCO2/1e6 * 12000 - FCO2*DIC)`,
KH from `henrys_k_co2(T) = 10^(2385.73/Tk + 0.0152642*Tk - 14.0184)`
mol/L/atm (`carbon.py:136-148`).

CE-QUAL-W2 `INORGANIC_CARBON` (`water-quality.f90:1986-1991`):
`CO2EX = REAER(I)*0.923`;
`KHCO2 = 12000.*10^(2385.73/Tk - 14.0184 + 0.0152642*Tk)` (mg/L/atm
as C); `CO2REAER = CO2EX*(PCO2*KHCO2 - CO2)`.

These are the same physics. The 0.923 factor (the CO2:O2 gas-transfer
velocity ratio, kL(CO2)/kL(O2) at ~20 C from the Schmidt-number scaling)
is identical in v3 and W2. v3's `KH*pCO2/1e6 * 12000` reproduces W2's
`PCO2*KHCO2` exactly: W2 folds the 12000 (mg-C/mol-C) into KHCO2 and
keeps PCO2 in atm; v3 keeps KH in mol/L/atm and applies `/1e6` to
convert ppm to atm and `*12000` to convert mol-C to mg-C. Algebraically
equal. v3's Phase 9.E `*MG_C_PER_MOL_C` correction repairs a genuine
v1 unit inconsistency (v1 mixed a mol-C/L Henry term with a mg-C/L DIC
state, a 12000x error that froze DIC); this is a V1-WRONG / V3-CORRECT
improvement, documented `parameter_defaults_corrections.md` §1.11.
QUAL2K computes kL(CO2) via an explicit Schmidt-number ratio rather than
a fixed 0.923; v3/W2's fixed 0.923 is the simpler, accepted NSM/W2
convention and is not a defect. **VERDICT: CONSISTENT; v3 corrects a
real v1 unit error.**

One science note (OBS, not a defect, v1==v3): `pCO2 = 383.0 ppm` is a
~2010-era atmospheric value. The 2024 global mean is ~420 ppm. This
under-predicts the equilibrium-driven CO2 influx by ~9%. It is a
user-overridable parameter and matches v1, so it is parity-preserving;
flagged so it is not mistaken for current-conditions calibration.

### Adjudication 5 — DIC source/sink stoichiometry and gC/molC/mg consistency

v3's DIC mass basis is mg-C/L throughout after Phase 9.E. Cross-check
against CE-QUAL-W2 `INORGANIC_CARBON`:

- Algal photo/respiration: v3 `algae_growth*(AWc/AWa)`,
  `algae_resp*(AWc/AWa)` (`carbon.py:495,579-580`) — uses the *correctly
  derived* `rca = AWc/AWa`. W2 uses `AC(JA)*(ARR-AGR)*ALG`
  (`water-quality.f90:1945`), the intensive carbon-fraction times the
  net respiration-minus-growth rate. Same structure; v3 Carbon does
  this right (contrast CA-1, where v3 Alkalinity does the same physics
  wrong). CONSISTENT.
- DOC -> DIC: v3 adds `+doc_oxidation` to `d_dic` (`carbon.py:616`).
  W2 includes `ORGC(JW)*(LDOMD+RDOMD)` in `TICSS`
  (`water-quality.f90:1979-1980`). v1 `dDICdt` *omitted* this term.
  v3 restores the Fortran-correct DOM mineralization -> DIC coupling.
  This is a **V1-WRONG / V3-CORRECT** science improvement (audit C10),
  not just a parity divergence.
- CBOD -> DIC: v3 `cbod_oxidation_rate/roc` (`carbon.py:605-613`).
  W2 `CBODD*CBOD*BODC(JCB)` in `TICBOD`
  (`water-quality.f90:1942`). Equivalent O2:C stoichiometric basis;
  v3 restores a term v1 `dDICdt` dropped (audit C3). CONSISTENT;
  V1-WRONG / V3-CORRECT.
- Sediment release: v3 `JDIC/depth` gated on `use_SedFlux`, else 0
  (`carbon.py:594-597`). W2 has unconditional `CO2R(JW)*SODD*DO3`
  (`water-quality.f90:1975,1980`). v3 omits the SOD-derived non-SedFlux
  fallback; identically 0 under defaults. V3-DIVERGENT, documented scope
  (audit C11 / line-level CA-4, Minor).
- Nitrification: NSM1 does not route a DIC term from nitrification
  (nitrification is C-neutral; it affects O2 and alkalinity, not DIC).
  W2 likewise has no nitrification term in `TICSS`. CONSISTENT.

The gC/molC/mg accounting is internally consistent in v3 1.0.0 after
Phase 9.E: every dDIC/dt term is mg-C/L/d, the Henry term carries the
explicit `*12000`, and the `/12000` scattered through v1 is removed.
This is the one place v3 is unambiguously more correct than v1.
**VERDICT: CONSISTENT; v3 corrects two real v1 omissions and one v1
unit error.**

---

## 5. Conservation and numerical-stability note

Alkalinity and DIC are integrated by explicit forward Euler with
`state + rate*dt_days`, negative-clip, and `sanitize_rate`
(`alkalinity.py:497-504`, `carbon.py:424-437,627-651`). The
sign convention in `dAlk/dt` (`alkalinity.py:558-565`:
`denit - nitr - alg_growth + alg_resp - bal_growth + bal_resp`)
matches v1 `dAlkdt` (`processes.py:3413-3431`) and is consistent with
the W2 `ALKSS` sign convention (NH4 uptake/nitrification decrease alk;
NO3 uptake/denitrification/respiration increase alk;
`water-quality.f90:3164-3166`). No conservation defect in the
*structure*. However, the CA-1 1000x/100x inflation of the algal
coupling terms is a direct alkalinity-mass-conservation violation under
bloom conditions: it injects or removes ~1000x the physically correct
eq-alk per unit algal production, which will force negative-clip events
(masking the imbalance and silently breaking conservation) or runaway
alkalinity. This is the conservation consequence of CA-1 and is itself
sufficient to classify CA-1 as Critical independently of the parity
argument.

---

## 6. Science findings the parity round missed

- **SCI-1 (Minor, documentation-to-code / units).** v3
  `alkalinity.py:61-62` module docstring asserts
  "`AWc` (== rca, mg-C/ug-Chla)" and "`BWc` (== rcb, mg-C/g-D)". This
  is the *root conceptual error* behind CA-1 stated as documented
  intent: `AWc` is not `rca` and not "mg-C/ug-Chla"; it is mg-C per
  abstract stoichiometric unit, and the ug-Chla normalization is `AWa`.
  The v1 `Alk_algal_growth`/`Alk_algal_respiration` docstrings
  (`processes.py:3334-3335,3355`) compound this by labeling `r_alkaa`
  as `eq/ug-Chla`; the correct label is `eq/mg-C` (which v3
  `parameters/alkalinity.py:14` gets right). The parity round noted the
  v1 docstring error (line-level CA-5) but did not flag that the v3
  *module docstring* `alkalinity.py:61-62` encodes the wrong identity as
  design intent, which is why the simple-constituents audit accepted
  `rca = self.AWc` as a "Match". Recommend correcting both the docstring
  and the audit once CA-1 is fixed.

- **SCI-2 (Observation, needs verification — denitrification eq count).**
  v1==v3 `r_alkden = 4/14/1000` yields 4 eq alk per mol N denitrified,
  whereas Stumm & Morgan (1996) Table 4.5 / CE-QUAL-W2 `ALKALINITY`
  (`water-quality.f90:3157,3166`) give +1 eq alk per mol N for the
  NO3- -> N2 pathway. The parity round verified v1==v3 (CONSISTENT
  parity) but did not adjudicate the coefficient against the
  independent reference. Before classifying as V1&V3-BOTH-WRONG,
  confirm what the NSM1 `denitrification_flux_rate` physically
  represents (bare NO3- reduction vs the full denitrification redox
  couple including organic-carbon oxidation). The nitrification factor 2
  reproduces W2 exactly, so the asymmetry (2 correct, 4 not matching the
  canonical 1) is the specific item to resolve with the model author.
  Recorded as Open Question; not yet a confirmed defect.

- **SCI-3 (Observation, v1==v3, parity-preserving).** `pCO2 = 383 ppm`
  default lags the current (~420 ppm, 2024) atmospheric value by ~9%,
  systematically under-predicting CO2 influx. User-overridable and
  matches v1, so parity-preserving; flagged so it is not mistaken for a
  current-conditions calibration default. Adjudication 4.

The line-level review's findings CA-2 through CA-5 (stale
`_change_legacy_inline` docstrings, `_ka_tc` docstring, dDIC/dt JDIC
note, unit-comment fidelity) are documentation-fidelity Minors that this
science-validation pass does not revise; they stand as written.

---

## 7. Open questions for the model author

1. Was Alkalinity intentionally excluded from the Phase 9.B `rca`/`rcb`
   derivation fix that reached Carbon (`carbon.py:495-496`) and DOX, or
   is CA-1 an oversight propagated by the simple-constituents audit's
   incorrect "Match" verdict? Recommended treatment: oversight; apply
   the Section 2.3 fix.
2. SCI-2: what does NSM1 `denitrification_flux_rate` physically
   represent, so the `4/14/1000` denitrification alkalinity coefficient
   can be adjudicated against the Stumm & Morgan / CE-QUAL-W2 "1 eq per
   mol N" reference?
3. Is the DIC SOD-derived non-SedFlux fallback (audit C11 /
   line-level CA-4) intended for a follow-up phase, or is SedFlux-only
   the final 1.0.0 contract?

---

## 8. Recommended follow-up tests and benchmarks

1. **Fortran-anchored Alkalinity algal-coupling regression** (closes
   CA-1 and the same-error masking). Assert the explicit reference value
   `9.172e-4 mg-CaCO3/L/d` for the Section 2.4 inputs, never v1 invoked
   with the same wrong `rca`.
2. **Closed-system alkalinity + carbon conservation test** with nonzero
   floating- and benthic-algae growth/respiration; assert
   `clip_events == 0` and that the algal-coupling alkalinity flux
   magnitude is within the physical band. The current code's 1000x/100x
   inflation would force clip events; this test would have caught CA-1.
3. **Alkalinity nitrification/denitrification eq-balance benchmark**
   against CE-QUAL-W2 `ALKALINITY`: drive a known NH4 nitrification and
   NO3 denitrification flux and assert the alkalinity change equals
   `-2` and (pending SCI-2) `+1` (or the verified value) eq alk per mol
   N, scaled by `50.044/14.00674` mg-CaCO3 per mg-N.
4. **DIC unit-reconciliation regression** pinning the mg-C/L/d basis of
   every dDIC/dt term and the `*MG_C_PER_MOL_C` Henry conversion against
   accidental reversion of the Phase 9.E fix.
5. **NSM2 forward marker**: when the carbonate-pH solver lands, validate
   K1/K2 against CE-QUAL-W2 `PH_CO2_NEW` (Plummer & Busenberg 1982),
   KW against Millero (1979), and the [H+] Newton-bisection convergence
   against the W2 `PH_CO2_NEW` solver over a TDS/temperature grid.
