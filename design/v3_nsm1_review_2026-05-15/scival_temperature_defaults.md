# v3 NSM1 Temperature-Correction Framework and Parameter-Default Library -- Science-Correctness Validation

**Review date:** 2026-05-15
**Reviewer:** water-quality-model-source-code-reviewer agent
**Branch:** streaming
**Commit:** 54f2b12
**Scope:** Arrhenius/van't Hoff temperature-correction primitive (`utils/conversions.py`),
the full v3 NSM1 parameter-default library (`parameters/*.py`), and the
team's own corrections log (`parameter_defaults_corrections.md`).
**References used:** CE-QUAL-W2 v2026.02 (LOCAL, authoritative:
`water-quality.f90` `TEMPERATURE_RATES`/`KINETIC_RATES`, `w2modules.F90`
`FR`/`FF`); v1 NSM1 `src/clearwater_modules/nsm1/constants.py` (LOCAL parity
reference); QUAL2K (Chapra, Pelletier & Tao 2008), QUAL2E (Brown & Barnwell
1987, EPA/600/3-87/007), WASP, Bowie et al. 1985 (EPA/600/3-85/040), Chapra
1997 (cited from domain knowledge, NOT read locally -- see confidence caveats).

---

## 1. Verdict

The v3 NSM1 temperature-correction framework is **scientifically correct in
its implementation**. The Arrhenius/van't Hoff primitive
`reaction_kinetics * theta ** (water_temperature - 20.0)` in
`utils/conversions.py:42` is the exact, canonical θ^(T-20) form with the
correct reference temperature (20 °C, in Celsius, not Kelvin-contaminated)
and the correct sign. There is no reference-temperature, sign, or
per-day/per-second θ-base error in the primitive itself.

The parameter-default library is in **good shape after the Phase 9.x
corrections**, with one genuine residual defect (CBOD settling θ, F-1, MAJOR
but currently dormant) and a small number of θ values that sit at the edge
of or just outside the conventional literature range but are individually
defensible. Of the roughly 70 rate constants and θ values audited, **zero
are off by a 10×/100×/1000× unit slip in the current code** (all four
historical unit-slip defects -- `lambdam`, `vb`, `pressure_mb`, the
sentinel-999 family -- are corrected). The most consequential
science-correctness issue this review surfaces that was *not* a parity
concern is the **CBOD settling formulation/θ mismatch (F-1)**, which is a
real form and θ error that is silent only because the default is zero.

The 16 claimed corrections in `parameter_defaults_corrections.md` were
audited individually. **Fifteen of sixteen move the value toward the
literature/W2-consistent value or are defensibly documented.** One (the
nitrogen θ transposition, item 1.10) is the highest-value finding in this
review: it is a real v1 bug that v3 inherited and then **correctly fixed**,
so v3 is now *more* correct than v1 -- a case where parity would have been
wrong.

**Framework-acceptability ruling:** θ^(T-20) is a scientifically acceptable
simplification of W2's 4-point Thornton & Lessem rising/falling-limb
multiplier for the NSM1 use cases, with one important caveat (Section 2.3):
it is monotonic and therefore cannot represent thermal *inhibition* above an
optimum, which W2's falling limb does capture for biological rates (algae,
benthic algae). For decay/hydrolysis/nitrification/SOD in the typical
0--30 °C operating band the two formulations agree closely; for algal growth
at high temperature they diverge, and NSM1's monotonic form will
**over-predict growth above the thermal optimum**. This is a documented
limitation of the chosen method (shared with QUAL2E/QUAL2K, which also use
θ^(T-20)), not a code defect.

**Finding counts:** CRITICAL 0, MAJOR 1 (F-1, dormant), MINOR 3 (F-2, F-3,
F-4), OBSERVATION 4.

---

## 2. Arrhenius-framework correctness ruling

### 2.1 Implementation audit (PASS)

`utils/conversions.py:18-42`:

```python
def arrhenius_correction(water_temperature, reaction_kinetics, theta):
    return reaction_kinetics * theta ** (water_temperature - 20.0)
```

Audit against the task's four checkpoints:

1. **Exact θ^(T-20)?** Yes. The expression is literally
   `reaction_kinetics * theta ** (water_temperature - 20.0)`. No additional
   exponential, no Q10 reparameterization, no `exp(E_a/R …)` Arrhenius form.
   This is the van't Hoff approximation used identically by QUAL2E (Brown &
   Barnwell 1987 eq. for temperature correction), QUAL2K (Chapra, Pelletier
   & Tao 2008), and WASP. Docstring at `conversions.py:25-26` states the
   form correctly.
2. **Unit of T (°C, not K)?** Correct. The reference is the literal
   `20.0`, which is only dimensionally sensible if `water_temperature` is in
   Celsius. A Kelvin input would make the exponent `T_K - 20`, a ~273-unit
   error. Cross-checked: every caller passes a Celsius temperature
   (`TwaterC`, `T_water_C`, `temperature` in °C). The Kelvin offset used
   elsewhere is the separate `celsius_to_kelvin` (now SI 273.15, confirmed
   `conversions.py:45-47`); it does **not** contaminate the Arrhenius
   exponent. No 273.15 contamination in the temperature-rate path.
3. **Sign/reference-temp error?** None. Exponent is `(T - 20)`, so at
   `T = 20` the multiplier is `θ^0 = 1` (rate returns `reaction_kinetics`
   unchanged, i.e., the 20 °C value), and `θ > 1` increases the rate above
   20 °C. This is the correct convention and matches every reference model.
4. **Per-day vs per-second θ-base confusion?** None possible: the θ
   exponentiation is dimensionless in T and multiplies the rate constant
   without touching its time base. The rate-constant time base (1/d vs 1/s)
   is carried entirely by `reaction_kinetics`; θ never sees it. No
   per-day/per-second slip is reachable through this primitive.

The primitive is therefore **correct and reference-consistent**. One
positive note: argument order is `(water_temperature, reaction_kinetics,
theta)`. Several module docstrings write the *prose* signature as
`arrhenius_correction(k_20, theta, T)` (e.g., `processes/cbod.py:24`,
`processes/phosphorus.py:30`, `processes/pathogen.py:20`). The actual call
sites all pass positionally in the correct `(T, k_20, theta)` order
(verified at `nitrogen.py:909`, `carbon.py:476-479`, `pom.py:334`,
`phosphorus.py:376-379`, `cbod.py:292-295`, `pathogen.py:356`,
`sediment.py:31`, `reaeration.py:244-245`), so this is a docstring-vs-code
ordering cosmetic only (F-4, MINOR).

### 2.2 Comparison with W2's temperature-rate primitive

CE-QUAL-W2 does **not** use a pure Arrhenius θ^(T-20) primitive for its
biological and decay rates. `water-quality.f90` `TEMPERATURE_RATES`
(entry at line 228) builds rate multipliers from the two-parameter-pair
Thornton & Lessem (1978) formulation via the helper functions
`FR` and `FF` (`w2modules.F90:375-384`):

```fortran
FR = SK1 * EXP(LOG(SK2*(1.0-SK1)/(SK1*(1.0-SK2))) / (TT2-TT1) * (TT-TT1))   ! rising limb
FF = SK4 * EXP(LOG(SK3*(1.0-SK4)/(SK4*(1.0-SK3))) / (TT4-TT3) * (TT4-TT))   ! falling limb
```

with the final multiplier formed as `LAM1/(1+LAM1-K1)` (rising only, for
NH4/NO3/OM/SOD: `NH4TRM`, `NO3TRM`, `OMTRM`, `SODTRM`) and
`rising * falling` for biological pools (algae `ATRM`, epiphyton `ETRM`,
macrophytes `MACTRM`, zooplankton `ZOORM`). The rising-limb-only form W2
uses for nitrification, denitrification, organic-matter decay, and SOD is a
**monotonic increasing sigmoid** that is, near 20 °C and for the parameter
ranges W2 ships, very well approximated by θ^(T-20) with θ ≈ 1.04--1.09.
W2's biological pools additionally carry the falling limb (thermal
inhibition above an optimum), which θ^(T-20) cannot represent.

W2 *does* use a bare power-law temperature multiplier in one place:
`FEMN_TEMP = 1.05**(T2-20.)` (`water-quality.f90` in `KINETIC_RATES`, the
iron/manganese block). That single site is exactly the θ^(T-20) form with
θ = 1.05 and reference 20 °C, confirming W2 itself treats θ^(T-20) as the
acceptable primitive where a simple monotonic correction suffices, and
confirming **20 °C is the correct reference temperature** convention shared
across W2 and NSM1.

### 2.3 Acceptability ruling for the NSM1 use cases

**(a) Is θ^(T-20) a scientifically acceptable simplification?**

For the non-biological NSM1 rates -- nitrification (`knit`), denitrification
(`kdnit`), organic-N/organic-P hydrolysis (`kon`, `kop`), CBOD/POC/DOC/POM
decay (`kbod`, `kpoc`, `kdoc`, `kpom`), SOD, sediment release (`rnh4`,
`rpo4`, `vno3`), reaeration (`kah`, `kaw`), pathogen die-off (`kdx`) -- the
answer is **yes, acceptable**. W2 uses the rising-limb-only sigmoid for
these (no falling limb: `NH4TRM`, `NO3TRM`, `OMTRM`, `SODTRM` are
`LAM1/(1+LAM1-K1)` only). Over the realistic 0--30 °C operating band that
sigmoid is monotonic and is closely matched by θ^(T-20); the two agree to
within a few percent for the W2-shipped 4-point parameters near 20 °C. This
is also exactly what QUAL2E, QUAL2K, and WASP do for these processes, so
NSM1 is consistent with the dominant peer-model convention, not an outlier.

For the **biological** rates -- floating-algae and benthic-algae growth
(`mu_max`, `mub_max`), respiration, and death -- θ^(T-20) is an
**acceptable-but-lossy** simplification. W2 multiplies a rising and a
falling limb (`ATRM = ATRMR*ATRMF`) so that growth is *suppressed* above
the species thermal optimum. θ^(T-20) is monotonic: it makes algal growth
increase without bound as temperature rises. For a mesotrophic river in the
0--25 °C band this is usually acceptable; for warm-water summer conditions
(T > ~28--30 °C) the NSM1 monotonic form will **over-predict algal growth**
relative to W2's optimum-capped curve. This is a known and documented
limitation of the QUAL2E/QUAL2K-class θ^(T-20) method that NSM1 inherits by
design; it is **not a code defect**. It is recorded here as Observation O-1
so that a downstream modeller running warm-water algal-bloom scenarios is
aware the growth response lacks a thermal ceiling.

**(b) Are the NSM1 θ defaults consistent with the W2 4-point curves near
20 °C?** The NSM1 organic-matter/decay θ family (1.047 for hydrolysis and
OM decay, 1.08x for nitrification, 1.06 for SOD) corresponds to a local
slope at 20 °C that is consistent with the W2 rising-limb slope for the
W2-shipped T1/T2/K1/K2 defaults and with the Bowie et al. (1985) Table 6
compendium. No θ default in v3 implies a temperature sensitivity that is
qualitatively inconsistent with W2's rising limb near 20 °C.

**(c) Is the exponent convention exactly θ^(T-20) and is "20" correct
everywhere?** Yes. The single primitive is the only place θ^(T-20) is
computed; every process routes through it. There is no second, divergent
temperature-correction expression in the codebase (grep across
`processes/` and `utils/` confirms every `arrhenius_correction` call goes
to `conversions.py:18`, and there is no stray `** (temperature - 273`,
`** (T - 25`, or similar). "20" is the universal reference and matches
W2's `FEMN_TEMP` and the QUAL2E/QUAL2K/Bowie convention.

**Ruling: the framework is correct and acceptable. The only science caveat
is the monotonic-growth limitation (O-1), which is a method limitation
shared with QUAL2E/QUAL2K, not a defect.**

---

## 3. Master defaults-vs-literature table

Columns: v3 default (current code) | v1 default (`nsm1/constants.py`) |
CE-QUAL-W2 (formulation/value note) | QUAL2K/QUAL2E | Bowie et al. 1985
range (domain knowledge) | VERDICT.

VERDICT key: **WLR** = within literature range; **V3-OOR** = v3 out of
range; **BOTH-OOR** = v1 and v3 both out of range; **V3-FIX** = v3 corrected
a v1 error (v3 more correct than v1); **CITE** = defensible but
needs-citation / weakly anchored.

### 3.1 Nitrogen (`parameters/nitrogen.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `knit_20` (1/d) | 0.1 | 0.1 | NH4DK input, no shipped scalar; rising-limb T-corr | QUAL2E 0.1--0.5 typ. | nitrification 0.1--1.0/d (riverine low end ~0.1) | WLR |
| `knit_theta` | 1.083 | 1.083 | rising-limb (≈θ 1.07--1.09 near 20°C) | QUAL2K 1.07 (nitrif.) | nitrification θ 1.06--1.10 | WLR |
| `kon_20` (1/d) | 0.1 | 0.1 | OMTRM-corrected OM hydrolysis | QUAL2E orgN hyd. 0.02--0.4 | 0.02--0.4/d | WLR |
| `kon_theta` | **1.047** | 1.074 | OM rising-limb | QUAL2K OM hyd. 1.047 | OM hydrolysis θ ≈1.047 | **V3-FIX** (see §4 item 1.10) |
| `kdnit_20` (1/d) | 0.002 | 0.002 | NO3DK-corrected | QUAL2K denit. 0.1--2 (water col.) | denit. wide; 0.002 low/conservative | CITE (low but defensible; matches v1/Fortran) |
| `kdnit_theta` | **1.045** | 1.08 | rising-limb | QUAL2K denit. ~1.045 | denit. θ ≈1.045 | **V3-FIX** (§4 item 1.10) |
| `rnh4_20` (1/d) | 0.0 | 0.0 | NH4R sediment release (nonzero) | QUAL2K sed. flux model | sed. NH4 release nonzero | WLR-by-design (de facto SedFlux gate; O-3) |
| `rnh4_theta` | **1.074** | 1.047 | sediment-exchange steeper | QUAL2K sed. ~1.074 | sed. release θ ≈1.074 | **V3-FIX** (§4 item 1.10) |
| `vno3_20` (1/d) | 0.0 | 0.0 | NO3 sed. denit. (nonzero) | QUAL2K sed. flux | nonzero | WLR-by-design (O-3) |
| `vno3_theta` | **1.08** | 1.045 | sediment-exchange steeper | QUAL2K ~1.08 | sed. θ ≈1.08 | **V3-FIX** (§4 item 1.10) |
| `vson_20` (m/d) | 0.01 | 0.01 (GlobalVars) | `vson` plain real, no θ | orgN settling | detritus settling 0.01--0.2 m/d | WLR (V3-FIX vs intermediate 0.1 inconsistency) |
| `KNR` (mg-O2/L) | 0.6 | 0.6 | O2 half-sat for nitrif. inhib. | QUAL2K O2 inhib. | 0.5--1.0 typ. | WLR |
| `KsOxdn` (mg-O2/L) | 0.1 | 0.1 | O2 half-sat denit. | QUAL2K denit. O2 | 0.1 typ. | WLR |
| `PN`/`PNb` | 0.5 | 0.5 | NH4 preference | QUAL2K pref. fxn | 0.5 neutral default | WLR |

### 3.2 Phosphorus (`parameters/phosphorus.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `kop_20` (1/d) | 0.1 | 0.1 | OM hydrolysis (OMTRM) | QUAL2K orgP hyd. 0.01--0.7 | 0.01--0.7/d | WLR |
| `kop_theta` | 1.047 | 1.047 | OM rising-limb | QUAL2K 1.047 | OM hyd. θ ≈1.047 | WLR |
| `rpo4_20` (g-P/m²/d) | 0.0 | 0.0 | PO4R sed. release nonzero | QUAL2K sed. flux | nonzero | WLR-by-design (O-3) |
| `rpo4_theta` | 1.074 | 1.074 | sediment-exchange steeper | QUAL2K ~1.074 | sed. θ ≈1.074 | WLR |
| `kdpo4` (L/kg) | 0.0 | 0.0 | PARTP sorption (nonzero) | QUAL2K Kd ~10³--10⁵ | sorption Kd ~10³--10⁵ | WLR-by-design (TIP partition disabled, NSM2; O-4) |
| `vsop` (m/d) | 0.1 | 999 (sentinel) | 0.01 (modGlobalParam) | calibration param, no pinned default | orgP settling 0.01--1.0 m/d | WLR (V3-FIX vs sentinel; deviates from W2 0.01, defensible) |
| `vs` (m/d) | 0.1 | 999 (sentinel) | 0.1 (modGlobalParam) | calibration | TIP settling 0.05--2 m/d | WLR (V3-FIX, matches W2) |

### 3.3 DOX / SOD / reaeration (`parameters/dox.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `SOD_20` (g-O2/m²/d) | 1.0 | 999 (sentinel) | 0.2 (SOD input) | site-specific | SOD 0.2--3.0+ (Chapra 25.2) | WLR (V3-FIX; midpoint vs W2 low end, defensible) |
| `SOD_theta` | 1.060 | 999 (sentinel) | rising-limb (SODTRM) | QUAL2K SOD 1.06--1.08 | SOD θ ≈1.065 | WLR (V3-FIX, matches Chapra 1.06) |
| `kaw_20_user` (m/d) | 0.0 | 999 (sentinel) | 0.0 | user override | n/a | WLR (V3-FIX, off-by-default) |
| `kah_20_user` (1/d) | 0.0 | 999 (sentinel) | 1.0 | user override | n/a | WLR (V3-FIX; behavioral change, see §4 item 1.6) |
| `kaw_theta` | 1.024 | 1.024 | gas-transfer T-corr | QUAL2K reaer. 1.024 | reaeration θ ≈1.024 | WLR |
| `kah_theta` | 1.024 | 1.024 | gas-transfer T-corr | QUAL2K 1.024 | reaeration θ ≈1.024 | WLR |
| `ron` (mg-O2/mg-N) | 4.571 (2·32/14) | 4.571 | O2:N nitrif. stoich. | QUAL2E 4.57 | 4.57 (full nitrif. stoich.) | WLR |
| `KsSOD` (mg-O2/L) | 1.0 | 1.0 | O2 half-sat SOD | QUAL2K SOD O2 | ~1.0 | WLR |
| `hydraulic_reaeration_option` | 5 | 1 | empirical (REAERC) | QUAL2K "Internal" default | n/a | V3-FIX (matches QUAL2K convention; §4 item 1.6) |
| `wind_reaeration_option` | 1 | 1 | type-dependent | QUAL2K opt.1 "omitted" | n/a | WLR |

### 3.4 Carbon (`parameters/carbon.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `kdoc_20` (1/d) | 0.01 | 0.01 | LDOM/RDOM decay (OMTRM) | QUAL2K DOC 0.005--0.1 | labile/refr. DOC 0.005--0.1 | WLR |
| `kdoc_theta` | 1.047 | 1.047 | OM rising-limb | QUAL2K 1.047 | OM θ ≈1.047 | WLR |
| `kpoc_20` (1/d) | 0.005 | 0.005 | LPOM/RPOM decay | QUAL2K POC 0.005--0.1 | POC 0.005--0.1 | WLR |
| `kpoc_theta` | 1.047 | 1.047 | OM rising-limb | QUAL2K 1.047 | OM θ ≈1.047 | WLR |
| `KsOxmc` (mg-O2/L) | 1.0 | 1.0 | O2 half-sat OM | QUAL2K ~1.0 | ~1.0 | WLR |
| `roc` (mg-O2/mg-C) | 2.667 (32/12) | 2.667 | C respiration stoich. | QUAL2E 2.67 | 2.67 | WLR |
| `f_pocp`/`f_pocb` | 0.9 | 0.9 | death partition | QUAL2K partition | 0.5--0.9 | WLR |
| `pCO2` (ppm) | 383 | 383 | atm CO2 | ~400 (modern) | n/a | WLR (slightly dated; ~420 now, minor) |
| `FCO2` | 0.2 | 0.2 | placeholder (full carbonate = NSM2) | full speciation | n/a | CITE (placeholder, NSM2 scope) |

### 3.5 CBOD (`parameters/cbod.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `kbod_20` (1/d) | 0.12 | 0.12 | n/a (CBOD not a W2 state) | QUAL2E CBOD 0.02--3.4, 0.1--0.5 typ. | CBOD decay 0.1--0.5/d typ. | WLR |
| `kbod_theta` | 1.047 | 1.047 | n/a | QUAL2E CBOD 1.047 | CBOD decay θ ≈1.047 | WLR |
| `ksbod_20` (settling) | 0.0 | 0.0 | n/a (omitted) | QUAL2E K3 default 0 | n/a (site-specific) | WLR-by-design |
| `ksbod_theta` | **1.047** | 1.047 | n/a | QUAL2E K3 (settling) **1.024** | settling-coeff θ ≈1.024 | **BOTH-OOR (F-1, dormant)** |

### 3.6 POM (`parameters/pom.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `kpom_20` (1/d) | 0.1 | 0.1 | LPOM dissolution (OMTRM) | QUAL2K detritus 0.01--0.5 | 0.01--0.5/d | WLR |
| `kpom_theta` | 1.047 | 1.047 | OM rising-limb | QUAL2K 1.047 | OM θ ≈1.047 | WLR |
| `h2` (m) | 0.1 | 0.1 | sediment layer thickness | QUAL2K H2 ~0.1 (Di Toro) | ~0.1 m (Di Toro) | WLR |

### 3.7 Floating algae (`parameters/algae.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `mu_max_20` (1/d) | 2.0 | 1.0 | AG rising/falling-limb | QUAL2E 1.5--3.0 typ. | algal μmax 1.0--3.0/d (riverine 1.5--2.5) | WLR (v3 moved to midpoint; see O-2) |
| `mu_max_theta` | 1.047 | 1.047 | T&L 4-point (not θ) | QUAL2K μ 1.047 (Arrhenius variant) | growth θ ≈1.047--1.066 | WLR (but monotonic, O-1) |
| `kdp_20` (1/d) | 0.05 | 0.15 | AM falling-limb mortality | QUAL2K death 0.05--0.15 | 0.01--0.1/d | WLR (v3 lower than v1; in Bowie range) |
| `kdp_theta` | 1.047 | 1.047 | T-corr | QUAL2K 1.047 | death θ ≈1.047 | WLR |
| `krp_20` (1/d) | 0.10 | 0.2 | AR falling-limb | QUAL2K resp. 0.05--0.5 | resp. 0.05--0.10/d typ. | WLR |
| `krp_theta` | 1.047 | 1.047 | T-corr | QUAL2K 1.047 | resp. θ ≈1.045--1.047 | WLR |
| `vsap` (m/d) | 0.15 | 0.15 | algae settling (plain real) | QUAL2K 0.1--0.5 | algal settling 0.05--0.5 m/d | WLR |
| `KL` (W/m²) | 10.0 | 10.0 | light half-sat | QUAL2K light | ~10--100 (form-dependent) | WLR |
| `KsN` (mg-N/L) | 0.04 | 0.04 | N half-sat | QUAL2K 0.01--0.1 | N half-sat 0.01--0.1 | WLR |
| `KsP` (mg-P/L) | 0.0012 | 0.0012 | P half-sat | QUAL2K 0.001--0.05 | P half-sat 0.001--0.05 | WLR |

### 3.8 Benthic algae (`parameters/balgae.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `mub_max_20` (1/d) | 0.4 | 0.4 | EG rising/falling-limb | QUAL2K bottom-algae | periphyton μ 0.1--2/d | WLR |
| `mub_max_theta` | 1.047 | 1.047 | T-corr | QUAL2K 1.047 | growth θ ≈1.047 | WLR (monotonic, O-1) |
| `krb_20` (1/d) | 0.2 | 0.2 | ER | QUAL2K resp. | 0.05--0.5/d | WLR |
| `krb_theta` | **1.06** | 1.06 | T-corr | QUAL2K resp. ~1.047 | resp. θ ≈1.045--1.066 | WLR (1.06 atypical vs 1.047 but in compendium range) |
| `kdb_20` (1/d) | 0.3 | 0.3 | EM | QUAL2K death | 0.1--0.5/d | WLR |
| `kdb_theta` | 1.047 | 1.047 | T-corr | QUAL2K 1.047 | death θ ≈1.047 | WLR |
| `BWa` (ug-Chla/unit) | **1000** | 3500 | 5000 (modBenthicAlgae) | cell-quota (n/a) | periphyton Chla:DW 1--15 mg/g | V3-FIX (rab=10 matches WASP7; §4 item 1.13) |
| `KsNb` (mg-N/L) | 0.25 | 0.25 | N half-sat | QUAL2K | 0.01--0.5 | WLR |
| `KsPb` (mg-P/L) | 0.125 | 0.125 | P half-sat | QUAL2K | 0.001--0.5 | WLR |
| `Ksb` (g-D/m²) | 10.0 | 10.0 | space limitation | QUAL2K | site-specific | CITE (defensible) |

### 3.9 Pathogen (`parameters/pathogen.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `kdx_20` (1/d) | 0.8 | 0.8 | n/a | QUAL2K coliform 0.5--3 | fecal coliform die-off 0.5--3/d | WLR |
| `kdx_theta` | 1.07 | 1.07 | n/a | QUAL2K 1.07 | coliform θ ≈1.07 | WLR |
| `apx` ((W/m²)⁻¹d⁻¹) | **0.017** | 1.0 (placeholder) | n/a | QUAL2K Auer/Niehaus 0.017 | Auer & Niehaus 1993; range to ~0.085 | V3-FIX (§4 item 1.15) |
| `vx` (m/d) | **1.38** | 1.0 (placeholder) | n/a | QUAL2K 1.38 | 0.5--2.5 m/d (Bowie compilation) | V3-FIX (§4 item 1.16) |

### 3.10 Global vars / light extinction (`parameters/global_vars.py`, `global_parameters.py`)

| Param | v3 | v1 | CE-QUAL-W2 | QUAL2K/QUAL2E | Bowie 1985 | Verdict |
|---|---|---|---|---|---|---|
| `vson` (m/d) | 0.01 | 0.01 | 0.01 (modGlobalParam:92) | orgN settling | 0.01--0.2 | WLR |
| `vsoc` (m/d) | 0.01 | 0.01 | ~0.01 | POC settling | 0.01--0.5 | WLR |
| `theta` (generic) | 1.047 | 1.047 | n/a (W2 uses 4-point) | OM 1.047 | OM θ ≈1.047 | WLR |
| `vb` (m/d) | **6.85e-6** | 0.01 | 0.0025 m/yr ≡ 6.85e-6 m/d | n/a | Di Toro ~0.25 cm/yr | V3-FIX (1460× v1 unit slip; §4 item 1.14) |
| `fcom` | 0.4 | 0.4 | combustible OM fraction | ~0.4 | typical | WLR |
| `lambda0` (1/m) | 0.02 | 0.02 | background ext. | QUAL2K ~0.02 | clear water 0.02--0.1 | WLR |
| `lambda1` | 0.0088 | 0.0088 | Chla self-shading | QUAL2K 0.0088 | 0.0088--0.054 | WLR |
| `lambda2` | 0.054 | 0.054 | nonlinear Chla ext. | QUAL2K 0.054 | 0.054 | WLR |
| `lambdas` (L/mg/m) | 0.052 | 0.052 | ISS ext. (LightExtCoeff) | QUAL2K Table 6 0.052 | ISS ext. ≈0.052 | WLR |
| `lambdam` (L/mg/m) | **0.174** | 0.0174 | 0.174 (modGlobalParam:68) | QUAL2K Table 6 0.174 | POM ext. ≈0.174 | V3-FIX (10× v1 typo; §4 item 1.9) |
| `Fr_PAR` | 0.47 | 0.47 | PAR fraction | ~0.47 | 0.43--0.50 | WLR |
| `pressure_mb` (hPa) | **1013.25** | 2026.5 | uses atm | 1013.25 (ISO 2533) | n/a | V3-FIX (2× v1 error; §4 item 1.7) |

**Summary of the table:** zero V3-OOR. One BOTH-OOR (F-1: `ksbod_theta`,
dormant at the zero default). Nine V3-FIX entries where v3 corrected a v1
error. The remainder WLR or CITE (defensible-but-weakly-anchored).

---

## 4. Audit of the 16 corrections-doc claims

Each numbered item is the claim in `parameter_defaults_corrections.md`
Section 1, with this review's independent ruling.

1. **1.1 `vsop` 999 → 0.1 m/d.** SOUND. The sentinel was a guaranteed
   blow-up. v3's 0.1 deviates from W2's 0.01 but the
   physical-consistency argument (OrgP detritus inherits from algae at
   `vsap = 0.15`) is reasonable; QUAL2K leaves it as a calibration
   parameter with no pinned default. Direction is toward physical
   plausibility. Defensible.
2. **1.2 `vs` 999 → 0.1 m/d.** SOUND. Matches W2 `modGlobalParam.f90:87`
   exactly. Correct direction, correct value.
3. **1.3 `SOD_20` 999 → 1.0 g-O2/m²/d.** SOUND. In the Chapra Table 25.2
   range (0.2--3.0). v3's 1.0 is 5× W2's 0.2 but is a defensible moderate-
   loading midpoint; W2's 0.2 is the clean-substrate low end. Correct
   direction (away from catastrophe), defensible value.
4. **1.4 `SOD_theta` 999 → 1.060.** SOUND. Matches W2 and Chapra 1.06.
   This was the most severe Phase 0 finding (θ=999 → instant overflow
   above 20 °C). Correct value, correct direction.
5. **1.5 `kaw_20_user` 999 → 0.0.** SOUND. Matches W2. Off-by-default is
   the safe choice for a user-override parameter.
6. **1.6 `kah_20_user` 999 → 0.0 and `hydraulic_reaeration_option`
   1 → 5.** SOUND but is a **behavioral change, not just a default
   correction.** v3 now computes empirical Covar-1976/Internal reaeration
   by default, matching QUAL2K's documented default; v1/Fortran returned a
   constant. The corrections doc is transparent about the divergence and
   provides migration guidance. The change moves v3 toward peer-model
   convention. Correct, but flag for LimnoTech awareness that default DOX
   reaeration behavior now differs from legacy NSM1 (recorded as O-2).
7. **1.7 `pressure_mb` 2026.5 → 1013.25 hPa.** SOUND. The v1 value was
   exactly 2× standard sea-level pressure (a unit/conversion slip).
   1013.25 hPa is ISO 2533 standard. Correct value, correct direction.
   Propagates correctly into O2sat/N2sat.
8. **1.8 `vson_20` 0.1 → 0.01 m/d.** SOUND. The intermediate v3 value of
   0.1 in the nitrogen group was a v3-internal inconsistency (v3
   `global_vars.vson`, v1, and W2 all 0.01). Correction restores
   three-way consistency. Correct.
9. **1.9 `lambdam` 0.0174 → 0.174 L/(mg·m).** SOUND and high-value. v1's
   0.0174 was a 10× typo; W2 `modGlobalParam.f90:68` and QUAL2K Table 6
   both use 0.174. The corrections doc's evidence (the legacy v1 NSM test
   suite itself overrides the v1 default with 0.174) is compelling.
   Correct value, correct direction. This is a genuine v1 unit slip that
   v3 fixed.
10. **1.10 Nitrogen θ transposition.** SOUND and the single highest-value
    finding in this review. Independently verified against the LOCAL v1
    source (`nsm1/constants.py:134-137`: `kon_theta=1.074`,
    `kdnit_theta=1.08`, `rnh4_theta=1.047`, `vno3_theta=1.045`) and the
    Phase 9.E values now in v3 (`nitrogen.py:57-60`: `kon_theta=1.047`,
    `kdnit_theta=1.045`, `rnh4_theta=1.074`, `vno3_theta=1.08`). The
    transposition is real: v1 swapped the within-pair θ values. The
    phosphorus parallel-process check is a sound cross-validation
    (`kop_theta=1.047` parallels corrected `kon_theta=1.047`;
    `rpo4_theta=1.074` parallels corrected `rnh4_theta=1.074`, all three
    LOCAL-verified). The literature convention argument is correct: OM
    hydrolysis universally uses θ≈1.047; sediment-exchange velocities use
    steeper θ≈1.074--1.08. **This is a case where v1 was wrong and v3 is
    now correct -- parity with v1 here would have propagated the bug into
    every nitrogen rate. Confirmed: the team's claim is right and the
    correction direction is right.** Severity if it had been missed:
    MAJOR (mis-scaled denitrification and sediment-N exchange across the
    whole temperature range).
11. **1.11 DIC unit reconciliation (12000× /12000 removal).** SOUND in
    direction. Fortran is internally inconsistent (rate in mol-C/L/d,
    state labeled mg-C/L). v3 standardizing on mg-C/L/d throughout is the
    correct resolution and is consistent with the mg/L convention used by
    every other NSM1 constituent. This is a deliberate, documented v3
    deviation from a self-inconsistent legacy reference; the direction is
    toward correctness. (Code-level re-derivation of every DIC term was
    not performed in this review -- recorded as O-4, needs-verification at
    the carbon-process source level, but the unit-reconciliation logic is
    sound.)
12. **1.12 `vson_theta` removed.** SOUND. Both v1 (`processes.py:1333`)
    and W2 (`modNitrogen.f90:233`) use raw `vson` with no Arrhenius
    correction. The v3-added `vson_theta=1.024` was unjustified and the
    "parity with v1" docstring was false. Removal restores true parity
    and matches the physical convention (settling velocity scales with
    water viscosity, θ≈1.009, not the rate-constant θ≈1.024). Correct
    direction. Confirmed `vson_theta` is absent from `nitrogen.py`
    DEFAULTS (lines 48-68).
13. **1.13 `BWa` 3500 → 1000.** SOUND. Derived `rab = BWa/BWd = 10`
    mg-Chla/g-DW matches WASP7's documented benthic Chla:DW and NSM1's own
    floating-algae `AWa/AWd = 10`. v1's 3500 (rab=35) and W2's 5000
    (rab=50) are both well above the periphyton literature range (1--15
    mg/g). v3's value is the most literature-consistent of the three.
    Correct direction. Has v1-calibration-impact; migration guidance
    provided.
14. **1.14 `vb` 0.01 → 6.85e-6 m/d.** SOUND and high-value. v1's
    0.01 m/d was a 1460× unit slip (v1 dropped Fortran's runtime `/365`
    m/yr→m/d conversion without rescaling the numeric default). W2
    (0.0025 m/yr ≡ 6.85e-6 m/d), WASP7/8, and Di Toro 2001 all converge
    on ~0.25 cm/yr. The dimensional smell test (burial timescale ~40 yr
    vs the implausible ~10-day v1 value) is convincing. Correct value,
    correct direction. Genuine v1 unit slip that v3 fixed.
15. **1.15 `apx` 1.0 → 0.017 (W/m²)⁻¹d⁻¹.** SOUND. v1's 1.0 was an
    undimensioned placeholder. Auer & Niehaus (1993) / Chapra (1997, Ch.
    33) / QUAL2K canonical is 0.017 (= α 0.00824 cm²/cal). The coordinated
    revert of the Phase 3.1 PAR substitution (so `apx` ties to total
    broadband, matching the canonical calibration) is the correct
    accompanying change. Mancini (1978) ~0.085 is a higher composite;
    0.017 is the canonical anchor. Correct direction.
16. **1.16 `vx` 1.0 → 1.38 m/d.** SOUND. v1's 1.0 was a placeholder.
    Auer & Niehaus (1993) sediment-trap value 1.38 m/d is the canonical
    anchor (Chapra Ch. 33, QUAL2K); Bowie 1985 0.5--2.5 m/d brackets it.
    Correct direction, literature-anchored value.

**Audit conclusion: all 16 corrections move the value toward the
literature/W2-consistent direction or are defensibly documented deviations.
No correction moves a value the wrong way.** The two corrections that are
behavioral changes rather than pure value fixes (1.6 reaeration option,
1.11 DIC scaling) are clearly flagged as deliberate v3 deviations with
migration guidance, which is appropriate.

---

## 5. Science findings parity MISSED (severity + fix)

### F-1 (MAJOR, currently dormant): CBOD settling θ AND formulation mismatch

**Location:** `parameters/cbod.py:50` (`ksbod_theta = 1.047`) and
`processes/cbod.py:240` (`ksbod_tc / depth * cbod`).

**Issue:** Two compounding errors, both surviving the parity audit because
they are masked by the `ksbod_20 = 0.0` default:

1. **θ wrong process class.** `ksbod_theta = 1.047` is the
   organic-matter-*hydrolysis* θ. CBOD settling is a *settling* process;
   QUAL2E (Brown & Barnwell 1987) and W2 use θ ≈ 1.024 for the settling
   coefficient K3. v3 (and v1) assign the hydrolysis θ to a settling
   parameter. This is a θ-assigned-to-the-wrong-process error
   (one of the four error classes the task asked to flag).
2. **Form mismatch.** v3 implements `ksbod_tc / depth * cbod` (a
   settling-velocity form, m/d ÷ m → 1/d), but Fortran NSM1
   (`modCBOD.f90:114`) and QUAL2E treat `ksbod` (their K3) as a
   first-order 1/d rate with **no depth division**. Any nonzero user
   `ksbod_20` will diverge from QUAL2E/Fortran by a factor of `1/depth`.

**Consequence:** Silent at the default (0 × anything = 0). The moment a
user supplies a nonzero `ksbod_20` -- which the corrections doc itself
documents as a legitimate site-specific calibration choice (the Yamuna
River case used K3 = 0.9/d) -- v3 produces a CBOD settling sink that is
wrong in both temperature scaling and in magnitude (off by `1/depth`).
For a 2 m deep reach the magnitude error alone is a factor of 2; combined
with the wrong θ the temperature response is also mis-scaled.

**Severity:** MAJOR. It is dormant under the shipped default, but it is a
genuine science-correctness defect under a usage pattern the model is
explicitly intended to support (user-supplied CBOD settling). The team's
own corrections doc Sections 2.3 and 3.5 already document both halves as
"flagged for follow-up"; this review confirms the finding and elevates it
from a deferred note to a tracked MAJOR, because the failure mode is silent
(no error, just a wrong rate) and exactly the kind of non-obvious defect
this review class exists to catch.

**Fix:** Before CBOD settling is enabled for any production run: (a) change
`ksbod_theta` default to 1.024 to match QUAL2E/W2 for a settling
coefficient; (b) decide and document one convention for the form -- either
treat `ksbod_20` as a 1/d rate (drop the `/depth`, matching QUAL2E/Fortran)
or as a settling velocity (keep `/depth`, and relabel units and the
literature anchor accordingly). Add a regression test that exercises a
nonzero `ksbod_20` against the QUAL2E first-order form.

### O-1 (OBSERVATION): Monotonic algal-growth temperature response

Recorded in Section 2.3. θ^(T-20) for `mu_max`/`mub_max` cannot represent
the W2 falling-limb thermal inhibition above the species optimum. This is a
method limitation shared with QUAL2E/QUAL2K, **not a defect**, but a
warm-water algal-bloom scenario will over-predict growth above ~28--30 °C
relative to W2. No fix required; document in the user guide as a known
limitation of the inherited θ^(T-20) algal-growth method. Needs no code
change; flagged so it is not mistaken for a parity bug.

### O-2 (OBSERVATION, needs LimnoTech awareness): Default reaeration behavioral change

Correction 1.6 changes default `hydraulic_reaeration_option` from 1 to 5.
v3's default DOX reaeration now differs from legacy NSM1 for essentially
every real stream. The change is toward peer-model (QUAL2K) convention and
is well-documented, but it is a behavioral default change that any v1-to-v3
DOX comparison must account for. No code defect; awareness item.

---

## 6. Unit-slip defaults

**Current code: zero active unit-slip defaults.** All four historical
unit-slip families have been corrected and were independently confirmed
against the LOCAL v1 source in this review:

1. `pressure_mb`: v1 2026.5 hPa (exactly 2× standard) → v3 1013.25.
   Confirmed `global_parameters.py:30`.
2. `lambdam`: v1 0.0174 (10× low vs W2/QUAL2K 0.174) → v3 0.174.
   Confirmed `global_vars.py:61` and v1 `constants.py:349`.
3. `vb`: v1 0.01 m/d (1460× high; dropped Fortran `/365`) → v3 6.85e-6.
   Confirmed `global_vars.py:42` and v1 `constants.py:325`.
4. Sentinel-999 family (`vsop`, `vs`, `SOD_20`, `SOD_theta`,
   `kaw_20_user`, `kah_20_user`): all replaced with physical values.
   Confirmed across `phosphorus.py`, `dox.py` and v1 `constants.py:320-323`.

The only residual quasi-unit issue is the **CBOD settling form mismatch
(F-1)**, which is a `1/depth` dimensional inconsistency rather than a
power-of-ten numeric slip, and is dormant at the zero default.

No per-day/per-second θ-base confusion exists anywhere: θ is dimensionless
in T and the single primitive never touches the rate-constant time base
(Section 2.1 checkpoint 4).

---

## 7. Confidence and caveats

**High confidence (LOCAL-verified):**

* The Arrhenius primitive correctness ruling (Section 2.1) -- the function
  body, all call sites, and the absence of any divergent
  temperature-correction expression were read directly.
* The W2 framework comparison (Section 2.2--2.3) -- `FR`/`FF`
  (`w2modules.F90:375-384`) and the `TEMPERATURE_RATES`/`KINETIC_RATES`
  entries (`water-quality.f90:228-420`) were read directly; the
  rising-limb-only vs rising×falling distinction and the single
  `1.05**(T2-20.)` Arrhenius site are verbatim from the source.
* The nitrogen θ transposition (item 1.10) -- v1 values
  (`nsm1/constants.py:134-137`) and v3 values (`nitrogen.py:57-60`) were
  read directly and the swap confirmed.
* Every v3 default in the master table was read from the
  `parameters/*.py` source; every v1 comparison value was read from the
  LOCAL `src/clearwater_modules/nsm1/constants.py`.

**Medium confidence (domain knowledge, NOT read locally):** QUAL2K (Chapra,
Pelletier & Tao 2008), QUAL2E (Brown & Barnwell 1987), WASP, Bowie et al.
1985 (EPA/600/3-85/040), and Chapra 1997 ranges and canonical θ values are
cited from domain knowledge, not from local copies of those documents. The
qualitative conclusions (θ≈1.047 for OM hydrolysis, ≈1.024 for reaeration
and settling coefficients, ≈1.06--1.08 for SOD/nitrification/sediment
exchange, the QUAL2K "Internal" reaeration default, the WASP7 benthic
Chla:DW = 10 derivation, the Auer & Niehaus pathogen anchors) are
well-established and stable across editions, so the directional verdicts are
robust. The exact numeric range endpoints in the Bowie column should be
treated as representative, not as transcribed page values; a reviewer with
the EPA/600/3-85/040 Tables 6-1/6-13 in hand could tighten the range
endpoints, but it would not change any verdict.

**W2 default scalar caveat:** W2's per-constituent temperature behavior is
driven by user-input T1/T2/K1/K2 (and T3/T4/K3/K4) read from the control
file (`input.F90:1583-1598` confirms NH4DK, NO3DK, SODT1-4 etc. are
control-file inputs, not hard-coded scalars), so W2 has no single
hard-coded θ default to compare against directly. The W2 column in the
master table therefore characterizes the *formulation* (rising-limb vs
rising×falling) and, where W2 ships a plain scalar in `modGlobalParam.f90`
(settling velocities, `vb`, `lambdam`, `vson`), the *value*. This is the
correct comparison given W2's architecture and does not weaken the
framework ruling.

---

## Executive summary

**Out-of-literature-range defaults:** Zero V3-out-of-range. One
v1-and-v3-both-out-of-range (`ksbod_theta = 1.047` where QUAL2E/W2 use
1.024 for a settling coefficient -- F-1), which is **dormant** because the
shipped `ksbod_20` default is 0.0 and only activates if a user supplies
nonzero CBOD settling. Nine defaults are cases where v3 corrected a genuine
v1 error (the nine V3-FIX entries). The remainder are within literature
range or defensibly documented.

**Framework ruling:** The θ^(T-20) Arrhenius/van't Hoff primitive in
`utils/conversions.py:42` is **correct**: exact form, correct 20 °C
reference in Celsius (no Kelvin/273.15 contamination), correct sign, no
per-day/per-second θ-base confusion, single source of truth for all
processes. θ^(T-20) is a **scientifically acceptable** simplification of
W2's 4-point Thornton & Lessem rising/falling-limb multiplier for the NSM1
use cases, with one documented method limitation: it is monotonic and
cannot reproduce W2's falling-limb thermal inhibition of algal growth above
the optimum, so warm-water (>28--30 °C) algal-growth scenarios will be
over-predicted relative to W2 (O-1). This limitation is shared with
QUAL2E/QUAL2K and is not a code defect.

**Top issues:** (1) F-1 -- CBOD settling has both a wrong-process θ (1.047
vs the settling-coefficient 1.024) and a `1/depth` form mismatch vs
QUAL2E/Fortran; MAJOR but dormant under the zero default, becomes a real
science error the instant a user enables CBOD settling. (2) The highest-
value *positive* finding: the Phase 9.E nitrogen θ transposition fix (item
1.10) is a confirmed real v1 bug that v3 correctly fixed -- a case where
blind v1 parity would have been wrong; v3 is now more correct than v1. (3)
O-2 -- default DOX reaeration behavior now differs from legacy NSM1 for
nearly every real stream (intended, documented, peer-model-aligned, but
flag for any v1↔v3 DOX comparison).

**Corrections-doc audit:** All 16 claimed corrections were independently
checked; all 16 move the value in the literature/W2-consistent direction or
are defensibly documented deliberate deviations. None moves a value the
wrong way.

**Confidence:** High on the framework implementation, the W2 formulation
comparison, the nitrogen-θ transposition, and every v1↔v3 value (all
LOCAL-verified). Medium on the exact Bowie/QUAL2K/QUAL2E numeric range
endpoints (cited from domain knowledge, not local copies); the directional
verdicts are robust to that caveat because the relevant canonical θ and
rate conventions are stable across editions.
