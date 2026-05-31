# NSM2 Alkalinity/pH — authoritative Fortran extraction + decision log (S4-3)

Authoritative transcription of NSM2 Fortran `modAlkalinity.f90`
(+ `modGlobalParam.f90` `pH_solution_option`) for v3 Step-4 phase
**S4-3 (Alkalinity/pH carbonate solver + NSM1-SCI-A1)**. Constants
verbatim from `/Users/todd/GitHub/ecohydrology/HEC-RAS-WQ/NSMIISource
files/modAlkalinity.f90`. CE-QUAL-W2 `water-quality.f90:3150-3171` is
the SCI-A1 / solver cross-check. Same discipline as the S4-1/S4-2
extraction docs.

## 1. Alkalinity source/sink kinetics (modAlkalinity `ComputeAlkalinityKinetics`)

`50000` converts eq-H⁺/L → mg-CaCO₃/L. Stoich (`InitializeAlkalinity`):

| Const | Value | Units | v3 1.0 name | Note |
|---|---|---|---|---|
| `ralkca` | 14/106/12/1000 | eq/mg-C | `r_alkaa`/`r_alkba` | algal photosynthesis (NH4 uptake) |
| `ralkcn` | 18/106/12/1000 | eq/mg-C | `r_alkan`/`r_alkbn` | algal NO3 uptake |
| `ralkn` | 2/14/1000 | eq/mg-N | `r_alkn` | nitrification |
| `ralkden` | **4**/14/1000 | eq/mg-N | `r_alkden`=**1**/14/1000 | **v3 already corrected (NSM1-SCI-N1: Fortran 4× too high; W2 + Stumm&Morgan = 1 eq/mol-N). Documented divergence shipped in v3 1.0.** |

```
Alk_AlgalUptake      = Σ_i (ralkca·PNp_i − ralkcn·(1−PNp_i)) · rca_i · AlgalGrowth_i · 50000
Alk_AlgalRespiration = Σ_i  ralkca · rca_i · AlgalRespiration_i · 50000
Alk_Nitrification    = ralkn   · NH4_NO3_Nitrification · 50000
Alk_Denitrification  = ralkden · NO3_Denitrification   · 50000
Alk_BenthicUptake      = Fb·(ralkca·PNb − ralkcn·(1−PNb))·rcb·BenthicGrowth/depth·50000   (use_BAlgae)
Alk_BenthicRespiration = Fb· ralkca · rcb · BenthicRespiration / depth · 50000
dAlkdt = −Alk_AlgalUptake + Alk_AlgalRespiration − Alk_Nitrification
         + Alk_Denitrification − Alk_BenthicUptake + Alk_BenthicRespiration
```

**The NSM2 Fortran algal/benthic term is carbon-routed**
(`rca·AlgalGrowth` = C flux × stoich × PN split) — *identical basis to
v3 1.0*. So **SCI-A1 (N-flux reformulation) is a deliberate divergence
from BOTH the NSM2 Fortran baseline AND v3 1.0** — the D2 / SCI-N1
pattern. Tier-5 Fortran-parity gets a documented algal/benthic-Alk
carve-out (the term that intentionally differs).

## 2. CE-QUAL-W2 N-flux basis (SCI-A1 reference; `water-quality.f90:3150-3171`)

Stumm & Morgan (1996) Table 4.5: NH4 uptake −14 eq / 16 mol-N;
NO3 uptake +18 eq / 16 mol-N; NH4 production (resp) +14/16;
nitrification −2 eq/mol; denitrification +1 eq/mol. Alk as
mg-CaCO₃/L, factor `50.044/14.00674`:

```
ALKSS = (50.044/14.00674) · [ 14/16·(NH4 produced by algae/resp − NH4 uptake)
                              + 18/16·(NO3 uptake) − 2·NH4_nitrif + NO3_denit ]
```

SCI-A1 = recompute the algal/benthic-Alk term from the **algal N flux**
(rna·AlgalGrowth split by NH4/NO3 preference), W2 stoichiometry,
replacing the carbon-routed `(ralkca·PN − ralkcn·(1−PN))·rca·AlgalGrowth`.

## 3. pH carbonate solver (modAlkalinity `ComputeAlkalinityDerivedVariables`)

`TwaterK = TwaterC + 273.15`. **Freshwater, T-only — NO ionic strength:**

```
Kw = 10^(−4787.3/TK − 7.1321·log10(TK) − 0.010365·TK + 22.80)
K1 = 10^(−356.3094 − 0.06091964·TK + 21834.37/TK + 126.8339·log10(TK) − 1684915/TK²)
K2 = 10^(−107.8871 − 0.03252849·TK +  5151.79/TK +  38.92561·log10(TK) −  563713.9/TK²)
```

Charge-balance residual (Chapra Eq. 3.58), `hh = 10^(−pH)`:

```
f(pH) = (K1·hh + 2·K1·K2)/(hh² + K1·hh + K1·K2)·DIC + Kw/hh − hh − Alk/50000
```

`pH_solution_option` (modGlobalParam, per region): **1 = Newton-Raphson**
(pH₀=7, analytic derivative), **2 = Bisection** (bracket [3,13]). Both:
`imax=13` iterations, `es=0.001` relative error. **On non-convergence
or pH∉[0,14] / bad bracket the Fortran executes `stop`** (hard process
abort) — unacceptable in v3.

Speciation (for [CO2*]/[HCO3⁻]/[CO3²⁻] and the Carbon DIC coupling):
`α0 = hh²/(hh²+K1·hh+K1·K2)`, `α1 = K1·hh/(…)`, `α2 = K1·K2/(…)`;
`[CO2*] = α0·DIC`. Carbon currently uses the `FCO2=0.2` placeholder in
`co2_reaeration = 0.923·ka_tc·(KH·pCO2/1e6·12000 − FCO2·DIC)`; the
solver supersedes `FCO2` with the real `α0`.

## 4. v3-side facts

- Existing `Alkalinity` process: state `alkalinity` (mg-CaCO₃/L);
  carbon-routed algal/benthic term (`processes/alkalinity.py:365-494`),
  `r_alkden=1/14/1000` (SCI-N1 shipped). Couples FloatingAlgae /
  BenthicAlgae / Nitrogen.
- Baseline parity = **all-variables `numpy.array_equal`** vs
  `baseline_coupled_trajectory_6c10f36.nc` (zero tolerance; no
  per-variable exclusion). A deliberate `alkalinity` change ⇒ the
  parity test fails unless the baseline is regenerated (§1.1 / §7
  re-baseline policy; NSM1-CA-1 / SCI-N1 precedent — each changed
  only `alkalinity`).
- **No Step-1 reserved CO2/pH growth hooks** (unlike silica's
  AWsi/KsSi/si_limitation_option). f_CO2 needs NEW per-group
  hooks. f_Si template: `floating_algae._compute_f_si` +
  the gated `rate_growth` Si block (S4-2).
- v3 1.0 pH stopgap = NSM1 design spec §14 "simple tracer; post-hoc
  pH worked example" — superseded by this live solver (plan S4-6).

## 5. Decision log — LOCKED 2026-05-17 (Todd)

- **D-A-1 — LOCKED: land SCI-A1 unconditionally + alkalinity-only
  re-baseline.** Reformulate the algal/benthic-Alk term to the W2 /
  Stumm&Morgan **N-flux basis** (14/16 NH4, 18/16 NO3, factor
  50.044/14.00674) from the algal N uptake flux (rna·AlgalGrowth split
  by the NH4/NO3 preference), replacing the carbon-routed
  `(ralkca·PN − ralkcn·(1−PN))·rca·AlgalGrowth`. Unconditional (not
  flagged). Gate: prove EVERY other state variable is still
  bit-identical vs `6c10f36` and **only `alkalinity` differs**, then
  regenerate the gold-standard baseline as a NEW
  `baseline_coupled_trajectory_<newhash>.nc` + a committed doc note of
  the alkalinity-only delta (NSM1-CA-1 / SCI-N1 precedent).
  **The new baseline artifact requires Todd's explicit sign-off at
  the S4-3 gate (a binding §1.1 action) — do NOT regenerate the
  baseline without it.** Tier-5 Fortran-parity gets a documented
  algal/benthic-Alk carve-out (this term intentionally ≠ NSM2 Fortran
  AND ≠ v3 1.0).
- **D-A-2 — LOCKED: global flag, default OFF.** New
  `global_parameters.py` flag (default False; name TBD, e.g.
  `use_carbonate` / `use_pH_solver`) gates the solver + f_CO2 + the
  Carbon `FCO2`→α0 supersession. OFF ⇒ no solver, Carbon keeps the
  `FCO2` placeholder ⇒ all non-`alkalinity` variables bit-identical
  (the §7 contract). SCI-A1 (D-A-1) is the separate
  unconditional-and-re-baselined piece — i.e. with the solver OFF the
  *only* trajectory change vs `6c10f36` is `alkalinity` (from SCI-A1),
  and that is exactly what the re-baseline captures.
- **D-A-3 — LOCKED: activate an ionic-strength correction now**
  (Todd chose the reservoir-grade path over replicate-NSM2-T-only).
  Implement I-dependent apparent constants K1(T,I), K2(T,I), Kw(T,I)
  (Millero / Davies or Plummer&Busenberg) with `I` an explicit solver
  input (from salinity / specific conductance; default 0).
  **MANDATORY carve-out design: the I-correction MUST vanish at I=0
  so the I=0 path collapses byte-exactly to the verbatim NSM2 Fortran
  Kw/K1/K2** — Tier-5 Fortran-parity is then asserted at I=0 (v3 ==
  NSM2 exactly), and the I>0 path is the documented reservoir
  extension cross-checked vs W2 / Millero (the D-Si-3 inert-default
  carve-out pattern). **`I` is sourced from the S4-2.5 `salinity`
  forcing via `utils.salinity.ionic_strength_from_salinity`**
  (`I(0)=0.0` exactly ⇒ the freshwater/no-salinity path is byte-exact
  to NSM2); the solver also accepts an explicit `ionic_strength`
  override.
- **D-A-4 — LOCKED: graceful failure, never abort.** Newton-Raphson
  primary → Bisection fallback (wide bracket) → hold previous-substep
  pH (or clamp) + emit a clip-style diagnostics event. Never `stop`.
  Documented robustness divergence (v3 clip_negative_state /
  sanitize_rate convention). Tier-5 parity asserts equality on the
  convergent path (where v3 == NSM2).
- **D-A-5 — LOCKED (recommended, proceeding):** new per-group
  `co2_limitation_option` (default 1 = unlimited → bit-identical) +
  `KsCO2` in algae/balgae GROUP_DEFAULTS; fold `f_CO2 =
  [CO2*]/(KsCO2+[CO2*])` into `rate_growth` symmetric with N/P
  (mirror S4-2 f_Si exactly); DSi-style lagged one sub-step; active
  only when the solver is ON and a group's option=2.
