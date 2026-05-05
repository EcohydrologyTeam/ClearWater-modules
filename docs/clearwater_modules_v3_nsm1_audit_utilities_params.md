# v3 NSM1 Utilities + Parameter Library — Three-way audit

**Date:** 2026-05-05
**References:** v3 (streaming branch), v1 (`src/clearwater_modules`), Fortran (`/Users/todd/Downloads/NSM_comparison/NSM1/Source Files`).

## Summary

1. Utilities: 6 functions audited. **5 of 6 match Fortran exactly** (`kah_20`, `kaw_20`, `ka_tc`, `L`, `fdp`). `SOD_tc` is a deliberate v3 architectural refactor (pure Arrhenius). `PAR` is a deliberate v3 refactor (toggle moved to consumer).
2. Parameters: ~145 distinct parameter entries audited across 13 v3 groups. **All 7 critical sentinel-999 corrections are confirmed genuinely needed** — Fortran has the correct physical defaults (e.g., 0.1 m/d for `vsop`, 1013.25 hPa is standard atm), so v1 alone introduced the 999/2026.5 sentinels; v3 restored Fortran-aligned values.
3. **Likely v1 flaw not corrected in v3**: `lambdam=0.0174` in v1/v3 vs Fortran's `0.174` — 10x discrepancy in POM contribution to Beer-Lambert light extinction. Needs reconciliation.
4. v3 deliberate improvements verified: `pressure_mb=1013.25`, sentinel rescue, `PAR` toggle inversion, `SOD_tc` pure-Arrhenius split. All consistent with `parameter_defaults_corrections.md` Sections 1-3.
5. Internal v3 inconsistency: `vson_20=0.1` in `parameters/nitrogen.py` vs `vson=0.01` in `parameters/global_vars.py` (Fortran uses 0.01). Either consolidate or document.

---

## Part 1 — Utility modules

### `reaeration.py`

#### `kah_20` (hydraulic, 9 options)

Fortran source: `modGlobalParam.f90:268-339`, subroutine `O2Reaeration`. v1 source: `shared/processes.py:65-98`. v3 source: `utils/reaeration.py:26-119`.

| Option | Fortran formula | v1 formula | v3 formula | Status |
|---|---|---|---|---|
| 1 | user `kah%rc20` | `kah_20_user` | `kah_20_user` | match |
| 2 | `(3.93 v^0.5)/h^1.5` (line 275) | identical | identical | match |
| 3 | `(5.32 v^0.67)/h^1.85` (line 281) | identical | identical | match |
| 4 | `5.026 v / h^1.67` (line 287) | identical | identical | match |
| 5 | depth-piecewise: `<0.61` Owens, `>0.61` O'Connor, `=0.61` Churchill (lines 289-306) | identical with `==0.61` Churchill branch | identical | match |
| 6 | flow-piecewise (lines 308-315), `<0.556` → `517 (vS)^0.524 Q^-0.242`, else `596 (vS)^0.528 Q^-0.136` | identical | identical | match |
| 7 | flow-piecewise (lines 317-324), `<0.556` → `88 (vS)^0.313 h^-0.353`, else `142 (vS)^0.333 h^-0.66 W^-0.243` | identical | identical | match |
| 8 | flow-piecewise (lines 326-333), `<0.425` → `31183 vS`, else `15308 vS` | identical | identical | match |
| 9 | Froude form: `2.16 (1 + 9 Fr^0.25) u*/h` with `Fr = v/sqrt(g h)` (lines 334-338) | identical | identical | match |

**Finding:** match across all three references. The docstring author attributions in v3 (Covar/Owens-Gibbs/Churchill/Tsivoglou-Wallace/Padden-Gloyna/USGS pool-and-riffle/Thackston-Krenkel/Langbien-Durum) disagree with Fortran's inline comments (Owens, O'Connor, Churchill, Cover, Melching-Flores, Tsivoglou-Neal, Thackson-Dawson). The author attributions differ but the formulas are identical. **Severity:** minor documentation inconsistency.

**v3 implementation detail (Phase 5.5):** `np.select` dim-stripping fix at `utils/reaeration.py:107-119` — match-preserving.

#### `kaw_20` (wind, 13 options)

Fortran source: `modGlobalParam.f90:341-414`. v1 source: `shared/processes.py:117-146`. v3 source: `utils/reaeration.py:122-204`.

Fortran labels its 13 options as: 1=user, 2=Broecker, 3=Gelda, 4=Banks-Herrera, 5=Wanninkhof, 6=Cole-Buchak, 7=Banks, 8=Smith, 9=Liss, 10=Downing-Truesdale, 11=Kanwisher, 12=Yu, 13=Weiler.

v1 and v3 match byte-for-byte numerically:

| Option | Fortran formula | v1 / v3 formula | Status |
|---|---|---|---|
| 1 | user `kaw%rc20` | `kaw_20_user` | match |
| 2 | `0.864 Uw10` | identical | match |
| 3 | `Uw10<=3.5: 0.2 Uw10; else 0.057 Uw10^2` | identical | match |
| 4 | `0.728 Uw10^0.5 - 0.317 Uw10 + 0.0372 Uw10^2` | identical | match |
| 5 | `0.0986 Uw10^1.64` | identical | match |
| 6 | `0.5 + 0.05 Uw10^2` | identical | match |
| 7 | `Uw10<=5.5: 0.362 sqrt(Uw10); else 0.0277 Uw10^2` | identical | match |
| 8 | `0.64 + 0.128 Uw10^2` | identical | match |
| 9 | `Uw10<=4.1: 0.156 Uw10^0.63; else 0.0269 Uw10^1.9` | identical | match |
| 10 | `0.0276 Uw10^2` | identical | match |
| 11 | `0.0432 Uw10^2` | identical | match |
| 12 | `0.319 Uw10` | identical | match |
| 13 | `Uw10<1.6: 0.398; else 0.155 Uw10^2` | identical | match |

`Uw10 = wind_speed * (10/2)^0.143` matches in all three.

**Finding:** match across all three on numerics. v3 docstring author attributions disagree with Fortran source comments (e.g., v3 calls option 2 "Banks 1975" while Fortran credits Broecker 1978). Recommend reconciling author attributions before LimnoTech review. **Severity:** observation.

#### `ka_tc` (combined, temperature-corrected)

Fortran: `modGlobalParam.f90:245-247`, `kah_tc + kaw_tc/depth`. v1: `shared/processes.py:165-178`. v3: `utils/reaeration.py:207-236`.

**Finding:** match. All three apply `Arrhenius_TempCorrection` to each component before summing. **Severity:** match.

### `sediment.py`

#### `SOD_tc`

Fortran: `modGlobalParam.f90:250-256`. v1: `shared/processes.py:180-200`. v3: `utils/sediment.py:16-31`.

- Fortran: `SOD_tc = Arrhenius(SOD%rc20, theta, T)`; if `use_DOX`, multiplies by `DOX/(DOX + KsSod)`.
- v1: same: `arrhenius_correction(...) * xr.where(use_DOX, DOX/(DOX+KsSOD), 1)`.
- v3: pure Arrhenius only. The DOX-Monod factor moved to the DOX Process call site.

**Finding:** v3 deliberate architectural refactor, documented in `parameter_defaults_corrections.md` Section 3.2. Numerically equivalent under matched fixtures (`tests/test_5_dox_calculations_v2.py::test_dox_sod_rate_matches_v1` passes `use_DOX=False` to v1 to obtain parity). **Severity:** intended improvement, not a deviation.

### `light.py`

#### `L` (Beer-Lambert extinction)

Fortran source: `modGlobalParam.f90:420-428`, `LightExtCoefficient`:

```
lambda = lambda0
do i=1,nGS:  lambda += lambdas * Solid(i)
if use_POC:  lambda += lambdam * POC / focm
if use_Algae: lambda += lambda1 * Ap + lambda2 * Ap^0.66667
```

v1 (`shared/processes.py:202-237`): identical, applies `lambdas * Solid` always (matching Fortran), then adds POC term inside `xr.where(use_POC,...)`, then adds algae term inside `xr.where(use_Algae,...)`.

v3 (`utils/light.py:13-53`): identical to v1.

**Finding:** match across all three. **Note:** `parameter_defaults_corrections.md` Section 2.8 incorrectly claims v1's `lambdas * Solid` is "commented out / defined but not used"; in fact `shared/processes.py:232` applies it unconditionally. v3 reproduces this. Recommend correcting the corrections doc Section 2.8.

Multi-solid difference: Fortran loops over `nGS` solid groups summing `lambdas * Solid(i)`. v3/v1 take a single scalar Solid concentration. For single-class use the formulas agree.

#### `PAR`

Fortran (`modGlobalParam.f90:222,234`):
```
real(R8) :: Fr_PAR = 0.47
if (use_Algae .or. use_BAlgae) PAR = q_solar * Fr_PAR
```

v1 (`shared/processes.py:240-253`): wraps with `xr.where(use_Algae or use_Balgae, q_solar*Fr_PAR)`. **Latent v1 bug:** `xr.where` with only two args returns NaN in the false branch.

v3 (`utils/light.py:56-70`): returns `q_solar * Fr_PAR` unconditionally; the `use_Algae/use_BAlgae` toggle moved to the consumer Process per Phase 1.1.

**Finding:** v3 deliberate refactor. Numerically equivalent inside the `use_Algae|use_BAlgae` branch. Avoids latent v1 NaN-propagation bug. **Severity:** v3 deliberate improvement.

### `partitioning.py`

#### `fdp`

Fortran (`modGlobalParam.f90:225-231`):
```
fdp = 1.0
do i=1,nGS:  fdp = fdp + kdpo4(i,r) * Solid(i) / 1.0E6
fdp = 1.0 / fdp
```

v1 (`shared/processes.py:256-271`): `xr.where(use_TIP, 1/(1 + kdpo4*Solid/0.000001), 0)` — single solid class.

v3 (`utils/partitioning.py:12-31`): identical to v1, with `0.000001` literal.

**Finding:** match (under single solid class). The `1.0E6` factor in Fortran (and `1e-6` denominator in v1/v3) is the unit conversion `(L/kg)(mg/L)(1 kg / 1e6 mg) = dimensionless`. The Phase 1.1 audit comment characterizing this as "suspicious" was incorrect — the formula is dimensionally consistent. All three references agree. **Severity:** match. Recommend retracting the Phase 1.1 "suspicious unit factor" flag.

Multi-solid difference: Fortran sums over `nGS` solid groups with possibly different `kdpo4(i,r)` per class. v3/v1 collapse to single scalar. If `nGS > 1` is ever activated, v3 would need a sum-over-classes form to match Fortran.

### `numerics.py` (v3-only)

`Diagnostics` dataclass + `clip_negative_state` at `utils/numerics.py:28-118`. No Fortran or v1 counterpart. Per Q7 in the design spec (Section 14), this is v3-only architectural infrastructure.

**Finding:** v3-only by design. Logic correct: clip target is exactly 0 (line 76), matching the Q7 contract. Detail-limit-per-call (default 10) properly rate-limits log records (line 84). Aggregate suppressed-count stub appended when n_clipped > limit (line 100). Returns DataArray with preserved coords/dims/attrs (line 112). **Severity:** intended addition, no issues.

---

## Part 2 — Parameter library

### Group: `algae` (`parameters/algae.py`, 17 entries)

All 17 v3 algae defaults match v1 exactly (`AWd=100`, `AWc=40`, `AWn=7.2`, `AWp=1`, `AWa=1000`, `KL=10`, `KsN=0.04`, `KsP=0.0012`, `mu_max_20=1`, `kdp_20=0.15`, `krp_20=0.2`, `mu_max_theta=1.047`, `kdp_theta=1.047`, `krp_theta=1.047`, `vsap=0.15`, `growth_rate_option=1`, `light_limitation_option=1`).

### Group: `balgae` (`parameters/balgae.py`, 19 entries)

All 19 entries in v3 match v1 byte-for-byte (`BWd=100`, `BWc=40`, `BWn=7.2`, `BWp=1`, `BWa=3500`, `KLb=10`, `KsNb=0.25`, `KsPb=0.125`, `Ksb=10`, `mub_max_20=0.4`, `krb_20=0.2`, `kdb_20=0.3`, `mub_max_theta=1.047`, `krb_theta=1.06`, `kdb_theta=1.047`, `b_growth_rate_option=1`, `b_light_limitation_option=1`, `Fw=0.9`, `Fb=0.9`).

### Group: `nitrogen` (`parameters/nitrogen.py`, 16 entries)

| Parameter | Fortran | v1 | v3 | Status |
|---|---|---|---|---|
| KNR | (modNitrogen) | 0.6 | 0.6 | match v1 |
| knit_20 | 0.1 | 0.1 | 0.1 | match v1 |
| kon_20 | 0.1 | 0.1 | 0.1 | match v1 |
| kdnit_20 | 0.002 | 0.002 | 0.002 | match v1 |
| rnh4_20 | 0 | 0 | 0.0 | match v1 — FIXME-flagged |
| vno3_20 | 0 | 0 | 0.0 | match v1 — FIXME-flagged |
| **vson_20** | **0.01** (in modGlobalParam.f90:92) | not in nitrogen group; v1 GlobalVars has `vson=0.01` | **0.1** | **structural relocation + 10x value change; undocumented** |
| knit_theta | 1.083 | 1.083 | 1.083 | match v1 |
| kon_theta | 1.074 | 1.074 | 1.074 | match v1 |
| kdnit_theta | 1.08 | 1.08 | 1.08 | match v1 |
| rnh4_theta | 1.047 | 1.047 | 1.047 | match v1 |
| vno3_theta | 1.045 | 1.045 | 1.045 | match v1 |
| **vson_theta** | n/a | not in v1 | **1.024** | **v3 addition, undocumented** |
| KsOxdn | 0.1 | 0.1 | 0.1 | match v1 |
| PN | 0.5 | 0.5 | 0.5 | match v1 |
| PNb | 0.5 | 0.5 | 0.5 | match v1 |
| use_OrgN | True | True (in GlobalParameters) | True | match (relocated) |

**Findings:**
- `vson_20=0.1` in v3 vs `vson=0.01` in v1 GlobalVars and `vson=0.01` in v3's own global_vars (`global_vars.py:26`). The v3 nitrogen group's `vson_20=0.1` is **10x v1's value, 10x v3's own global_vars value, and 10x Fortran's value**. Categorize as **undocumented v3 deviation needing review**.
- `vson_theta=1.024`: also new in v3, not in v1.
- Fortran `vson` (modGlobalParam.f90:92) initializes to `0.01`, matching v1's GlobalVars value.

### Group: `phosphorus` (`parameters/phosphorus.py`, 7 entries)

| Parameter | Fortran modGlobalParam | v1 | v3 | Status |
|---|---|---|---|---|
| kop_20 | (modPhosphorus) | 0.1 | 0.1 | match v1 |
| rpo4_20 | 0 | 0 | 0.0 | match v1 — FIXME-flagged |
| kop_theta | 1.047 | 1.047 | 1.047 | match v1 |
| rpo4_theta | 1.074 | 1.074 | 1.074 | match v1 |
| kdpo4 | line 82: `kdpo4 = 0.0` | 0.0 | 0.0 | match all three — FIXME-flagged |
| **vsop** | line 98: **0.01** | **999** | **0.1** | **likely v1 flaw; v3 corrects but to 10x Fortran's 0.01** |
| **vs** | line 87: **0.1** | **999** | **0.1** | **likely v1 flaw; v3 matches Fortran exactly** |

### Group: `carbon` (`parameters/carbon.py`, 10 entries)

All 10 v3 carbon parameters match v1 exactly: `f_pocp=0.9`, `kdoc_20=0.01`, `kdoc_theta=1.047`, `f_pocb=0.9`, `kpoc_20=0.005`, `kpoc_theta=1.047`, `KsOxmc=1.0`, `pCO2=383.0`, `FCO2=0.2`, `roc=32/12`.

### Group: `cbod` (`parameters/cbod.py`, 5 entries)

5 of 5 match v1 (`KsOxbod=0.5`, `kbod_20=0.12`, `ksbod_20=0.0`, `kbod_theta=1.047`, `ksbod_theta=1.047`).

### Group: `dox` (`parameters/dox.py`, 10 entries)

| Parameter | Fortran modGlobalParam | v1 GlobalVars | v3 dox | Status |
|---|---|---|---|---|
| ron | n/a | 2*32/14 = 4.5714 | 2*32/14 | match v1 |
| KsSOD | line 127: `KsSod = 1.0` | 1 | 1.0 | match all three |
| **SOD_20** | line 122: **0.2** | **999** | **1.0** | **all-three disagreement: Fortran=0.2, v1=999, v3=1.0** |
| **SOD_theta** | line 122: **1.06** | **999** | **1.060** | **likely v1 flaw; v3 matches Fortran** |
| **kaw_20_user** | line 117: **0.0** | **999** | **0.0** | **likely v1 flaw; v3 matches Fortran** |
| **kah_20_user** | line 113: **1.0** | **999** | **0.0** | **all-three disagreement: Fortran=1.0, v1=999, v3=0.0** |
| kaw_theta | line 117: 1.024 | 1.024 | 1.024 | match all three |
| kah_theta | line 113: 1.024 | 1.024 | 1.024 | match all three |
| hydraulic_reaeration_option | line 143: 1 | 1 | 1 | match all three |
| wind_reaeration_option | line 147: 1 | 1 | 1 | match all three |

**Critical findings:**
- **SOD_theta**: v3's correction to 1.060 vindicated by Fortran (Chapra 1997 standard). **Likely v1 flaw**.
- **SOD_20**: Fortran 0.2, v3 chose 1.0 (5x Fortran's value). Both defensible from literature; v3 deviates from Fortran by 5x. Recommend documenting.
- **kaw_20_user**: v3 matches Fortran exactly. **Likely v1 flaw, v3 correctly restored.**
- **kah_20_user**: Fortran 1.0, v3 0.0. **Behavioral consequence:** at default `hydraulic_reaeration_option=1`, v3 reaeration = 0, Fortran reaeration = 1.0 1/d. v3 makes user-override path explicit; Fortran has a non-zero hidden default. v3 + Fortran will produce different DOX trajectories at default settings.

### Group: `pathogen` (`parameters/pathogen.py`, 4 entries)

`kdx_20=0.8`, `kdx_theta=1.07`, `apx=1`, `vx=1`. All 4 match v1. `apx` and `vx` FIXME-flagged for unknown literature basis.

### Group: `alkalinity` (`parameters/alkalinity.py`, 6 entries)

All 6 stoichiometric ratios match v1 exactly.

### Group: `n2` (`parameters/n2.py`, 0 entries)

Empty in v3 and v1. Match.

### Group: `pom` (`parameters/pom.py`, 3 entries)

| Parameter | Fortran modGlobalParam | v1 | v3 | Status |
|---|---|---|---|---|
| kpom_20 | n/a | 0.1 | 0.1 | match v1 |
| h2 | line 134: `h2 = 0.1` | 0.1 | 0.1 | match all three |
| kpom_theta | n/a | 1.047 | 1.047 | match v1 |

### Group: `global_parameters` (`parameters/global_parameters.py`, 17 entries)

All 16 boolean `use_*` flags match v1 (use_NH4, use_NO3, use_OrgN, use_OrgP, use_TIP, use_SedFlux=False, use_POC, use_DOC, use_DOX, use_DIC, use_Algae, use_Balgae, use_N2, use_Pathogen, use_Alk, use_POM all True except SedFlux).

Fortran defaults differ: `use_BAlgae=.false.`, `use_POC=.false.`, `use_DOC=.false.`, `use_DIC=.false.`, `use_N2=.false.`, `use_Pathogen=.false.`, `use_Alk=.false.`, `use_POM2=.false.`. v3 enables all by default (matching v1's "all on" stance), differing from Fortran's selective-enable approach. **Inherited v1 convention**, not a deviation in v3 from v1.

| Parameter | Fortran | v1 GlobalVars | v3 | Status |
|---|---|---|---|---|
| **pressure_mb** | n/a in NSM1 (Fortran uses pressure_atm) | **2026.5** | **1013.25** | **likely v1 flaw, v3 correctly restored to ISO standard** |

### Group: `global_vars` (`parameters/global_vars.py`, 21 entries)

| Parameter | Fortran modGlobalParam | v1 GlobalVars | v3 | Status |
|---|---|---|---|---|
| vson | line 92: 0.01 | 0.01 | 0.01 | match all three |
| vsoc | line 104: 0.01 | 0.01 | 0.01 | match all three |
| theta | n/a | 1.047 | 1.047 | match v1 |
| **vb** | line 138: **0.0025 m/yr** (line 201 divides by 365) | **0.01 m/d** | **0.01 m/d** | **Fortran/v1 unit reconciliation issue**; FIXME-flagged in v3 |
| fcom (Fortran focm) | line 108: `focm = 0.4` | 0.4 | 0.4 | match all three |
| dt | n/a | 1 | 1.0 | match v1 |
| depth | n/a (runtime) | 1.5 | 1.5 | match v1 |
| TwaterC | n/a | 20 | 20.0 | match v1 |
| velocity | n/a | 1 | 1.0 | match v1 |
| flow | n/a | 2 | 2.0 | match v1 |
| topwidth | n/a | 1 | 1.0 | match v1 |
| slope | n/a | 2 | 2.0 | match v1 |
| shear_velocity | n/a | 4 | 4.0 | match v1 |
| wind_speed | n/a | 4 | 4.0 | match v1 |
| q_solar | n/a | 500 | 500.0 | match v1; FIXME for unit-doc mismatch |
| Solid | n/a | 1 | 1 | match v1 |
| lambda0 | line 60: 0.02 | 0.02 | 0.02 | match all three |
| lambda1 | line 72: 0.0088 | 0.0088 | 0.0088 | match all three |
| lambda2 | line 76: 0.054 | 0.054 | 0.054 | match all three |
| lambdas | line 64: 0.052 | 0.052 | 0.052 | match all three |
| **lambdam** | line 68: **0.174** | **0.0174** | **0.0174** | **likely v1 flaw — Fortran has 0.174, v1/v3 have 0.0174 (10x lower)** |
| Fr_PAR | line 222: 0.47 | 0.47 | 0.47 | match all three |

**Critical finding (lambdam):** Fortran `modGlobalParam.f90:68` initializes `lambdam = 0.174` while v1 and v3 use `0.0174` — 10x discrepancy in POM contribution to Beer-Lambert light extinction. Fortran's 0.174 is consistent with QUAL2K Table 6 references; v1's 0.0174 may be a typo. **Likely v1 flaw, propagated to v3. Requires reconciliation with LimnoTech before review.**

---

## Critical-correction verification

For each of the 7 corrections in `parameter_defaults_corrections.md` Section 1:

| Correction | v1 sentinel | v3 fix | Fortran value | Verdict |
|---|---|---|---|---|
| 1.1 vsop | 999 | 0.1 | **0.01** | v1 flaw confirmed; v3 corrects but to 10x Fortran's |
| 1.2 vs | 999 | 0.1 | **0.1** | v1 flaw confirmed; v3 matches Fortran exactly |
| 1.3 SOD_20 | 999 | 1.0 | **0.2** | v1 flaw confirmed; v3 corrects but to 5x Fortran's |
| 1.4 SOD_theta | 999 | 1.060 | **1.060** | v1 flaw confirmed; v3 matches Fortran exactly |
| 1.5 kaw_20_user | 999 | 0.0 | **0.0** | v1 flaw confirmed; v3 matches Fortran exactly |
| 1.6 kah_20_user | 999 | 0.0 | **1.0** | v1 flaw confirmed; v3 zeroes user-override branch but disagrees with Fortran's 1.0 default |
| 1.7 pressure_mb | 2026.5 | 1013.25 | **n/a** (Fortran uses pressure_atm) | v1 flaw confirmed (2x); v3 restored to ISO 2533 standard |

**Summary:** All 7 v3 corrections are vindicated as genuine v1 flaws. For 4 of 7 (`vs`, `SOD_theta`, `kaw_20_user`, `pressure_mb`) v3 chose a value matching or equivalent to Fortran. For 3 of 7 (`vsop`, `SOD_20`, `kah_20_user`) v3 chose a value differing from Fortran by O(1)-O(10) but defensible from literature.

---

## Conclusions

### v3 deliberate improvements (confirmed correct)

1. `pressure_mb=1013.25` (Section 1.7): v1 had 2026.5 (2x error). v3 restored to ISO standard.
2. `SOD_theta=1.060` (Section 1.4): v1 had sentinel 999. v3 matches Fortran's 1.06 and Chapra (1997).
3. `vs=0.1`, `kaw_20_user=0.0` (Sections 1.2, 1.5): v1 had 999. v3 restored Fortran-aligned defaults.
4. `vsop=0.1` (Section 1.1): v1 had 999. v3 chose 0.1 (10x Fortran's 0.01); literature-defensible.
5. `SOD_tc` pure Arrhenius split (Section 3.2): cleaner architecture; DOX-Monod moved to consumer Process.
6. `PAR` toggle inversion (Section 3.4 / Phase 1.1): avoids latent v1 NaN-on-disable bug from `xr.where(cond, value)` two-arg form.
7. `clip_negative_state` + `Diagnostics` (Q7): v3-only safety net, not a deviation.
8. `np.select` dim-stripping fix in `kah_20`/`kaw_20`: v3 implementation detail correctly preserves xarray broadcasting.

### Likely v1 flaws (v3 correctly bypassed; not v3 issues)

1. The 7 sentinel-999 / 2026.5 defaults — all confirmed by Fortran-coded defaults to be flaws in v1.
2. v1 `PAR` two-arg `xr.where` returns NaN when both algae modules disabled (Fortran/v3 avoid).

### Likely v1 flaw NOT corrected in v3 (needs action)

1. **`lambdam=0.0174`** in v1/v3 vs Fortran's `0.174` — 10x discrepancy in POM contribution to Beer-Lambert light extinction. Fortran's value consistent with QUAL2K Table 6; v1 likely has a typo that propagated to v3. **Reconcile with LimnoTech before review.**

### Undocumented v3 deviations (flag for review)

1. **`vson_20=0.1`** in `parameters/nitrogen.py:16`: v1 GlobalVars has `vson=0.01`, v3 global_vars also has `vson=0.01`, Fortran has 0.01. The new `vson_20=0.1` in nitrogen group is 10x v1's value and 10x v3's own global_vars value. **Internal v3 inconsistency.** Recommend either consolidate to 0.01 (matching Fortran/v1) or document the 0.1 choice with rationale.
2. **`vson_theta=1.024`** in `parameters/nitrogen.py:22`: not present in v1. New v3 parameter; document.
3. **`vsop=0.1` (Section 1.1)**: v3 chose 0.1 over Fortran's 0.01. Acknowledged in corrections doc but Fortran value not noted.
4. **`SOD_20=1.0` (Section 1.3)**: v3 chose 1.0 over Fortran's 0.2. Document.
5. **`kah_20_user=0.0` (Section 1.6)**: v3 chose 0.0 over Fortran's 1.0. Behavioral consequence: at default `hydraulic_reaeration_option=1`, v3 reaeration = 0, Fortran reaeration = 1.0 1/d. **Document this divergence prominently** — could surface as "v3 has no DOX recovery" puzzle in side-by-side runs.

### Documentation defects in `parameter_defaults_corrections.md`

1. Section 2.8 claims v1's `lambdas * Solid` term is "commented out / defined but not used"; the actual v1 source (`shared/processes.py:232`) applies it unconditionally. Recommend correcting.
2. Phase 1.1's "suspicious unit factor `1/(1 + kdpo4 * Solid / 0.000001)`" comment for `fdp`: the formula is dimensionally consistent. Recommend retracting the suspicion.
3. Section 3.6/3.7 (Kelvin offset, mb-to-atm scaling): v3's `utils/conversions.py` re-exports from v2, and v2's `celsius_to_kelvin` returns `T_C + 273.16` (with comment "for testing consistency with v1"), **not** 273.15 as Section 3.6 claims. Recommend updating Section 3.6 to reflect actual v3 behavior.

### Required actions before LimnoTech review

1. Resolve the **lambdam 0.174 vs 0.0174** discrepancy (Fortran vs v1/v3 10x mismatch).
2. Resolve the **vson_20 nitrogen group 0.1 vs global_vars 0.01** internal v3 inconsistency.
3. Update `parameter_defaults_corrections.md` Sections 1.1, 1.3, 1.6 to record the Fortran-coded default and the rationale for v3's chosen value where it differs from Fortran.
4. Correct `parameter_defaults_corrections.md` Section 2.8 (lambdas not commented out) and Section 3.6 (Kelvin offset).
5. Reconcile docstring author attributions in `utils/reaeration.py` (option labels disagree with Fortran source comments).

---

## Finding count

- **Utility findings:** 5 matches, 2 deliberate v3 improvements (PAR, SOD_tc), 1 v3-only addition (numerics), 0 undocumented deviations. Documentation issues: 3.
- **Parameter findings across ~145 entries:** ~135 match v1, 7 critical sentinel corrections (all vindicated), 1 likely v1 flaw not corrected (`lambdam`), 5 undocumented v3 deviations or value choices needing rationale (`vson_20`, `vson_theta`, `vsop` value vs Fortran, `SOD_20` value vs Fortran, `kah_20_user` value vs Fortran).
