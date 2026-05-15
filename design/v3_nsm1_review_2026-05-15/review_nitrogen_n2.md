# v3 NSM1 Nitrogen + N2 — Line-Level Source / Science / Documentation Review

Review date: 2026-05-15
Reviewer: water-quality model source-code reviewer (Claude)
Branch: `streaming` @ HEAD (54f2b12)
Scope (read line-by-line):

- `src/clearwater_modules_v3/processes/nitrogen.py` (995 lines)
- `src/clearwater_modules_v3/processes/n2.py` (444 lines)
- `src/clearwater_modules_v3/parameters/nitrogen.py` (69 lines)
- `src/clearwater_modules_v3/parameters/n2.py` (15 lines)

Cross-read for parity: `src/clearwater_modules/nsm1/processes.py` (v1 reference algorithm),
`src/clearwater_modules/nsm1/constants.py`, `static_variables.py`, `dynamic_variables.py`;
Fortran NSM1 `modNitrogen.f90` (`/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/`);
v3 `utils/conversions.py`, `utils/numerics.py`, `utils/reaeration.py`,
`processes/floating_algae.py`, `processes/benthic_algae.py`, `parameters/balgae.py`.
Classification authority: `design/clearwater_modules_v3_nsm1_audit_n_p.md`,
`clearwater_modules_v3_nsm1_audit_summary.md`,
`clearwater_modules_v3_nsm1_design_specification.md` (Section 6 bug list),
`parameter_defaults_corrections.md`.

---

## 1. Summary verdict

The v3 Nitrogen and N2 Processes are scientifically faithful to the v1
reference algorithm and, where they deviate, the deviations are documented
intentional improvements that align v3 with the authoritative Fortran NSM1.

The 2026-05-05 three-way audit (`audit_n_p.md`, `audit_summary.md`)
enumerated five critical Nitrogen defects in the **v2 overlay**
(`clearwater_modules_v2/processes/nitrogen.py`): the phantom
`ammonium_decay_nitrate` source, the legacy-kwarg wiring defect, the
static NO3 floating-algae uptake fraction, the broken benthic-algae NO3
uptake stoichiometry, and the sediment-release default-value wiring.
**All five are fixed in the scoped v3 file.** I traced each fix:

- The phantom NH4 source is removed from the `change_ammonium` rate sum
  (nitrogen.py:547--562); `ammonium_decay_nitrate` survives only as a
  dead, uncalled back-compat shim.
- Every kinetic method reads the `NITROGEN_DEFAULTS` attribute
  (`self.knit_20`, `self.kdnit_20`, `self.rnh4_20`, `self.vno3_20`,
  `self.KNR`, and their thetas), not the legacy `1.0`-valued kwargs.
- `nitrate_uptake_floating_algae` reads the dynamic
  `1 - floating_algae_process.algal_nh4_uptake_fraction`
  (nitrogen.py:810--817), so NH4 + NO3 algal uptake sum to
  `rna * algal_growth_rate`.
- `nitrate_uptake_benthic_algae` is reconstructed with `rnb = BWn/BWd`,
  `Fb`, dynamic `1 - balgae_nh4_uptake_fraction`, and the `/depth`
  divisor (nitrogen.py:861--875).
- `ammonium_from_bed` / `nitrate_bed_denitrification` read the v1/Fortran
  `rnh4_20=0` / `vno3_20=0` defaults.

The xarray refactor is complete on the live integration path. The
Forward-Euler-in-days integrator matches v1's `X + dXdt * dt` (v1 `dt`
units = days, confirmed at `static_variables.py:405`). NaN handling uses
`.isnull()` / `xr.where` with a correct scalar fallback branch
throughout. The `use_SedFlux=True` `NotImplementedError` guard is a
correctly-deferred NSM2 boundary, not a defect.

Findings are concentrated in three areas: (a) a class of stale or
now-misleading comments and docstrings that describe pre-fix behavior or
reference v1/v2 line numbers that no longer correspond to anything in v3;
(b) a residual mass-balance asymmetry in the `change_ammonium` NH4-uptake
path (it sums `ammonium_floating_growth + ammonium_benthic_growth` via the
algae caches while `change_nitrate` uses the Nitrogen-side
`nitrate_uptake_*` helpers — these are consistent with each other and with
v1, so this is an OBSERVATION, not a defect); and (c) a small number of
documentation-unit nits.

Severity counts: CRITICAL 0, MAJOR 0, MINOR 7, OBSERVATION 8.

Overall confidence in scientific correctness of the scoped Nitrogen/N2
code: **high**. No correctness-affecting defect was found on the live
`run` -> `_change_with_components` -> `change_*` path. The remaining work
is comment hygiene and documentation precision, plus retiring dead
back-compat shims to reduce the documentation-trap surface the audit
warned about.

---

## 2. Findings table

| ID | Severity | file:line | Category | Description | Recommended fix |
|----|----------|-----------|----------|-------------|-----------------|
| F1 | MINOR | nitrogen.py:23, 32, 38 | stale-comment | Class docstring says "v2 NSM1 Nitrogen Process (Phase 2.B fixes applied)" and "the 16 known v2 Nitrogen bugs are now fixed". The scoped file is the post-Phase-9 v3 refactor, not a v2 overlay; the Phase 9.A.2 wiring/algal fixes (the substantive corrections) are not mentioned in the class docstring, only in method-level comments. A reader of the class docstring would conclude Phase 2.B was the terminal fix state. | Update the class docstring to state this is the v3 Process and that the Phase 9.A.2 audit findings N1/N2/N4/N10/N11/N12/N13 and Phase 9.C/9.E parameter corrections are applied on top of Phase 2.B. |
| F2 | MINOR | nitrogen.py:62--82, 314--317, 347--351, 412--418 | stale-comment | Multiple comments reference "pre-Phase-4", "Phase 4 keeps the names exact", "DOX, Alkalinity, and N2 already consume them via getattr(nitrogen_process, ...)". These are accurate for the cache-name contract but the repeated "pre-Phase-4 / verbatim / bit-identical" framing describes a code-motion refactor history that is no longer load-bearing and obscures the current contract. | Compress the historical phase narrative into a single short note; keep only the load-bearing statement: the names `nitrification_flux_rate` / `denitrification_flux_rate` are a public contract consumed by DOX/Alkalinity/N2 and must not be renamed. |
| F3 | MINOR | nitrogen.py:727--738 (`ammonium_decay_nitrate`), 672--696 (`ammonium_uptake_floating_algae`), 698--699 (`ammonium_rate_settling`), 701--704 (`ammonium_rate_death`) | open-issue | Four dead methods retained "for back-compat". None is called on any live path (verified by grep across `src/clearwater_modules_v3` and `tests/`; only attribute-inspection tests touch the related attributes). `ammonium_decay_nitrate` in particular is the phantom source the audit flagged as critical; leaving a callable method that still computes `ammonium * arrhenius(T, ammonium_decay_rate=0.0, ...)` is harmless at the zeroed default but is exactly the "documentation trap" the audit summary warns about (Escalation item 2). | Either delete the four dead methods, or add a single explicit deprecation marker (e.g. raise on call, or `warnings.warn(DeprecationWarning)`) and a one-line note that they are not on the integration path. Do not leave them silently callable with the old semantics. |
| F4 | MINOR | nitrogen.py:790--805 docstring | documentation | The `nitrate_uptake_floating_algae` docstring narrates the *old broken* implementation ("previously used the static `float_algea_faction_uptake_from_nitrate` ... the two paths therefore did not sum") at length before describing current behavior. The bug narrative is longer than the behavior description and the first sentence a maintainer reads is about a defect that no longer exists. | Lead with the current behavior and equation; move the pre-fix narrative to a short "History (Phase 9.A.2 N12)" trailing note. Same pattern in `nitrate_uptake_benthic_algae` docstring (nitrogen.py:826--855). |
| F5 | MINOR | nitrogen.py:251--253 | documentation | Comments label `ALGAE_DEFAULTS["AWn"]` as "mg-N per stoichiometric unit" and `ALGAE_DEFAULTS["AWa"]` as "ug-Chla per stoichiometric unit", and `rna = AWn/AWa` is annotated "mg-N/ug-Chla" at line 816. This is internally consistent and correct, but the inline comment at 243--250 asserting the historical "1000x conservation bug (Phase 9.G commit ee31218)" is a commit-archaeology note that belongs in the corrections doc, not in the constructor body where it adds no actionable information. | Replace the commit-archaeology paragraph with a one-line pointer to `parameter_defaults_corrections.md`; keep the unit annotations. |
| F6 | MINOR | n2.py:133 (class docstring) | documentation | The N2 class docstring says the denitrification source is read "from `nitrogen_process.denitrification_rate`". The actual attribute read in `_change_with_components` (n2.py:424) is `denitrification_flux_rate`. `denitrification_rate` does not exist on the Nitrogen Process; only `denitrification_flux_rate` is set (nitrogen.py:271, 352). The docstring names a non-existent attribute. | Correct the class docstring to `denitrification_flux_rate`. The module-level docstring (n2.py:16) already uses the correct name; only the class docstring is wrong. |
| F7 | MINOR | n2.py:11, 99; nitrogen.py none | documentation | `pwv` docstring cites "v1 source: processes.py:2878-2886" and the formula `exp(11.8571 - 3840.70/T - 216961/T^2)`. This matches v1 `processes.py:2886` exactly (verified). However the n2.py module docstring (line 11) writes the constant as `216961/T^2` while the literature/v1 form is unitful in K; the value is correct but the docstring does not state that `T` must be Kelvin (the function `pwv(t_water_k)` does take Kelvin, so the code is correct; only the prose omits the unit). | Add "(T in K)" to the module docstring's `p_wv` formula line for unit clarity. Code is correct; this is prose-only. |
| O1 | OBSERVATION | nitrogen.py:547--562 | algorithm-parity | NH4 budget sign/term audit vs v1 `dNH4dt` (processes.py:1584): v3 = `-nitrification + from_bed + float_resp - float_growth + benthic_resp - benthic_growth + orgn_to_nh4`. v1 = `OrgN_NH4_Decay - NH4_Nitrification + NH4fromBed + NH4_ApRespiration - NH4_ApGrowth + NH4_AbRespiration - NH4_AbGrowth`. Term-by-term and sign-by-sign MATCH. Phantom `ammonium_decay_nitrate` correctly absent. | None — confirms the audit's critical N2 finding is fixed. |
| O2 | OBSERVATION | nitrogen.py:621--651 | algorithm-parity | NO3 budget vs v1 `dNO3dt` (processes.py:1729): v3 = `+nitrification - denitrification - bed_denitrification - float_uptake - benthic_uptake`. v1 = `NH4_Nitrification - NO3_Denit - NO3_BedDenit - NO3_ApGrowth - NO3_AbGrowth`. MATCH. | None — confirms audit findings N9/N12/N13 fixed. |
| O3 | OBSERVATION | parameters/nitrogen.py:57--60 | algorithm-parity | `kon_theta=1.047`, `kdnit_theta=1.045`, `rnh4_theta=1.074`, `vno3_theta=1.08`. These DIVERGE from v1 (`constants.py:134-137`: 1.074/1.08/1.047/1.045 — transposed in pairs) but MATCH Fortran NSM1 `modNitrogen.f90:77,82,89,95,100` exactly (verified: knit 1.083, rnh4 1.074, kon 1.047, kdnit 1.045, vno3 1.08). This is the documented Phase 9.E transposition correction (`parameter_defaults_corrections.md` Sections 1.10, 4.1). IMPROVED, not a discrepancy. | None — intentional v3 improvement, correctly documented and Fortran-anchored. |
| O4 | OBSERVATION | n2.py:419--432 | algorithm-parity | v3 `rate = atm_exchange + denit_source` where v1 `dN2dt = 1.034 * ka_tc * (N2sat - N2)` only (processes.py:3504; confirmed v1 has no N2 denitrification source). The `denit_source` term is the documented Category 2 #12 intentional improvement (closes NO3 -> N2 mass balance). Collapses exactly to the v1 form when `use_nitrogen=False`. IMPROVED(doc: audit_summary.md Cat 2 #12). | None. |
| O5 | OBSERVATION | n2.py:113--117 vs v1 processes.py:3484 | algorithm-parity | v3 `n2sat_henry` uses `MB_TO_ATM = 1/1013.25 = 9.86923e-4`; v1 `N2sat` uses literal `pressure_mb * 0.000986923`. `1/1013.25 = 9.869232667e-4`; v1's truncated literal `0.000986923`. Relative difference ~3e-7, far below any physical or numerical tolerance. The negative-saturation trap (`< 0 -> 1e-6`) matches v1. MATCH (within float rounding). | None; optionally note in a comment that `MB_TO_ATM` is the exact reciprocal vs v1's truncated literal. |
| O6 | OBSERVATION | nitrogen.py:877--885 (`nitrification_inhibition`) | algorithm-parity | v3 returns `1 - exp(-self.KNR * DOX)` and `1.0` when `not self.use_nitrate`. v1 `NitrificationInhibition` (processes.py:1437) returns `xr.where(use_DOX, 1 - exp(-KNR*DOX), 1.0)`. v3 gates on `use_nitrate` rather than `use_DOX`; in v3 1.0.0 NSM1 the DOX state is always present when nitrate is active, so the behavioral difference is nil under supported configurations. KNR reads `self.KNR=0.6` (NITROGEN_DEFAULTS), matching v1/Fortran (audit N1 fix confirmed). MATCH under supported configs. | None; if a `use_DOX=False` configuration is ever supported, revisit the gating predicate. |
| O7 | OBSERVATION | nitrogen.py:457--472 | xarray | `_change_with_components` sets `orgn_settling = 0.0` (Python scalar) in the `use_OrgN=False` else-branch while other components may be `xr.DataArray`. This scalar is only stored in the `components` dict for diagnostics and is never arithmetically combined with arrays here (the integrator gets `orgn_rate` from `change_organic_nitrogen`, which early-returns `0.0` consistently). No broadcasting hazard on the live path. The opportunistic registry write (`set_at_time`) of a scalar `0.0` to a multi-cell registry variable would broadcast correctly under xarray assignment. OBSERVATION (needs verification only if a future consumer does array arithmetic on `components["orgn_settling_rate"]`). | None now; note for future consumers. |
| O8 | OBSERVATION | nitrogen.py:480--503 | algorithm-parity | The `components` dict computes `nh4_algal_growth = ammonium_floating_growth() + ammonium_benthic_growth()` and `no3_algal_growth = nitrate_uptake_floating_algae(...) + nitrate_uptake_benthic_algae(...)`. These are diagnostic recomputes that mirror the same helper calls inside `change_ammonium`/`change_nitrate`; arguments are identical so values are bit-identical. The NH4-uptake path routes through the FloatingAlgae/BenthicAlgae `ammonium_growth()` caches (dynamic `algal_nh4_uptake_fraction`) and the NO3 path through the Nitrogen-side helpers (dynamic `1 - algal_nh4_uptake_fraction`); the two therefore sum to `rna * growth` (mass balance closes), matching v1. MATCH. | None — confirms audit N6/N12 closure. |

---

## 3. Algorithm parity matrix

Reference: v1 `src/clearwater_modules/nsm1/processes.py`. Verdicts:
MATCH = v3 math equals v1 math; IMPROVED = documented intentional v3
deviation aligned to Fortran/spec; DISCREPANCY = unexplained divergence.

| v3 term (nitrogen.py) | v1 reference | Verdict |
|---|---|---|
| `ammonium_nitrification` (706--725): `NH4 * arrhenius(T,knit_20,knit_theta) * (1-exp(-KNR*DOX))` | `NH4_Nitrification` (1454) = `inhibition * knit_tc * NH4` | MATCH (KNR now `self.KNR=0.6`; audit N1 fix verified) |
| `nitrification_inhibition` (877--885): `1 - exp(-KNR*DOX)` | `NitrificationInhibition` (1437) | MATCH (gating predicate differs but nil under supported configs; see O6) |
| `nitrate_denitrification` (740--767): `NO3 * arrhenius(T,kdnit_20,kdnit_theta) * (1 - DOX/(DOX+KsOxdn))`, NaN->0 | `NO3_Denit` (1623--1637) `np.select` with NaN fallback `kdnit_tc*NO3` | MATCH in formula; NaN-fallback differs (v3 -> 0, v1/Fortran -> `kdnit_tc*NO3`). Fortran reachable only at `DOX = -KsOxdn` (unphysical). Documented audit N10 minor; no practical consequence. MATCH (practical) |
| `nitrate_bed_denitrification` (769--785): `NO3 * arrhenius(T,vno3_20,vno3_theta) / depth` | `NO3_BedDenit` (1655) = `vno3_tc * NO3 / depth` | MATCH; `vno3_20=0` v1/Fortran default (audit N11 fix verified) |
| `ammonium_from_bed` (660--670): `arrhenius(T,rnh4_20,rnh4_theta) / depth` | `NH4fromBed` (1470) = `rnh4_tc / depth` | MATCH; `rnh4_20=0` v1/Fortran default (audit N4 fix verified) |
| `organic_nitrogen_to_ammonium_hydrolysis` (896--910): `arrhenius(T,kon_20,kon_theta) * OrgN`, `use_OrgN` gate | `OrgN_NH4_Decay` (1330) = `where(use_OrgN, kon_tc*OrgN, 0)` | MATCH |
| `organic_nitrogen_settling` (912--934): `vson_20 / depth * OrgN` (raw vson, no Arrhenius) | `OrgN_Settling` (1345) = `vson / depth * OrgN` | MATCH (Phase 9.E removed the erroneous `vson_theta`; now exactly v1/Fortran) |
| `nitrate_uptake_floating_algae` (787--817): `rna * algal_growth_rate * (1 - algal_nh4_uptake_fraction)`, `rna=AWn/AWa` | `NO3_ApGrowth` (1675) = `(1-ApUptakeFr_NH4) * rna * ApGrowth` | MATCH (dynamic fraction; audit N12 fix verified) |
| `nitrate_uptake_benthic_algae` (819--875): `(1-balgae_nh4_uptake_fraction) * (BWn/BWd) * Fb * balgae_growth_rate / depth` | `NO3_AbGrowth` (1697) = `(1-AbUptakeFr_NH4) * rnb * Fb * AbGrowth / depth`, `rnb=BWn/BWd` | MATCH (audit N13 three-defect fix verified: BWd not AWa; `/depth` present; Fb not fraction_bottom_area; dynamic fraction) |
| `ammonium_floating_respiration` -> `floating_algae.ammonium_respiration` = `rna*algal_respiration_rate` | `NH4_ApRespiration` (1486) = `rna * ApRespiration` | MATCH |
| `ammonium_floating_growth` -> `floating_algae.ammonium_growth` = `algal_nh4_uptake_fraction*rna*algal_growth_rate` | `NH4_ApGrowth` (1504) = `ApUptakeFr_NH4 * rna * ApGrowth` | MATCH (dynamic fraction recomputed in FloatingAlgae.run) |
| `ammonium_benthic_respiration` -> `benthic_algae.ammonium_respiration` = `rnb*balgae_respiration_rate*Fb/depth` | `NH4_AbRespiration` (1525) = `(rnb*AbRespiration*Fb)/depth` | MATCH |
| `ammonium_benthic_growth` -> `benthic_algae.ammonium_growth` = `balgae_nh4_uptake_fraction*rnb*Fb*balgae_growth_rate/depth` | `NH4_AbGrowth` (1547) = `(AbUptakeFr_NH4*rnb*Fb*AbGrowth)/depth` | MATCH |
| `organic_nitrogen_from_floating_algae_mortality` -> `algal_orgn_from_mortality_rate` = `rna*ap_death` | `ApDeath_OrgN` (1360) = `rna * ApDeath` | MATCH |
| `organic_nitrogen_from_benthic_algae_mortality` -> `balgae_orgn_from_mortality_rate` = `rnb*Fw*Fb*ab_death/depth` | `AbDeath_OrgN` (1381) = `rnb*Fw*Fb*AbDeath/depth` | MATCH |
| `change_ammonium` (519--570) NH4 budget | `dNH4dt` (1584) | MATCH (term/sign verified, O1) |
| `change_nitrate` (592--658) NO3 budget | `dNO3dt` (1729) | MATCH (term/sign verified, O2) |
| `change_organic_nitrogen` (960--994): `ap_death+ab_death - hydrolysis - settling`, `use_OrgN` gate | `dOrgNdt` (1402) = `where(use_OrgN, ApDeath_OrgN+AbDeath_OrgN-OrgN_NH4_Decay-OrgN_Settling, 0)` | MATCH |
| `run` integrator (356--359): `X_new = X + rate * dt_days`, `dt_days = total_seconds()/86400` | `NH4/NO3/OrgN` (1602/1747/1420) = `X + dXdt * dt`, `dt` units = days | MATCH (additive Forward Euler; v1 `dt` confirmed days at static_variables.py:405) |
| `ammonium_uptake_floating_algae` (672--696) | `ApUptakeFr_NH4` (1226--1247) | MATCH in formula but DEAD CODE (not on live path; FloatingAlgae owns the dynamic fraction). See F3. |
| `khn2_tc` (n2.py:83--91) | `KHN2_tc` (3467) = `0.00065*exp(1300*(1/T-1/298.15))` | MATCH |
| `pwv` (n2.py:94--99) | `pwv` (2886) = `exp(11.8571 - 3840.70/T - 216961/T^2)` | MATCH |
| `n2sat_henry` (n2.py:102--117) | `N2sat` (3484) = `2.8e4*KHN2_tc*0.79*(P_mb*0.000986923 - pwv)`, neg-trap | MATCH (mb->atm constant differs by <3e-7; O5) |
| N2 `atm_exchange` (n2.py:419) = `1.034*ka_tc*(N2sat-N2)` | `dN2dt` (3504) = `1.034*ka_tc*(N2sat-N2)` | MATCH |
| N2 `rate = atm_exchange + denit_source` (n2.py:432) | `dN2dt` (3504), no denit source | IMPROVED (Cat 2 #12; collapses to v1 when use_nitrogen=False; O4) |
| N2 `tdg = n2_new / n2_sat` (n2.py:330) | `TDG` (3541) = `where(use_DOX, 79*N2/N2sat + 21*DOX/DOX_sat, N2/N2sat)` | DEVIATION — v3 1.0.0 implements only the simple `N2/N2sat` form; O2-weighted form deferred to Phase 5. Documented (n2.py:23--24; design spec). Correctly-deferred, not a discrepancy. |

No DISCREPANCY rows. Every divergence from v1 is either a documented
intentional improvement aligned to Fortran (O3, O4) or a deferred-to-NSM2
/ deferred-to-Phase-5 scope boundary (TDG O2-weighted form).

---

## 4. Stale-comment list

Each item was checked against the current code; classification per the
review brief (a comment that still says "broken/returns 0/not implemented"
when the code is fixed is a FINDING).

1. **nitrogen.py:23, 32, 38 (class docstring)** — STALE (F1). Describes
   the Process as the "v2 NSM1 Nitrogen Process (Phase 2.B fixes
   applied)". The substantive Phase 9.A.2 correctness fixes (the ones the
   audit graded critical) are not reflected in the class docstring. The
   docstring is not wrong about Phase 2.B but is incomplete and frames
   Phase 2.B as terminal.

2. **nitrogen.py:62--82, 314--317, 347--351, 412--418** — STALE-LEANING
   (F2). Heavy "pre-Phase-4 / verbatim / bit-identical / Phase 4 keeps
   the names exact" narrative. The load-bearing fact (cache names are a
   public contract) is buried in refactor-history prose. Not incorrect,
   but the phase archaeology is no longer actionable and obscures the
   contract.

3. **nitrogen.py:790--805 (`nitrate_uptake_floating_algae` docstring)**
   and **826--855 (`nitrate_uptake_benthic_algae` docstring)** —
   MISLEADING-ON-FIRST-READ (F4). The docstrings lead with a long
   description of the *old broken* implementation. The behavior is
   correct in code; the docstring structure inverts importance.

4. **nitrogen.py:243--250** — ARCHAEOLOGY-IN-CODE (F5). The inline
   "Phase 9.G commit ee31218 ... 1000x ... closed-system N conservation
   bug" paragraph in the constructor body. The accompanying unit
   annotations are correct and should stay; the commit-archaeology
   sentence belongs in the corrections doc.

5. **n2.py:133 (class docstring)** — WRONG (F6). States the
   denitrification source is read from
   `nitrogen_process.denitrification_rate`. The code reads
   `denitrification_flux_rate` (n2.py:424); `denitrification_rate` does
   not exist on the Nitrogen Process. This is a factual error in the
   docstring, not merely stale phrasing.

6. **nitrogen.py:117--134 (`use_SedFlux` guard comment)** — ACCURATE,
   NOT STALE. The comment correctly describes the current behavior: v3
   1.0.0 does not gate `ammonium_from_bed` / `nitrate_bed_denitrification`
   by `use_SedFlux`; the `rnh4_20=0` / `vno3_20=0` defaults are the de
   facto gate; `use_SedFlux=True` raises `NotImplementedError`. Verified
   against code (nitrogen.py:126--134) and `parameters/nitrogen.py:53-54`.
   No finding. Listed here to record that this flagged marker line was
   checked and is correct.

7. **nitrogen.py:531--546, 603--607, 629--631, 715--718, 750--753,
   775--778, 881--884 ("Phase 9.A.2 audit finding Nx" comments)** —
   ACCURATE, NOT STALE. Each of these comments describes a fix that is in
   fact applied in the adjacent code. I verified every one: N1 (KNR/knit
   wiring), N2 (phantom decay dropped), N4 (rnh4 default), N10 (kdnit
   wiring), N11 (vno3 wiring), N12 (dynamic NO3 float fraction), N13
   (benthic reconstruction). These comments correctly describe corrected
   status and are appropriate to keep. No finding.

8. **nitrogen.py:308--310 ("Bug #1 / #2 / #16 ... preserved")** —
   ACCURATE. Additive Forward Euler (356--359), `set_at_time` persistence
   (369--372), clip-with-log (362--366) all present and correct. Bug #1,
   #2, #16 from design spec Section 6 are genuinely fixed. No finding.

Net: F1, F2, F4, F5, F6 are the actionable stale/misleading-comment
findings. Items 6, 7, 8 are explicitly recorded as checked-and-correct so
the marker lines the brief flagged are not left ambiguous.

---

## 5. Correctly-deferred-to-NSM2 / later-phase list

These are explicit, documented scope deferrals. They are NOT defects.

1. **`use_SedFlux=True` -> `NotImplementedError`** (nitrogen.py:126--134).
   The full sediment-flux feature requires the NSM2 diagenesis path.
   Refusing explicitly (rather than silently producing partial behavior)
   is the correct design. Backed by `parameters/nitrogen.py:53-54` and
   `parameter_defaults_corrections.md` Section 2.1. Correctly deferred.

2. **Constant sediment release via `rnh4_20` / `vno3_20`**
   (nitrogen.py:660--670, 769--785; parameters/nitrogen.py:53--54). The
   ungated `ammonium_from_bed` / `nitrate_bed_denitrification` consumers
   are intentionally held silent by the `0.0` defaults; site-specific
   constant release is available by setting `rnh4_20` / `vno3_20`
   directly without `use_SedFlux`. This matches v1/Fortran defaults
   (`rnh4=0`, `vno3=0`). Correctly deferred / correctly defaulted.

3. **N2 O2-weighted TDG form** (n2.py:23--24, 325--336). v3 1.0.0
   implements `TDG = N2 / N2sat` (the v1 `use_DOX=False` branch). The
   O2-weighted form `0.79*N2/N2sat + 0.21*DOX/DOX_sat` requires the DOX
   Process and is deferred to Phase 5. Documented in the module and class
   docstrings and the design spec. Correctly deferred.

4. **`ammonium_decay_rate` / `ammonium_decay_theta` attributes retained
   but dropped from the NH4 budget** (nitrogen.py:192--203, 727--738).
   These are zeroed-default back-compat attributes with no v1/Fortran
   analogue; the audit's critical N2 finding is resolved by removing the
   term from `change_ammonium`. Retaining the inert attributes is a
   deliberate back-compat choice (not NSM2 scope, but a deliberate
   non-defect). See F3 for the recommendation to make the dead
   `ammonium_decay_nitrate` method non-silently-callable; the attributes
   themselves are correctly inert.

---

## 6. Cross-checks performed and their outcomes

- **Integrator parity.** v1 `dt` units confirmed = days
  (`static_variables.py:405`); v1 state update `X + dXdt*dt` (additive).
  v3 `X + rate*dt_days`, `dt_days = total_seconds()/86400`. Equivalent.
  No Forward-Euler-vs-multiplicative defect (design spec Bug #1/#2 fixed).
- **xarray completeness on the live path.** All NaN guards use
  `xr.where(rate.isnull(), ...)` for `xr.DataArray` and a `np.where`/
  scalar fallback otherwise (nitrogen.py:566--569, 654--657, 694--696,
  765--767, 990--993; n2.py via `sanitize_rate`). No `== np.nan`, no
  truthiness test on arrays, no Python loop over cells on the integration
  path. `clip_negative_state` and `sanitize_rate` are container-type-aware
  and preserve coords/dims. Bug #5--#8 fixed.
- **Algae cache contract.** `floating_algae.ammonium_growth/respiration`,
  `algal_growth_rate`, `algal_nh4_uptake_fraction`,
  `algal_orgn_from_mortality_rate` and the benthic equivalents
  (`BWn/BWd`, `Fb`, `balgae_nh4_uptake_fraction`,
  `balgae_orgn_from_mortality_rate`, `balgae_growth_rate`) exist and carry
  the v1-matching formulas (verified in floating_algae.py / benthic_algae.py
  and parameters/balgae.py). The Nitrogen Process consumes them via
  `getattr` with safe numeric fallbacks.
- **Fortran anchor for the Phase 9.E theta correction.** Read
  `modNitrogen.f90:77,82,89,95,100` directly: `knit%theta=1.083`,
  `rnh4%theta=1.074`, `kon%theta=1.047`, `kdnit%theta=1.045`,
  `vno3%theta=1.08`. v3 `parameters/nitrogen.py` matches Fortran exactly;
  v1 has the pairs transposed. The v3 deviation from v1 is the correct
  action.

---

## 7. Recommended follow-up tests / hygiene

These strengthen confidence but are not blockers for the scoped code's
correctness.

1. **Mass-balance closure test** for the algal-N split: over a
   multi-step integration with non-trivial NH4, NO3, and active algae,
   assert `nh4_algal_growth_rate + no3_algal_growth_rate ==
   rna * algal_growth_rate` (floating) and the benthic analogue, to
   within tolerance. The audit (Pattern D) notes existing parity tests do
   not exercise mass-balance closure; the v3 code is correct here but the
   property is unguarded by a test.
2. **Default-instantiation regression**: assert `Nitrogen()` (no kwargs)
   produces v1-matching `change_ammonium` / `change_nitrate` /
   `change_organic_nitrogen` for a fixed fixture, specifically to lock in
   the Phase 9.A.2 wiring (the audit's root-cause Pattern A was that
   tests passed explicit kwargs that masked wiring).
3. **N2 collapse test**: assert that with `use_nitrogen=False`, v3 N2
   reproduces the v1 `dN2dt = 1.034*ka_tc*(N2sat-N2)` trajectory exactly,
   confirming the denit-source extension is additive-only.
4. **Stale-comment lint**: a doc-fidelity check that the cache names in
   the N2 class docstring (F6) match the attributes actually read.

---

## 8. Open questions

1. **Dead back-compat shims (F3).** Should `ammonium_decay_nitrate`,
   `ammonium_uptake_floating_algae`, `ammonium_rate_settling`,
   `ammonium_rate_death` be deleted now, or kept inert until the v1/v2
   retirement in v3 1.1.0? The audit summary (Escalation 2) recommends
   full retirement; the in-code comments commit to back-compat. This is a
   product decision, not a correctness defect.
2. **`use_DOX=False` configurability (O6).** `nitrification_inhibition`
   gates on `use_nitrate` where v1 gated on `use_DOX`. Under all v3 1.0.0
   supported configurations DOX is present, so this is nil. Confirm
   whether a `use_DOX=False` NSM1 configuration is in scope for v3 1.x; if
   so, the gating predicate should be revisited.
