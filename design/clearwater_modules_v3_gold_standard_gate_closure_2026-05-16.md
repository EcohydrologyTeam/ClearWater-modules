# ClearWater Modules v3 — Gold-Standard Gate Closure Record

**Date:** 2026-05-16
**Status:** as-executed closure of `design/clearwater_modules_v3_gold_standard_specification.md` (the plan of record, committed `216ab6a`). This document records what was done; it does not change scope.
**Branch:** `streaming`

This is the Workstream-E capstone: the E1–E5 author/PI decisions resolved with commit pointers, the NSM2-deferral guard audit, and the §1/§10 exit-criteria checklist.

## 1. Workstream commit lineage (A → D)

| Item | Fix commit | Re-baseline commit | Trajectory? |
|---|---|---|---|
| A1 · NSM1-CA-1 (CRITICAL) — alkalinity intensive C ratio | `624ed7c` | `d16a2ea` (→`624ed7c` baseline) | alkalinity only |
| A2 · NSM1-SCI-N1 — `r_alkden` 1 eq/mol-N | `b51df71` | `8df4f37` (→`b51df71`) | alkalinity only |
| B1 · NSM1-SCI-A3 — algae PAR (Fr_PAR=0.47) | `3a8c188` | `ac8265f` (→`3a8c188`) | broad cascade |
| C1 · NSM1-SCI-A2 — `f_pocp/f_pocb`=0.8 | `6c10f36` | `7c9e798` (→`6c10f36`) | broad (carbon) |
| C2 · NSM1-SCI-CB1 — CBOD 1/d, θ=1.024 | `d2cdb4c` | none | dormant (`ksbod_20=0`) |
| C3 · NSM1-DOX-F1 — freshwater DO-sat doc + guard | `bc1af16` | none | non-perturbing |
| C4 · NSM1-DOX-F2 — silent-zero warn + opt-in floor | `cb2d8c2` | none | non-perturbing (floor OFF) |
| C5 · NSM1-SCI-A1 — documented NSM2 deferral | `3ff50cb` + `6d6de18` | none | doc-only |
| D1 · NSM1 doc-staleness cluster | `3eb1d7f` | none | doc-only |
| D2 · TSM doc hygiene + wind_c callout | `331c2dc` | none | doc-only |
| D3 · MMS test + joint bench + zero-defaults note | `d59da59` | none | test/doc-only |

**Terminal Phase-0 baseline:** `6c10f36`
(`baseline_coupled_trajectory_6c10f36.nc`). A1/A2/B1/C1 are the only
trajectory-perturbing changes; C2–C5/D do not perturb the coupled-demo
trajectory (re-verified per item — coupled-demo parity stayed green
through C2–D3). NSM2 plan committed `bb362c0`.

## 2. Workstream-E author/PI decisions (E1–E5) — resolved

| ID | Decision | Resolution | Where |
|---|---|---|---|
| **E1** | SCI-A2 `f_pocp/f_pocb` | **0.8** (CE-QUAL-W2 `APOM`; cited deliberate value, `mu_max_20` precedent) | C1 `6c10f36`; `parameter_defaults_corrections.md` §3.9 |
| **E2** | SCI-A1 timing | **Stage with NSM2**; documented known limitation now | C5 `3ff50cb`; NSM2 plan §4.B + Phase S4-3 `6d6de18`; corrections §3.11 |
| **E3** | SCI-CB1 form | **Match Fortran 1/d; θ=1.024** (per research doc) | C2 `d2cdb4c`; corrections §2.3/§3.5 |
| **E4** | TSM 273.15 vs 273.16 | **Hold 273.15 (SI).** No v3 1.0 action — v3 `utils/constants.py:29` `KELVIN_OFFSET=273.15`; the 273.16 reconciliation is a *conditional* LimnoTech item, only if v2 numerical parity is renegotiated (TSM review §5.3). Not a gold-standard gate item. | `utils/constants.py:29` (unchanged, already SI) |
| **E5** | MMS test | **Land before the gold-standard claim** | D3 `d59da59` — `test_mms_end_to_end_water_sediment_energy_conservation_under_ramp` |

## 3. NSM2-deferral guard audit (spec §8 / review §7) — all guarded, none mislabeled

| Deferral | Guard / disposition (verified) |
|---|---|
| `Nitrogen use_SedFlux=True` | **`raise NotImplementedError`** at `nitrogen.py:136-144` (explicit refuse, message + corrections §2.1 pointer). Verified in code. |
| Sediment-flux release (`rnh4_20`, `vno3_20`, `rpo4_20` = 0) | Intentionally-zero inherited defaults; term identically zero. Documented: corrections §2.1 + the new consolidated reader-note box (D3). |
| `kdpo4=0` TIP partitioning | Intentionally-zero; `fdp≡1`. Corrections §2.2 + consolidated box. v3 `fdp` unit form is correct (D1 corrected the stale audit_n_p claim). |
| Full carbonate/pH solver | Documented NSM2 scope (alkalinity is a simple tracer in v3 1.0); NSM2 plan §4.B / Phase S4-3. |
| SCI-A1 N-flux algal-alkalinity reformulation | Documented known limitation (alkalinity docstring + corrections §3.11); scheduled NSM2 plan §4.B / Phase S4-3 (`6d6de18`). |
| SOD-derived DIC sediment release; alkalinity DOX-Monod single-source; POM→DOC source; multi-group algae | Documented deferrals (review §7; NSM2 plan / multi-group-algae design spec). Not defects. |

## 4. Exit-criteria checklist (spec §1 (1)–(5) + §10)

- **§1(1) Correctness** — ✅ No CRITICAL or science-MAJOR open. CA-1 (CRITICAL) fixed; SCI-N1/A3/A2/CB1 (MAJOR) fixed; DOX-F1/F2 addressed; every divergence reference-anchored (CE-QUAL-W2 / Stumm & Morgan / QUAL2K / Bowie) and documented inline + in `parameter_defaults_corrections.md`.
- **§1(2) v3 ⊇ v1** — ✅ Fixes either restore v1 (CA-1: v1 was correct; SCI-A3: restored v1 Fr_PAR) or are deliberate documented divergences from v1/Fortran where v1 was itself wrong (SCI-N1, SCI-CB1) or sub-optimal (wind_c). Parity-≠-correctness payoffs retained as evidence in the audit/review docs.
- **§1(3) Documentation** — ✅ D1 cleared the NSM1 stale-doc cluster (README, headers, n2 false attr, carbon dead-refs, audit_n_p stale-fdp); D2 cleared TSM hygiene + the prominent wind_c=2.0≠v1-3.0 reviewer callout. Intentional divergences flagged at point of use.
- **§1(4) Tests** — ✅ Every gate fix has a **non-shared-path** regression test (hardcoded literals, not the process's own symbols): CA-1, SCI-N1, SCI-A3, SCI-A2, CB1, DOX-F1, DOX-F2, plus the joint CA-1+SCI-N1 closed-system benchmark. The highest-value deferred conservation test (TSM MMS, §5.2) is **landed** (E5). Frozen parity references that embedded defects (CA-1 algal refs, SCI-N1 denit ref, CBOD `/depth` ref) were re-derived/reframed, not left masking the bug.
- **§1(5) Deferrals clean** — ✅ See §3 above; all guarded, documented as deferrals, not mislabeled as defects.
- **§10 exit list** — ✅ CA-1 & SCI-N1 fixed w/ non-shared-path tests; SCI-A3 resolved (confirmed by trace, fixed); SCI-A2/CB1/DOX-F1/DOX-F2 fixed; SCI-A1 documented + NSM2-scheduled; NSM1 doc cluster + TSM doc hygiene cleared; MMS test landed.

## 5. The one remaining exit item — v1 archival (author action, NOT done here)

Spec §1(2)/§10: "v1 is deprecated and archived **only after this gate passes**." The correctness/documentation/test gate is now satisfied (§4). **Archiving/removing v1 (`src/clearwater_modules/`) is a separate, outward-facing repository action requiring explicit author authorization** — it is intentionally **not** performed as part of this gate work (and v1 source is still needed in-tree by the faithful v1-parity test re-derivations until those references are frozen independently). Recommended as a distinct, signed-off follow-up commit once the gate closure is accepted.

## 6. Status

All Workstream A–E correctness, documentation, and test gates are met. Full v3 suite green at `d59da59` (595 passed, 0 failed; coupled-demo bit-identical to the terminal `6c10f36` baseline). Pending author actions: (a) accept this closure; (b) authorize v1 deprecation/archival as a separate commit.
