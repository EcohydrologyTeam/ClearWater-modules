# v1 Retirement Plan

**Status:** Tracked plan for retiring `src/clearwater_modules/` (v1) and making `src/clearwater_modules_v3/` the only supported version.

**Target window:** Next few weeks.

---

## 1. Current state (as of 2026-05-10)

### v3 coverage

`tests/v3/` has **284 passing tests, 0 xfailed**, covering:

- v3 TSM (temperature) — `test_5_tsm_calculations_v3.py`, `test_tsm_*_v3.py` (10 modules)
- v3 NSM1 kinetics and integrator (Nitrogen, Phosphorus, Carbon, DOX, CBOD, POM, N2, Algae, Pathogen, Alkalinity) — `tests/v3/nsm1/test_*_tier1.py` and adjacent
- v3 Model orchestration, hotstart, wet/dry handling, scheduling
- Tier 1 + Tier 1.5 closed-system mass conservation (helper + active-kinetics)
- Cross-cutting numerics (NaN propagation, sanitize_rate, mixing-ratio)

### What's still in `tests/` (top-level)

| Category | Files | Status |
|---|---|---|
| Pure v1 framework | `test_2_equation_sort.py`, `test_3_model.py`, `test_4_tsm_module.py`, `test_5_tsm_calculations.py`, `test_6_nsm_module.py` | Delete with v1 src |
| Pure v1 NSM1 kinetics | `test_7..17_nsm_*_calculations.py` (11 files) | Delete with v1 src |
| Pure v1 hotstart/TSM | `test_hotstart_roundtrip.py`, `test_tsm_latent_heat.py`, `test_tsm_stability_ramp.py` | v3 equivalents exist; delete with v1 src |
| Parity tests (v3 vs v1) | `test_5_*_v2.py` (11 files, despite name) | **Needs decision — see Section 3** |
| Pure v3 (misplaced) | `test_phase9fc_documentation.py`, `test_sanitize_rate.py` | Move into `tests/v3/` |
| Empty/no clearwater imports | `test_1_imports.py` | Keep or delete |

Total non-v3 tests today: **521 passing + 2 xfailed**.

### Source layout

```
src/clearwater_modules/         <- v1 (to retire)
   nsm1/                        <- v1 NSM kinetics + Model
   tsm/                         <- v1 TSM kinetics + Model
   simlab/, csm/, msm/, nsm2/   <- legacy/empty stubs, also retire
   base.py, sorter.py, etc.     <- v1 framework
src/clearwater_modules_v3/      <- v3 (canonical going forward)
```

---

## 2. The parity tests are the load-bearing question

`tests/test_5_*_v2.py` (11 modules, ~3,400 lines) assert that v3 Process implementations produce numerically equivalent outputs to the legacy v1 kinetics functions (`from clearwater_modules.nsm1 import processes as v1`). They are valuable regression coverage but **depend on `clearwater_modules.nsm1` source remaining live**.

If v1 source is deleted, every parity test goes red unless they are migrated first.

---

## 3. Three options for the parity tests

### Option A — Snapshot v1 outputs to frozen fixtures (recommended)

Run each parity test under the current v1, capture the v1-side reference array, and pickle/json/npz it into a `tests/v3/fixtures/v1_reference/` directory. Rewrite the parity test to load the snapshot instead of calling `v1.<kinetic_func>` live.

**Pros:**
- Preserves the regression coverage that v3 numerically matches v1.
- v1 source can be deleted as soon as the snapshot is captured.
- Snapshots are version-controlled, reviewable.

**Cons:**
- Initial work to capture and load fixtures (~1-2 days for 11 modules).
- Snapshots are static; if v3 intentionally diverges from v1 in some Phase 9 correction, the snapshot has to be updated by hand (which is appropriate — divergence should be deliberate and documented).

### Option B — Delete parity tests outright

The Phase 1-9 v3 audit (`design/clearwater_modules_v3_nsm1_audit_*.md`) plus the existing test_5_*_v2 PRs already established parity at landing. The v3 Tier 1 / Tier 1.5 / TSM unit tests in `tests/v3/` cover all the kinetics independently.

**Pros:**
- Fast — delete 11 files, done.
- Removes the only constraint that ties v1 retirement to test migration.

**Cons:**
- Loses the v3-vs-v1 numerical-parity guard. If a future v3 change accidentally diverges from v1, only the closed-system / conservation tests would catch it (which is coarser).

### Option C — Migrate parity tests to v3-internal regression (least recommended)

Convert each `test_5_*_v2.py` to assert v3's own kinetic helpers against analytic expected values (e.g., for first-order decay, expected = `state * exp(-k * dt)`). No reference to v1 source.

**Pros:**
- Tests are conceptually cleaner (against physics, not against legacy code).
- No frozen fixtures.

**Cons:**
- Largest work effort — each test needs a hand-derived analytic expected value.
- Loses the bit-exact-match property that parity tests provide.

**Recommendation: Option A.** Mid effort, preserves regression coverage, fully decouples from v1 source.

---

## 4. Proposed sequence

### Phase 1 — clean the deck (this PR or a follow-up)

- [x] Fix the 3 immediate POM test failures from PR #17's `pom_doc_source_rate` rename
- [x] Delete the broken `test_5_tsm_calculations_v2.py` orphan (collection error, no test functions)
- [x] Write this plan

### Phase 2 — migrate parity tests (1-2 days)

- [ ] Create `tests/v3/fixtures/v1_reference/` directory
- [ ] For each `test_5_*_v2.py`: run it under v1, capture v1 reference arrays to npz, rewrite the test to load from npz
- [ ] Move the rewritten files into `tests/v3/nsm1/` with the `_v3` suffix convention
- [ ] Delete the originals

### Phase 3 — move pure-v3 tests

- [ ] Move `tests/test_phase9fc_documentation.py` and `tests/test_sanitize_rate.py` into `tests/v3/`
- [ ] Move `tests/test_1_imports.py` if kept, or delete

### Phase 4 — retire v1 source and v1-only tests

- [ ] Delete `src/clearwater_modules/` entirely (v1 src + nsm1, tsm, simlab, csm, msm, nsm2, base, sorter, shared, etc.)
- [ ] Delete `tests/test_2..17_*.py` and `tests/test_5_tsm_calculations.py`, `tests/test_hotstart_roundtrip.py`, `tests/test_tsm_latent_heat.py`, `tests/test_tsm_stability_ramp.py`
- [ ] Update `pyproject.toml` to drop v1 entrypoints / packages if any
- [ ] Update `README.md` to remove v1 references
- [ ] Tag a v3-only release

### Phase 5 — verify and announce

- [ ] Confirm `pixi run -e dev pytest tests/` is green and only references v3
- [ ] Announce v1 retirement to downstream consumers (ESM, riverine, anyone else)

---

## 5. Risk and rollback

- **Risk:** A downstream consumer (e.g. ESM, riverine, an example notebook) imports from `clearwater_modules` (v1). Phase 4 of this plan should grep the ecosystem before deleting.
- **Rollback:** v1 source is preserved in git history; if a downstream consumer surfaces, restoring is `git checkout <pre-deletion-sha> -- src/clearwater_modules/`.
- **Coordination:** The merge from `streaming` → `main` should happen at the end of Phase 4 (so `main` lands v3-only in a clean state).
