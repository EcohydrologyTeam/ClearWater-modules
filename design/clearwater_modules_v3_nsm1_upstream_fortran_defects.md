# Upstream NSM1 Fortran Defects — Report of Record

**Purpose:** track scientific defects found in the upstream NSM1 Fortran
source that ClearWater v3 (Python) **deliberately diverges from**, with
the reference anchor for each correction. These are to be filed with the
NSM1 Fortran maintainers. A defect listed here is *not* a v3 bug — v3 is
correct by reference; the divergence is intentional and documented at the
point of use.

**Status:** open report. Filing to upstream maintainers is an external
action tracked separately.

---

## UF-1 · `r_alkden` denitrification alkalinity coefficient is 4× too high

- **Where:** `modAlkalinity.f90:54` — `ralkden = 4.0/14.0/1000.0` eq/mg-N.
- **Propagation:** inherited verbatim by ClearWater v1
  (`clearwater_modules/nsm1/constants.py`) and, pre-fix, v3. The error is
  shared **Fortran = v1 = v3**, so it is invisible to any Fortran↔v1 or
  v1↔v3 regression/parity check (the canonical "wrong at all stages"
  case).
- **Correct value:** `ralkden = 1.0/14.0/1000.0` eq/mg-N.
- **Basis:** denitrification reduces NO₃⁻ to N₂ and produces **1
  equivalent of alkalinity per mole of NO₃-N reduced**:
  - CE-QUAL-W2 `water-quality.f90:3157` routes denitrification alkalinity
    at 1 eq per mol N.
  - Stumm & Morgan, *Aquatic Chemistry*, 3rd ed., Ch. 4 — alkalinity
    changes from N redox: nitrification consumes 2 eq/mol-N,
    denitrification produces 1 eq/mol-N.
  - Internal cross-check: NSM1's own `ralkn` (nitrification) =
    `2/14/1000` eq/mg-N is the textbook 2 eq/mol-N on the identical
    basis, confirming the denitrification term should be `1/14/1000`,
    not `4/14/1000`. There is no stoichiometric pathway that yields a
    factor of 4.
- **Magnitude of effect:** every alkalinity budget with active
  water-column denitrification overstates the denitrification alkalinity
  *source* by exactly 4×. Dormant only where denitrification is inactive
  (anoxic-zone / high-NO₃ systems are most affected).
- **v3 disposition:** corrected to `1/14/1000` in
  `src/clearwater_modules_v3/parameters/alkalinity.py` with an inline
  reference-anchored divergence note (NSM1-SCI-N1, gold-standard spec
  A2). Audit entry corrected in
  `clearwater_modules_v3_nsm1_audit_simple_constituents.md` (Alkalinity
  parameter table, `r_alkden` row). Regression-guarded by a
  closed-system test asserting Δalkalinity = 1 eq per mol-N denitrified
  (`tests/v3/nsm1/test_alkalinity_scin1_regression.py`).
- **Recommended upstream fix:** change `modAlkalinity.f90:54` to
  `ralkden = 1.0/14.0/1000.0` and add a citation comment
  (CE-QUAL-W2 `water-quality.f90:3157`; Stumm & Morgan 3rd ed., Ch. 4).
