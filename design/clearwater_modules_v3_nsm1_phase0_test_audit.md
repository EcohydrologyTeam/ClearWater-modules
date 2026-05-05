# NSM1 v3 Phase 0 Test Audit: Inventory of v1 Calculation Tests & Expected-Value Fixtures

**Phase:** 0.3 — Test fixture audit  
**Date:** 2026-05-04  
**Scope:** v1 NSM1 calculation tests, v2 NSM1 parity tests, and expected-value tables in `tests/NSM Manual Calcs/`  
**Purpose:** Complete inventory to support Phase 7 test plan mapping onto v3's Process-based dispatch

---

## 1. v1 NSM1 Calculation Test Files (tests/test_7–test_17)

All v1 tests follow a common structure:
- **Framework:** Import and instantiate full `NutrientBudget` v1 Model
- **Fixture setup:** Pytest fixtures for initial state, all parameter groups (Algae, Alkalinity, Balgae, Nitrogen, Carbon, CBOD, DOX, N2, POM, Pathogen, Phosphorus, GlobalParameters, GlobalVars)
- **Test pattern:** Modify one or two parameters from defaults, run 1 time step, extract final state value, assert against hard-coded expected value
- **Expected values:** Hard-coded in assertion statements (`assert pytest.approx(result, tolerance) == expected_value`)
- **Tolerance:** Two regimes: 0.000001 (high precision: Algae, Balgae, Nitrogen, PX, POM, N2) and 0.01 (relaxed: Carbon, DOX, Alkalinity, CBOD)

### Per-Test-File Inventory

| Test File | Constituent(s) | Test Count | Expected Values | Tolerance | Tolerance Regime | Notes |
|---|---|---|---|---|---|---|
| test_7_nsm_algae_calculations.py | Ap (floating algae) | 32 | Hard-coded in test assertions (31 expected values) | 0.000001 | High precision | Tests growth rate options, light limitation, nutrient limitation (KsN, KsP), settling, PAR, POC effects. All single-parameter perturbations from defaults. |
| test_8_nsm_balgae_calculations.py | Ab (benthic algae) | 33 | Hard-coded in test assertions | 0.000001 | High precision | Parallel to test_7 but for benthic. Tests density limitation (FSb, Ksb), respiration/death rates, hydraulic factors. |
| test_9_nsm_nitrogen_calculations.py | NH4, NO3, OrgN | 27 | Hard-coded in test assertions (81 expected values across 27 tests) | 0.000001 | High precision | Tests nitrification, denitrification, NH4 from sediment, nitrification inhibition by DOX, OrgN hydrolysis. |
| test_10_nsm_carbon_calculations.py | POC, DOC, DIC | 73 | Hard-coded in test assertions (215 expected values) | 0.01 | Relaxed | Largest v1 test file. Tests POC/DOC hydrolysis, DIC reaeration, CO2 equilibrium, light-dependent processes, settling. |
| test_11_nsm_DOX_calculations.py | DOX | 61 | Hard-coded in test assertions (61 expected values) | 0.01 | Relaxed | Tests reaeration (9 hydraulic + wind options), nitrification/denitrification O2 demand, respiration, SOD, reaeration temperature correction. |
| test_12_nsm_alkalinity_calculations.py | Alk | 52 | Hard-coded in test assertions (52 expected values) | 0.01 | Relaxed | Tests nitrification/denitrification stoichiometry, algae/CBOD effects on alkalinity. |
| test_13_nsm_CBOD_calculations.py | CBOD | 10 | Hard-coded in test assertions (10 expected values) | 0.01 | Relaxed | Small file. Tests CBOD oxidation kinetics, DOX inhibition, temperature correction. |
| test_14_nsm_phosphrous_calculations.py | TIP, OrgP | 15 | Hard-coded in test assertions (30 expected values) | 0.000001 | High precision | Note: filename typo (phosphrous not phosphorus). Tests OrgP hydrolysis, TIP settling, sediment flux. |
| test_15_nsm_POM_calculations.py | POM | 16 | Hard-coded in test assertions (16 expected values) | 0.000001 | High precision | Tests POM hydrolysis/dissolution, settling, sediment delivery. |
| test_16_nsm_PX_calculations.py | PX (pathogen) | 9 | Hard-coded in test assertions (9 expected values) | 0.000001 | High precision | Small file. Tests pathogen decay kinetics, temperature correction. |
| test_17_nsm_N2_calculations.py | N2 | 30 | Hard-coded in test assertions (30 expected values) | 0.000001 | High precision | Tests N2 production from denitrification, TDG (total dissolved gas) interactions. |

**Totals for v1 NSM1:**
- **358 test functions** across 11 files
- **485 hard-coded expected values**
- **Framework:** Full NutrientBudget model instantiation; not isolated process functions
- **State under test:** Full 16-constituent coupled system with one or two parameters varied per test
- **Dispatch:** Single time step (1 day default)

### Test Structure Details

#### Fixture Architecture

Every v1 test file includes:

```python
@pytest.fixture(scope='function')
def initial_nsm1_state() -> dict[str, float]:
    """Return initial state values for the model."""
    return {
        'Ap': 40, 'Ab': 24, 'NH4': 0.05, 'NO3': 5, 'OrgN': 1.726, 'N2': 1,
        'TIP': 0.07, 'OrgP': 0.24, 'POC': 4, 'DOC': 1, 'DIC': 1, 
        'POM': 10, 'CBOD': 5, 'DOX': 8, 'PX': 1, 'Alk': 1
    }

@pytest.fixture(scope='function')
def default_[constituent]_params() -> [TypedDict]:
    """Returns default static variable values."""
    return [Constituent]StaticVariables(...)
    # 13 such fixtures covering: Algae, Alkalinity, Balgae, Nitrogen, Carbon,
    # CBOD, DOX, N2, POM, Pathogen, Phosphorus, GlobalParameters, GlobalVars

def get_nutrient_budget_instance(...) -> NutrientBudget:
    """Factory helper; constructs full model from fixtures."""
    return NutrientBudget(
        time_steps=time_steps,
        initial_state_values=initial_nsm1_state,
        algae_parameters=default_algae_params,
        # ... all 13 parameter groups
        time_dim='nsm1_time_step',
    )
```

#### Test Function Pattern

```python
def test_changed_[PARAMETER](
    time_steps, initial_nsm1_state, default_algae_params, ..., tolerance
) -> None:
    """Test the model with [PARAMETER] changed from default."""
    
    # 1. Modify one or two parameters
    initial_state_dict['NH4'] = 0.5  # or
    default_algae_params['KL'] = 20
    
    # 2. Instantiate model
    nsm1: NutrientBudget = get_nutrient_budget_instance(...)
    
    # 3. Run single time step
    nsm1.increment_timestep()
    
    # 4. Extract final state
    Ap = nsm1.dataset.isel(nsm1_time_step=-1).Ap.values.item()
    
    # 5. Assert hard-coded expected value
    assert isinstance(Ap, float)
    assert pytest.approx(Ap, tolerance) == 52.668069
```

#### Tolerance Semantics

- `pytest.approx(actual, tolerance)` accepts values within `abs(actual - expected) <= tolerance * max(abs(actual), abs(expected))`
- **0.000001:** Algae/Balgae/Nitrogen/PX/POM/N2 — expect ~1 ppm relative error
- **0.01:** Carbon/DOX/Alkalinity/CBOD — expect ~1% relative error

### No External Data Files in v1 Tests

**Key finding:** v1 test files do NOT import or read from `tests/NSM Manual Calcs/`. All expected values are hard-coded inside test assertions. The Excel workbooks are present but unlinked to the test suite.

---

## 2. tests/NSM Manual Calcs/ Directory — Framework-Independent Expected Values

Directory location: `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/tests/NSM Manual Calcs/`

### File Inventory

| Filename | Size | Constituent(s) | Format | Sheet Count | Status |
|---|---|---|---|---|---|
| Manual Calculations for Alkalinity.xlsx | 36 KB | Alk | Excel 2007+ | Unknown (not extracted) | Present, not linked to tests |
| Manual Calculations for Carbon.xlsx | 61 KB | POC, DOC, DIC | Excel 2007+ | Unknown | Present, not linked to tests |
| Manual Calculations for CBOD.xlsx | 11 KB | CBOD | Excel 2007+ | Unknown | Present, not linked to tests |
| Manual Calculations for DOX.xlsx | 47 KB | DOX | Excel 2007+ | Unknown | Present, not linked to tests |

**Coverage:** 4 of 11 constituents have manual-calc files (Alk, Carbon, CBOD, DOX).

**Missing workbooks:** No files for Ap, Ab, NH4, NO3, OrgN, TIP, OrgP, POM, PX, N2.

### Excel Format Assessment

- **Format:** Microsoft Excel 2007+ (`.xlsx` binary format)
- **Extractability:** Standard; readable with `openpyxl`, `pandas.read_excel()`, or equivalent
- **Current v3 reuse:** Would require ETL to convert to CSV or structured test-data format; not currently compatible with pytest parameterization without intermediate tooling
- **Likely content:** Expected-value tables with parameter combinations and computed outputs for validation

---

## 3. v2 NSM1 Parity Tests

Three test files validate v2 sub-rate methods against equivalent v1 helper functions (not full model parity). These tests are **narrow in scope and high in precision** — they test individual kinetic sub-routines in isolation, not coupled whole-model behavior.

### Per-Test-File Inventory

| Test File | Constituent Process | Test Count | Expected Values | Test Pattern | Notes |
|---|---|---|---|---|---|
| test_5_nitrogen_calculations_v2.py | Nitrogen (NH4, NO3, OrgN) | 7 | Hard-coded assertions via `np.testing.assert_allclose(rtol=1e-6)` | Synthetic 5-cell DataArrays; compare v2 sub-rate to v1 helper-function output | Tests: nitrification_inhibition, ammonium_nitrification, ammonium_from_bed, nitrate_bed_denitrification, nitrate_denitrification, change_ammonium_no_algae, change_nitrate_no_algae. Includes NaN-guard verification (bug fix #4). |
| test_5_floating_algae_calculations_v2.py | FloatingAlgae (Ap) | 6 | Hard-coded assertions via `np.testing.assert_allclose(rtol=1e-6)` | 5-cell synthetic xarray DataArrays; compare v2 method to v1 helper | Tests: rate_respiration, rate_death, rate_settling, limit_phosphorus, limit_nitrogen, rate_growth (multiplicative). Includes in-place mutation bug (v2 limit_nitrogen aliases input). |
| test_5_benthic_algae_calculations_v2.py | BenthicAlgae (Ab) | 4 | Hard-coded assertions via `np.testing.assert_allclose(rtol=1e-6)` | 5-cell synthetic xarray DataArrays; compare v2 to v1 | Tests: rate_respiration, rate_death, limit_density, rate_growth. |

**Totals for v2 NSM1:**
- **17 test functions**
- **Hard-coded tolerance:** rtol=1e-6 (1 ppm) for all; no relative or absolute tolerance loosening
- **Framework:** xarray DataArray arrays (5-cell synthetic data); v2 Process instances; v1 helper-function calls for reference
- **Scope:** Sub-rate kinetic validation; NOT full model runs
- **Known issues tested:** Nitrification inhibition NaN check (bug #4), in-place mutation in limit_nitrogen (bug not yet fixed)

### Why v2 Tests Matter for v3

v2 tests lock in the intended behavior of 4 constituents (NH4/NO3 in Nitrogen, Ap in FloatingAlgae, Ab in BenthicAlgae) despite known bugs in v2's full integrator. They demonstrate that:
- v1 helper functions produce correct sub-rate values
- v2 Process methods *call* those helpers correctly
- v2's integrator bugs (multiplicative instead of additive) are orthogonal to the kinetics
- v3 must reproduce these sub-rate results exactly (modulo floating-point) to claim correctness

---

## 4. Assessment: Porting Fixtures to v3

### What Reuses Cleanly

1. **Hard-coded expected values (485 v1 + 17 v2):** Can be migrated directly to v3 tests as long as matching Process instantiation and parameter passing is possible. v3 tests need not change tolerance values if v3's numerical behavior is identical to v1.

2. **Parameter fixtures:** All 13 parameter groups (AlgaeStaticVariables, etc.) in v1 tests have clear v3 equivalents via `Process.DEFAULTS` dicts + YAML overrides. Direct copy-paste of fixture definitions is feasible.

3. **Initial state dict:** The 16-constituent initial state dict transfers unchanged; v3's registry will read/write the same variable names (Ap, Ab, NH4, etc.).

4. **Test structure:** Individual single-parameter perturbations → isolated expected values is a clean pattern for v3. Each test can become a per-Process test:
   - `test_floating_algae_growth_rate_variation` (moved to FloatingAlgae process test file)
   - `test_nitrogen_nitrification_rate_variation` (moved to Nitrogen process test file)
   - ...and so on for all 11 processes.

### What Needs Adaptation

1. **Framework target:** v1 tests run the full NutrientBudget coupled model (16 constituents, all rates coupled in a single step). v3 tests must either:
   - Run the equivalent v3 Model with all 11 Processes declared and dispatched (integration tests), OR
   - Run individual Process.run() calls in isolation with mocked registry reads/writes (unit tests).
   
   **Recommendation:** Both. Integration tests reuse expected values for the coupled behavior (requires v3 Model functional); unit tests reuse sub-rate comparisons from v2 parity suite.

2. **Expected values for new constituents:** 8 v1 test files have no Excel workbooks (Ap, Ab, NH4, NO3, OrgN, TIP, OrgP, N2, PX, POM—10 total). Expected values exist only as hard-coded assertions. v3 can reuse these, but new tests for Carbon/DOX/Alkalinity integration must decide whether to:
   - Port the hard-coded values as-is (assumes v3 duplicates v1 kinetics exactly), OR
   - Regenerate expected values from v1 via standalone calculations, then update hard-coded values.

3. **Excel workbooks (4 files, 4 constituents):** Format is not compatible with pytest parameterization without preprocessing. Options:
   - Extract to CSV and write parameterized pytest fixtures (Phase 1 work)
   - Leave as documentation; reuse hard-coded values only
   - Implement Excel→Python converter if they contain richer expected-value tables (TBD on workbook content)

4. **Test execution model:** v1 tests run `nsm1.increment_timestep()` and extract xarray dataset results. v3 tests must adapt to v3's registry-based model. Likely patterns:
   - `registry.get_at_time(t_next, constituent_name)` for state reads
   - `registry.get_rate_variable(rate_name)` for rate variable validation
   - Equivalence: "does the Process dispatch produce the expected state transition?"

### What Is Missing

1. **Isolated sub-rate test coverage for 8 constituents:** v2 parity tests exist only for Nitrogen, FloatingAlgae, BenthicAlgae. New v3 tests need equivalent sub-rate coverage for:
   - Phosphorus (TIP, OrgP)
   - Carbon (POC, DOC, DIC)
   - POM
   - CBOD
   - DOX
   - Pathogen (PX)
   - N2
   - Alkalinity

   These can be built in Phase 7 by following the v2 parity pattern (5-cell arrays, v1 helper-function reference).

2. **Integration test coverage for v3 coupled behavior:** No v1 tests are structured as "run all 11 processes, check closure for a closed system" or "vary one parameter across all dependent processes." Such tests would require Phase 7 design.

3. **Regression test coverage for reaeration menu (22 options):** test_11_nsm_DOX_calculations.py has ~61 tests, but coverage of all 9 hydraulic + 13 wind reaeration options is incomplete. v3 can inherit and extend this via Phase 1's shared `reaeration.py` utility, but v1 tests may not exercise all 22 combinations.

---

## 5. Recommendations for Phase 7 Test Plan

1. **Tier 1 (reuse as-is):** All 485 v1 hard-coded expected values → v3 integration tests, one test per constituent Process. Tolerance values unchanged (0.000001 or 0.01 per regime).

2. **Tier 2 (new, modeled on v2 parity):** Create isolated sub-rate tests for the 8 missing constituents (Phosphorus, Carbon, POM, CBOD, DOX, Pathogen, N2, Alkalinity) using 5-cell synthetic arrays and v1 helper-function references. Tolerance: rtol=1e-6.

3. **Tier 3 (design-dependent):** Whole-system conservation and coupling tests (TBD on Phase 7 spec). These will NOT reuse v1 expected values; they verify new v3 invariants.

4. **Excel workbooks:** Extract to CSV during Phase 1 or early Phase 7. Decide whether to parameterize (make broad Coverage) or document-only (save effort). Recommend document-only for v3.0.0; parameterization deferred to 1.1.

5. **Tolerance tightening:** v3's Jacobi semantics for state reads (vs v2/v1's mix of Jacobi and Gauss-Seidel) may require small tolerance loosening (e.g., 0.000001 → 0.00001) for a few coupled multi-process tests. Measure empirically in Phase 7; do not assume v1 tolerance applies.

---

## 6. Summary Table: Audit Findings

| Category | Finding | Implication for v3 |
|---|---|---|
| **v1 test count** | 358 tests, 485 hard-coded expected values | Reusable as integration test expected values |
| **v2 test count** | 17 tests, sub-rate kinetics only | Lock in v3 Process.rate_* method behavior |
| **Manual Calc files** | 4 Excel workbooks (Alk, Carbon, CBOD, DOX); 4 constituents (Ap, Ab, NH4/NO3, OrgN/TIP/OrgP, N2, PX) missing | Extract to CSV (Phase 1); leave as documentation otherwise |
| **Expected value source** | All hard-coded in assertions; no programmatic link to Excel | v3 must decide: regenerate or reuse hard-coded values |
| **Test framework** | v1 tests run full NutrientBudget; v2 tests run sub-rate methods | v3 must support both integration and unit test execution |
| **Parameter fixtures** | All 13 groups present in v1 test files; directly portable | Copy fixture definitions; adapt to v3 Process.DEFAULTS pattern |
| **Missing sub-rate tests** | 8 constituents lack parity-style isolated kinetics tests | Phase 7 task: create 8 new sub-rate test files (modeled on v2) |
| **Tolerance regimes** | 0.000001 (11 tests) and 0.01 (3 tests) | Measure actual v3/v1 differences; tighten/loosen as needed |

---

## Appendix A: Excel Workbook Content (Not Extracted)

The four Excel workbooks in `tests/NSM Manual Calcs/` are present but their sheet structure, cell layout, and test-case organization have not been examined. Recommendations:

1. **Phase 1 discovery:** Assign one team member to open each workbook and document:
   - Sheet names and count
   - Column headers (parameter names, expected outputs)
   - Row count (number of test cases per workbook)
   - Whether cell values are formulae (recalculate on parameter change) or static (one-time calculations)

2. **Phase 1 decision:** Given the discovery, decide whether to:
   - **Option A (low effort):** Leave workbooks as documentation; test suite uses hard-coded values only.
   - **Option B (medium effort):** Extract to CSV and write parameterized pytest fixtures for v3 tests.
   - **Option C (high effort):** Implement VBA-to-Python conversion (if workbooks have custom logic) or link workbooks programmatically (if they compute live).

3. **Recommendation:** Option A for v3.0.0. Workbooks serve as human-readable kinetics validation; test suite reuses hard-coded v1 expected values. Parameterization from Excel deferred to v3.1.0 if needed.

---

**End of Phase 0.3 Audit Report**
