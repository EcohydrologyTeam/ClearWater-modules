# Verification and Validation Status of the Legacy Water Quality Modules and Implications for ClearWater Modules v2

**Author:** Todd E. Steissberg
**Status:** Internal working memo
**Scope:** Reviews the verification and validation (V&V) record for the ERDC-developed nutrient, temperature, contaminant, and mercury simulation modules (TSM, NSM I, NSM II, CSM, MSM/HgSM, GC) as documented in the four authoritative ERDC reports, identifies the V&V evidence that exists and the gaps that remain, and lays out how the ClearWater Modules v2 (`ClearWater-modules-streaming`) effort can build on this foundation in a way that is defensible to sponsors and reviewers.

This memo is intended as a constructive working document, not a critique of any individual or team. The legacy modules represent a substantial body of work — clear, well-documented theory, careful integration with HEC-RAS, and a foundation that the v2 effort would not be possible without.

---

## 1. Background

The legacy water quality module suite consists of:

- **TSM** — Water Temperature Simulation Module (full and simplified energy balance).
- **NSM I** — Nutrient Simulation Module I (16 state variables; simplified N, P, C, DO, algae, BOD, pathogen).
- **NSM II** — Nutrient Simulation Module II (24 state variables; multiple algal groups, full N/P/C cycles, methane, sulfides, silica, optional benthic sediment diagenesis sub-module).
- **GC/GCSM** — General Constituent Simulation Module (multi-class solids and user-defined constituents).
- **CSM** — Contaminant Simulation Module (ionization, multi-phase partitioning, degradation, photolysis, hydrolysis, volatilization, second-order transformations).
- **MSM / HgSM** — Mercury Simulation Module (Hg⁰, Hg(II), MeHg speciation, methylation, demethylation, partitioning).

These modules were developed by the U.S. Army Engineer Research and Development Center (ERDC) Environmental Laboratory and are coupled to HEC-RAS 1D and other hydrologic and hydraulic platforms as plug-in dynamic-link libraries. The kinetic formulations in NSM I and NSM II draw on a deep heritage of established surface water quality models — QUAL2E, QUAL2K, CE-QUAL-RIV1, RCA, CE-QUAL-ICM, and SedFlux — and the sediment diagenesis sub-module follows the Di Toro / Cerco lineage.

The four reports that constitute the public V&V record are:

| Report | Citation | Type |
|---|---|---|
| ERDC/EL TR-16-1 | Zhang & Johnson (2016), *Aquatic Nutrient Simulation Modules (NSMs) Developed for Hydrologic and Hydraulic Models* | Theory and formulations |
| ERDC/EL TR-16-8 | Zhang & Johnson (2016), *Aquatic Contaminant and Mercury Simulation Modules Developed for Hydrologic and Hydraulic Models* | Theory and formulations |
| ERDC/EL TR-16-11 | Johnson & Zhang (2016), *Testing and Validation Studies of the NSMII-Benthic Sediment Diagenesis Module* | Code intercomparison and sensitivity analysis |
| ERDC TN-EMRRP-SR-47 | Zhang & Johnson (2014), *Application and Evaluation of the HEC-RAS-Nutrient Simulation Module (NSM I)* | Field application and statistical evaluation |

ClearWater Modules v2 is a Python reimplementation of these modules, organized for streaming and reactive-transport coupling and intended to enable a level of unit testing, regression testing, and analytical benchmarking that the original DLL-based delivery did not require.

---

## 2. What the Reports Establish

### 2.1 Theory and formulation documentation (TR-16-1 and TR-16-8)

Both 2016 theory reports are substantial technical documents that fully derive the governing equations, parameter definitions, default rate coefficients, temperature-correction relationships, stoichiometric ratios, and mass-balance closure for every state variable in the suite. They tie each formulation back to the relevant primary literature (QUAL2K, ICM, SedFlux, Di Toro, Chapra, and others) and provide the mathematical record needed for an independent reimplementation. TR-16-8 acknowledges expert review of the formulations by K. Farley (Manhattan College), J. Martin (Mississippi State), J. Wang (UC Davis), and M. Dortch (ERDC).

This formulation record is the foundation that makes a v2 port feasible at all and is one of the strongest aspects of the legacy delivery.

### 2.2 NSM II benthic sediment diagenesis (TR-16-11)

TR-16-11 documents a structured testing program for the sediment diagenesis sub-module that includes:

- **Three-way code intercomparison** of NSM II against SedFlux and CE-QUAL-ICM, using a consistent set of kinetic coefficients, initial conditions, and water column depositional fluxes drawn from a calibrated ICM application to Chesapeake Bay (1991–2000).
- **Steady-state vs. unsteady-state internal consistency checks** for sediment organic matter (POC, PON, POP) across all three reactivity classes (G1, G2, G3).
- **Analytical vs. numerical comparison** for sediment methane under saturating conditions.
- **Sensitivity analyses** for water-column depositional rates, deep burial velocities, particle mixing coefficients, and benthic stress.

The reported agreement among the three codes for sediment diagenesis fluxes (carbon, nitrogen, phosphorus), particle mixing coefficient, diffusion coefficient, and sediment-water transfer coefficient is good. The authors clearly note in the report that a full validation against contemporary field measurements of all sediment state variables was outside the scope of the program. Within that explicit scope, the testing demonstrates that the NSM II implementation produces results consistent with two well-established peer-reviewed sediment diagenesis frameworks under shared forcing.

### 2.3 NSM I Lower Minnesota River application (TN-EMRRP-SR-47)

TN-EMRRP-SR-47 documents the application of the integrated HEC-RAS / NSM I model to a 90-segment representation of the Lower Minnesota River (LMNR) over the 2001–2006 period, with calibration and verification against four MCES long-term monitoring stations. Coefficients of determination (R²), Nash–Sutcliffe efficiency (NSE), and percent bias (PBIAS) were reported for hydraulics, temperature, total dissolved solids, inorganic suspended solids, organic and inorganic N and P species, algal biomass, dissolved oxygen, and CBOD.

Performance was strong across most constituents and stations. Hydraulics and water temperature were reproduced with R² and NSE values approaching 0.97–0.98 and PBIAS within a few percent. Most nutrient species and DO had PBIAS values within ±10% at most stations. Some constituents — notably CBOD at the most downstream station and inorganic suspended solids — had weaker fits, which the authors attribute to known structural simplifications (CBOD as a lumped first-order variable, the absence of a true sediment resuspension mechanism in NSM I as opposed to the dedicated HEC-RAS sediment transport module). These are honestly identified in the report and are appropriate scoping choices for a simplified-kinetics module.

### 2.4 Summary table of the public V&V record

| Module | Formulation documented | Code-to-code testing | Field application reported |
|---|---|---|---|
| TSM | Yes (TR-16-1 §2) | — | Exercised inside the LMNR NSM I run; not as a standalone study |
| NSM I | Yes (TR-16-1 §3) | — | Yes — LMNR, 4 stations, 6 years (TN-EMRRP-SR-47) |
| NSM II — water column | Yes (TR-16-1 §4) | — | Application study deferred to future work in TN-EMRRP-SR-47 |
| NSM II — benthic sediment diagenesis | Yes (TR-16-1 §5) | Yes — vs. SedFlux and ICM (TR-16-11) | Driven by Chesapeake Bay forcing through the ICM parameter set; not a standalone field study |
| GC | Yes (TR-16-8 §2) | — | Application studies not reported in this set |
| CSM | Yes (TR-16-8 §3) | — | Application studies not reported in this set |
| MSM / HgSM | Yes (TR-16-8 §4) | — | Application studies not reported in this set |

---

## 3. What ClearWater Modules v2 Inherits, and What It Does Not

The kinetic structures embedded in NSM I and NSM II are not novel inventions — they are careful recompositions of formulations with decades of literature support. The same is true of the sediment diagenesis lineage (Di Toro–Cerco–SedFlux), the temperature energy-balance approach in TSM (Deas, Lowney, Chapra), the contaminant kinetics in CSM (standard ionization, sorption, and degradation forms), and the mercury speciation framework in MSM (well-established Hg cycling). The v2 effort therefore inherits a kinetic structure that is well-grounded in the peer-reviewed literature.

What the v2 effort does *not* inherit is a comprehensive multi-site field validation record or a programmatic test suite (unit tests, analytical-solution benchmarks, regression tests, continuous integration). The legacy delivery format — compiled DLLs documented through technical reports — was not structured around those artifacts, and producing them was outside the scope of the funded program. This is a normal characteristic of the era and class of model delivery and is not unique to this lineage.

A clean-room reimplementation accordingly produces a code base whose:

- **Kinetic correctness against the literature** can be argued from the formulation reports.
- **Implementation parity with the legacy Fortran** can be argued from v1-vs-v2 comparison tests (the `tests/v1_parity` family on the streaming branch).
- **Predictive skill against observational data** must be re-established through new application studies, since only NSM I has a published field application and that application is single-site.
- **Numerical correctness** must be re-established through analytical and method-of-manufactured-solutions (MMS) tests, since only the sediment-methane analytical comparison in TR-16-11 covers this category.

---

## 4. Recommended Framing for Sponsors, Reviewers, and Stakeholders

A defensible characterization of ClearWater Modules v2 is the following:

> The v2 modules are a clean-room Python reimplementation of the ERDC nutrient, temperature, contaminant, and mercury simulation modules. The kinetic formulations follow the same well-documented mathematical framework given in the ERDC technical reports and rest on decades of peer-reviewed surface water quality modeling literature (QUAL2K, ICM, SedFlux, Di Toro, Chapra, and related sources). The reimplementation is paired with a programmatic testing infrastructure — unit tests for individual kinetic terms, analytical and method-of-manufactured-solutions benchmarks, parity tests against the v1 Python translations of the original Fortran kernels, and conservation-of-mass checks — that complements the formulation documentation and field application work already in the public record. The v2 effort is positioned to extend the published evidence base by enabling new site applications and reactive-transport benchmarks (including PhreeqcRM coupling) that were not feasible under the original DLL-based delivery format.

This framing is accurate, transparent, and respectful of the legacy contribution. It does not overstate the inherited V&V record (e.g., it does not claim that NSM II water column or MSM, CSM, GC have published field validation), and it does not understate the kinetic and theoretical foundation, which is genuinely strong.

What to avoid in stakeholder communications:

- Claiming that the v2 modules are "validated" because the legacy modules are. The legacy formulations are documented and partially demonstrated; that is not the same as validated against multi-site field data.
- Conversely, framing the legacy modules as inadequate or untested. The published record establishes formulation correctness, code consistency with peer models for sediment diagenesis, and a credible single-site application of NSM I. That is meaningful.
- Suggesting that closing remaining V&V gaps falls solely on the v2 effort. Field validation of an 1D water quality model is properly a community activity and benefits from collaboration with HEC, ERDC, district offices, and external users.

---

## 5. Recommended V&V Additions in v2

The most productive additions, ordered roughly by leverage:

1. **Unit tests for individual kinetic terms.** These are inexpensive, give regression coverage, and document expected behavior at the rate-equation level. The TSM thin-water stability test and v1-parity tests on the streaming branch are good models to extend to all kinetic source/sink terms in NSM I, NSM II, MSM, and CSM.

2. **Analytical and MMS benchmarks.** Examples include single-cell decay tests with known half-life, analytical advection–reaction solutions for first-order constituents, the sediment-methane analytical solution noted in TR-16-11, and conservation-of-mass tests over closed-system simulations. These provide numerical correctness evidence that is not strongly represented in the published record.

3. **A new NSM II water-column field application.** This is the single most valuable addition to the published evidence base, because the NSM II water column is the most ambitious of the kinetic modules and currently has only formulation documentation behind it. Suitable candidate sites are systems with co-located W2 or RAS-WQ models, multi-year monitoring, and active stakeholder interest.

4. **Field benchmarks for CSM and MSM.** Even single-site demonstration studies would substantially strengthen the published record for these modules. Mercury and contaminant applications often have rich monitoring data attached to regulatory drivers, which makes them well-suited to demonstration cases.

5. **Re-execution of the NSM II sediment diagenesis intercomparison from TR-16-11.** Reproducing the SedFlux / ICM / NSM II comparison with the v2 Python implementation would provide a direct continuity check between the legacy and v2 lineages and would surface any implementation drift introduced by the port. The known divergence between analytical and numerical sediment methane at saturation is worth re-checking specifically.

6. **Reactive-transport benchmarks via PhreeqcRM.** These tests are available to v2 in a way that they were not available to the legacy modules and offer a fundamentally new line of evidence for the speciation, partitioning, and equilibrium components of NSM II, MSM, and CSM.

---

## 6. Closing Note

The picture painted by the four reports is of a substantial theoretical and integration effort delivered under the constraints of an applied research program: comprehensive formulation documentation, careful HEC-RAS integration, code-level consistency with established peer models for sediment diagenesis, and one demonstrated field application for the simplest of the kinetic modules. The remaining V&V gaps — multi-site field validation for NSM II, demonstration applications for CSM and MSM, a programmatic numerical test suite — were a function of program scope rather than oversight, and they are exactly the gaps that the v2 streaming effort is structured to close.

Treating the v2 work as the next chapter in a continuing lineage, rather than as a replacement, is both accurate and the strongest position from which to engage with HEC, ERDC, and external users on what comes next.
