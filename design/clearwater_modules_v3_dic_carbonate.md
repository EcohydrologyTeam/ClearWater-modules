# Design Spec: DIC carbonate chemistry and reactive carbon (NSM1-I port completion)

**Date:** 2026-05-30
**Component:** `src/clearwater_modules_v3/processes/carbon.py`,
`src/clearwater_modules_v3/processes/alkalinity.py`
**Severity:** Incomplete Fortran -> Python port. DIC and pH are currently inert.
**Status:** RESOLVED 2026-05-30. Findings on review:

1. **Reactive DIC (`carbon.py`) was already complete** — DOC/CBOD oxidation,
   algal/benthic photosynthesis-respiration exchange, and CO2 reaeration are all
   implemented. The CO2 reaeration term uses a constant free-CO2 fraction
   `FCO2 = 0.2`, which **matches NSM1-I** (the Fortran `nsmi_carbon.f90` likewise
   uses a static `Fco2 = 0.2` parameter, not a pH-derived value). So the spec's
   "DIC is inert / simple-tracer" premise did not hold for current v3.
2. **pH is a diagnostic, not a coupled state.** Because `Fco2` is constant, pH
   feeds no kinetics (DIC, alkalinity, or otherwise), so computing it does not
   change the state trajectory — no re-baseline.
3. **The carbonate solver already existed on the `nsm2-and-hab` line**
   (`utils/carbonate.py`, NSM2 step S4-3). It was adopted here byte-for-byte
   (Newton -> bisection[3,13] -> hold-previous graceful failure; Davies
   ionic-strength correction that vanishes at I=0 -> exact NSM1-I constants).
   `Alkalinity.run` computes pH from (alkalinity, dic, water_temperature) at
   freshwater I=0 and writes it only when `"pH"` is pre-registered (opportunistic
   diagnostic). The `DICfromBed` sediment term remains NSM-II / `use_SedFlux`.

Tests: `tests/v3/nsm1/test_carbonate_ph_v3.py` (shared-solver tests are
byte-consistent with the `nsm2-and-hab` line; plus this branch's `Alkalinity`
pH-diagnostic wiring). Provenance: `design/nsm2_alkalinity_ph_fortran_extraction.md`.

---

**Original proposal:** Proposed. Confirmed against the Fortran NSM1-I source.

## Summary

v3 carries DIC as a "simple-tracer placeholder until the carbonate solver lands"
(`carbon.py:66`). The Fortran NSM1-I implements a **reactive** dissolved
inorganic carbon pool plus a carbonate-equilibrium pH solver. This is a port
gap: DIC, pH, and the CO2 exchange are part of NSM1-I water-column chemistry,
not a future release.

## Fortran NSM1-I reference

DIC source/sink terms (`fortran/NSM1/03_dissolved_oxygen/nsmi_carbon.f90`):

- `DIC_Reaeration` — atmospheric CO2 exchange.
- `DIC_ApGrowth` / `ApRespiration_DIC` — DIC removed by photosynthesis, returned
  by algal respiration (and the benthic-algae equivalents `DIC_AbGrowth` /
  `AbRespiration_DIC`).
- `DOC_DIC_Oxidation` and `CBOD_DIC_Oxidation` — DIC produced by oxidation of
  dissolved organic carbon and CBOD.
- `DICfromBed` — sediment DIC flux (only with `use_SedFlux`; that path is NSM-II
  and remains future work).

Carbonate-equilibrium pH solver (`fortran/NSM1/05_additional_variables/nsmi_alkalinity.f90`),
Newton iteration on:

```
f(pH) = (K1*h + 2*K1*K2)/(h^2 + K1*h + K1*K2) * DIC + Kw/h - h - Alk/50000
```

with temperature-dependent K1, K2, Kw. Alkalinity and DIC are coupled through pH.

## Required change

1. Implement the DIC reactive terms in `carbon.py`: CO2 reaeration, the
   algal/benthic photosynthesis-respiration DIC exchange, and DIC production from
   DOC and CBOD oxidation, matching the Fortran term-by-term. Exclude the
   `DICfromBed` sediment term (NSM-II / `use_SedFlux`, out of scope).
2. Implement the carbonate pH solver in `alkalinity.py` (Newton iteration over
   the carbonate system) so pH, DIC, and alkalinity are coupled as in NSM1-I.
3. Default the carbonate constants (K1, K2, Kw temperature dependence, the
   `Alk/50000` equivalence) to the NSM1-I values; verify against the Fortran.
4. Keep a switch to run DIC as an inert tracer for didactic/benchmark runs, but
   make the reactive carbonate behavior the default.

## Verification

- Unit: pH solver reproduces the Fortran pH for known (DIC, Alk, T) triples;
  each DIC term matches the Fortran for known inputs.
- Conservation: carbon closes across DOC/POC/DIC/CBOD exchanges (no spurious
  source/sink) in a no-flux box test.
- Integration: the coupled Willowbend run shows DIC responding to photosynthesis,
  respiration, and oxidation rather than holding its tracer value.

## Report follow-up (Report 2 repo)

After the fix: re-run the demonstration (DIC, pH, and alkalinity become dynamic)
and document the carbonate system and reactive DIC in the Methods carbon /
alkalinity sections, replacing the "transported tracer" description. The sediment
DIC flux remains explicitly future work (NSM-II). Tracked as item A4.
