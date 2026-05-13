# Appendix A Diff — NSM1 1.0.0 vs Pattern-Alignment Spec §4

**Date:** 2026-05-13
**Purpose:** Phase 0.5 of `clearwater_modules_v3_nsm1_pattern_alignment_specification.md`. Diff the per-Process component inventory in the pattern-alignment spec §4 against the existing `clearwater_modules_v3_nsm1_design_specification.md` Appendix A registry-coupling cheatsheet. The amended Appendix A is in this file's §3 below and is applied to the design spec in a follow-up commit.

---

## 1. The original Appendix A scope

The 1.0.0 Appendix A is a **coupling cheatsheet** — it enumerates only the rate variables that one Process writes to the registry for another Process to read. It is scoped to *inter-process communication*, not to all rate/flux diagnostics. The original 11 entries are:

| Producer | Registry variable | Consumed by |
|---|---|---|
| FloatingAlgae | `algal_growth_rate` | DOX, Carbon, Nitrogen, Phosphorus, Alkalinity |
| FloatingAlgae | `algal_respiration_rate` | DOX, Carbon, Nitrogen |
| FloatingAlgae | `algal_death_rate` | OrgN, OrgP, POC, DOC, POM |
| FloatingAlgae | `algal_nh4_uptake_fraction` | Nitrogen |
| BenthicAlgae | analogous to FloatingAlgae | same |
| Nitrogen | `nitrification_rate` | DOX, Alkalinity |
| Nitrogen | `denitrification_rate` | Alkalinity, N2 |
| Carbon | `doc_dic_oxidation_rate` | DOX |
| CBOD | `cbod_oxidation_rate` | DOX |
| Sediment-globals | `sod_rate`, `nh4_from_bed`, `dip_from_bed`, `no3_from_bed_denit`, `dic_from_bed` | DOX, Nitrogen, Phosphorus, Carbon |

## 2. The pattern-alignment spec §4 scope

The §4 per-Process component inventory broadens the surface to *every* component flux/rate that each Process computes inside its `_change_with_components` helper. This is the diagnostics surface a calibration user wants to subscribe to — not just the inter-process coupling subset. The §4 inventory thus *includes* the original Appendix A entries as a subset and *adds* a larger set of intra-Process diagnostics.

## 3. Amended Appendix A — full registry-diagnostics surface per Process

The amendment expands Appendix A from a coupling cheatsheet into a complete registry-diagnostics catalog. Each name is unique across the catalog. Each name maps to exactly one (Process, cached `self.<name>`, opportunistic-write target).

| Process | Registry diagnostic | Sibling consumers (if any) | Notes |
|---|---|---|---|
| **FloatingAlgae** | `algal_growth_rate` | DOX, Carbon, Nitrogen, Phosphorus, Alkalinity | Inherited from 1.0.0 Appendix A |
| FloatingAlgae | `algal_respiration_rate` | DOX, Carbon, Nitrogen, Alkalinity | Inherited; Alkalinity added |
| FloatingAlgae | `algal_death_rate` | OrgN, OrgP, POC, DOC, POM | Inherited |
| FloatingAlgae | `algal_settling_rate` | POM | New |
| FloatingAlgae | `algal_orgn_from_mortality_rate` | Nitrogen | New (already cached on `self` post Phase 2.A; now exposed) |
| FloatingAlgae | `algal_orgp_from_mortality_rate` | Phosphorus | New (same) |
| FloatingAlgae | `algal_poc_from_mortality_rate` | Carbon | New (same) |
| FloatingAlgae | `algal_doc_from_mortality_rate` | Carbon | New (same) |
| FloatingAlgae | `algal_pom_from_settling_rate` | POM | New (same) |
| FloatingAlgae | `algal_nh4_uptake_fraction` | Nitrogen, DOX | Inherited; DOX consumer added |
| FloatingAlgae | `algal_light_limitation` | (diagnostic) | New |
| FloatingAlgae | `algal_nutrient_limitation_n` | (diagnostic) | New |
| FloatingAlgae | `algal_nutrient_limitation_p` | (diagnostic) | New |
| **BenthicAlgae** | `balgae_growth_rate` | DOX, Carbon, Nitrogen, Phosphorus, Alkalinity | Inherited |
| BenthicAlgae | `balgae_respiration_rate` | DOX, Carbon, Nitrogen, Alkalinity | Inherited |
| BenthicAlgae | `balgae_death_rate` | OrgN, OrgP, POC, DOC, POM | Inherited |
| BenthicAlgae | `balgae_orgn_from_mortality_rate` | Nitrogen | New |
| BenthicAlgae | `balgae_orgp_from_mortality_rate` | Phosphorus | New |
| BenthicAlgae | `balgae_poc_from_mortality_rate` | Carbon | New |
| BenthicAlgae | `balgae_doc_from_mortality_rate` | Carbon | New |
| BenthicAlgae | `balgae_nh4_uptake_fraction` | Nitrogen, DOX | New |
| BenthicAlgae | `balgae_light_limitation` | (diagnostic) | New |
| BenthicAlgae | `balgae_nutrient_limitation_n` | (diagnostic) | New |
| BenthicAlgae | `balgae_nutrient_limitation_p` | (diagnostic) | New |
| **Nitrogen** | `nitrification_flux_rate` | DOX, Alkalinity | Renamed from `nitrification_rate` (now distinguishes the *flux* from the kinetic *rate constant* `knit_20`). Migration note in spec §11. |
| Nitrogen | `denitrification_flux_rate` | Alkalinity, N2 | Renamed from `denitrification_rate`. Same rationale. |
| Nitrogen | `nh4_from_bed` | (Nitrogen-internal; was a sediment global) | Moved from Sediment-globals row to Nitrogen producer |
| Nitrogen | `no3_from_bed_denit` | (Nitrogen-internal) | Moved as above |
| Nitrogen | `orgn_hydrolysis_rate` | (diagnostic) | New |
| Nitrogen | `orgn_settling_rate` | (diagnostic) | New |
| Nitrogen | `nh4_algal_growth_rate` | (diagnostic) | New |
| Nitrogen | `no3_algal_growth_rate` | (diagnostic) | New |
| Nitrogen | `nh4_algal_resp_rate` | (diagnostic) | New |
| Nitrogen | `nh4_balgae_resp_rate` | (diagnostic) | New |
| **Phosphorus** | `orgp_hydrolysis_rate` | (diagnostic) | New |
| Phosphorus | `orgp_settling_rate` | (diagnostic) | Already on `self` |
| Phosphorus | `tip_settling_rate` | (diagnostic) | Already on `self` |
| Phosphorus | `dip_from_bed` | (Phosphorus-internal; was a sediment global) | Moved |
| Phosphorus | `orgp_algal_mortality_rate` | (diagnostic) | New |
| Phosphorus | `tip_algal_growth_rate` | (diagnostic) | New |
| Phosphorus | `tip_balgae_growth_rate` | (diagnostic) | New |
| **Carbon** | `doc_dic_oxidation_rate` | DOX | Inherited |
| Carbon | `poc_hydrolysis_rate` | (diagnostic) | New |
| Carbon | `dic_atm_exchange_rate` | (diagnostic) | New |
| Carbon | `dic_sed_release_rate` | (Carbon-internal; was `dic_from_bed` sediment global) | Renamed and moved |
| Carbon | `carbon_algal_resp_rate` | (diagnostic) | New |
| Carbon | `carbon_balgae_resp_rate` | (diagnostic) | New |
| Carbon | `carbon_algal_photo_rate` | (diagnostic) | New |
| Carbon | `carbon_balgae_photo_rate` | (diagnostic) | New |
| Carbon | `carbon_cbod_oxidation_rate` | (diagnostic) | New |
| **DOX** | `dox_sat` | (diagnostic) | Already on `self` |
| DOX | `atm_reaeration_rate` | (diagnostic) | Already on `self` |
| DOX | `dox_nitrification_rate` | (diagnostic) | Already on `self` |
| DOX | `dox_sod_rate` | (diagnostic) | Already on `self` |
| DOX | `dox_doc_oxidation_rate` | (diagnostic) | New |
| DOX | `dox_cbod_oxidation_rate` | (diagnostic) | New |
| DOX | `dox_algal_photo_rate` | (diagnostic) | New |
| DOX | `dox_algal_resp_rate` | (diagnostic) | New |
| DOX | `dox_balgae_photo_rate` | (diagnostic) | New |
| DOX | `dox_balgae_resp_rate` | (diagnostic) | New |
| DOX | `sod_rate` | (DOX-internal; was a sediment global) | Moved |
| **POM** | `pom_hydrolysis_rate` | (diagnostic) | Already on `self` (cache currently lives in `rate()`; moves to `run` in Phase 7) |
| POM | `pom_settling_rate` | (diagnostic) | New |
| POM | `pom_algal_mortality_rate` | (diagnostic) | New |
| POM | `pom_balgae_mortality_rate` | (diagnostic) | New |
| **CBOD** | `cbod_oxidation_rate` | DOX | Inherited (sum over groups) |
| CBOD | `cbod_settling_rate` | (diagnostic) | New |
| **N2** | `n2_atm_exchange_rate` | (diagnostic) | Already on `self` |
| N2 | `n2_sat` | (diagnostic) | Already on `self` |
| N2 | `total_dissolved_gas` | (diagnostic) | Already on `self`, already opportunistically exposed (the sole pre-existing example) |
| N2 | `n2_denit_source_rate` | (diagnostic) | New |
| **Pathogen** | `pathogen_natural_death_rate` | (diagnostic) | New |
| Pathogen | `pathogen_light_death_rate` | (diagnostic) | New |
| Pathogen | `pathogen_settling_rate` | (diagnostic) | New |
| **Alkalinity** | `alk_nitrification_sink_rate` | (diagnostic) | New |
| Alkalinity | `alk_denitrification_source_rate` | (diagnostic) | New |
| Alkalinity | `alk_algal_growth_rate` | (diagnostic) | New |
| Alkalinity | `alk_algal_respiration_rate` | (diagnostic) | New |
| Alkalinity | `alk_balgae_growth_rate` | (diagnostic) | New |
| Alkalinity | `alk_balgae_respiration_rate` | (diagnostic) | New |

**Total names:** 70 (vs 11 in the 1.0.0 Appendix A coupling cheatsheet).

## 4. Renames and migrations

Two renames in §3 deserve explicit migration documentation in the LimnoTech-facing review packet (Phase 11):

| Old name (1.0.0) | New name (1.0.1) | Reason |
|---|---|---|
| `nitrification_rate` (registry / Appendix A) | `nitrification_flux_rate` (registry); `nitrification_rate` stays as a *kinetic rate constant* attribute name | The 1.0.0 entry conflated the kinetic rate constant (`knit_20`, units 1/d) with the realized flux (mg-N/L/d). The instance-attribute names `self.nitrification_flux_rate` and `self.denitrification_flux_rate` already exist in code post-Phase 9; the registry name is harmonized to match. |
| `denitrification_rate` | `denitrification_flux_rate` | Same. |

In-code `self.nitrification_flux_rate` and `self.denitrification_flux_rate` already exist (no code change); only the registry name conventions are updated. DOX, Alkalinity, and N2 already read from the sibling-cache attributes (`getattr(nitrogen_process, "nitrification_flux_rate")`) so no consumer-side change is needed either.

Three sediment-global names move from a single "Sediment-globals" producer row to per-Process ownership:

| Old producer | New producer | Name |
|---|---|---|
| Sediment-globals | Nitrogen | `nh4_from_bed`, `no3_from_bed_denit` |
| Sediment-globals | Phosphorus | `dip_from_bed` |
| Sediment-globals | Carbon | `dic_sed_release_rate` (renamed from `dic_from_bed`) |
| Sediment-globals | DOX | `sod_rate` |

This is a documentation migration, not a code migration — the sediment-flux scalars are already implemented inside the corresponding Processes per 1.0.0 design spec §14 (the "Sediment-globals" attribution in the original Appendix A was an artifact of the early-design SedFlux-Process plan that did not materialize for 1.0.0).

## 5. Naming-convention check

All 70 names in §3 conform to the §14 conventions:
- snake_case.
- `_rate` suffix for time-derivative quantities (62 of 70).
- `_fraction` suffix for dimensionless ratios (2 of 70: `algal_nh4_uptake_fraction`, `balgae_nh4_uptake_fraction`).
- Source-named prefixes for sediment fluxes (`*_from_bed`, `*_sed_release_rate`, `sod_rate`).
- Process-name disambiguation where the same physics appears in multiple Processes' diagnostics: `dox_nitrification_rate` (DOX's contribution from nitrification) vs `nitrification_flux_rate` (Nitrogen's actual nitrification flux).

No collisions across Processes. No conflicts with existing v3 state-variable names.

## 6. Code-side implications

- **Class attribute `REGISTRY_DIAGNOSTICS: tuple[str, ...]`** must be added to every Process in Phase 1. Phase 1 lands the empty tuple (no opportunistic writes yet); each per-Process phase populates its Process's tuple.
- The `self.<name>` cache attributes for the 70 names are already a mix of "already exists" (~12 names cached today on the various Process instances) and "to be added by Phase 2–9". The component inventory in pattern-alignment spec §4 will be the implementation reference.
- **No registry-side schema change.** The registry is implemented in `clearwater_data` and accepts any `set_at_time(name, time, value)` for an already-registered name. The opportunistic-write pattern depends only on registry membership (`name in registry`), not on schema declaration.
