"""Phase 10 Appendix A completeness — every name in the
``clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3 catalog maps to
exactly one Process's ``REGISTRY_DIAGNOSTICS`` tuple. No orphans, no
duplicates.

Retained indefinitely. A future Appendix A amendment must also update
this catalog (or the test fails, surfacing the divergence).
"""

from __future__ import annotations

import pytest

from clearwater_modules_v3.processes import (
    Alkalinity,
    BenthicAlgae,
    CBOD,
    Carbon,
    DOX,
    FloatingAlgae,
    N2,
    Nitrogen,
    POM,
    Pathogen,
    Phosphorus,
)


# The full Appendix A diff §3 catalog grouped by Process. This is the
# load-bearing contract: each entry below MUST appear in the
# corresponding ``REGISTRY_DIAGNOSTICS`` tuple.
APPENDIX_A_CATALOG: dict[type, set[str]] = {
    FloatingAlgae: {
        "algal_growth_rate",
        "algal_respiration_rate",
        "algal_death_rate",
        "algal_settling_rate",
        "algal_orgn_from_mortality_rate",
        "algal_orgp_from_mortality_rate",
        "algal_poc_from_mortality_rate",
        "algal_doc_from_mortality_rate",
        "algal_pom_from_settling_rate",
        "algal_nh4_uptake_fraction",
        "algal_light_limitation",
        "algal_nutrient_limitation_n",
        "algal_nutrient_limitation_p",
    },
    BenthicAlgae: {
        "balgae_growth_rate",
        "balgae_respiration_rate",
        "balgae_death_rate",
        "balgae_orgn_from_mortality_rate",
        "balgae_orgp_from_mortality_rate",
        "balgae_poc_from_mortality_rate",
        "balgae_doc_from_mortality_rate",
        "balgae_nh4_uptake_fraction",
        "balgae_light_limitation",
        "balgae_nutrient_limitation_n",
        "balgae_nutrient_limitation_p",
    },
    Nitrogen: {
        "nitrification_flux_rate",
        "denitrification_flux_rate",
        "nh4_from_bed",
        "no3_from_bed_denit",
        "orgn_hydrolysis_rate",
        "orgn_settling_rate",
        "nh4_algal_growth_rate",
        "no3_algal_growth_rate",
        "nh4_algal_resp_rate",
        "nh4_balgae_resp_rate",
    },
    Phosphorus: {
        "orgp_hydrolysis_rate",
        "orgp_settling_rate",
        "tip_settling_rate",
        "dip_from_bed",
        "orgp_algal_mortality_rate",
        "tip_algal_growth_rate",
        "tip_balgae_growth_rate",
    },
    Carbon: {
        "poc_hydrolysis_rate",
        "doc_dic_oxidation_rate",
        "dic_atm_exchange_rate",
        "dic_sed_release_rate",
        "carbon_algal_resp_rate",
        "carbon_balgae_resp_rate",
        "carbon_algal_photo_rate",
        "carbon_balgae_photo_rate",
        "carbon_cbod_oxidation_rate",
    },
    DOX: {
        "dox_sat",
        "atm_reaeration_rate",
        "dox_nitrification_rate",
        "dox_sod_rate",
        "dox_doc_oxidation_rate",
        "dox_cbod_oxidation_rate",
        "dox_algal_photo_rate",
        "dox_algal_resp_rate",
        "dox_balgae_photo_rate",
        "dox_balgae_resp_rate",
        "sod_rate",
    },
    POM: {
        "pom_hydrolysis_rate",
        "pom_settling_rate",
        "pom_algal_mortality_rate",
        "pom_balgae_mortality_rate",
    },
    CBOD: {
        "cbod_oxidation_rate",
        "cbod_settling_rate",
    },
    N2: {
        "n2_atm_exchange_rate",
        "n2_sat",
        "total_dissolved_gas",
        "n2_denit_source_rate",
    },
    Pathogen: {
        "pathogen_natural_death_rate",
        "pathogen_light_death_rate",
        "pathogen_settling_rate",
    },
    Alkalinity: {
        "alk_nitrification_sink_rate",
        "alk_denitrification_source_rate",
        "alk_algal_growth_rate",
        "alk_algal_respiration_rate",
        "alk_balgae_growth_rate",
        "alk_balgae_respiration_rate",
    },
}


# ---------------------------------------------------------------------------
# Per-Process completeness
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cls", list(APPENDIX_A_CATALOG.keys()), ids=lambda c: c.__name__)
def test_process_registry_diagnostics_matches_appendix_a(cls) -> None:
    """Each Process's REGISTRY_DIAGNOSTICS equals the Appendix A catalog
    entry exactly (same set of names)."""
    expected = APPENDIX_A_CATALOG[cls]
    actual = set(cls.REGISTRY_DIAGNOSTICS)
    missing = expected - actual
    extra = actual - expected
    assert not missing, (
        f"{cls.__name__}: REGISTRY_DIAGNOSTICS missing Appendix A names: "
        f"{sorted(missing)}"
    )
    assert not extra, (
        f"{cls.__name__}: REGISTRY_DIAGNOSTICS has names not in Appendix A: "
        f"{sorted(extra)}"
    )


# ---------------------------------------------------------------------------
# Cross-Process global uniqueness
# ---------------------------------------------------------------------------


def test_no_duplicate_names_across_processes() -> None:
    """Each Appendix A name belongs to exactly one Process. A duplicate
    would mean two Processes both write the same registry variable on
    pre-registration — ambiguous semantics, registered explicitly as a
    spec rule."""
    name_to_process: dict[str, str] = {}
    duplicates: dict[str, list[str]] = {}

    for cls, names in APPENDIX_A_CATALOG.items():
        for name in names:
            if name in name_to_process:
                duplicates.setdefault(name, [name_to_process[name]]).append(
                    cls.__name__
                )
            else:
                name_to_process[name] = cls.__name__

    assert not duplicates, (
        f"Appendix A names appearing in multiple Processes' "
        f"REGISTRY_DIAGNOSTICS: {duplicates}"
    )


def test_appendix_a_catalog_covers_all_pattern_aligned_processes() -> None:
    """The catalog covers all 11 pattern-aligned NSM1 Processes."""
    expected = {
        FloatingAlgae, BenthicAlgae, Nitrogen, Phosphorus, Carbon, DOX,
        POM, CBOD, N2, Pathogen, Alkalinity,
    }
    assert set(APPENDIX_A_CATALOG.keys()) == expected


def test_catalog_total_name_count_is_expected() -> None:
    """Sanity check: the total number of Appendix A names across all
    Processes is exactly 80 (per Appendix A diff §3 row count). A
    silent addition or removal of a name surfaces here."""
    total = sum(len(names) for names in APPENDIX_A_CATALOG.values())
    expected_total = 80
    assert total == expected_total, (
        f"Total Appendix A names: {total}, expected {expected_total}. "
        "If this is an intentional Appendix A amendment, update the "
        "expected_total and document the change in the closeout."
    )
