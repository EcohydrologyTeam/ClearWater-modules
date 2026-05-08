"""Tests for the headless NSM1 demo builder (Phase 8.A).

Covers ``clearwater_modules_v3.examples.nsm1_demo_setup``:

* :func:`build_nsm1_demo` returns a runnable handle whose 11 Process
  instances all advance state without raising. This validates the
  Phase 8.A acceptance criterion that ``Nitrogen()`` instantiates
  and runs without ``AttributeError`` on the five legacy
  algal-uptake attributes (``floating_algae_nitrogen_weight``,
  ``benthic_algae_nitrogen_weight``, ``algal_chlorophyll``,
  ``benthic_algea_faction_uptake_from_nitrate``,
  ``fraction_bottom_area``).
* The v2 algae overlays read the canonical v3 inorganic-P state
  name ``tip`` (no ``phosphorus_total_inorganic`` mirror needed).
* The default initial-conditions and parameter dicts cover every
  state variable and every Process class declared by the v3 NSM1
  module.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.examples import (
    InMemoryRegistry,
    build_nsm1_demo,
    default_initial_conditions,
    default_process_parameters,
)
from clearwater_modules_v3.examples.nsm1_demo_setup import N_CELLS_DEFAULT


# ---------------------------------------------------------------------------
# Default initial conditions / parameters
# ---------------------------------------------------------------------------


def test_default_initial_conditions_covers_all_nsm1_states():
    """Every state variable read by the 11 NSM1 Processes must be
    present in the default IC dict, otherwise the demo build raises
    ``KeyError`` at the first ``run`` call."""
    ic = default_initial_conditions()
    expected = {
        "ammonium",
        "nitrate",
        "organic_nitrogen",
        "n2",
        "tip",
        "organic_phosphorus",
        "poc",
        "doc",
        "dic",
        "pom",
        "cbod",
        "oxygen_dissolved",
        "pathogen",
        "alkalinity",
        "algae_floating",
        "benthic_algae",
    }
    forcings = {
        "water_temperature",
        "depth",
        "solar_radiation",
        "atmospheric_pressure",
    }
    assert expected.issubset(ic.keys())
    assert forcings.issubset(ic.keys())


def test_default_initial_conditions_no_legacy_phosphorus_name():
    """Phase 8.A Fix 2: the canonical v3 inorganic-P state name is
    ``tip``. The legacy v2 name ``phosphorus_total_inorganic`` must
    not appear in the default IC dict."""
    ic = default_initial_conditions()
    assert "tip" in ic
    assert "phosphorus_total_inorganic" not in ic


def test_default_initial_conditions_n_cells_param():
    ic = default_initial_conditions(n_cells=7)
    assert ic["ammonium"].sizes == {"cell": 7}


def test_default_process_parameters_covers_all_11_processes():
    params = default_process_parameters()
    expected = {
        "FloatingAlgae",
        "BenthicAlgae",
        "Nitrogen",
        "Phosphorus",
        "Carbon",
        "POM",
        "CBOD",
        "DOX",
        "N2",
        "Pathogen",
        "Alkalinity",
    }
    assert expected == set(params.keys())


# ---------------------------------------------------------------------------
# build_nsm1_demo
# ---------------------------------------------------------------------------


def test_build_nsm1_demo_returns_runnable_handle():
    demo = build_nsm1_demo()
    assert isinstance(demo.registry, InMemoryRegistry)
    assert len(demo.processes) == 11
    assert demo.time_step == timedelta(minutes=5)


def test_build_nsm1_demo_step_advances_without_attribute_error():
    """Phase 8.A Fix 1 acceptance: a fresh ``Nitrogen()`` runs without
    ``AttributeError`` on the five legacy algal-uptake attributes."""
    demo = build_nsm1_demo()
    t = datetime(2026, 1, 1, 0, 0, 0)
    # A single step exercises every Process's ``run`` and the
    # rate-cache propagation between producers and consumers.
    demo.step(t)


def test_build_nsm1_demo_run_60_substeps_stable():
    """A short 5-hour simulation should run end-to-end without
    raising and leave every state variable finite (no NaN/Inf)."""
    demo = build_nsm1_demo()
    start = datetime(2026, 1, 1, 0, 0, 0)
    demo.run(start, n_steps=60)  # 5 minutes * 60 = 5 hours
    for name in demo.registry.keys():
        da = demo.registry.get(name)
        # Every variable must be finite. NaN propagation through any
        # rate calc would show up here.
        assert np.isfinite(da.values).all(), (
            f"Variable {name!r} contains non-finite values after 60 substeps"
        )


def test_build_nsm1_demo_phosphorus_persists_to_tip():
    """Phase 8.A Fix 2 acceptance: after running, the inorganic-P
    state lives at the registry key ``tip`` (not
    ``phosphorus_total_inorganic``)."""
    demo = build_nsm1_demo()
    start = datetime(2026, 1, 1, 0, 0, 0)
    demo.step(start)
    assert "tip" in demo.registry
    # The legacy name must NOT have been mirrored by any v3 Process.
    assert "phosphorus_total_inorganic" not in demo.registry


def test_nitrogen_legacy_attributes_seeded_from_v3_defaults():
    """Phase 8.A Fix 1: ``Nitrogen()`` seeds the five legacy
    algal-uptake attributes from the v3 ``ALGAE_DEFAULTS`` /
    ``BALGAE_DEFAULTS`` so the uptake paths run with sensible
    physics-aware defaults rather than raising ``AttributeError``."""
    from clearwater_modules_v3.parameters.algae import DEFAULTS as ALGAE
    from clearwater_modules_v3.parameters.balgae import DEFAULTS as BALGAE
    from clearwater_modules_v3.processes.nitrogen import Nitrogen

    nitrogen = Nitrogen()
    assert nitrogen.floating_algae_nitrogen_weight == ALGAE["AWn"]
    assert nitrogen.benthic_algae_nitrogen_weight == BALGAE["BWn"]
    assert nitrogen.algal_chlorophyll == ALGAE["AWa"]
    # Note the legacy "algea" typo on the attribute name (preserved
    # for back-compat). 0.5 matches the floating-algae default.
    assert nitrogen.benthic_algea_faction_uptake_from_nitrate == 0.5
    assert nitrogen.fraction_bottom_area == 1.0


def test_build_nsm1_demo_custom_initial_conditions():
    """Caller can override initial conditions wholesale."""
    custom_ic = default_initial_conditions(n_cells=3)
    custom_ic["ammonium"] = xr.DataArray(np.array([1.0, 2.0, 3.0]), dims="cell")
    demo = build_nsm1_demo(initial_conditions=custom_ic)
    assert demo.registry.get("ammonium").sizes == {"cell": 3}
    np.testing.assert_array_equal(
        demo.registry.get("ammonium").values, np.array([1.0, 2.0, 3.0])
    )


def test_build_nsm1_demo_custom_process_parameters():
    """Caller can override per-Process parameter dicts."""
    overrides = {"DOX": {"pressure_mb": 950.0}}
    demo = build_nsm1_demo(process_parameters=overrides)
    dox = next(p for p in demo.processes if type(p).__name__ == "DOX")
    assert dox.pressure_mb == 950.0


# ---------------------------------------------------------------------------
# YAML config exists
# ---------------------------------------------------------------------------


def test_nsm1_default_yaml_exists():
    """Phase 8.A Fix 3: ``nsm1_default.yml`` ships with the v3
    package as the canonical reference for "what NSM1 looks like with
    v3 defaults"."""
    from pathlib import Path

    import clearwater_modules_v3

    yaml_path = (
        Path(clearwater_modules_v3.__file__).parent / "config" / "nsm1_default.yml"
    )
    assert yaml_path.is_file(), (
        f"Expected NSM1 default YAML at {yaml_path}; not found"
    )


def test_nsm1_default_yaml_declares_all_11_processes():
    """The YAML processes list must enumerate all 11 NSM1 Process
    names so the file remains a faithful template even without I/O
    plumbing."""
    from pathlib import Path

    import clearwater_modules_v3
    from clearwater_modules_v3.config.read import read_config

    yaml_path = (
        Path(clearwater_modules_v3.__file__).parent / "config" / "nsm1_default.yml"
    )
    config = read_config(yaml_path)
    process_names = {next(iter(p.keys())) for p in config["processes"]}
    expected = {
        "floating_algae",
        "benthic_algae",
        "nitrogen",
        "phosphorus",
        "carbon",
        "pom",
        "cbod",
        "dox",
        "n2",
        "pathogen",
        "alkalinity",
    }
    assert process_names == expected
