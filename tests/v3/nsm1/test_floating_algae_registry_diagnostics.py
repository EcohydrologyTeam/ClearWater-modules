"""Phase 5 opportunistic-registry-write contract for FloatingAlgae.

Asserts pattern G (pattern-alignment spec §3): when a name in
``FloatingAlgae.REGISTRY_DIAGNOSTICS`` is pre-registered, the value
cached on ``FloatingAlgae`` is written to the registry each substep;
when no name is pre-registered, the algae trajectory is bit-identical
to the no-diagnostics baseline.

Mirrors ``test_carbon_registry_diagnostics.py`` /
``test_dox_registry_diagnostics.py`` /
``test_nitrogen_registry_diagnostics.py``; retained through Phase 10.
"""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.examples import (
    InMemoryRegistry,
    build_nsm1_demo,
    default_initial_conditions,
)


N_STEPS = 60
START = datetime(2026, 1, 1, 0, 0, 0)


def _floating_algae_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "FloatingAlgae")


def _as_array(value) -> np.ndarray:
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    from clearwater_modules_v3.processes.floating_algae import FloatingAlgae

    ic = default_initial_conditions()
    reference = ic["algae_floating"]
    for name in FloatingAlgae.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    fa = _floating_algae_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in fa.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all(), (
            f"{name} has non-finite values after {N_STEPS} substeps"
        )


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    a = demo_a.registry.get("algae_floating").values.copy()

    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    b = demo_b.registry.get("algae_floating").values.copy()

    np.testing.assert_array_equal(a, b)


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    """Pattern G zero-cost-when-unused: subscribing to FloatingAlgae
    diagnostics does NOT change algae_floating trajectory. Also
    verifies that downstream sibling consumers (Carbon / DOX / Nitrogen
    / POM read FloatingAlgae caches via getattr) see identical values
    in both runs."""
    from clearwater_modules_v3.processes.floating_algae import FloatingAlgae

    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    ic_b = default_initial_conditions()
    reference = ic_b["algae_floating"]
    for name in FloatingAlgae.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    # The downstream-consumer state variables (the whole closed-system
    # set of NSM1 states) must be bit-identical between subscribed and
    # unsubscribed runs. Otherwise pattern G has leaked side effects
    # into kinetics.
    for state_name in (
        "algae_floating", "ammonium", "nitrate", "organic_nitrogen",
        "tip", "organic_phosphorus", "poc", "doc", "dic", "pom",
        "oxygen_dissolved",
    ):
        a = demo_a.registry.get(state_name).values
        b = demo_b.registry.get(state_name).values
        np.testing.assert_array_equal(
            a, b,
            err_msg=(
                f"{state_name} differs between subscribed and unsubscribed "
                "runs; pattern G zero-cost-when-unused invariant violated"
            ),
        )


def test_g_diagnostics_attribute_caches_match_registry_writes() -> None:
    from clearwater_modules_v3.processes.floating_algae import FloatingAlgae

    ic = default_initial_conditions()
    reference = ic["algae_floating"]
    for name in FloatingAlgae.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    fa = _floating_algae_process(demo)

    demo.step(START)

    for name in fa.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(fa, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(
            cached, from_registry,
            err_msg=(
                f"{name}: registry value differs from self.{name}; "
                "pattern F/G single-source-of-truth violated"
            ),
        )


def test_g_partial_subscription_writes_only_requested_names() -> None:
    from clearwater_modules_v3.processes.floating_algae import FloatingAlgae

    subscribed = tuple(FloatingAlgae.REGISTRY_DIAGNOSTICS[::2])
    unsubscribed = tuple(FloatingAlgae.REGISTRY_DIAGNOSTICS[1::2])
    assert subscribed and unsubscribed

    ic = default_initial_conditions()
    reference = ic["algae_floating"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry
    for name in unsubscribed:
        assert name not in demo.registry


def test_g_preserved_names_consumed_by_siblings_remain_visible() -> None:
    """The seven preserved-name caches that Carbon / DOX / Nitrogen /
    Phosphorus / POM read via ``getattr`` must remain populated on
    self after run, regardless of registry subscription."""
    demo = build_nsm1_demo()
    fa = _floating_algae_process(demo)
    demo.step(START)

    for name in (
        "algal_growth_rate",
        "algal_respiration_rate",
        "algal_death_rate",
        "algal_orgn_from_mortality_rate",
        "algal_orgp_from_mortality_rate",
        "algal_poc_from_mortality_rate",
        "algal_doc_from_mortality_rate",
        "algal_pom_from_settling_rate",
        "algal_nh4_uptake_fraction",
    ):
        value = _as_array(getattr(fa, name))
        assert np.isfinite(value).all(), (
            f"sibling-consumer attribute {name} has non-finite values"
        )
