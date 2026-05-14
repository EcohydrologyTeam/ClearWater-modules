"""Phase 5 opportunistic-registry-write contract for BenthicAlgae.

Asserts pattern G (pattern-alignment spec §3): when a name in
``BenthicAlgae.REGISTRY_DIAGNOSTICS`` is pre-registered, the value
cached on ``BenthicAlgae`` is written to the registry each substep;
when no name is pre-registered, the benthic_algae trajectory is
bit-identical to the no-diagnostics baseline.

Mirrors the Phase 2/3/4/5 registry-diagnostics templates; retained
through Phase 10.
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


def _benthic_algae_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "BenthicAlgae")


def _as_array(value) -> np.ndarray:
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae

    ic = default_initial_conditions()
    reference = ic["benthic_algae"]
    for name in BenthicAlgae.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    ba = _benthic_algae_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in ba.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all(), (
            f"{name} has non-finite values after {N_STEPS} substeps"
        )


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    a = demo_a.registry.get("benthic_algae").values.copy()

    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    b = demo_b.registry.get("benthic_algae").values.copy()

    np.testing.assert_array_equal(a, b)


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    """Pattern G zero-cost-when-unused: subscribing to BenthicAlgae
    diagnostics does NOT change the closed-system NSM1 state. This
    also implicitly verifies the Phase 5 rate_death dedup: if the
    dedup had broken bit-identicality, this test would fail."""
    from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae

    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    ic_b = default_initial_conditions()
    reference = ic_b["benthic_algae"]
    for name in BenthicAlgae.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    for state_name in (
        "benthic_algae", "ammonium", "nitrate", "organic_nitrogen",
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
    from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae

    ic = default_initial_conditions()
    reference = ic["benthic_algae"]
    for name in BenthicAlgae.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    ba = _benthic_algae_process(demo)

    demo.step(START)

    for name in ba.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(ba, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(
            cached, from_registry,
            err_msg=(
                f"{name}: registry value differs from self.{name}; "
                "pattern F/G single-source-of-truth violated"
            ),
        )


def test_g_partial_subscription_writes_only_requested_names() -> None:
    from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae

    subscribed = tuple(BenthicAlgae.REGISTRY_DIAGNOSTICS[::2])
    unsubscribed = tuple(BenthicAlgae.REGISTRY_DIAGNOSTICS[1::2])
    assert subscribed and unsubscribed

    ic = default_initial_conditions()
    reference = ic["benthic_algae"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry
    for name in unsubscribed:
        assert name not in demo.registry


def test_g_pom_routing_cache_remains_set() -> None:
    """``balgae_pom_from_mortality_rate`` is consumed by POM but is
    not in REGISTRY_DIAGNOSTICS. Pattern G does not write it to the
    registry, but pattern F still sets it on self for sibling reads."""
    demo = build_nsm1_demo()
    ba = _benthic_algae_process(demo)
    demo.step(START)

    pom_cache = _as_array(getattr(ba, "balgae_pom_from_mortality_rate"))
    assert np.isfinite(pom_cache).all()
    assert "balgae_pom_from_mortality_rate" not in demo.registry, (
        "balgae_pom_from_mortality_rate is not in REGISTRY_DIAGNOSTICS; "
        "pattern G should not write it to the registry"
    )
