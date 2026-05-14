"""Phase 9 opportunistic-registry-write contract for Alkalinity."""

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


def _alkalinity_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "Alkalinity")


def _as_array(value) -> np.ndarray:
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    from clearwater_modules_v3.processes.alkalinity import Alkalinity

    ic = default_initial_conditions()
    reference = ic["alkalinity"]
    for name in Alkalinity.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    alk = _alkalinity_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in alk.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all()


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    a = demo_a.registry.get("alkalinity").values.copy()

    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    b = demo_b.registry.get("alkalinity").values.copy()

    np.testing.assert_array_equal(a, b)


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    from clearwater_modules_v3.processes.alkalinity import Alkalinity

    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    ic_b = default_initial_conditions()
    reference = ic_b["alkalinity"]
    for name in Alkalinity.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    a = demo_a.registry.get("alkalinity").values
    b = demo_b.registry.get("alkalinity").values
    np.testing.assert_array_equal(a, b)


def test_g_diagnostics_attribute_caches_match_registry_writes() -> None:
    from clearwater_modules_v3.processes.alkalinity import Alkalinity

    ic = default_initial_conditions()
    reference = ic["alkalinity"]
    for name in Alkalinity.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    alk = _alkalinity_process(demo)

    demo.step(START)

    for name in alk.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(alk, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(cached, from_registry)


def test_g_partial_subscription_writes_only_requested_names() -> None:
    from clearwater_modules_v3.processes.alkalinity import Alkalinity

    subscribed = tuple(Alkalinity.REGISTRY_DIAGNOSTICS[::2])
    unsubscribed = tuple(Alkalinity.REGISTRY_DIAGNOSTICS[1::2])
    assert subscribed and unsubscribed

    ic = default_initial_conditions()
    reference = ic["alkalinity"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry
    for name in unsubscribed:
        assert name not in demo.registry


def test_g_legacy_attribute_aliases_remain_set() -> None:
    """The four legacy attribute names
    (``alk_nitrification_rate``, ``alk_denitrification_rate``,
    ``alk_benthic_algae_growth_rate``, ``alk_benthic_algae_respiration_rate``)
    consumed by ``test_alkalinity_v1_parity_v3.py`` and
    ``test_alkalinity_tier1.py`` must remain populated on self after
    run, regardless of registry subscription. NOT in REGISTRY_DIAGNOSTICS
    so pattern G does not write them to the registry."""
    demo = build_nsm1_demo()
    alk = _alkalinity_process(demo)
    demo.step(START)

    for legacy_name in (
        "alk_nitrification_rate",
        "alk_denitrification_rate",
        "alk_benthic_algae_growth_rate",
        "alk_benthic_algae_respiration_rate",
    ):
        cached = _as_array(getattr(alk, legacy_name))
        assert np.isfinite(cached).all(), (
            f"legacy attribute {legacy_name} has non-finite values"
        )
        assert legacy_name not in demo.registry, (
            f"{legacy_name} is a legacy alias; not in REGISTRY_DIAGNOSTICS"
        )
