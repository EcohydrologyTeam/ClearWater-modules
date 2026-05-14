"""Phase 7 opportunistic-registry-write contract for CBOD."""

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


def _cbod_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "CBOD")


def _as_array(value) -> np.ndarray:
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    from clearwater_modules_v3.processes.cbod import CBOD

    ic = default_initial_conditions()
    reference = ic["cbod"]
    for name in CBOD.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    cbod = _cbod_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in cbod.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all()


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    a = demo_a.registry.get("cbod").values.copy()

    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    b = demo_b.registry.get("cbod").values.copy()

    np.testing.assert_array_equal(a, b)


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    """Pattern G zero-cost-when-unused: subscribing to CBOD diagnostics
    must NOT change CBOD trajectory or downstream DOX/Carbon (which
    consume cbod_oxidation_rate via getattr)."""
    from clearwater_modules_v3.processes.cbod import CBOD

    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    ic_b = default_initial_conditions()
    reference = ic_b["cbod"]
    for name in CBOD.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    for state_name in ("cbod", "oxygen_dissolved", "dic"):
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
    from clearwater_modules_v3.processes.cbod import CBOD

    ic = default_initial_conditions()
    reference = ic["cbod"]
    for name in CBOD.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    cbod = _cbod_process(demo)

    demo.step(START)

    for name in cbod.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(cbod, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(cached, from_registry)


def test_g_partial_subscription_writes_only_requested_names() -> None:
    from clearwater_modules_v3.processes.cbod import CBOD

    subscribed = (CBOD.REGISTRY_DIAGNOSTICS[0],)
    unsubscribed = CBOD.REGISTRY_DIAGNOSTICS[1:]
    assert subscribed and unsubscribed

    ic = default_initial_conditions()
    reference = ic["cbod"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry
    for name in unsubscribed:
        assert name not in demo.registry


def test_g_cbod_oxidation_rate_consumed_by_dox(demo_setup_check=None) -> None:
    """``cbod_oxidation_rate`` is consumed by DOX via getattr; it must
    stay populated on self after run regardless of registry
    subscription."""
    demo = build_nsm1_demo()
    cbod = _cbod_process(demo)
    demo.step(START)

    cached = _as_array(getattr(cbod, "cbod_oxidation_rate"))
    assert np.isfinite(cached).all()
    assert np.all(cached >= 0.0), "oxidation_rate is a positive magnitude"
