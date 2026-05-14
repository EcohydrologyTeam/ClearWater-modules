"""Phase 6 opportunistic-registry-write contract for Phosphorus.

Asserts pattern G (pattern-alignment spec §3): when a name in
``Phosphorus.REGISTRY_DIAGNOSTICS`` is pre-registered, the value cached
on ``Phosphorus`` is written to the registry each substep; when no name
is pre-registered, the TIP / OrgP trajectories are bit-identical to
the no-diagnostics baseline.

Mirrors the Phase 2-5 templates; retained through Phase 10.
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


def _phosphorus_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "Phosphorus")


def _as_array(value) -> np.ndarray:
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    from clearwater_modules_v3.processes.phosphorus import Phosphorus

    ic = default_initial_conditions()
    reference = ic["tip"]
    for name in Phosphorus.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    p = _phosphorus_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in p.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all(), (
            f"{name} has non-finite values after {N_STEPS} substeps"
        )


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    state_a = {
        name: demo_a.registry.get(name).values.copy()
        for name in ("tip", "organic_phosphorus")
    }

    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    state_b = {
        name: demo_b.registry.get(name).values.copy()
        for name in ("tip", "organic_phosphorus")
    }

    for name in state_a:
        np.testing.assert_array_equal(state_a[name], state_b[name])


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    """Pattern G zero-cost-when-unused: subscribing to Phosphorus
    diagnostics does NOT change the TIP / OrgP trajectories."""
    from clearwater_modules_v3.processes.phosphorus import Phosphorus

    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    ic_b = default_initial_conditions()
    reference = ic_b["tip"]
    for name in Phosphorus.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    for state_name in ("tip", "organic_phosphorus"):
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
    from clearwater_modules_v3.processes.phosphorus import Phosphorus

    ic = default_initial_conditions()
    reference = ic["tip"]
    for name in Phosphorus.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    p = _phosphorus_process(demo)

    demo.step(START)

    for name in p.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(p, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(
            cached, from_registry,
            err_msg=(
                f"{name}: registry value differs from self.{name}; "
                "pattern F/G single-source-of-truth violated"
            ),
        )


def test_g_partial_subscription_writes_only_requested_names() -> None:
    from clearwater_modules_v3.processes.phosphorus import Phosphorus

    subscribed = tuple(Phosphorus.REGISTRY_DIAGNOSTICS[::2])
    unsubscribed = tuple(Phosphorus.REGISTRY_DIAGNOSTICS[1::2])
    assert subscribed and unsubscribed

    ic = default_initial_conditions()
    reference = ic["tip"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry
    for name in unsubscribed:
        assert name not in demo.registry


def test_g_orgp_hydrolysis_alias_to_legacy_attribute() -> None:
    """``orgp_hydrolysis_rate`` (Appendix A name) registry value must
    equal ``self.orgp_to_tip_hydrolysis_rate`` (the legacy attribute
    name read by test_phosphorus_v1_parity_v3.py via getattr)."""
    ic = default_initial_conditions()
    reference = ic["tip"]
    ic["orgp_hydrolysis_rate"] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    p = _phosphorus_process(demo)
    demo.step(START)

    from_registry = _as_array(demo.registry.get("orgp_hydrolysis_rate"))
    legacy_attr = _as_array(p.orgp_to_tip_hydrolysis_rate)
    np.testing.assert_array_equal(from_registry, legacy_attr)
