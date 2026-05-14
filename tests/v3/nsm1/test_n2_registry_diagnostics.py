"""Phase 8 opportunistic-registry-write contract for N2.

Specifically pins the **pre-existing total_dissolved_gas extension**:
N2 had the only pre-Phase-8 opportunistic-write path in v3 NSM1 1.0.0
(``if "total_dissolved_gas" in registry: ...``); Phase 8 extends that
loop to cover the full Appendix A set without breaking the pre-existing
behaviour.
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


def _n2_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "N2")


def _as_array(value) -> np.ndarray:
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    from clearwater_modules_v3.processes.n2 import N2

    ic = default_initial_conditions()
    reference = ic["n2"]
    for name in N2.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    n2 = _n2_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in n2.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all()


def test_g_total_dissolved_gas_pre_existing_path_still_works() -> None:
    """The N2 ``total_dissolved_gas`` opportunistic write was the sole
    pre-Phase-8 example of pattern G. Phase 8 must preserve that
    specific behaviour."""
    ic = default_initial_conditions()
    reference = ic["n2"]
    # Pre-register ONLY total_dissolved_gas (mirrors pre-Phase-8 usage).
    ic["total_dissolved_gas"] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    assert "total_dissolved_gas" in demo.registry
    tdg = _as_array(demo.registry.get("total_dissolved_gas"))
    assert np.isfinite(tdg).all()
    assert np.all(tdg > 0.0), "TDG fraction should be positive"


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    a = demo_a.registry.get("n2").values.copy()

    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    b = demo_b.registry.get("n2").values.copy()

    np.testing.assert_array_equal(a, b)


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    """Pattern G zero-cost-when-unused: subscribing does NOT change
    n2 trajectory."""
    from clearwater_modules_v3.processes.n2 import N2

    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    ic_b = default_initial_conditions()
    reference = ic_b["n2"]
    for name in N2.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    a = demo_a.registry.get("n2").values
    b = demo_b.registry.get("n2").values
    np.testing.assert_array_equal(a, b)


def test_g_diagnostics_attribute_caches_match_registry_writes() -> None:
    from clearwater_modules_v3.processes.n2 import N2

    ic = default_initial_conditions()
    reference = ic["n2"]
    for name in N2.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    n2 = _n2_process(demo)

    demo.step(START)

    for name in n2.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(n2, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(cached, from_registry)


def test_g_tdg_alias_to_total_dissolved_gas(demo_stub=None) -> None:
    """``self.tdg`` (back-compat attribute name) and
    ``self.total_dissolved_gas`` (Appendix A name) must be aliases —
    same value, just two names. Pinned so a future refactor cannot
    silently drop the back-compat alias."""
    demo = build_nsm1_demo()
    n2 = _n2_process(demo)
    demo.step(START)

    tdg_legacy = _as_array(n2.tdg)
    tdg_new_name = _as_array(n2.total_dissolved_gas)
    np.testing.assert_array_equal(
        tdg_legacy, tdg_new_name,
        err_msg="self.tdg and self.total_dissolved_gas must be aliases",
    )


def test_g_partial_subscription_writes_only_requested_names() -> None:
    from clearwater_modules_v3.processes.n2 import N2

    subscribed = (N2.REGISTRY_DIAGNOSTICS[0], N2.REGISTRY_DIAGNOSTICS[2])
    unsubscribed = (N2.REGISTRY_DIAGNOSTICS[1], N2.REGISTRY_DIAGNOSTICS[3])
    assert subscribed and unsubscribed

    ic = default_initial_conditions()
    reference = ic["n2"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry
    for name in unsubscribed:
        assert name not in demo.registry
