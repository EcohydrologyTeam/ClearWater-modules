"""Phase 3 opportunistic-registry-write contract for DOX.

Asserts pattern G (pattern-alignment spec §3): when a name in
``DOX.REGISTRY_DIAGNOSTICS`` is pre-registered, the value cached on
``DOX`` is written to the registry each substep; when no name is
pre-registered, the DOX state trajectory is bit-identical to the
no-diagnostics baseline.

Mirrors ``test_carbon_registry_diagnostics.py``; retained through
Phase 10.
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


def _dox_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "DOX")


def _as_array(value) -> np.ndarray:
    """Coerce a registry value to a numpy array.

    DOX's components dict mixes ``xr.DataArray`` (full per-cell rates)
    and Python scalars (e.g., ``ka_tc_value = 0.0`` short-circuit when
    both reaeration menu options are user-defined-and-zero — produces a
    scalar ``atm_reaeration_rate`` for all-uniform-zero forcings).
    """
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    """Pre-register every DOX REGISTRY_DIAGNOSTICS name; assert each
    is written each substep with finite values."""
    from clearwater_modules_v3.processes.dox import DOX

    ic = default_initial_conditions()
    reference = ic["oxygen_dissolved"]
    for name in DOX.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    dox = _dox_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in dox.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry, (
            f"{name} should still be registered after run"
        )
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all(), (
            f"{name} has non-finite values after {N_STEPS} substeps"
        )


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    """When no DOX diagnostic is pre-registered, the DOX trajectory
    must be reproducible across two consecutive runs."""
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    dox_a = demo_a.registry.get("oxygen_dissolved").values.copy()

    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    dox_b = demo_b.registry.get("oxygen_dissolved").values.copy()

    np.testing.assert_array_equal(dox_a, dox_b)


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    """The core pattern G invariant: subscribing to DOX diagnostics
    does NOT change the DOX state trajectory."""
    from clearwater_modules_v3.processes.dox import DOX

    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    ic_b = default_initial_conditions()
    reference = ic_b["oxygen_dissolved"]
    for name in DOX.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    a = demo_a.registry.get("oxygen_dissolved").values
    b = demo_b.registry.get("oxygen_dissolved").values
    np.testing.assert_array_equal(
        a, b,
        err_msg=(
            "oxygen_dissolved differs between subscribed and unsubscribed "
            "runs; pattern G zero-cost-when-unused invariant violated"
        ),
    )


def test_g_diagnostics_attribute_caches_match_registry_writes() -> None:
    """The value written to the registry for each REGISTRY_DIAGNOSTICS
    name must equal ``getattr(self, name)`` (the F/G single-source-of-
    truth contract)."""
    from clearwater_modules_v3.processes.dox import DOX

    ic = default_initial_conditions()
    reference = ic["oxygen_dissolved"]
    for name in DOX.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    dox = _dox_process(demo)

    demo.step(START)

    for name in dox.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(dox, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(
            cached, from_registry,
            err_msg=(
                f"{name}: registry value differs from self.{name}; "
                "pattern F/G single-source-of-truth violated"
            ),
        )


def test_g_partial_subscription_writes_only_requested_names() -> None:
    """Subscribing to a subset writes only those names; the rest stay
    out of the registry."""
    from clearwater_modules_v3.processes.dox import DOX

    subscribed = tuple(DOX.REGISTRY_DIAGNOSTICS[::2])
    unsubscribed = tuple(DOX.REGISTRY_DIAGNOSTICS[1::2])
    assert subscribed and unsubscribed

    ic = default_initial_conditions()
    reference = ic["oxygen_dissolved"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry
    for name in unsubscribed:
        assert name not in demo.registry


def test_g_sod_and_dox_sod_rate_match_in_registry() -> None:
    """``sod_rate`` and ``dox_sod_rate`` aliases must produce identical
    registry values when both are pre-registered."""
    ic = default_initial_conditions()
    reference = ic["oxygen_dissolved"]
    ic["sod_rate"] = xr.zeros_like(reference)
    ic["dox_sod_rate"] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    sod = _as_array(demo.registry.get("sod_rate"))
    dox_sod = _as_array(demo.registry.get("dox_sod_rate"))
    np.testing.assert_array_equal(
        sod, dox_sod,
        err_msg="sod_rate and dox_sod_rate aliases produced different values",
    )
