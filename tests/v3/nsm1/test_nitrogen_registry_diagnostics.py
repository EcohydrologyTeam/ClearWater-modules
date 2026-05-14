"""Phase 4 opportunistic-registry-write contract for Nitrogen.

Asserts pattern G (pattern-alignment spec §3): when a name in
``Nitrogen.REGISTRY_DIAGNOSTICS`` is pre-registered, the value cached
on ``Nitrogen`` is written to the registry each substep; when no name
is pre-registered, the NH4 / NO3 / OrgN trajectories are bit-identical
to the no-diagnostics baseline.

Mirrors ``test_carbon_registry_diagnostics.py`` /
``test_dox_registry_diagnostics.py``; retained through Phase 10.
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


def _nitrogen_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "Nitrogen")


def _as_array(value) -> np.ndarray:
    """Coerce a registry value to a numpy array."""
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    """Pre-register every Nitrogen REGISTRY_DIAGNOSTICS name; assert
    each is written each substep with finite values."""
    from clearwater_modules_v3.processes.nitrogen import Nitrogen

    ic = default_initial_conditions()
    reference = ic["ammonium"]
    for name in Nitrogen.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    nitrogen = _nitrogen_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in nitrogen.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry, (
            f"{name} should still be registered after run"
        )
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all(), (
            f"{name} has non-finite values after {N_STEPS} substeps"
        )


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    """Reproducibility under no-subscription."""
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    state_a = {
        name: demo_a.registry.get(name).values.copy()
        for name in ("ammonium", "nitrate", "organic_nitrogen")
    }

    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    state_b = {
        name: demo_b.registry.get(name).values.copy()
        for name in ("ammonium", "nitrate", "organic_nitrogen")
    }

    for name in state_a:
        np.testing.assert_array_equal(state_a[name], state_b[name])


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    """Pattern G zero-cost-when-unused: subscribing to Nitrogen
    diagnostics does NOT change NH4 / NO3 / OrgN trajectories."""
    from clearwater_modules_v3.processes.nitrogen import Nitrogen

    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    ic_b = default_initial_conditions()
    reference = ic_b["ammonium"]
    for name in Nitrogen.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    for state_name in ("ammonium", "nitrate", "organic_nitrogen"):
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
    """Cached ``self.<name>`` matches ``registry.get(name)`` after each
    substep. Specifically pins the ``nitrification_flux_rate`` /
    ``denitrification_flux_rate`` preserved-name contract: the values
    DOX / Alkalinity / N2 read from the sibling Nitrogen Process via
    ``getattr`` are the same values exposed to the registry."""
    from clearwater_modules_v3.processes.nitrogen import Nitrogen

    ic = default_initial_conditions()
    reference = ic["ammonium"]
    for name in Nitrogen.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    nitrogen = _nitrogen_process(demo)

    demo.step(START)

    for name in nitrogen.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(nitrogen, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(
            cached, from_registry,
            err_msg=(
                f"{name}: registry value differs from self.{name}; "
                "pattern F/G single-source-of-truth violated"
            ),
        )


def test_g_partial_subscription_writes_only_requested_names() -> None:
    """Subscribing to a subset writes only those names."""
    from clearwater_modules_v3.processes.nitrogen import Nitrogen

    subscribed = tuple(Nitrogen.REGISTRY_DIAGNOSTICS[::2])
    unsubscribed = tuple(Nitrogen.REGISTRY_DIAGNOSTICS[1::2])
    assert subscribed and unsubscribed

    ic = default_initial_conditions()
    reference = ic["ammonium"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry
    for name in unsubscribed:
        assert name not in demo.registry


def test_g_preserved_names_are_consumer_visible() -> None:
    """``nitrification_flux_rate`` and ``denitrification_flux_rate`` are
    the consumer-facing attribute names. After ``run``, both must be
    populated on ``self`` with finite, non-negative values regardless
    of registry pre-registration (pattern F is unconditional)."""
    demo = build_nsm1_demo()
    nitrogen = _nitrogen_process(demo)

    demo.step(START)

    nitr = _as_array(getattr(nitrogen, "nitrification_flux_rate"))
    denit = _as_array(getattr(nitrogen, "denitrification_flux_rate"))
    assert np.all(nitr >= 0.0), "nitrification_flux_rate has negative values"
    assert np.all(denit >= 0.0), "denitrification_flux_rate has negative values"
    assert np.isfinite(nitr).all()
    assert np.isfinite(denit).all()
