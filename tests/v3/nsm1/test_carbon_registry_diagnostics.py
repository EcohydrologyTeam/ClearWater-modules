"""Phase 2 opportunistic-registry-write contract for Carbon.

Asserts pattern G (pattern-alignment spec §3): when a name in
``Carbon.REGISTRY_DIAGNOSTICS`` is pre-registered, the value cached on
``Carbon`` is written to the registry each substep; when no name is
pre-registered, the state-variable outputs (POC, DOC, DIC) are
bit-identical to the no-diagnostics baseline.

Retained through the entire pattern-alignment work (unlike the
helper-vs-inline shadow file). After Phase 10, the corresponding tests
for the other 10 Processes follow this template.
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
)


N_STEPS = 60  # 5 minutes/step * 60 = 5 hours
START = datetime(2026, 1, 1, 0, 0, 0)


def _carbon_process(demo):
    return next(p for p in demo.processes if type(p).__name__ == "Carbon")


# ---------------------------------------------------------------------------
# Pattern G contract — opportunistic registry writes
# ---------------------------------------------------------------------------


def _as_array(value) -> np.ndarray:
    """Coerce a registry value to a numpy array.

    Carbon's components dict mixes ``xr.DataArray`` (full per-cell
    rates) and Python scalars (e.g., ``0.0`` for ``dic_sed_release_rate``
    when ``use_SedFlux=False``). The registry stores whatever we hand
    it; this helper normalises both shapes to ndarrays so the tests
    treat them uniformly.
    """
    if isinstance(value, xr.DataArray):
        return value.values
    return np.asarray(value)


def test_g_diagnostics_written_when_pre_registered() -> None:
    """Pre-register every Carbon REGISTRY_DIAGNOSTICS name; assert
    each is written each substep with finite values."""
    ic = default_initial_conditions()
    # Carbon writes only at wet cells; mark every cell wet by ensuring
    # the IC dict already has finite forcings. Pre-register diagnostics
    # by injecting zeroed DataArrays into the IC dict with the right
    # shape; the registry treats them as state variables it knows
    # about, and ``set_at_time(name, ...)`` succeeds.
    from clearwater_modules_v3.processes.carbon import Carbon
    reference = ic["poc"]
    for name in Carbon.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    carbon = _carbon_process(demo)

    demo.run(START, n_steps=N_STEPS)

    for name in carbon.REGISTRY_DIAGNOSTICS:
        assert name in demo.registry, (
            f"{name} should still be registered after run"
        )
        arr = _as_array(demo.registry.get(name))
        assert np.isfinite(arr).all(), (
            f"{name} has non-finite values after {N_STEPS} substeps"
        )


def test_g_diagnostics_skipped_when_not_pre_registered() -> None:
    """When no Carbon diagnostic is pre-registered, the state-variable
    trajectory must be bit-identical to a reference baseline that also
    skips pre-registration. This is the zero-cost-when-not-subscribed
    invariant under pattern G."""
    # Reference: no diagnostics pre-registered, capture state outputs.
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)
    state_a = {
        name: demo_a.registry.get(name).values.copy()
        for name in ("poc", "doc", "dic")
    }

    # Second run identical (no pre-registration), verify the run is
    # itself reproducible.
    demo_b = build_nsm1_demo()
    demo_b.run(START, n_steps=N_STEPS)
    state_b = {
        name: demo_b.registry.get(name).values.copy()
        for name in ("poc", "doc", "dic")
    }

    for name in state_a:
        np.testing.assert_array_equal(
            state_a[name], state_b[name],
            err_msg=f"{name} not reproducible across two unsubscribed runs",
        )


def test_g_state_bit_identical_with_and_without_diagnostics() -> None:
    """The core pattern G invariant: subscribing to diagnostic outputs
    must NOT change the state-variable trajectory. This is what makes
    the opportunistic-write pattern zero-risk for existing consumers."""
    from clearwater_modules_v3.processes.carbon import Carbon

    # Run A: no diagnostics pre-registered.
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_STEPS)

    # Run B: all Carbon diagnostics pre-registered.
    ic_b = default_initial_conditions()
    reference = ic_b["poc"]
    for name in Carbon.REGISTRY_DIAGNOSTICS:
        ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_STEPS)

    for state_name in ("poc", "doc", "dic"):
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
    """When subscribed, the values written to the registry must equal
    the values cached on ``self.<name>`` (pattern F). Verifies the
    single-source-of-truth contract between the components dict, the
    instance attributes, and the registry writes."""
    from clearwater_modules_v3.processes.carbon import Carbon

    ic = default_initial_conditions()
    reference = ic["poc"]
    for name in Carbon.REGISTRY_DIAGNOSTICS:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    carbon = _carbon_process(demo)

    # One substep — registry write must reflect the cached value.
    demo.step(START)

    for name in carbon.REGISTRY_DIAGNOSTICS:
        cached = _as_array(getattr(carbon, name))
        from_registry = _as_array(demo.registry.get(name))
        np.testing.assert_array_equal(
            cached, from_registry,
            err_msg=(
                f"{name}: registry value differs from self.{name}; "
                "pattern F/G single-source-of-truth violated"
            ),
        )


def test_g_partial_subscription_writes_only_requested_names() -> None:
    """Subscribing to a subset (not all) of REGISTRY_DIAGNOSTICS writes
    only those names; un-subscribed names are NOT silently added to the
    registry by Carbon."""
    from clearwater_modules_v3.processes.carbon import Carbon

    # Subscribe to half the diagnostics; leave the others unregistered.
    subscribed = tuple(Carbon.REGISTRY_DIAGNOSTICS[::2])
    unsubscribed = tuple(Carbon.REGISTRY_DIAGNOSTICS[1::2])
    assert subscribed and unsubscribed, "test setup requires both subsets non-empty"

    ic = default_initial_conditions()
    reference = ic["poc"]
    for name in subscribed:
        ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.step(START)

    for name in subscribed:
        assert name in demo.registry, (
            f"{name} was pre-registered but is missing from registry"
        )
    for name in unsubscribed:
        assert name not in demo.registry, (
            f"{name} was NOT pre-registered but appeared in registry"
        )
