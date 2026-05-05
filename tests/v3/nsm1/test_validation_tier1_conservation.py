"""Tier 1 closed-system mass-conservation tests for v3 NSM1.

The Tier 1 contract (design spec Section 9, Section 14 resolved Q7):

    Closed system + no boundaries + no settling + balanced source/sink
    pairs --> total mass of N, P, C, O2-equivalents, Alk constant to
    floating-point roundoff AND ``diagnostics.clip_events == {}``.

This module is the harness scaffold for Phase 1.4. Per-constituent
conservation tests are added in Phases 2-6 as the corresponding Process
classes land:

* Phase 2 (Nitrogen integrator fix + OrgN)        -- total-N test
* Phase 3 (POM, CBOD, N2/TDG, Pathogen)           -- N2 conservation,
                                                    CBOD-as-O2 contribution
* Phase 4 (Phosphorus)                            -- total-P test
* Phase 5 (Carbon, DOX)                           -- total-C, total-O2
                                                    Streeter-Phelps Tier 2
* Phase 6 (Alkalinity)                            -- total-Alk test

The single test in this Phase 1.4 file demonstrates the harness end-to-
end against the existing v2-overlay Nitrogen process. Because v2
``Nitrogen`` carries the multiplicative-integrator bug (NH4 update at
``nitrogen.py:101`` and NO3 update at ``nitrogen.py:112``), the test
is marked ``xfail`` and is expected to flip to passing the moment
Phase 2 lands. That makes this test a regression-test entry point
that signals when the bug is fixed.

Why xfail and not skip: an xfail that *unexpectedly passes* (XPASS)
shows up in the pytest summary. If a future change accidentally
papers over the integrator bug without correctly fixing it, the
XPASS will be the first thing the next test run advertises.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.nitrogen import Nitrogen
from clearwater_modules_v3.utils.numerics import Diagnostics, clip_negative_state

# These helpers live in ``conftest.py`` and are auto-discovered as
# fixtures, but they are also imported directly so the test module can
# call them as plain functions.
from .conftest import (
    InMemoryRegistry,
    total_n,
    total_p,
    total_c,
    total_o2_equivalents,
    total_alkalinity,
)


# ---------------------------------------------------------------------------
# Tier 1 stub test (Phase 1.4)
# ---------------------------------------------------------------------------


def test_tier1_total_n_conservation_closed_system_nitrogen_only(
    in_memory_registry: InMemoryRegistry,
    closed_system_config,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Total nitrogen is conserved across 100 substeps of v3 Nitrogen.

    Setup:
    * 5-cell mesh, ``initial_state_5cell`` initial conditions
    * Floating algae and benthic algae are NOT instantiated, so the
      Nitrogen process's ``use_floating_algae`` / ``use_benthic_algae``
      branches return zero contributions. Mass-changing pathways are
      NH4 <-> NO3 (nitrification + denitrification), OrgN -> NH4
      hydrolysis (mass-conserving within the N pool), and sediment
      fluxes (zero in the closed-system config).
    * ``vson_20=0`` disables OrgN settling (mass-removing) so the
      closed-system invariant holds.
    * ``ammonium_decay_rate=0`` removes the legacy ``ammonium_decay_nitrate``
      conversion that is a non-conservative artifact of the v2
      sub-rate methods (mirrors the v1 behavior under
      ``ammonium_decay_rate==0``).

    Expected (Phase 2.B and beyond):
    * ``total_n(state_final) == total_n(state_initial)`` to roundoff
      (rtol=1e-12)
    * ``diagnostics.clip_events == {}``
    """
    start, end, time_step = closed_system_time_window

    # Build a v3 Nitrogen with v1-aligned defaults plus closed-system
    # overrides. ``parameters`` overrides drive the v3 NITROGEN_DEFAULTS
    # entries; ``vson_20=0`` makes the OrgN reservoir closed (no settling
    # to bed). ``kon_20`` is left at its default so the OrgN -> NH4
    # hydrolysis pathway exercises end-to-end, but it is mass-conserving
    # within the total-N pool so the Tier 1 invariant still holds.
    nitrogen = Nitrogen(
        parameters={
            "vson_20": 0.0,        # OrgN settling disabled (closed system)
            "rnh4_20": closed_system_config.rnh4_20,
            "vno3_20": closed_system_config.vno3_20,
        },
        time_step=time_step,
        # Process kinetics rates -- v1-aligned defaults from
        # clearwater_modules_v3.parameters.nitrogen.
        # nitrification_rate is non-zero (NH4 <-> NO3 is mass-conserving
        # within total-N). denitrification_rate is set to 0 because
        # denitrification routes N out of NO3 and into N2; the N2
        # Process is Phase 3, so for Phase 2.B's Tier 1 we close the
        # NO3-loss pathway by zeroing the rate.
        nitrification_rate=0.1,
        nitrification_theta=1.083,
        denitrification_rate=0.0,
        denitrification_theta=1.045,
        # Sediment exchange disabled by closed-system config
        sediment_ammonium_release_rate=closed_system_config.rnh4_20,
        sediment_ammonium_release_theta=1.0,
        sediment_denitrification_rate=closed_system_config.vno3_20,
        sediment_denitrification_theta=1.0,
        # ammonium_decay_nitrate is a v2-side artifact; v1 has no such
        # term. Set to 0 to keep the available kinetics pure NH4 <-> NO3
        # plus the new OrgN <-> NH4 hydrolysis.
        ammonium_decay_rate=0.0,
        ammonium_decay_theta=1.0,
        nitrification_oxygen_inhibition_factor=0.6,
        floating_algae_preference_factor=0.5,
        settling_velocity=closed_system_config.vsop,
        death_rate=0.0,
        float_algea_faction_uptake_from_nitrate=0.0,
    )
    # The v2 Nitrogen.run() depends on these flags being set; normally
    # they are populated by ``init_process(model, registry)`` when a
    # Model is wired up. For the harness we set them directly so we can
    # exercise ``run`` in a Model-free unit-test mode.
    nitrogen.use_nitrate = True
    nitrogen.use_ammonium = True
    nitrogen.use_floating_algae = False
    nitrogen.use_benthic_algae = False

    # Snapshot the closed-system initial total-N pool.
    n_initial = float(total_n(in_memory_registry).values)

    # Run 100 substeps. v2 Nitrogen.run mutates the registry in place.
    diagnostics = Diagnostics()
    current_time = start
    step_index = 0
    while current_time < end:
        nitrogen.run(current_time, in_memory_registry)
        # Apply the v3 clip-with-log contract to the two state variables
        # the process touches. The v2 Nitrogen.run already does an
        # ``xr.where(X < 0, 0, X)`` clamp inline; in a v3-native
        # Nitrogen this would be replaced by ``clip_negative_state``,
        # at which point ``diagnostics`` would be the canonical record
        # of clip events. We invoke ``clip_negative_state`` here too so
        # the harness covers the diagnostic-checking branch even
        # against the v2 overlay.
        for state_name in ("ammonium", "nitrate"):
            current = in_memory_registry.get_at_time(state_name, current_time)
            clipped = clip_negative_state(
                current, state_name, diagnostics, step=step_index
            )
            in_memory_registry.set_at_time(state_name, current_time, clipped)
        current_time += time_step
        step_index += 1

    # Tier 1 invariant 1: total-N constant to roundoff.
    n_final = float(total_n(in_memory_registry).values)
    np.testing.assert_allclose(
        n_final,
        n_initial,
        rtol=1e-12,
        err_msg=(
            "Closed-system total-N mass conservation failed. "
            f"initial={n_initial!r}, final={n_final!r}, "
            f"absolute drift={(n_final - n_initial)!r}"
        ),
    )

    # Tier 1 invariant 2: no clipping under closed-system + physically
    # reasonable initial conditions. A clip event is the canonical
    # diagnostic that either the parameters or the test case is
    # malformed (design spec Section 14 Q7).
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Tier 1 conditions: "
        f"{diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


# ---------------------------------------------------------------------------
# Helper-function smoke tests (Phase 1.4)
# ---------------------------------------------------------------------------
# These tests verify the conservation helpers themselves return
# sensible totals against a known initial state. They are not
# Tier 1 conservation assertions per se; they are sanity checks that
# the helpers don't have arithmetic typos that would mask a real
# conservation failure when Phase 2-6 tests are added.


def test_total_n_helper_returns_finite_positive(in_memory_registry: InMemoryRegistry) -> None:
    """``total_n`` over the 5-cell initial state returns a finite positive 0-d DataArray."""
    result = total_n(in_memory_registry)
    assert result.ndim == 0
    assert np.isfinite(float(result.values))
    assert float(result.values) > 0.0


def test_total_p_helper_returns_finite_positive(in_memory_registry: InMemoryRegistry) -> None:
    result = total_p(in_memory_registry)
    assert result.ndim == 0
    assert np.isfinite(float(result.values))
    assert float(result.values) > 0.0


def test_total_c_helper_returns_finite_positive(in_memory_registry: InMemoryRegistry) -> None:
    result = total_c(in_memory_registry)
    assert result.ndim == 0
    assert np.isfinite(float(result.values))
    assert float(result.values) > 0.0


def test_total_o2_equivalents_helper_returns_finite_positive(
    in_memory_registry: InMemoryRegistry,
) -> None:
    result = total_o2_equivalents(in_memory_registry)
    assert result.ndim == 0
    assert np.isfinite(float(result.values))
    assert float(result.values) > 0.0


def test_total_alkalinity_helper_returns_finite_positive(
    in_memory_registry: InMemoryRegistry,
) -> None:
    result = total_alkalinity(in_memory_registry)
    assert result.ndim == 0
    assert np.isfinite(float(result.values))
    assert float(result.values) > 0.0


def test_helpers_silently_skip_missing_reservoirs() -> None:
    """A registry with only NH4 should still produce a finite total-N
    (algal contribution drops out, but the NH4 pool sums correctly)."""
    registry = InMemoryRegistry()
    registry.register("ammonium", xr.DataArray(np.array([1.0, 2.0, 3.0]), dims="cell"))
    result = total_n(registry)
    np.testing.assert_allclose(float(result.values), 6.0, rtol=1e-12)


def test_diagnostics_import_smoke() -> None:
    """Smoke: ``Diagnostics`` constructs and starts empty.

    The Tier 1 invariant ``diagnostics.clip_events == {}`` depends on
    a freshly-constructed ``Diagnostics`` defaulting to an empty
    ``clip_events`` dict; this test pins that contract.
    """
    diag = Diagnostics()
    assert diag.clip_events == {}
    assert diag.clip_log == []
    assert diag.detail_limit_per_call == 10
