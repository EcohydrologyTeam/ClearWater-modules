"""Tier 1 closed-system mass-conservation test for v3 NSM1 Alkalinity.

Phase 6 (v3 NSM1 design spec, Section 11 Phase 6, Section 5 Alkalinity
design notes, Section 14 resolved Q "Alkalinity simple-tracer"):
asserts that the v3-native ``Alkalinity`` Process is conservative when
all source/sink coupling (Nitrogen + algae) is disabled. Under these
conditions:

* The ``alkalinity`` state must equal its initial value at every cell
  to roundoff (``rtol=1e-12``).
* The clip-with-log diagnostics must remain empty
  (``diagnostics.clip_events == {}``).

How the closed system is constructed:

* No FloatingAlgae / BenthicAlgae Processes are wired up, so the algae
  coupling helpers all return 0 (``self.use_floating_algae == False``,
  ``self.use_benthic_algae == False``).
* No Nitrogen process is wired up, so the nitrification sink and
  denitrification source are both identically zero
  (``self.use_nitrogen == False``).
* Per Section 14 resolved Q, v3 1.0.0 has no carbonate equilibrium /
  pH solver. Alk is a pure source/sink tracer; with all sources/sinks
  zeroed it is invariant.

Optional second test (``test_tier1_alk_with_nitrogen_zeroed_kinetics``):
Wires Nitrogen into Alkalinity with all kinetic-rate constants zeroed.
This verifies the production rate-DAG contract: when Nitrogen.run caches
zero ``nitrification_flux_rate`` and zero ``denitrification_flux_rate``,
the Alkalinity coupling helpers produce zero rate contributions and Alk
remains invariant.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.alkalinity import Alkalinity
from clearwater_modules_v3.processes.nitrogen import Nitrogen
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


def test_tier1_alkalinity_conservation_closed_system_loss_disabled(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system Alkalinity conservation when all source/sink
    coupling is disabled. With no algae, no nitrogen process, Alk is
    invariant under the simple-tracer model (resolved Q in Section 14
    of the design spec).

    Setup:
    * 5-cell mesh; ``initial_state_5cell`` initial conditions
      (``alkalinity`` initialized to ``[100.0, 105.0, 110.0, 115.0, 120.0]``
      mg-CaCO3/L).
    * No Nitrogen, no algae Processes registered with the model, so
      ``self.use_nitrogen == False``, ``self.use_floating_algae == False``,
      ``self.use_benthic_algae == False`` after ``__init__`` and every
      coupling helper returns 0.

    Expected:
    * ``alkalinity_final == alkalinity_initial`` per cell (rtol=1e-12).
    * ``diagnostics.clip_events == {}``.
    """
    start, end, time_step = closed_system_time_window

    alkalinity = Alkalinity(time_step=time_step)
    diagnostics = Diagnostics()
    alkalinity.diagnostics = diagnostics

    # No siblings wired up in this Tier 1 harness; verify the coupling
    # flags default to False so all source/sink helpers return 0.
    assert alkalinity.use_nitrogen is False
    assert alkalinity.use_floating_algae is False
    assert alkalinity.use_benthic_algae is False
    assert alkalinity.nitrogen_process is None
    assert alkalinity.floating_algae_process is None
    assert alkalinity.benthic_algae_process is None

    # Snapshot the initial Alk state for the per-cell equality check.
    alk_initial = in_memory_registry.get_at_time("alkalinity", start).copy()

    # Run 100 substeps. ``Alkalinity.run`` mutates the registry in place.
    current_time = start
    while current_time < end:
        alkalinity.run(current_time, in_memory_registry)
        current_time += time_step

    # Tier 1 invariant 1: Alk state per-cell invariance under
    # zeroed-coupling closed-system conditions.
    alk_final = in_memory_registry.get_at_time("alkalinity", end)
    np.testing.assert_allclose(
        alk_final.values,
        alk_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system Alkalinity invariance failed under "
            "no-Nitrogen, no-algae conditions. "
            f"initial={alk_initial.values!r}, "
            f"final={alk_final.values!r}, "
            f"absolute drift={(alk_final.values - alk_initial.values)!r}"
        ),
    )

    # Tier 1 invariant 2: no clipping under closed-system conditions.
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Tier 1 Alkalinity "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


def test_tier1_alk_with_nitrogen_zeroed_kinetics(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system Alk invariance with Nitrogen wired but kinetics zeroed.

    Mirrors the ``test_tier1_n2_with_nitrogen_zeroed_kinetics`` pattern
    (Item 1 contract): when Nitrogen is wired into Alkalinity (so Alk
    reads ``nitrogen_process.nitrification_flux_rate`` and
    ``denitrification_flux_rate``) but every Nitrogen kinetic-rate
    constant is zero, the cached fluxes are identically zero and Alk
    remains invariant.

    Setup:
    * 5-cell mesh, ``initial_state_5cell`` initial conditions.
    * Nitrogen instance with all kinetic-rate constants set to 0.0.
    * Manually wire ``alkalinity.use_nitrogen = True`` and
      ``alkalinity.nitrogen_process = nitrogen`` (no Model harness).
    * Run order per substep: Nitrogen.run -> Alkalinity.run.

    Expected:
    * ``alk_final == alk_initial`` per cell (rtol=1e-12).
    * ``nitrogen.nitrification_flux_rate`` and
      ``nitrogen.denitrification_flux_rate`` are zero everywhere after
      ``Nitrogen.run``.
    * Alkalinity ``diagnostics.clip_events == {}``.
    """
    start, end, time_step = closed_system_time_window

    # Build a Nitrogen with all kinetic-rate constants zeroed so the
    # cached fluxes are identically zero. Mirrors the v3 N2 Tier 1 test.
    nitrogen = Nitrogen(
        parameters={
            "vson_20": 0.0,
            "rnh4_20": 0.0,
            "vno3_20": 0.0,
            "knit_20": 0.0,
            "kdnit_20": 0.0,
            "kon_20": 0.0,
        },
        time_step=time_step,
        nitrification_rate=0.0,
        denitrification_rate=0.0,
        sediment_ammonium_release_rate=0.0,
        sediment_denitrification_rate=0.0,
        ammonium_decay_rate=0.0,
        death_rate=0.0,
        settling_velocity=0.0,
    )
    nitrogen.use_nitrate = True
    nitrogen.use_ammonium = True
    nitrogen.use_floating_algae = False
    nitrogen.use_benthic_algae = False

    alkalinity = Alkalinity(time_step=time_step)
    diagnostics = Diagnostics()
    alkalinity.diagnostics = diagnostics

    # Manually wire Nitrogen onto Alkalinity (what ``init_process`` would
    # do under a v3 Model). Confirm both halves of the contract.
    alkalinity.use_nitrogen = True
    alkalinity.nitrogen_process = nitrogen
    assert alkalinity.use_nitrogen is True
    assert alkalinity.nitrogen_process is nitrogen

    alk_initial = in_memory_registry.get_at_time("alkalinity", start).copy()

    current_time = start
    while current_time < end:
        nitrogen.run(current_time, in_memory_registry)
        # Sanity-check both flux caches: with zeroed kinetic constants,
        # both fluxes must be zero everywhere.
        for cache_name in ("nitrification_flux_rate",
                           "denitrification_flux_rate"):
            flux = getattr(nitrogen, cache_name)
            if hasattr(flux, "values"):
                np.testing.assert_allclose(
                    flux.values,
                    np.zeros_like(flux.values),
                    atol=0.0,
                    err_msg=(
                        f"Nitrogen.{cache_name} non-zero under zeroed "
                        f"kinetics: {flux.values!r}"
                    ),
                )
            else:
                assert flux == 0.0, (
                    f"Nitrogen.{cache_name} non-zero under zeroed "
                    f"kinetics: {flux!r}"
                )
        alkalinity.run(current_time, in_memory_registry)
        current_time += time_step

    alk_final = in_memory_registry.get_at_time("alkalinity", end)
    np.testing.assert_allclose(
        alk_final.values,
        alk_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system Alkalinity invariance failed with Nitrogen "
            "wired and zeroed kinetics. "
            f"initial={alk_initial.values!r}, "
            f"final={alk_final.values!r}, "
            f"absolute drift={(alk_final.values - alk_initial.values)!r}"
        ),
    )

    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Alkalinity-with-Nitrogen "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


# ---------------------------------------------------------------------------
# Alkalinity instantiation smoke tests
# ---------------------------------------------------------------------------


def test_alkalinity_instantiates_with_defaults() -> None:
    """``Alkalinity()`` constructs cleanly with no arguments and pulls a
    composed DEFAULTS dict (from alkalinity / global_parameters / algae /
    balgae) onto the instance.
    """
    alk = Alkalinity()
    # Pulled from alkalinity (six stoichiometric ratios from Phase 1.2).
    assert alk.r_alkaa == 14.0 / 106.0 / 12.0 / 1000.0
    assert alk.r_alkan == 18.0 / 106.0 / 12.0 / 1000.0
    assert alk.r_alkn == 2.0 / 14.0 / 1000.0
    assert alk.r_alkden == 4.0 / 14.0 / 1000.0
    assert alk.r_alkba == 14.0 / 106.0 / 12.0 / 1000.0
    assert alk.r_alkbn == 18.0 / 106.0 / 12.0 / 1000.0
    # Pulled from global_parameters (feature flags).
    assert alk.use_NH4 is True
    assert alk.use_NO3 is True
    assert alk.use_Algae is True
    assert alk.use_Balgae is True
    # Pulled from algae / balgae.
    assert alk.AWc > 0.0
    assert alk.BWc > 0.0
    assert alk.Fb > 0.0
    # Coupling defaults to disconnected.
    assert alk.use_nitrogen is False
    assert alk.use_floating_algae is False
    assert alk.use_benthic_algae is False
    assert alk.nitrogen_process is None
    assert alk.floating_algae_process is None
    assert alk.benthic_algae_process is None
    # Step-scoped caches start zeroed.
    assert alk.alk_nitrification_rate == 0.0
    assert alk.alk_denitrification_rate == 0.0
    assert alk.alk_algal_growth_rate == 0.0
    assert alk.alk_algal_respiration_rate == 0.0
    assert alk.alk_benthic_algae_growth_rate == 0.0
    assert alk.alk_benthic_algae_respiration_rate == 0.0
    assert alk.alk_rate == 0.0


def test_alkalinity_accepts_parameter_override() -> None:
    """``Alkalinity(parameters={...})`` overrides composed defaults and
    warns on unknown keys (no exception raised)."""
    alk = Alkalinity(
        parameters={
            "r_alkn": 0.0,            # disable nitrification stoich
            "r_alkden": 0.0,          # disable denitrification stoich
            "use_Algae": False,
        },
    )
    assert alk.r_alkn == 0.0
    assert alk.r_alkden == 0.0
    assert alk.use_Algae is False
    # Unchanged composed defaults.
    assert alk.r_alkaa == 14.0 / 106.0 / 12.0 / 1000.0
    assert alk.use_NH4 is True
