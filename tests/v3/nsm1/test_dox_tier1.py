"""Tier 1 closed-system mass-conservation test for v3 NSM1 DOX.

Phase 5.B (v3 NSM1 design spec, Section 11 Phase 5, Section 5 DOX
design notes): asserts that the v3-native ``DOX`` Process is conservative
(DOX invariant) when atmospheric exchange, SOD, and every coupled
sibling source/sink (algae, nitrogen, carbon, CBOD) are disabled. Under
these conditions:

* The DOX state must equal its initial value at every cell to roundoff
  (``rtol=1e-12``).
* The clip-with-log diagnostics must remain empty
  (``diagnostics.clip_events == {}``); a clip event under closed-system
  conditions signals either a malformed test or an integrator bug.

How the closed system is constructed:

* Atmospheric exchange is disabled by setting ``kah_20_user=0`` and
  ``kaw_20_user=0`` AND selecting the user-defined menu options
  (``hydraulic_reaeration_option=1``, ``wind_reaeration_option=1``) so
  the menu picks zero. This makes ``ka_tc == 0`` and zeroes the
  ``ka_tc * (O2sat - DOX)`` reaeration term.
* SOD is disabled by setting ``SOD_20=0``.
* The nitrification O2 sink is disabled by setting ``use_NH4=False``
  (which short-circuits the nitrification flux to zero regardless of
  NH4 / DOX state).
* No FloatingAlgae, BenthicAlgae, Nitrogen, Carbon, or CBOD Processes
  are wired into the model — DOX's coupling flags default to False, so
  every algal / DOC / CBOD source-or-sink path is identically zero.

DOX-coupled conservation tests (DOX + CBOD == invariant under closed
system, DOX + DOC == invariant, etc.) are the responsibility of Phase
5.5 integration once all 7 Processes are registered.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

from clearwater_modules_v3.processes.dox import DOX
from clearwater_modules_v3.processes.nitrogen import Nitrogen
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


def test_tier1_dox_conservation_closed_system_loss_disabled(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system DOX conservation when all sources and sinks are
    disabled. With no coupled Processes, no atmospheric exchange, and
    no SOD, DOX is invariant.

    Setup:
    * 5-cell mesh; ``initial_state_5cell`` initial conditions
      (``oxygen_dissolved`` initialized to ``[8.0, 8.5, 9.0, 9.5, 10.0]``
      mg-O2/L).
    * ``kah_20_user=0``, ``kaw_20_user=0`` AND
      ``hydraulic_reaeration_option=1``, ``wind_reaeration_option=1``
      together force ``ka_tc == 0``.
    * ``SOD_20=0`` zeroes the SOD sink.
    * ``use_NH4=False`` zeroes the nitrification sink.
    * No coupled Processes are wired up; ``self.use_floating_algae``,
      ``self.use_benthic_algae``, ``self.use_nitrogen``,
      ``self.use_carbon``, ``self.use_cbod`` all default to False, so
      every coupled source/sink term is identically zero.

    Expected:
    * ``dox_final == dox_initial`` per cell (rtol=1e-12).
    * ``diagnostics.clip_events == {}``.
    """
    start, end, time_step = closed_system_time_window

    dox = DOX(
        parameters={
            # Disable atmospheric exchange.
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            # Disable SOD.
            "SOD_20": 0.0,
            "SOD_theta": 1.0,
            # Disable nitrification O2 sink at the gate (also zeroes
            # the term even if NH4 happens to be present in the
            # registry fixture).
            "use_NH4": False,
            # Standard sea-level pressure (the corrected v3 default).
            "pressure_mb": 1013.25,
        },
        time_step=time_step,
    )
    # Make the test own the diagnostics reference.
    diagnostics = Diagnostics()
    dox.diagnostics = diagnostics

    # Verify the coupling flags default to False so every coupled
    # source/sink path is identically zero in ``run``.
    assert dox.use_floating_algae is False
    assert dox.use_benthic_algae is False
    assert dox.use_nitrogen is False
    assert dox.use_carbon is False
    assert dox.use_cbod is False
    assert dox.floating_algae_process is None
    assert dox.benthic_algae_process is None
    assert dox.nitrogen_process is None
    assert dox.carbon_process is None
    assert dox.cbod_process is None

    # Snapshot the initial DOX state for the per-cell equality check.
    dox_initial = in_memory_registry.get_at_time("oxygen_dissolved", start).copy()

    # Run 100 substeps. ``DOX.run`` mutates the registry in place.
    current_time = start
    while current_time < end:
        dox.run(current_time, in_memory_registry)
        current_time += time_step

    # Tier 1 invariant 1: DOX state per-cell invariance under
    # zeroed-flux closed-system conditions.
    dox_final = in_memory_registry.get_at_time("oxygen_dissolved", end)
    np.testing.assert_allclose(
        dox_final.values,
        dox_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system DOX invariance failed under "
            "kah_20_user=kaw_20_user=0, SOD_20=0, use_NH4=False, and "
            "no coupled Processes. "
            f"initial={dox_initial.values!r}, "
            f"final={dox_final.values!r}, "
            f"absolute drift={(dox_final.values - dox_initial.values)!r}"
        ),
    )

    # Tier 1 invariant 2: no clipping under closed-system conditions.
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Tier 1 DOX "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


# ---------------------------------------------------------------------------
# DOX instantiation smoke tests
# ---------------------------------------------------------------------------


def test_dox_instantiates_with_defaults() -> None:
    """``DOX()`` constructs cleanly with no arguments and pulls a composed
    DEFAULTS dict (from dox / global_parameters / global_vars / algae /
    balgae / nitrogen / carbon) onto the instance.
    """
    dox = DOX()
    # Pulled from DOX_DEFAULTS
    assert hasattr(dox, "ron")
    assert hasattr(dox, "KsSOD")
    assert hasattr(dox, "SOD_20")
    assert hasattr(dox, "SOD_theta")
    assert dox.kah_20_user == 0.0
    assert dox.kaw_20_user == 0.0
    assert dox.kah_theta == 1.024
    assert dox.kaw_theta == 1.024
    assert dox.hydraulic_reaeration_option == 1
    assert dox.wind_reaeration_option == 1
    # Pulled from global_parameters
    assert dox.pressure_mb == 1013.25
    assert dox.use_NH4 is True
    assert dox.use_DOC is True
    assert dox.use_Algae is True
    assert dox.use_Balgae is True
    # Pulled from global_vars (toy hydraulic forcings)
    assert dox.velocity == 1.0
    assert dox.flow == 2.0
    # Stoichiometric ratios composed from algae / balgae / carbon
    assert dox.AWc == 40.0           # rca
    assert dox.BWc == 40.0           # rcb
    assert dox.roc == 32.0 / 12.0
    # Nitrogen kinetics composed for the local nitrification flux.
    assert dox.KNR == 0.6
    assert dox.knit_20 == 0.1
    assert dox.knit_theta == 1.083
    # Coupling defaults to disconnected.
    assert dox.use_floating_algae is False
    assert dox.use_benthic_algae is False
    assert dox.use_nitrogen is False
    assert dox.use_carbon is False
    assert dox.use_cbod is False
    # Step-scoped caches start zeroed.
    assert dox.dox_sat == 0.0
    assert dox.atm_reaeration_rate == 0.0
    assert dox.dox_nitrification_rate == 0.0
    assert dox.dox_sod_rate == 0.0
    assert dox.dox_rate == 0.0
    # Diagnostics handle should be live.
    assert dox.diagnostics is not None
    assert dox.diagnostics.clip_events == {}


def test_dox_accepts_parameter_override() -> None:
    """``DOX(parameters={...})`` overrides composed defaults and warns on
    unknown keys (no exception raised)."""
    dox = DOX(
        parameters={
            "kah_20_user": 0.5,
            "SOD_20": 2.0,
            "pressure_mb": 1000.0,
        },
    )
    assert dox.kah_20_user == 0.5
    assert dox.SOD_20 == 2.0
    assert dox.pressure_mb == 1000.0
    # Unchanged composed defaults
    assert dox.kaw_20_user == 0.0
    assert dox.SOD_theta == 1.060


def test_tier1_dox_with_nitrogen_zeroed_kinetics(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system DOX invariance with Nitrogen wired but kinetics zeroed.

    Integration Item 1: verifies that when Nitrogen is wired into DOX
    (so DOX reads ``nitrogen_process.nitrification_flux_rate`` for its
    O2 sink) but every Nitrogen kinetic-rate constant is zero, the
    cached nitrification flux is identically zero and DOX remains
    invariant.

    This is the contract the production rate-DAG relies on: the DOX
    nitrification O2 sink is ``ron * nitrification_flux_rate``, where
    the flux is the *step-scoped* mg-N/L/d value cached by Nitrogen
    (NOT recomputed locally from NH4 / T / DOX). With zeroed rate
    constants the flux is zero, so DOX must be invariant to roundoff
    when all other O2 sources / sinks are also disabled.

    Setup:
    * 5-cell mesh, ``initial_state_5cell`` initial conditions.
    * Atmospheric exchange disabled (zeroed user reaeration menu).
    * SOD disabled (``SOD_20=0``).
    * Nitrogen instance with all kinetic-rate constants set to 0.0.
    * Manually wire ``dox.use_nitrogen = True`` and
      ``dox.nitrogen_process = nitrogen``.
    * Run order per substep: Nitrogen.run -> DOX.run.

    Expected:
    * ``dox_final == dox_initial`` per cell (rtol=1e-12).
    * ``nitrogen.nitrification_flux_rate`` is zero everywhere after
      ``run`` (sanity check on the new cache attribute).
    * DOX ``diagnostics.clip_events == {}``.
    """
    start, end, time_step = closed_system_time_window

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

    dox = DOX(
        parameters={
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "SOD_20": 0.0,
            "SOD_theta": 1.0,
            # NH4 enabled here; the gate that zeroes the sink is
            # the absent-Nitrogen / zeroed-flux path, not use_NH4.
            "use_NH4": True,
            "pressure_mb": 1013.25,
        },
        time_step=time_step,
    )
    diagnostics = Diagnostics()
    dox.diagnostics = diagnostics

    # Manually wire Nitrogen onto DOX (what ``init_process`` would do
    # under a v3 Model).
    dox.use_nitrogen = True
    dox.nitrogen_process = nitrogen
    assert dox.use_nitrogen is True
    assert dox.nitrogen_process is nitrogen

    dox_initial = in_memory_registry.get_at_time("oxygen_dissolved", start).copy()

    current_time = start
    while current_time < end:
        nitrogen.run(current_time, in_memory_registry)
        # Sanity-check the new flux cache.
        nitr_flux = nitrogen.nitrification_flux_rate
        if hasattr(nitr_flux, "values"):
            np.testing.assert_allclose(
                nitr_flux.values,
                np.zeros_like(nitr_flux.values),
                atol=0.0,
            )
        else:
            assert nitr_flux == 0.0
        dox.run(current_time, in_memory_registry)
        current_time += time_step

    dox_final = in_memory_registry.get_at_time("oxygen_dissolved", end)
    np.testing.assert_allclose(
        dox_final.values,
        dox_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system DOX invariance failed with Nitrogen wired and "
            "zeroed kinetics. "
            f"initial={dox_initial.values!r}, "
            f"final={dox_final.values!r}, "
            f"absolute drift={(dox_final.values - dox_initial.values)!r}"
        ),
    )

    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system DOX-with-Nitrogen "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )
