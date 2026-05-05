"""Tier 1 closed-system mass-conservation test for v3 NSM1 N2 / TDG.

Phase 3.4 (v3 NSM1 design spec, Section 11 Phase 3, Section 5 N2/TDG
design notes): asserts that the v3-native ``N2`` Process is conservative
(N2 invariant) when atmospheric exchange and the denitrification source
are both disabled. Under these conditions:

* The N2 state must equal its initial value at every cell to roundoff
  (``rtol=1e-12``).
* The clip-with-log diagnostics must remain empty
  (``diagnostics.clip_events == {}``); a clip event under closed-system
  conditions signals either a malformed test or an integrator bug.

How the closed system is constructed:

* Atmospheric exchange is disabled by setting both ``kah_20_user=0`` and
  ``kaw_20_user=0`` AND selecting the user-defined menu options
  (``hydraulic_reaeration_option=1``, ``wind_reaeration_option=1``) so
  the menu picks zero. This makes ``ka_tc == 0`` and zeroes the
  ``1.034 * ka_tc * (N2sat - N2)`` exchange term regardless of the
  saturation deficit.
* No Nitrogen process is wired into the model, so the
  ``denitrification_rate`` source is identically zero
  (``self.use_nitrogen == False`` after ``run`` without ``init_process``).
* No floating algae, no benthic algae, no DOX -- this Process is pure
  state ``n2`` only.

The TDG derived variable is computed but not asserted on for this Tier 1
test (TDG is a diagnostic ratio, not a conserved mass; closed-system
TDG behavior is exercised by Phase 5 once DOX lands).
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.n2 import N2
from clearwater_modules_v3.processes.nitrogen import Nitrogen
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


def test_tier1_n2_conservation_closed_system_loss_disabled(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system N2 conservation when atmospheric exchange and
    denitrification source are both disabled. N2 should be invariant.

    Setup:
    * 5-cell mesh; ``initial_state_5cell`` initial conditions
      (``n2`` initialized to ``[10.0, 10.5, 11.0, 11.5, 12.0]`` mg-N/L).
    * ``kah_20_user=0``, ``kaw_20_user=0`` AND
      ``hydraulic_reaeration_option=1``, ``wind_reaeration_option=1``
      together force the reaeration menu to pick zero ->
      ``ka_tc == 0`` -> atmospheric exchange flux == 0.
    * No Nitrogen process is registered with the model, so
      ``self.use_nitrogen`` is False after ``__init__`` and the
      denitrification source is identically zero.

    Expected:
    * ``n2_final == n2_initial`` per cell (rtol=1e-12).
    * ``diagnostics.clip_events == {}``.
    """
    start, end, time_step = closed_system_time_window

    n2 = N2(
        parameters={
            # Disable atmospheric exchange: zero user-defined rates AND
            # select the user-defined menu options so ka_tc == 0.
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            # Standard sea-level pressure (the corrected v3 default).
            "pressure_mb": 1013.25,
        },
        time_step=time_step,
    )
    # Make the test own the diagnostics reference (replaces the
    # locally-instantiated one in __init__).
    diagnostics = Diagnostics()
    n2.diagnostics = diagnostics

    # No Nitrogen process is wired up in this Tier 1 harness; verify
    # that ``self.use_nitrogen`` defaults to False so the denit source
    # is identically zero in ``run``.
    assert n2.use_nitrogen is False
    assert n2.nitrogen_process is None

    # Snapshot the initial N2 state for the per-cell equality check.
    n2_initial = in_memory_registry.get_at_time("n2", start).copy()

    # Run 100 substeps. ``N2.run`` mutates the registry in place.
    current_time = start
    while current_time < end:
        n2.run(current_time, in_memory_registry)
        current_time += time_step

    # Tier 1 invariant 1: N2 state per-cell invariance under
    # zeroed-flux closed-system conditions.
    n2_final = in_memory_registry.get_at_time("n2", end)
    np.testing.assert_allclose(
        n2_final.values,
        n2_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system N2 invariance failed under "
            "kah_20_user=kaw_20_user=0 and no Nitrogen denit source. "
            f"initial={n2_initial.values!r}, "
            f"final={n2_final.values!r}, "
            f"absolute drift={(n2_final.values - n2_initial.values)!r}"
        ),
    )

    # Tier 1 invariant 2: no clipping under closed-system conditions.
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Tier 1 N2 "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


# ---------------------------------------------------------------------------
# N2 instantiation smoke tests
# ---------------------------------------------------------------------------


def test_n2_instantiates_with_defaults() -> None:
    """``N2()`` constructs cleanly with no arguments and pulls a composed
    DEFAULTS dict (from global_parameters / dox / global_vars) onto the
    instance, since v3 ``parameters.n2`` is empty per the Phase 1.2 audit.
    """
    n2 = N2()
    # Pulled from global_parameters
    assert n2.pressure_mb == 1013.25
    # Pulled from dox (reaeration menu defaults)
    assert n2.kah_20_user == 0.0
    assert n2.kaw_20_user == 0.0
    assert n2.kah_theta == 1.024
    assert n2.kaw_theta == 1.024
    assert n2.hydraulic_reaeration_option == 1
    assert n2.wind_reaeration_option == 1
    # Pulled from global_vars (toy hydraulic forcings)
    assert n2.velocity == 1.0
    assert n2.flow == 2.0
    # Coupling defaults to disconnected
    assert n2.use_nitrogen is False
    assert n2.nitrogen_process is None
    # Step-scoped caches start zeroed.
    assert n2.n2_sat == 0.0
    assert n2.n2_atm_exchange_rate == 0.0
    assert n2.tdg == 0.0


def test_n2_accepts_parameter_override() -> None:
    """``N2(parameters={...})`` overrides composed defaults and warns on
    unknown keys (no exception raised)."""
    n2 = N2(
        parameters={
            "kah_20_user": 0.5,
            "pressure_mb": 1000.0,
        },
    )
    assert n2.kah_20_user == 0.5
    assert n2.pressure_mb == 1000.0
    # Unchanged composed defaults
    assert n2.kaw_20_user == 0.0
    assert n2.kah_theta == 1.024


def test_n2_run_computes_tdg_diagnostic(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """``run`` computes the TDG diagnostic and exposes it on ``self.tdg``.

    Under closed-system zero-exchange settings the N2 state is invariant
    and TDG = N2 / N2sat is a finite, positive ratio.
    """
    start, _end, time_step = closed_system_time_window

    n2 = N2(
        parameters={
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        time_step=time_step,
    )
    n2.run(start, in_memory_registry)
    # TDG is exposed as an attribute (per spec acceptance criterion).
    assert hasattr(n2, "tdg")
    assert isinstance(n2.tdg, xr.DataArray)
    assert np.all(np.isfinite(n2.tdg.values))
    # N2sat at 20 deg C and 1013.25 mb is ~14-16 mg-N/L; with N2 around
    # 10-12 mg-N/L the TDG ratio should be < 1 and > 0.
    assert np.all(n2.tdg.values > 0.0)
    assert np.all(n2.tdg.values < 1.5)
    # n2_sat cached too
    assert hasattr(n2, "n2_sat")
    assert np.all(np.isfinite(n2.n2_sat.values))


def test_tier1_n2_with_nitrogen_zeroed_kinetics(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system N2 invariance with Nitrogen wired but kinetics zeroed.

    Integration Item 1: verifies that when Nitrogen is wired into N2
    (so N2 reads ``nitrogen_process.denitrification_flux_rate``) but
    every Nitrogen kinetic-rate constant is zero, the cached
    denitrification flux is identically zero and N2 remains invariant.

    This is the contract the production rate-DAG relies on: the
    denitrification source is the *step-scoped* flux (mg-N/L/d),
    not the kinetic-rate constant (1/d). With zeroed rate constants,
    the flux is zero too, and N2 conservation must hold to floating-
    point roundoff.

    Setup:
    * 5-cell mesh, ``initial_state_5cell`` initial conditions.
    * Atmospheric exchange disabled (zeroed user reaeration menu).
    * Nitrogen instance with all kinetic-rate constants set to 0.0.
    * Manually wire ``n2.use_nitrogen = True`` and
      ``n2.nitrogen_process = nitrogen`` (no Model harness).
    * Run order per substep: Nitrogen.run -> N2.run.

    Expected:
    * ``n2_final == n2_initial`` per cell (rtol=1e-12).
    * ``nitrogen.denitrification_flux_rate`` is zero everywhere after
      ``run`` (sanity check on the new cache attribute).
    * N2 ``diagnostics.clip_events == {}``.
    """
    start, end, time_step = closed_system_time_window

    # Build a Nitrogen with all kinetic-rate constants zeroed so the
    # cached fluxes are identically zero.
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
    # ``init_process`` would set these; for the harness we pin them.
    nitrogen.use_nitrate = True
    nitrogen.use_ammonium = True
    nitrogen.use_floating_algae = False
    nitrogen.use_benthic_algae = False

    n2 = N2(
        parameters={
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "pressure_mb": 1013.25,
        },
        time_step=time_step,
    )
    diagnostics = Diagnostics()
    n2.diagnostics = diagnostics

    # Manually wire Nitrogen onto N2 (what ``init_process`` would do
    # under a v3 Model). Confirm both halves of the contract.
    n2.use_nitrogen = True
    n2.nitrogen_process = nitrogen
    assert n2.use_nitrogen is True
    assert n2.nitrogen_process is nitrogen

    # Snapshot the initial N2 state.
    n2_initial = in_memory_registry.get_at_time("n2", start).copy()

    current_time = start
    while current_time < end:
        nitrogen.run(current_time, in_memory_registry)
        # Sanity-check the new flux cache: with zeroed kinetic-rate
        # constants, the denitrification flux must be zero everywhere.
        denit_flux = nitrogen.denitrification_flux_rate
        if hasattr(denit_flux, "values"):
            np.testing.assert_allclose(
                denit_flux.values,
                np.zeros_like(denit_flux.values),
                atol=0.0,
            )
        else:
            assert denit_flux == 0.0
        n2.run(current_time, in_memory_registry)
        current_time += time_step

    # Tier 1 invariant: N2 invariant under zero-flux closed system.
    n2_final = in_memory_registry.get_at_time("n2", end)
    np.testing.assert_allclose(
        n2_final.values,
        n2_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system N2 invariance failed with Nitrogen wired and "
            "zeroed kinetics. "
            f"initial={n2_initial.values!r}, "
            f"final={n2_final.values!r}, "
            f"absolute drift={(n2_final.values - n2_initial.values)!r}"
        ),
    )

    # No clipping under closed-system conditions.
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system N2-with-Nitrogen "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )
