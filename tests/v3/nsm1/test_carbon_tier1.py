"""Tier 1 closed-system conservation tests for v3 NSM1 ``Carbon`` Process.

The Tier 1 contract (design spec Section 9, Section 14 resolved Q7):

    Closed system + no boundaries + no settling + balanced source/sink
    pairs --> total mass of N, P, C, O2-equivalents, Alk constant to
    floating-point roundoff AND ``diagnostics.clip_events == {}``.

This module covers the Phase 5.A Carbon Process. The closed-system test
asserts that total carbon (POC + DOC + DIC) is *invariant* when:

* POC hydrolysis disabled  (``kpoc_20 = 0``) -- no POC -> DOC flux
* DOC oxidation  disabled  (``kdoc_20 = 0``) -- no DOC -> DIC flux
* POC settling   disabled  (``vsoc = 0``)
* CO2 reaeration disabled  (``kah_20_user = 0``, ``kaw_20_user = 0`` AND
  user-defined menu options so ``ka_tc == 0``)
* Sediment release disabled (``JDIC = 0``, ``use_SedFlux = False``)
* No FloatingAlgae, BenthicAlgae, POM, or DOX Processes wired up, so all
  algal coupling and POM-hydrolysis source/sink terms drop out.

Under these conditions the per-cell rates ``dPOC/dt``, ``dDOC/dt``,
``dDIC/dt`` are exactly zero and the Forward Euler integrator should
leave POC, DOC, and DIC unchanged to roundoff. The total carbon (sum of
the three pools across cells) is therefore invariant.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.carbon import Carbon
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


def test_tier1_carbon_conservation_closed_system_loss_disabled(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system total-C conservation when settling, atmospheric
    exchange, and sediment release are disabled. With no algae, POM,
    or DOX in the model, the POC<->DOC<->DIC chain is internally
    balanced and total C is invariant.
    """
    start, end, time_step = closed_system_time_window

    # Snapshot the initial POC / DOC / DIC pools (5-cell mesh from
    # ``initial_state_5cell``).
    poc_initial = in_memory_registry.get("poc").copy()
    doc_initial = in_memory_registry.get("doc").copy()
    dic_initial = in_memory_registry.get("dic").copy()
    total_c_initial = (
        poc_initial.sum() + doc_initial.sum() + dic_initial.sum()
    )

    # Construct Carbon with all internal kinetics zeroed, all algal /
    # POM / DOX coupling defaulted off (no Model wired up), and
    # atmospheric / sediment fluxes disabled.
    carbon_process = Carbon(
        parameters={
            # POC <-> DOC <-> DIC chain rate constants disabled.
            "kpoc_20": 0.0,
            "kdoc_20": 0.0,
            "kpoc_theta": 1.0,    # benign even at kpoc_20 == 0
            "kdoc_theta": 1.0,    # benign even at kdoc_20 == 0
            # POC settling disabled.
            "vsoc": 0.0,
            # CO2 atmospheric reaeration disabled: zero user-defined rates
            # AND select the user-defined menu options so ka_tc == 0
            # regardless of the saturation deficit.
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            # Sediment release disabled.
            "JDIC": 0.0,
            "use_SedFlux": False,
            # Algal / POM coupling explicitly disabled (also enforced
            # below by setting the use_* flags on the instance).
            "use_Algae": False,
            "use_Balgae": False,
            "use_POM": False,
            "use_POC": False,    # disable POM->POC settling source path
            "use_DOC": False,    # disable POM->DOC source path
        },
        time_step=time_step,
    )
    # No Model is wired up; ensure coupling flags are False so ``run``
    # skips the floating/benthic algae / POM / DOX branches.
    carbon_process.use_floating_algae = False
    carbon_process.use_benthic_algae = False
    carbon_process.use_pom = False
    carbon_process.use_dox = False

    # Make the test own the diagnostics reference (replaces the
    # locally-instantiated one in __init__).
    diagnostics = Diagnostics()
    carbon_process.diagnostics = diagnostics

    # Drive 100 substeps.
    current_time = start
    while current_time < end:
        carbon_process.run(current_time, in_memory_registry)
        current_time += time_step

    # Tier 1 invariant 1: per-cell POC / DOC / DIC equality to roundoff.
    poc_final = in_memory_registry.get("poc")
    doc_final = in_memory_registry.get("doc")
    dic_final = in_memory_registry.get("dic")

    np.testing.assert_allclose(
        poc_final.values,
        poc_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system POC invariance failed. "
            f"initial={poc_initial.values!r}, "
            f"final={poc_final.values!r}, "
            f"absolute drift={(poc_final.values - poc_initial.values)!r}"
        ),
    )
    np.testing.assert_allclose(
        doc_final.values,
        doc_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system DOC invariance failed. "
            f"initial={doc_initial.values!r}, "
            f"final={doc_final.values!r}, "
            f"absolute drift={(doc_final.values - doc_initial.values)!r}"
        ),
    )
    np.testing.assert_allclose(
        dic_final.values,
        dic_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system DIC invariance failed. "
            f"initial={dic_initial.values!r}, "
            f"final={dic_final.values!r}, "
            f"absolute drift={(dic_final.values - dic_initial.values)!r}"
        ),
    )

    # Tier 1 invariant 2: total carbon (sum across all three pools and
    # all cells) is constant.
    total_c_final = poc_final.sum() + doc_final.sum() + dic_final.sum()
    np.testing.assert_allclose(
        float(total_c_final),
        float(total_c_initial),
        rtol=1e-12,
        err_msg=(
            "Closed-system total-C conservation failed. "
            f"initial={float(total_c_initial)!r}, "
            f"final={float(total_c_final)!r}, "
            f"absolute drift={float(total_c_final - total_c_initial)!r}"
        ),
    )

    # Tier 1 invariant 3: no clipping under closed-system + physically
    # reasonable initial conditions.
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Tier 1 Carbon "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


# ---------------------------------------------------------------------------
# Carbon instantiation smoke tests
# ---------------------------------------------------------------------------


def test_carbon_instantiates_from_defaults() -> None:
    """``Carbon()`` constructs from the composed CARBON_DEFAULTS without
    errors.

    Acceptance check that the lazy-loaded composed-DEFAULTS pattern is
    wired up correctly and that the carbon-specific kinetics and
    coupling defaults are populated on the instance.
    """
    carbon_process = Carbon()
    # Required parameters from v3 CARBON_DEFAULTS
    assert hasattr(carbon_process, "kpoc_20")
    assert hasattr(carbon_process, "kdoc_20")
    assert hasattr(carbon_process, "kpoc_theta")
    assert hasattr(carbon_process, "kdoc_theta")
    assert hasattr(carbon_process, "KsOxmc")
    assert hasattr(carbon_process, "pCO2")
    assert hasattr(carbon_process, "FCO2")
    assert hasattr(carbon_process, "roc")
    # Reaeration menu (composed from DOX_DEFAULTS)
    assert hasattr(carbon_process, "kah_20_user")
    assert hasattr(carbon_process, "kaw_20_user")
    assert hasattr(carbon_process, "kah_theta")
    assert hasattr(carbon_process, "kaw_theta")
    # POC settling velocity (composed from GLOBAL_VAR_DEFAULTS)
    assert hasattr(carbon_process, "vsoc")
    # Algal stoichiometric ratios (composed from algae / balgae)
    assert hasattr(carbon_process, "AWc")
    assert hasattr(carbon_process, "BWc")
    assert hasattr(carbon_process, "Fb")
    # Sediment release knob (Phase 5.A standalone)
    assert hasattr(carbon_process, "JDIC")
    assert carbon_process.JDIC == 0.0
    # Coupling defaults to disconnected
    assert carbon_process.use_floating_algae is False
    assert carbon_process.use_benthic_algae is False
    assert carbon_process.use_pom is False
    assert carbon_process.use_dox is False
    # Step-scoped caches start zeroed (Q10 GS-rates contract)
    assert carbon_process.doc_dic_oxidation_rate == 0.0
    assert carbon_process.poc_hydrolysis_rate == 0.0
    # Diagnostics handle should be live
    assert carbon_process.diagnostics is not None
    assert carbon_process.diagnostics.clip_events == {}


def test_carbon_run_caches_doc_dic_oxidation_rate(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """``run`` caches ``self.doc_dic_oxidation_rate`` (mg-C/L/d) for the
    Phase 5.B DOX Process consumer.

    Under the Tier 1 closed-system parameters (``kdoc_20 == 0``) this
    cache should be zero post-run; the load-bearing assertion is that
    the attribute exists and is consistent with the kinetics dict.
    """
    start, _end, time_step = closed_system_time_window

    carbon_process = Carbon(
        parameters={
            "kpoc_20": 0.0,
            "kdoc_20": 0.0,
            "kpoc_theta": 1.0,
            "kdoc_theta": 1.0,
            "vsoc": 0.0,
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "JDIC": 0.0,
            "use_SedFlux": False,
            "use_Algae": False,
            "use_Balgae": False,
            "use_POM": False,
        },
        time_step=time_step,
    )
    carbon_process.use_floating_algae = False
    carbon_process.use_benthic_algae = False
    carbon_process.use_pom = False
    carbon_process.use_dox = False

    carbon_process.run(start, in_memory_registry)

    # Cache exists and is consistent: kdoc_20 == 0 implies zero flux.
    assert hasattr(carbon_process, "doc_dic_oxidation_rate")
    rate = carbon_process.doc_dic_oxidation_rate
    if isinstance(rate, xr.DataArray):
        assert np.all(np.asarray(rate.values) == 0.0)
    else:
        assert rate == 0.0
