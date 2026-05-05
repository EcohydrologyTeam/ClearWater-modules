"""Tier 1 closed-system conservation tests for v3 NSM1 ``POM`` Process.

The Tier 1 contract (design spec Section 9, Section 14 resolved Q7):

    Closed system + no boundaries + no settling + balanced source/sink
    pairs --> total mass of N, P, C, O2-equivalents, Alk constant to
    floating-point roundoff AND ``diagnostics.clip_events == {}``.

This module covers the Phase 3.2 POM Process. The single closed-system
test asserts that POM is *invariant* when:

* Dissolution rate ``kpom_20 = 0``      (no POM -> DOC sink)
* Burial velocity ``vb = 0``            (no sediment loss)
* POC settling velocity ``vsoc = 0``    (no POC -> POM source)
* No FloatingAlgae or BenthicAlgae Processes instantiated, so the
  algal mortality / settling source terms drop out.

Under these conditions the POM rate is exactly zero and the Forward
Euler integrator should leave POM unchanged to roundoff.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.pom import POM
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


def test_tier1_pom_conservation_closed_system_loss_disabled(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system POM conservation when dissolution and settling are
    disabled. With no algal sources (FloatingAlgae/BenthicAlgae not in
    Model), POM should be invariant.
    """
    start, end, time_step = closed_system_time_window

    # Snapshot initial POM (5-cell mesh from ``initial_state_5cell``).
    pom_initial = in_memory_registry.get("pom").copy()

    # Construct POM with all loss/source pathways disabled. We pass the
    # ``use_*`` flags as kwargs to suppress algal/POC coupling even if
    # the registry happens to have ``poc`` registered (it does in the
    # Tier 1 fixture, but with vsoc=0 it would not contribute mass
    # anyway).
    pom_process = POM(
        parameters={
            "kpom_20": 0.0,         # POM -> DOC dissolution disabled
            "kpom_theta": 1.0,      # benign even at kpom_20 == 0
            "vb": 0.0,              # burial disabled
            "vsoc": 0.0,            # POC settling source disabled
            "use_POC": False,       # gate POC source explicitly
            "use_Algae": False,     # gate floating algae source explicitly
            "use_Balgae": False,    # gate benthic algae source explicitly
        },
        time_step=time_step,
    )
    # No Model is wired up; ensure coupling flags are False so
    # ``run`` skips the floating/benthic algae branches.
    pom_process.use_floating_algae = False
    pom_process.use_benthic_algae = False

    # Drive 100 substeps.
    diagnostics = pom_process.diagnostics
    current_time = start
    while current_time < end:
        pom_process.run(current_time, in_memory_registry)
        current_time += time_step

    # Tier 1 invariant 1: POM cell-wise constant to roundoff.
    pom_final = in_memory_registry.get("pom")
    np.testing.assert_allclose(
        pom_final.values,
        pom_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system POM conservation failed. "
            f"initial={pom_initial.values!r}, "
            f"final={pom_final.values!r}, "
            f"absolute drift={(pom_final.values - pom_initial.values)!r}"
        ),
    )

    # Tier 1 invariant 2: no clipping under closed-system + physically
    # reasonable initial conditions.
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Tier 1 conditions: "
        f"{diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


def test_pom_instantiates_from_defaults() -> None:
    """``POM()`` constructs from POM_DEFAULTS without errors.

    Acceptance check that the lazy-loaded DEFAULTS pattern is wired up
    correctly and that the inline coupling defaults are populated.
    """
    pom_process = POM()
    # Required parameters from v3 POM_DEFAULTS
    assert hasattr(pom_process, "kpom_20")
    assert hasattr(pom_process, "kpom_theta")
    assert hasattr(pom_process, "h2")
    # Required coupling defaults from POM_GLOBAL_DEFAULTS
    assert hasattr(pom_process, "vsoc")
    assert hasattr(pom_process, "fcom")
    assert hasattr(pom_process, "vb")
    assert hasattr(pom_process, "use_POC")
    assert hasattr(pom_process, "use_Algae")
    assert hasattr(pom_process, "use_Balgae")
    # Diagnostics handle should be live
    assert pom_process.diagnostics is not None
    assert pom_process.diagnostics.clip_events == {}
