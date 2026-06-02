"""Every NSM1 reactive process declares its written state as output_variables.

Guards the residual-2 (dry-cell kinetics inflation) fix. ``Model.__apply_wet_mask``
masks a process's OUTPUT variables to NaN on dry cells when a wet-mask is
configured; for a process that does NOT declare ``output_variables`` it falls
back to masking every name in ``process.variables`` -- including scalar input
forcings, which (a) corrupts forcings for the next substep and (b) raises a
``set_at_time`` dimension error on 0-d parameters. The coupled HEC-RAS-2D run
needs the wet-mask to exclude thin/dry cells from the kinetics output, so each
reactive process must declare exactly the per-cell state(s) it writes via
``registry.set_at_time`` in ``run()``.

This test pins those declarations to the canonical state names each process
persists, so a future edit that drops or renames a declaration is caught.
"""
from __future__ import annotations

import pytest

from clearwater_modules_v3.processes.nitrogen import Nitrogen
from clearwater_modules_v3.processes.carbon import Carbon
from clearwater_modules_v3.processes.dox import DOX
from clearwater_modules_v3.processes.phosphorus import Phosphorus
from clearwater_modules_v3.processes.cbod import CBOD
from clearwater_modules_v3.processes.pom import POM
from clearwater_modules_v3.processes.alkalinity import Alkalinity
from clearwater_modules_v3.processes.pathogen import Pathogen
from clearwater_modules_v3.processes.n2 import N2
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae


# (process class, expected output_variables) -- the per-cell states each
# process persists via registry.set_at_time in run(). Order-independent.
EXPECTED = [
    (Nitrogen, {"ammonium", "nitrate", "organic_nitrogen"}),
    (Carbon, {"poc", "doc", "dic"}),
    (DOX, {"oxygen_dissolved"}),
    (Phosphorus, {"tip", "organic_phosphorus"}),
    (CBOD, {"cbod"}),
    (POM, {"pom"}),
    (Alkalinity, {"alkalinity"}),
    (Pathogen, {"pathogen"}),
    # Already declared before the residual-2 work; pinned here for completeness.
    (N2, {"n2"}),
    (FloatingAlgae, {"algae_floating"}),
    (BenthicAlgae, {"benthic_algae"}),
]


@pytest.mark.parametrize("cls,expected", EXPECTED, ids=[c.__name__ for c, _ in EXPECTED])
def test_process_declares_output_variables(cls, expected):
    declared = getattr(cls, "output_variables", None)
    assert declared is not None, (
        f"{cls.__name__} must declare class-level output_variables so "
        f"Model.__apply_wet_mask masks only outputs, not input forcings."
    )
    assert set(declared) == expected, (
        f"{cls.__name__}.output_variables={list(declared)} != expected {sorted(expected)}"
    )


@pytest.mark.parametrize("cls,expected", EXPECTED, ids=[c.__name__ for c, _ in EXPECTED])
def test_output_variables_exclude_common_forcings(cls, expected):
    """Sanity: a process must never list input forcings as outputs (masking
    those on dry cells corrupts the next substep -- the C5 failure mode)."""
    forcings = {
        "depth", "water_temperature", "wind_speed", "solar_radiation",
        "air_temperature", "pressure_mb", "velocity", "flow", "wetted_surface_area",
    }
    declared = set(getattr(cls, "output_variables", []) or [])
    leaked = declared & forcings
    assert not leaked, f"{cls.__name__}.output_variables leaks input forcings: {leaked}"
