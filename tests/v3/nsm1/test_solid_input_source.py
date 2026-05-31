"""Shared `Solid` (suspended-solids) input source.

Guards `design/clearwater_modules_v3_solid_input_source.md`: `Solid`
(mg/L) is a single canonical, OPTIONAL registry input read
registry-first / scalar-fallback (the `wind_shelter_coefficient`
pattern). `phosphorus`, `floating_algae`, and `benthic_algae` previously
used a constructor-scalar `self.Solid`; their `run()` now prefers a
registered per-cell `Solid` field and falls back to `self.Solid` when
none is provided, threading the resolved value into
`_change_with_components` -> `utils.partitioning.fdp`. (`pathogen`
already read `Solid` from the registry and is unchanged.)

Two layers of coverage:

1. ``test_run_resolves_solid_*`` (all three processes) monkeypatch
   ``_change_with_components`` to capture the ``solid`` value ``run()``
   resolves and forwards — proving registry-first + scalar-fallback at
   the ``run()`` boundary. This is independent of the downstream
   kinetics (FloatingAlgae/BenthicAlgae are not runnable bare end-to-end
   without the configured model wiring ``use_phosphate``; that is a
   separate pre-existing harness limitation).
2. ``test_phosphorus_solid_affects_output`` runs Phosphorus end-to-end
   (it runs standalone cleanly) with ``kdpo4 > 0`` so the resolved
   ``Solid`` actually moves ``fdp = 1/(1 + kdpo4*Solid*1e-6)`` and the
   TIP-settling output — proving the resolved value is really used, not
   just forwarded.
"""
from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.examples import InMemoryRegistry
from clearwater_modules_v3.processes.phosphorus import Phosphorus
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae


N = 5
START = datetime(2026, 1, 1, 0, 0, 0)


def _da(value) -> xr.DataArray:
    return xr.DataArray(np.full(N, value, dtype=float), dims=["cell"])


def _phosphorus_inputs() -> dict:
    return {
        "tip": _da(0.10),
        "organic_phosphorus": _da(0.05),
        "water_temperature": _da(20.0),
        "depth": _da(1.0),
    }


def _floating_inputs() -> dict:
    return {
        "algae_floating": _da(40.0),
        "ammonium": _da(0.10),
        "nitrate": _da(0.20),
        "tip": _da(0.10),
        "depth": _da(1.0),
        "water_temperature": _da(20.0),
        "solar_radiation": _da(300.0),
    }


def _benthic_inputs() -> dict:
    return {
        "benthic_algae": _da(5.0),
        "ammonium": _da(0.10),
        "nitrate": _da(0.20),
        "tip": _da(0.10),
        "depth": _da(1.0),
        "water_temperature": _da(20.0),
        "solar_radiation": _da(300.0),
    }


# (label, class, inputs-builder)
CASES = [
    ("phosphorus", Phosphorus, _phosphorus_inputs),
    ("floating_algae", FloatingAlgae, _floating_inputs),
    ("benthic_algae", BenthicAlgae, _benthic_inputs),
]


class _StopAfterResolve(Exception):
    """Raised by the spy helper to stop run() right after Solid resolution."""


def _solid_forwarded_by_run(
    monkeypatch, cls, base_inputs, *, solid_param, registry_solid
):
    """Capture the ``solid`` value ``run()`` resolves and forwards to
    ``_change_with_components``. Stops execution before the (kinetics)
    helper body runs, so this is immune to downstream wiring."""
    seen = {}

    def _spy(self, **kwargs):
        seen["solid"] = kwargs.get("solid")
        raise _StopAfterResolve

    monkeypatch.setattr(cls, "_change_with_components", _spy)

    proc = cls(parameters={"Solid": solid_param})
    reg = InMemoryRegistry()
    for key, value in base_inputs.items():
        reg.register(key, value)
    if registry_solid is not None:
        reg.register("Solid", _da(registry_solid))

    with pytest.raises(_StopAfterResolve):
        proc.run(START, reg)
    return seen["solid"]


@pytest.mark.parametrize("name, cls, inputs_fn", CASES)
def test_run_resolves_registered_solid(monkeypatch, name, cls, inputs_fn):
    """When a per-cell ``Solid`` is registered, ``run()`` forwards THAT
    field (not the constructor scalar) to the kinetics helper."""
    seen = _solid_forwarded_by_run(
        monkeypatch, cls, inputs_fn(), solid_param=99.0, registry_solid=7.0
    )
    seen_arr = np.asarray(seen)
    assert seen_arr.shape == (N,), f"{name}: expected the per-cell registry field"
    assert np.all(seen_arr == 7.0), (
        f"{name}: run() forwarded {seen!r}, not the registered Solid=7.0 "
        "(registry-first violated; the constructor scalar 99.0 leaked through)"
    )


@pytest.mark.parametrize("name, cls, inputs_fn", CASES)
def test_run_falls_back_to_param_when_absent(monkeypatch, name, cls, inputs_fn):
    """With no ``Solid`` registered, ``run()`` forwards the constructor
    scalar ``self.Solid``."""
    seen = _solid_forwarded_by_run(
        monkeypatch, cls, inputs_fn(), solid_param=3.0, registry_solid=None
    )
    assert np.all(np.asarray(seen) == 3.0), (
        f"{name}: run() forwarded {seen!r}, not the self.Solid fallback 3.0"
    )


def test_phosphorus_solid_affects_output():
    """End-to-end: with kdpo4>0 the resolved Solid moves fdp and the TIP
    output. Two registered Solid values must give different tip; and a
    registered field must override the constructor scalar."""

    def _run(*, solid_param, registry_solid):
        proc = Phosphorus(parameters={"kdpo4": 1.0e5, "use_TIP": True, "Solid": solid_param})
        reg = InMemoryRegistry()
        for key, value in _phosphorus_inputs().items():
            reg.register(key, value)
        if registry_solid is not None:
            reg.register("Solid", _da(registry_solid))
        proc.run(START, reg)
        return np.asarray(reg.get("tip").values)

    # Registered Solid is actually used: 2 vs 50 mg/L -> different fdp -> different tip.
    out_lo = _run(solid_param=2.0, registry_solid=2.0)
    out_hi = _run(solid_param=2.0, registry_solid=50.0)
    assert not np.array_equal(out_lo, out_hi), (
        "registered Solid (2 vs 50 mg/L, kdpo4>0) did not move the TIP output"
    )

    # Registry-first AND fallback both resolve to a 2.0 field -> identical.
    out_registry = _run(solid_param=99.0, registry_solid=2.0)  # registry wins over param 99
    out_fallback = _run(solid_param=2.0, registry_solid=None)  # absent -> param 2.0
    np.testing.assert_array_equal(
        out_registry,
        out_fallback,
        err_msg="registry(field=2.0,param=99) must equal absent(param=2.0)",
    )


def test_demo_registers_no_solid():
    """The coupled demo registers no ``Solid``, so every consumer uses its
    per-process fallback — which is why this change is byte-identical to the
    `d530a3a` baseline (no re-baseline). The bit-identical trajectory itself
    is asserted by `test_coupled_demo_parity.py`; this pins the premise.
    """
    from clearwater_modules_v3.examples import build_nsm1_demo

    demo = build_nsm1_demo()
    assert "Solid" not in demo.registry, (
        "demo registry unexpectedly carries 'Solid'; the no-re-baseline "
        "rationale assumes no Solid provider in the demo"
    )
