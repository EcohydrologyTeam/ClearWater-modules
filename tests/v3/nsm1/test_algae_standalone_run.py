"""FloatingAlgae / BenthicAlgae are runnable when bare-constructed.

Regression for a latent bug: the nutrient-availability flags
``use_nitrate`` / ``use_ammonium`` / ``use_phosphate`` (read by
``run()``'s ``limit_nitrogen`` / ``limit_phosphorus`` terms) were set
only in ``init_process(model, ...)``. A process constructed directly and
run without going through the Model — as unit tests and isolated
harnesses do — never had them set, so ``run()`` raised
``AttributeError: 'FloatingAlgae' object has no attribute
'use_phosphate'`` (and likewise for BenthicAlgae, which inherits
FloatingAlgae.__init__).

The fix sets the three flags to ``True`` in ``FloatingAlgae.__init__``
(BenthicAlgae inherits it). ``init_process`` still overrides them when
the process is wired to a Model, so configured / coupled runs are
unchanged (verified by ``test_coupled_demo_parity.py``).
"""
from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.examples import InMemoryRegistry
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae


N = 5
START = datetime(2026, 1, 1, 0, 0, 0)


def _da(value) -> xr.DataArray:
    return xr.DataArray(np.full(N, value, dtype=float), dims=["cell"])


def _registry(state_name, state_value) -> InMemoryRegistry:
    reg = InMemoryRegistry()
    reg.register(state_name, _da(state_value))
    reg.register("ammonium", _da(0.10))
    reg.register("nitrate", _da(0.20))
    reg.register("tip", _da(0.10))
    reg.register("depth", _da(1.0))
    reg.register("water_temperature", _da(20.0))
    reg.register("solar_radiation", _da(300.0))
    return reg


# (class, primary-state name, initial value)
CASES = [
    (FloatingAlgae, "algae_floating", 40.0),
    (BenthicAlgae, "benthic_algae", 5.0),
]


@pytest.mark.parametrize("cls, state_name, state_value", CASES)
def test_bare_construction_sets_nutrient_flags(cls, state_name, state_value):
    """The three flags are present (True) immediately after __init__,
    before any init_process call."""
    proc = cls()
    assert proc.use_nitrate is True
    assert proc.use_ammonium is True
    assert proc.use_phosphate is True


@pytest.mark.parametrize("cls, state_name, state_value", CASES)
def test_bare_construction_runs_one_substep(cls, state_name, state_value):
    """A bare-constructed process runs one substep without AttributeError
    and updates its primary state to a finite value."""
    proc = cls()
    reg = _registry(state_name, state_value)
    proc.run(START, reg)  # previously raised AttributeError on use_phosphate
    out = np.asarray(reg.get(state_name).values)
    assert out.shape == (N,)
    assert np.isfinite(out).all()
