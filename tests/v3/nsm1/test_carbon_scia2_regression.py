"""NSM1-SCI-A2 (MAJOR) regression: algal/benthic mortality carbon is
routed predominantly to **POC** (f_pocp = f_pocb = 0.8), not split
~half to DOC.

Gold-standard spec Workstream C1; E1 author decision 2026-05-16
(f_pocp = f_pocb = 0.8, CE-QUAL-W2 ``APOM``).

The operative routing fraction (``floating_algae.py`` ``_FDP_DEFAULTS`` /
``benthic_algae.py`` ``_BENTHIC_FDP_DEFAULTS``; Carbon consumes the
cached mortality rates) was 0.5 pre-fix, mis-routing ~40% of mortality
C from POC to DOC and biasing DOC->DIC / DO demand. v1 used 0.9;
CE-QUAL-W2 ``APOM`` ~0.8.

Non-shared-path contract (spec Section 1(4)): the expected split is
built from an **independently hardcoded** literal (0.8), NOT by reading
the process ``f_pocp`` / ``f_pocb`` or importing the DEFAULTS.
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae

from .conftest import InMemoryRegistry


F_POC_LITERAL = 0.8          # independently hardcoded; the constant under test
F_POC_PREFIX_DEFECT = 0.5    # the pre-fix operative value
DT = timedelta(minutes=5)


def _registry(state_name: str) -> InMemoryRegistry:
    reg = InMemoryRegistry()
    one = lambda v: xr.DataArray(np.array([v]), dims="cell")
    reg.register(state_name, one(30.0))
    reg.register("ammonium", one(0.15))
    reg.register("nitrate", one(3.0))
    reg.register("phosphorus_total_inorganic", one(0.10))
    reg.register("depth", one(1.5))
    reg.register("water_temperature", one(20.0))
    reg.register("solar_radiation", one(300.0))
    return reg


def test_scia2_floating_mortality_split_is_80_20_not_50_50():
    fa = FloatingAlgae(time_step=DT, death_rate=0.15)
    fa.use_nitrate = True
    fa.use_ammonium = True
    fa.use_phosphate = True
    reg = _registry("algae_floating")
    t = datetime(2026, 5, 16)
    fa.run(t, reg)

    poc = np.asarray(fa.algal_poc_from_mortality_rate)
    doc = np.asarray(fa.algal_doc_from_mortality_rate)
    total = poc + doc
    assert np.all(total > 0.0), "need nonzero algal mortality to test the split"

    poc_frac = poc / total
    # Independent expectation (hardcoded literal, not fa.f_pocp).
    np.testing.assert_allclose(poc_frac, F_POC_LITERAL, rtol=1e-12)
    np.testing.assert_allclose(doc / total, 1.0 - F_POC_LITERAL, rtol=1e-12)
    # Hard anti-regression: must NOT be the pre-fix 50/50 mis-routing.
    assert np.all(np.abs(poc_frac - F_POC_PREFIX_DEFECT) > 0.25)
    # DOC sensitivity: the corrected split routes strictly less mortality
    # C to DOC than the pre-fix 0.5 would have (0.2*total vs 0.5*total).
    assert np.all(doc < F_POC_PREFIX_DEFECT * total)


def test_scia2_benthic_mortality_split_is_80_20():
    ba = BenthicAlgae(time_step=DT)
    ba.use_nitrate = True
    ba.use_ammonium = True
    ba.use_phosphate = True
    reg = _registry("benthic_algae")
    t = datetime(2026, 5, 16)
    ba.run(t, reg)

    poc = np.asarray(ba.balgae_poc_from_mortality_rate)
    doc = np.asarray(ba.balgae_doc_from_mortality_rate)
    total = poc + doc
    assert np.all(total > 0.0), "need nonzero benthic mortality"

    np.testing.assert_allclose(poc / total, F_POC_LITERAL, rtol=1e-12)
    np.testing.assert_allclose(doc / total, 1.0 - F_POC_LITERAL, rtol=1e-12)
    assert np.all(doc < F_POC_PREFIX_DEFECT * total)
