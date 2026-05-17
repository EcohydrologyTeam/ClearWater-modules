"""Joint closed-system alkalinity benchmark — NSM1-CA-1 + NSM1-SCI-N1
together (gold-standard spec D3).

CA-1 and SCI-N1 each have a dedicated non-shared-path regression
(``test_alkalinity_ca1_regression.py``,
``test_alkalinity_scin1_regression.py``). Per spec D3 this adds a
**joint** closed-system benchmark exercising the algal-coupling
(CA-1: intensive ``rca = AWc/AWa``) and the denitrification term
(SCI-N1: ``r_alkden = 1/14/1000``) *simultaneously* in one box, so a
regression in either — or a sign/wiring error in how they combine in
``dAlk/dt`` — is caught.

Non-shared-path contract (spec Section 1(4)): every expected quantity
is built from **independently hardcoded literals** (CA-1 ``rca`` and
SCI-N1 ``r_alkden`` supplied by hand), never read from the Alkalinity
process or the parameter DEFAULTS.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.alkalinity import Alkalinity
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


# --- Independently hardcoded reference constants (no shared symbol). ---
R_ALKAA = 14.0 / 106.0 / 12.0 / 1000.0   # eq/mg-C, NH4 photosynthesis
R_ALKAN = 18.0 / 106.0 / 12.0 / 1000.0   # eq/mg-C, NO3 photosynthesis
RCA_INTENSIVE = 40.0 / 1000.0            # CA-1: AWc/AWa = 0.04 mg-C/ug-Chla
R_ALKDEN = 1.0 / 14.0 / 1000.0           # SCI-N1: 1 eq alkalinity / mol NO3-N
EQ_TO_MG_CACO3 = 50000.0
DT = timedelta(minutes=5)
DT_DAYS = DT.total_seconds() / 86400.0


@dataclass
class _MockNitrogen:
    nitrification_flux_rate: xr.DataArray
    denitrification_flux_rate: xr.DataArray


@dataclass
class _MockFloatingAlgae:
    algal_growth_rate: xr.DataArray
    algal_respiration_rate: xr.DataArray
    algal_nh4_uptake_fraction: xr.DataArray | float


def test_ca1_scin1_joint_closed_system_alkalinity_balance():
    """Closed box: denitrification (SCI-N1) + floating-algae growth and
    respiration (CA-1) active together. The net alkalinity change over
    one substep equals the independently-computed sum of the three
    terms with the correct CA-1 intensive ratio and SCI-N1 1-eq/mol-N
    coefficient — and is NOT the pre-fix raw-weight / 4-eq value.
    """
    denit_flux = 0.5      # mg-N/L/d (water-column NO3 denitrification)
    ap_growth = 8.0       # ug-Chla/L/d
    ap_resp = 3.0         # ug-Chla/L/d
    fnh4 = 1.0            # all-NH4 uptake (growth is a net Alk sink)

    # --- Independent (v1-mirror) expectation, hardcoded literals ---
    # dAlk/dt = denit_source - algal_growth_sink + algal_resp_source
    denit_source = R_ALKDEN * denit_flux * EQ_TO_MG_CACO3
    algal_growth_sink = (
        (R_ALKAA * fnh4 - R_ALKAN * (1.0 - fnh4))
        * ap_growth * RCA_INTENSIVE * EQ_TO_MG_CACO3
    )
    algal_resp_source = ap_resp * R_ALKAA * RCA_INTENSIVE * EQ_TO_MG_CACO3
    expected_rate = denit_source - algal_growth_sink + algal_resp_source
    expected_delta = expected_rate * DT_DAYS

    # Pre-fix counterfactuals (must NOT match): CA-1 raw weight (rca=40,
    # 1000x) and SCI-N1 4-eq (r_alkden=4/14/1000).
    raw_growth_sink = algal_growth_sink * 1000.0
    raw_resp_source = algal_resp_source * 1000.0
    fourx_denit = denit_source * 4.0
    prefix_rate = fourx_denit - raw_growth_sink + raw_resp_source

    alk = Alkalinity(time_step=DT)
    alk.diagnostics = Diagnostics()
    alk.use_nitrogen = True
    alk.use_floating_algae = True
    alk.nitrogen_process = _MockNitrogen(
        nitrification_flux_rate=xr.DataArray(np.array([0.0]), dims="cell"),
        denitrification_flux_rate=xr.DataArray(
            np.array([denit_flux]), dims="cell"
        ),
    )
    alk.floating_algae_process = _MockFloatingAlgae(
        algal_growth_rate=xr.DataArray(np.array([ap_growth]), dims="cell"),
        algal_respiration_rate=xr.DataArray(np.array([ap_resp]), dims="cell"),
        algal_nh4_uptake_fraction=xr.DataArray(np.array([fnh4]), dims="cell"),
    )

    reg = InMemoryRegistry()
    one = lambda v: xr.DataArray(np.array([v]), dims="cell")
    reg.register("alkalinity", one(120.0))
    reg.register("water_temperature", one(20.0))
    reg.register("depth", one(1.0))

    t = datetime(2026, 5, 16)
    alk.run(t, reg)

    # Net rate (cached) matches the independently composed expectation.
    np.testing.assert_allclose(
        np.asarray(alk.alk_rate), expected_rate, rtol=1e-12
    )
    # State delta over the substep matches.
    alk_final = np.asarray(reg.get_at_time("alkalinity", t))
    np.testing.assert_allclose(
        alk_final - 120.0, expected_delta, rtol=1e-12
    )
    # Hard joint anti-regression: not the pre-fix CA-1/SCI-N1 value.
    assert np.all(np.abs(np.asarray(alk.alk_rate) - prefix_rate) > 1.0)

    # Per-term cross-checks (each independently hardcoded).
    np.testing.assert_allclose(
        np.asarray(alk.alk_denitrification_rate), denit_source, rtol=1e-12
    )
    np.testing.assert_allclose(
        np.asarray(alk.alk_algal_growth_rate), algal_growth_sink, rtol=1e-12
    )
    np.testing.assert_allclose(
        np.asarray(alk.alk_algal_respiration_rate),
        algal_resp_source, rtol=1e-12,
    )
    assert alk.diagnostics.clip_events == {}, (
        "closed-system joint CA-1+SCI-N1 box fired a negative-clip"
    )
