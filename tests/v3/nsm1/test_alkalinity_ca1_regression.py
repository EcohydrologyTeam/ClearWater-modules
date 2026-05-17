"""NSM1-CA-1 (CRITICAL) regression: alkalinity algal/benthic coupling
must use the *intensive* carbon ratio (``rca = AWc/AWa``,
``rcb = BWc/BWd``), not the raw stoichiometric weight (``AWc``/``BWc``).

Gold-standard spec Section 3, Workstream A1.

Non-shared-path contract (spec Section 1(4)): the expected values here
are computed from **independently hardcoded literal constants**, NOT by
importing ``parameters.alkalinity`` / ``parameters.algae`` DEFAULTS or by
reading ``alk_proc.r_alkaa`` / ``alk_proc.AWc`` / ``alk_proc.AWa``. The
v3 Process and this "v1-mirror" reference therefore do not share the
constant or code path under test -- the structural defect that let
CA-1 escape the existing parity suite.

Acceptance (spec A1): floating-algae alkalinity flux at default
stoichiometry is the intensive-ratio value (~9.172e-4 mg-CaCO3/L per
5-min step under the synthetic bloom below), NOT the raw-weight value
(~9.172e-1, exactly 1000x larger); benthic terms 100x; and a synthetic
bloom does not trigger silent negative-clip events.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.alkalinity import Alkalinity
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


# --- Independently hardcoded reference constants (the "v1-mirror" side).
# These literals are NOT imported from the v3 parameter library; they
# reproduce the v1 / Fortran / CE-QUAL-W2 convention by hand. ---
R_ALKAA = 14.0 / 106.0 / 12.0 / 1000.0   # eq/mg-C, NH4 photosynthesis path
R_ALKAN = 18.0 / 106.0 / 12.0 / 1000.0   # eq/mg-C, NO3 photosynthesis path
R_ALKBA = 14.0 / 106.0 / 12.0 / 1000.0   # eq/mg-C, benthic NH4 path
AWC_RAW = 40.0      # mg-C per stoichiometric unit
AWA_RAW = 1000.0    # ug-Chla per stoichiometric unit
BWC_RAW = 40.0      # mg-C per stoichiometric unit
BWD_RAW = 100.0     # mg-D per stoichiometric unit
RCA_INTENSIVE = AWC_RAW / AWA_RAW   # 0.04 mg-C/ug-Chla  (correct)
RCB_INTENSIVE = BWC_RAW / BWD_RAW   # 0.40 mg-C/mg-D     (correct)
EQ_TO_MG_CACO3 = 50000.0

DT = timedelta(minutes=5)
DT_DAYS = DT.total_seconds() / 86400.0


@dataclass
class _MockFloatingAlgae:
    algal_growth_rate: xr.DataArray
    algal_respiration_rate: xr.DataArray
    algal_nh4_uptake_fraction: xr.DataArray | float


@dataclass
class _MockBenthicAlgae:
    balgae_growth_rate: xr.DataArray
    balgae_respiration_rate: xr.DataArray
    balgae_nh4_uptake_fraction: xr.DataArray | float


def _registry(alk0: float) -> InMemoryRegistry:
    reg = InMemoryRegistry()
    reg.register("alkalinity", xr.DataArray(np.array([alk0]), dims="cell"))
    reg.register(
        "water_temperature", xr.DataArray(np.array([20.0]), dims="cell")
    )
    reg.register("depth", xr.DataArray(np.array([1.0]), dims="cell"))
    return reg


def test_ca1_floating_growth_uses_intensive_ratio_not_raw_weight():
    """Synthetic bloom: ApGrowth = 12 ug-Chla/L/d, all-NH4 uptake.

    Intensive ratio  -> flux ~ 9.172e-4 mg-CaCO3/L per 5-min step.
    Raw weight (bug)  -> flux ~ 9.172e-1 (exactly 1000x larger).
    """
    ap_growth = 12.0
    fnh4 = 1.0

    # Independently computed expectation (v1-mirror, hardcoded literals).
    expected_rate = (
        (R_ALKAA * fnh4 - R_ALKAN * (1.0 - fnh4))
        * ap_growth
        * RCA_INTENSIVE
        * EQ_TO_MG_CACO3
    )
    raw_weight_rate = (
        (R_ALKAA * fnh4 - R_ALKAN * (1.0 - fnh4))
        * ap_growth
        * AWC_RAW            # the pre-fix defect: raw weight, not AWc/AWa
        * EQ_TO_MG_CACO3
    )
    # Sanity on the spec's stated magnitudes (spec A1: "~9.172e-4 ...
    # not 9.172e-1"; the precise check is the rtol=1e-12 assert below).
    np.testing.assert_allclose(
        abs(expected_rate * DT_DAYS), 9.172e-4, rtol=1e-3
    )
    np.testing.assert_allclose(raw_weight_rate, 1000.0 * expected_rate,
                               rtol=1e-12)

    alk = Alkalinity(time_step=DT)
    alk.diagnostics = Diagnostics()
    alk.use_floating_algae = True
    alk.floating_algae_process = _MockFloatingAlgae(
        algal_growth_rate=xr.DataArray(np.array([ap_growth]), dims="cell"),
        algal_respiration_rate=xr.DataArray(np.array([0.0]), dims="cell"),
        algal_nh4_uptake_fraction=xr.DataArray(np.array([fnh4]), dims="cell"),
    )

    reg = _registry(100.0)
    alk.run(datetime(2026, 5, 16), reg)

    flux = np.asarray(alk.alk_algal_growth_rate)
    np.testing.assert_allclose(flux, expected_rate, rtol=1e-12)
    # Hard anti-regression: must NOT be the raw-weight (1000x) value.
    assert np.all(np.abs(flux - raw_weight_rate) > 1.0)

    # Growth is a net sink under all-NH4 uptake: alk drops by ~9.172e-4
    # over one 5-min step, NOT ~9.172e-1.
    alk_final = np.asarray(reg.get_at_time("alkalinity", datetime(2026, 5, 16)))
    delta = alk_final - 100.0
    np.testing.assert_allclose(delta, -expected_rate * DT_DAYS, rtol=1e-12)
    assert np.all(np.abs(delta) < 1.0e-2)  # pre-fix would be ~0.92


def test_ca1_benthic_growth_uses_intensive_ratio_not_raw_weight():
    """Benthic coupling uses rcb = BWc/BWd (100x anti-regression)."""
    ab_growth = 5.0  # g-D/m^2/d
    fnh4 = 1.0
    fb = 0.9
    depth = 1.0

    expected_rate = (
        (1.0 / depth)
        * (R_ALKBA * fnh4 - R_ALKAN * (1.0 - fnh4))
        * ab_growth
        * fb
        * RCB_INTENSIVE
        * EQ_TO_MG_CACO3
    )
    raw_weight_rate = (
        (1.0 / depth)
        * (R_ALKBA * fnh4 - R_ALKAN * (1.0 - fnh4))
        * ab_growth
        * fb
        * BWC_RAW            # pre-fix defect
        * EQ_TO_MG_CACO3
    )
    np.testing.assert_allclose(raw_weight_rate, 100.0 * expected_rate,
                               rtol=1e-12)

    alk = Alkalinity(time_step=DT)
    alk.diagnostics = Diagnostics()
    alk.use_benthic_algae = True
    alk.benthic_algae_process = _MockBenthicAlgae(
        balgae_growth_rate=xr.DataArray(np.array([ab_growth]), dims="cell"),
        balgae_respiration_rate=xr.DataArray(np.array([0.0]), dims="cell"),
        balgae_nh4_uptake_fraction=xr.DataArray(np.array([fnh4]), dims="cell"),
    )

    reg = _registry(100.0)
    alk.run(datetime(2026, 5, 16), reg)

    flux = np.asarray(alk.alk_benthic_algae_growth_rate)
    np.testing.assert_allclose(flux, expected_rate, rtol=1e-12)
    assert np.all(np.abs(flux - raw_weight_rate) > 1.0)


def test_ca1_synthetic_bloom_mass_conservation_no_silent_clip():
    """Sustained synthetic bloom over a 1-day window.

    With the intensive-ratio fix, growth and respiration alkalinity
    fluxes are O(1e-1)/d, so a day of integration perturbs a 100
    mg-CaCO3/L pool by O(1) and never drives it negative -> no
    ``clip_events``. The pre-fix raw-weight code (1000x) would crash
    alkalinity through zero within a few substeps and fire silent
    negative-clips; asserting an empty clip log is the regression guard.
    """
    ap_growth = 10.0
    ap_resp = 8.0
    fnh4 = 0.6

    alk = Alkalinity(time_step=DT)
    diagnostics = Diagnostics()
    alk.diagnostics = diagnostics
    alk.use_floating_algae = True
    alk.floating_algae_process = _MockFloatingAlgae(
        algal_growth_rate=xr.DataArray(np.array([ap_growth]), dims="cell"),
        algal_respiration_rate=xr.DataArray(np.array([ap_resp]), dims="cell"),
        algal_nh4_uptake_fraction=xr.DataArray(np.array([fnh4]), dims="cell"),
    )

    reg = _registry(100.0)
    start = datetime(2026, 5, 16)
    end = start + timedelta(days=1)
    t = start
    while t < end:
        alk.run(t, reg)
        t += DT

    alk_final = np.asarray(reg.get_at_time("alkalinity", t))
    # No silent negative-clip under the corrected magnitude.
    assert diagnostics.clip_events == {}, (
        f"Synthetic bloom fired clip events with the CA-1 fix in place "
        f"(would indicate the raw-weight 1000x defect): "
        f"{diagnostics.clip_events!r}"
    )
    # Physically plausible: O(1) perturbation, not the O(1e3) runaway.
    assert np.all(alk_final > 90.0)
    assert np.all(alk_final < 110.0)
    assert np.all(np.isfinite(alk_final))
