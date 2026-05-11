"""v3 POM kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_pom_calculations_v2.py``.

Phase 9.G note: ``pom_hydrolysis_rate`` was renamed to
``pom_doc_source_rate`` and now carries units of mg-C/L_water/d
(``fcom * kpom_tc * pom * h2 / depth``). The original raw v1
dissolution rate (``kpom_tc * pom``, mg-D/L_sed/d) still appears as
the in-process sink term in dPOM/dt; the dissolution-parity test
below asserts both forms (the lumped negative rate matches the raw
v1 form, and the doc-source cache matches the post-9.G corrected
form).
"""
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.pom import POM


V1_DISSOLUTION_REFERENCE = np.array([
    0.2384447948053778,
    0.31928218067904934,
    0.4,
    0.49329404999999993,
    0.6290764288750033,
])

V1_BURIAL_REFERENCE = np.array([
    1.5000000000000002,
    1.75,
    2.0,
    2.25,
    2.5,
])

V1_POC_SETTLING_REFERENCE = np.array([
    0.49999999999999994,
    0.6,
    0.6999999999999998,
    0.7999999999999999,
    0.9,
])

V1_ALGAL_SETTLING_REFERENCE = np.array([
    1.5000000000000002,
    3.0000000000000004,
    4.5,
    6.000000000000001,
    7.5,
])

V1_BENTHIC_MORTALITY_REFERENCE = np.array([
    1.0730015766241996,
    2.46303396523838,
    4.049999999999999,
    5.919528599999997,
    8.492531789812544,
])


@pytest.fixture(scope="module")
def pom_5cell():
    return xr.DataArray(np.array([3.0, 3.5, 4.0, 4.5, 5.0]))


@pytest.fixture(scope="module")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="module")
def poc_5cell():
    return xr.DataArray(np.array([1.0, 1.2, 1.4, 1.6, 1.8]))


@pytest.fixture(scope="module")
def algae_5cell():
    return xr.DataArray(np.array([10.0, 20.0, 30.0, 40.0, 50.0]))


@pytest.fixture(scope="module")
def benthic_algae_5cell():
    return xr.DataArray(np.array([5.0, 10.0, 15.0, 20.0, 25.0]))


@pytest.fixture(scope="module")
def dummy_time() -> datetime:
    return datetime(2026, 1, 1)


@pytest.fixture(scope="module")
def dummy_registry():
    """Stub registry supplying only ``depth`` (the only registry read
    inside POM.rate after Phase 9.G)."""
    class _Stub:
        def __contains__(self, name: str) -> bool:
            return name == "depth"

        def get_at_time(self, name, time):
            if name == "depth":
                return xr.DataArray(np.full(5, 1.0, dtype=float))
            raise KeyError(name)
    return _Stub()


def _make_pom(
    *, kpom_20=0.0, kpom_theta=1.0, vb=0.0, vsoc=0.0, fcom=0.4, h2=0.1,
    use_POC=False, use_Algae=False, use_Balgae=False,
) -> POM:
    pom = POM(
        parameters={
            "kpom_20": kpom_20, "kpom_theta": kpom_theta,
            "vb": vb, "vsoc": vsoc, "fcom": fcom, "h2": h2,
            "use_POC": use_POC, "use_Algae": use_Algae, "use_Balgae": use_Balgae,
        },
        time_step=timedelta(minutes=5),
    )
    pom.use_floating_algae = False
    pom.use_benthic_algae = False
    return pom


def test_pom_dissolution_matches_v1_POM_dissolution(
    pom_5cell, water_temp_5cell, poc_5cell, dummy_time, dummy_registry
):
    """v3 POM dissolution: lumped rate == -v1 ``POM_dissolution``
    (mg-D/L_sed/d); cached ``pom_doc_source_rate`` ==
    ``fcom * v1_dissolution * h2 / depth`` (mg-C/L_water/d, post-Phase-9.G)."""
    pom = _make_pom(kpom_20=0.1, kpom_theta=1.047)
    v3_rate = pom.rate(
        pom=pom_5cell,
        water_temperature=water_temp_5cell,
        poc=poc_5cell,
        time=dummy_time,
        registry=dummy_registry,
    )

    np.testing.assert_allclose(
        np.asarray(v3_rate), -V1_DISSOLUTION_REFERENCE, rtol=1e-6
    )
    fcom = float(pom.fcom)
    h2 = float(pom.h2)
    depth = 1.0  # from dummy_registry
    expected_doc_source = V1_DISSOLUTION_REFERENCE * fcom * h2 / depth
    np.testing.assert_allclose(
        np.asarray(pom.pom_doc_source_rate),
        expected_doc_source,
        rtol=1e-6,
    )


def test_pom_burial_matches_v1_POM_burial(
    pom_5cell, water_temp_5cell, poc_5cell, dummy_time, dummy_registry
):
    """v3 POM burial: lumped rate == -v1 ``POM_burial``."""
    pom = _make_pom(vb=0.05, h2=0.1)
    v3_rate = pom.rate(
        pom=pom_5cell,
        water_temperature=water_temp_5cell,
        poc=poc_5cell,
        time=dummy_time,
        registry=dummy_registry,
    )
    np.testing.assert_allclose(
        np.asarray(v3_rate), -V1_BURIAL_REFERENCE, rtol=1e-6
    )


def test_pom_poc_settling_matches_v1_POM_POC_settling(
    pom_5cell, water_temp_5cell, poc_5cell, dummy_time, dummy_registry
):
    """v3 POC-settling source: lumped rate == +v1 ``POM_POC_settling``."""
    pom = _make_pom(vsoc=0.02, fcom=0.4, h2=0.1, use_POC=True)
    v3_rate = pom.rate(
        pom=pom_5cell,
        water_temperature=water_temp_5cell,
        poc=poc_5cell,
        time=dummy_time,
        registry=dummy_registry,
    )
    np.testing.assert_allclose(
        np.asarray(v3_rate), V1_POC_SETTLING_REFERENCE, rtol=1e-6
    )


def test_pom_algal_settling_input_matches_v1_POM_algal_settling(
    algae_5cell, water_temp_5cell
):
    """v3 FloatingAlgae ``algal_pom_from_settling_rate`` matches frozen
    v1 ``POM_algal_settling = vsap * Ap * (AWd/AWa) / h2`` reference."""
    fa = FloatingAlgae(
        time_step=timedelta(minutes=5),
        death_rate=0.15,
        death_rate_correction_factor=1.047,
    )
    fa.use_nitrate = True
    fa.use_ammonium = True
    fa.use_phosphate = True

    fa._cache_mortality_rates(algae_5cell, water_temp_5cell)

    np.testing.assert_allclose(
        np.asarray(fa.algal_pom_from_settling_rate),
        V1_ALGAL_SETTLING_REFERENCE,
        rtol=1e-6,
    )


def test_pom_benthic_mortality_input_matches_v1_POM_benthic_algae_mortality(
    benthic_algae_5cell, water_temp_5cell
):
    """v3 BenthicAlgae ``balgae_pom_from_mortality_rate`` matches frozen
    v1 ``POM_benthic_algae_mortality = Ab * kdb_tc * Fb * (1 - Fw) / h2``."""
    ba = BenthicAlgae(
        time_step=timedelta(minutes=5),
        death_rate=0.3,
        death_rate_correction_factor=1.047,
    )
    ba.use_nitrate = True
    ba.use_ammonium = True
    ba.use_phosphate = True

    depth = xr.DataArray(np.array([1.0, 1.5, 2.0, 2.5, 3.0]))
    ba._cache_benthic_mortality_rates(benthic_algae_5cell, water_temp_5cell, depth)

    np.testing.assert_allclose(
        np.asarray(ba.balgae_pom_from_mortality_rate),
        V1_BENTHIC_MORTALITY_REFERENCE,
        rtol=1e-6,
    )


def test_phase9fa_vb_value_pinned():
    """Phase 9.F.A: vb = 6.85e-6 m/d (= 0.0025 m/yr = 0.25 cm/yr)."""
    from clearwater_modules_v3.parameters.global_vars import (
        DEFAULTS as GLOBAL_VAR_DEFAULTS,
    )

    vb = GLOBAL_VAR_DEFAULTS["vb"]

    canonical = 0.0025 / 365.0
    np.testing.assert_allclose(vb, canonical, rtol=1e-3)
    assert vb < 0.01 / 100


def test_phase9fa_vb_dimensional_smell_test():
    """Phase 9.F.A: burial timescale vb/h2 is physically reasonable."""
    from clearwater_modules_v3.parameters.global_vars import (
        DEFAULTS as GLOBAL_VAR_DEFAULTS,
    )
    from clearwater_modules_v3.parameters.pom import (
        DEFAULTS as POM_DEFAULTS,
    )

    vb = GLOBAL_VAR_DEFAULTS["vb"]
    h2 = POM_DEFAULTS["h2"]

    burial_rate_per_day = vb / h2  # 1/d
    burial_timescale_days = 1.0 / burial_rate_per_day
    burial_timescale_years = burial_timescale_days / 365.0

    # ~ 40 years e-folding timescale (h2 = 0.1 m, vb = 6.85e-6 m/d)
    assert 30.0 < burial_timescale_years < 50.0, (
        f"Expected ~40 yr burial timescale; got {burial_timescale_years} yr"
    )
