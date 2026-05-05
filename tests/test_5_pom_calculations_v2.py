"""Parity tests: v3 POM sub-rate terms vs v1 nsm1.processes helpers.

Each test constructs a v3 ``POM`` instance, drives its lumped ``rate``
method (or in two cases the FloatingAlgae / BenthicAlgae cached
attribute that POM consumes), with all-but-one source/sink terms
zeroed out, and compares the isolated v3 term to the equivalent v1
helper-function output computed with the same inputs.

v1 reference: ``clearwater_modules.nsm1.processes`` ``kpom_tc``,
``POM_dissolution``, ``POM_burial``, ``POM_POC_settling``,
``POM_algal_settling``, ``POM_benthic_algae_mortality``, ``dPOMdt``.

The integrator branch in ``POM.run`` (Forward Euler + clip-with-log +
registry write) is exercised in ``tests/v3/nsm1/test_pom_tier1.py``;
this module covers the kinetic forms only.

Synthetic mesh: 5-cell numpy arrays.
"""
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v2.processes.benthic_algae import BenthicAlgae
from clearwater_modules_v2.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.pom import POM


@pytest.fixture(scope="module")
def pom_5cell():
    return xr.DataArray(np.array([3.0, 3.5, 4.0, 4.5, 5.0]))


@pytest.fixture(scope="module")
def poc_5cell():
    return xr.DataArray(np.array([1.0, 1.2, 1.4, 1.6, 1.8]))


@pytest.fixture(scope="module")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="module")
def algae_5cell():
    return xr.DataArray(np.array([5.0, 6.0, 7.0, 8.0, 10.0]))


@pytest.fixture(scope="module")
def benthic_algae_5cell():
    return xr.DataArray(np.array([1.0, 1.5, 2.0, 2.5, 3.0]))


# A dummy time / registry; POM.rate doesn't read them in the no-coupling
# branches we exercise here, but the signature requires them.
@pytest.fixture(scope="module")
def dummy_time() -> datetime:
    return datetime(2026, 1, 1)


@pytest.fixture(scope="module")
def dummy_registry():
    """A trivial dummy registry; POM.rate does not consult it when the
    floating/benthic algae coupling is disabled and POC is supplied
    directly as a kwarg."""
    class _Stub:
        def __contains__(self, name: str) -> bool:
            return False

        def get_at_time(self, name, time):
            raise KeyError(name)
    return _Stub()


def _make_pom(
    *,
    kpom_20: float = 0.0,
    kpom_theta: float = 1.0,
    vb: float = 0.0,
    vsoc: float = 0.0,
    fcom: float = 0.4,
    h2: float = 0.1,
    use_POC: bool = False,
    use_Algae: bool = False,
    use_Balgae: bool = False,
) -> POM:
    """Helper to build a POM instance with arbitrary parameter overrides."""
    pom = POM(
        parameters={
            "kpom_20": kpom_20,
            "kpom_theta": kpom_theta,
            "vb": vb,
            "vsoc": vsoc,
            "fcom": fcom,
            "h2": h2,
            "use_POC": use_POC,
            "use_Algae": use_Algae,
            "use_Balgae": use_Balgae,
        },
        time_step=timedelta(minutes=5),
    )
    # No model wired -- ensure flags align with the test scenario.
    pom.use_floating_algae = False
    pom.use_benthic_algae = False
    return pom


def test_pom_dissolution_matches_v1_POM_dissolution(
    pom_5cell, water_temp_5cell, poc_5cell, dummy_time, dummy_registry
):
    """v3 POM dissolution term matches v1 ``POM_dissolution = POM * kpom_tc``.

    Setup: enable dissolution only; vb=vsoc=0, no algal coupling. The
    lumped ``rate`` then equals ``-POM_dissolution`` (negative, since
    dissolution is a sink). The cached ``pom_hydrolysis_rate`` is the
    *positive* dissolution flux (mg/L/d), which is what we assert
    against v1.
    """
    pom = _make_pom(kpom_20=0.1, kpom_theta=1.047)
    v3_rate = pom.rate(
        pom=pom_5cell,
        water_temperature=water_temp_5cell,
        poc=poc_5cell,
        time=dummy_time,
        registry=dummy_registry,
    )

    kpom_tc = v1.kpom_tc(water_temp_5cell, 0.1, 1.047)
    v1_dissolution = v1.POM_dissolution(POM=pom_5cell, kpom_tc=kpom_tc)

    # The lumped v3 rate is -POM_dissolution under this configuration.
    np.testing.assert_allclose(
        np.asarray(v3_rate), -np.asarray(v1_dissolution), rtol=1e-6
    )
    # And the cached hydrolysis rate is the positive flux.
    np.testing.assert_allclose(
        np.asarray(pom.pom_hydrolysis_rate),
        np.asarray(v1_dissolution),
        rtol=1e-6,
    )


def test_pom_burial_matches_v1_POM_burial(
    pom_5cell, water_temp_5cell, poc_5cell, dummy_time, dummy_registry
):
    """v3 POM burial term matches v1 ``POM_burial = vb * POM / h2``.

    Setup: enable burial only; kpom_20=0, vsoc=0, no algal coupling.
    The lumped ``rate`` then equals ``-POM_burial`` (sink).
    """
    vb = 0.05
    h2 = 0.1
    pom = _make_pom(vb=vb, h2=h2)
    v3_rate = pom.rate(
        pom=pom_5cell,
        water_temperature=water_temp_5cell,
        poc=poc_5cell,
        time=dummy_time,
        registry=dummy_registry,
    )

    v1_burial = v1.POM_burial(vb=vb, POM=pom_5cell, h2=h2)

    np.testing.assert_allclose(
        np.asarray(v3_rate), -np.asarray(v1_burial), rtol=1e-6
    )


def test_pom_poc_settling_matches_v1_POM_POC_settling(
    pom_5cell, water_temp_5cell, poc_5cell, dummy_time, dummy_registry
):
    """v3 POC-settling source term matches v1
    ``POM_POC_settling = vsoc * POC / h2 / fcom``.

    Setup: enable POC settling only; kpom_20=vb=0, no algal coupling.
    The lumped ``rate`` then equals ``+POM_POC_settling`` (source).
    """
    vsoc = 0.02
    fcom = 0.4
    h2 = 0.1
    pom = _make_pom(vsoc=vsoc, fcom=fcom, h2=h2, use_POC=True)
    v3_rate = pom.rate(
        pom=pom_5cell,
        water_temperature=water_temp_5cell,
        poc=poc_5cell,
        time=dummy_time,
        registry=dummy_registry,
    )

    v1_poc_settling = v1.POM_POC_settling(
        POC=poc_5cell, vsoc=vsoc, h2=h2, fcom=fcom, use_POC=True
    )

    np.testing.assert_allclose(
        np.asarray(v3_rate), np.asarray(v1_poc_settling), rtol=1e-6
    )


def test_pom_algal_settling_input_matches_v1_POM_algal_settling(
    algae_5cell, water_temp_5cell
):
    """v3 FloatingAlgae cached ``algal_pom_from_settling_rate`` matches
    v1 ``POM_algal_settling = vsap * Ap * (AWd/AWa) / h2``.

    The Phase 3.5 inter-process coupling stores the consumer-ready POM
    source flux on the FloatingAlgae instance after each ``run`` (and on
    ``_cache_mortality_rates``). POM.rate then reads that attribute via
    ``getattr``. This test calls ``_cache_mortality_rates`` directly to
    avoid wiring a full Model and registry.
    """
    fa = FloatingAlgae(
        time_step=timedelta(minutes=5),
        death_rate=0.15,
        death_rate_correction_factor=1.047,
    )
    fa.use_nitrate = True
    fa.use_ammonium = True
    fa.use_phosphate = True

    # Trigger the cache.
    fa._cache_mortality_rates(algae_5cell, water_temp_5cell)

    # v1 reference: vsap * Ap * (AWd/AWa) / h2.
    rda = fa.AWd / fa.AWa
    v1_settling = v1.POM_algal_settling(
        Ap=algae_5cell, vsap=fa.vsap, rda=rda, h2=fa.h2, use_Algae=True
    )

    np.testing.assert_allclose(
        np.asarray(fa.algal_pom_from_settling_rate),
        np.asarray(v1_settling),
        rtol=1e-6,
    )


def test_pom_benthic_mortality_input_matches_v1_POM_benthic_algae_mortality(
    benthic_algae_5cell, water_temp_5cell
):
    """v3 BenthicAlgae cached ``balgae_pom_from_mortality_rate`` matches
    v1 ``POM_benthic_algae_mortality = Ab * kdb_tc * Fb * (1 - Fw) / h2``.

    Mirrors the Phase 3.5 inter-process coupling test for floating
    algae. ``_cache_benthic_mortality_rates`` is called directly with
    arbitrary depth (the POM cache term does not consume depth — only
    the OrgN/OrgP/POC/DOC routings do — but the helper signature takes
    it).
    """
    # BenthicAlgae's inherited ``rate_death`` uses legacy v2 kwargs
    # death_rate / death_rate_correction_factor, so we set those to the
    # v1-equivalent kdb_20 / kdb_theta values.
    ba = BenthicAlgae(
        time_step=timedelta(minutes=5),
        death_rate=0.3,
        death_rate_correction_factor=1.047,
    )
    ba.use_nitrate = True
    ba.use_ammonium = True
    ba.use_phosphate = True

    # Depth dummy (only used by sister mortality routings, not the POM term).
    depth = xr.DataArray(np.array([1.0, 1.5, 2.0, 2.5, 3.0]))
    ba._cache_benthic_mortality_rates(benthic_algae_5cell, water_temp_5cell, depth)

    # v1 reference. kdb_tc uses water-temp Arrhenius correction with
    # kdb_20=0.3, kdb_theta=1.047 (the values set via the legacy kwargs
    # above). NOTE: v1's ``kdb_tc`` signature is
    # ``(kdb_20, TwaterC, kdb_theta)`` -- positional order differs from
    # the analogous ``kdp_tc(TwaterC, kdp_20, kdp_theta)`` signature.
    kdb_tc = v1.kdb_tc(0.3, water_temp_5cell, 1.047)
    v1_mortality = v1.POM_benthic_algae_mortality(
        Ab=benthic_algae_5cell,
        kdb_tc=kdb_tc,
        Fb=ba.Fb,
        Fw=ba.Fw,
        h2=ba.h2,
        use_Balgae=True,
    )

    np.testing.assert_allclose(
        np.asarray(ba.balgae_pom_from_mortality_rate),
        np.asarray(v1_mortality),
        rtol=1e-6,
    )
