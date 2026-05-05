"""Parity tests: v3 CBOD sub-rate cached attributes vs v1 nsm1.processes helpers.

Each test constructs a v3 ``CBOD`` instance, drives ``run`` against an
in-memory registry, and compares the cached step-scoped rate
attributes (``cbod_oxidation_rate``, ``cbod_settling_rate``) to the
equivalent v1 helper-function output computed with the same inputs.

v1 reference: ``clearwater_modules.nsm1.processes`` ``kbod_tc``,
``ksbod_tc``, ``CBOD_oxidation``, ``CBOD_sedimentation``, ``dCBODdt``.

v3 deviation note (settling): v1 ``CBOD_sedimentation = CBOD * ksbod_tc``
(units treat ``ksbod_tc`` as 1/d). v3 ``cbod_settling_rate =
ksbod_tc / depth * cbod`` (treats ``ksbod_tc`` as m/d and divides by
depth to get 1/d). The two forms differ by the factor ``1/depth``.
The Phase 3.3 spec documents this as the intentional v3 convention.
The settling parity test below divides the v1 reference by depth so
the assertion holds. The default ``ksbod_20=0.0`` makes this term
identically zero in production use anyway.

The integrator branch (Forward Euler + clip-with-log + registry
set_at_time) is exercised in the Phase 3 Tier 1 conservation tests.
This module covers only the kinetic forms.

Synthetic mesh: 5-cell numpy arrays in an in-memory registry.
"""
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.cbod import CBOD


class _StubRegistry:
    """Minimal stand-in for VariableRegistry; supports get_at_time /
    set_at_time / __contains__. Mirrors the InMemoryRegistry pattern in
    ``tests/v3/nsm1/conftest.py`` but lives inline here so this file
    does not depend on the v3 conftest fixtures.
    """

    def __init__(self) -> None:
        self._data: dict[str, xr.DataArray] = {}

    def register(self, name: str, value: xr.DataArray) -> None:
        self._data[name] = value

    def get_at_time(self, name: str, time: datetime) -> xr.DataArray:
        if name not in self._data:
            raise KeyError(name)
        return self._data[name]

    def set_at_time(self, name: str, time: datetime, value: xr.DataArray) -> None:
        self._data[name] = value

    def __contains__(self, name: str) -> bool:
        return name in self._data


@pytest.fixture(scope="function")
def cbod_5cell():
    return xr.DataArray(np.array([2.0, 2.5, 3.0, 3.5, 4.0]))


@pytest.fixture(scope="function")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="function")
def depth_5cell():
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))


@pytest.fixture(scope="function")
def dox_5cell():
    return xr.DataArray(np.array([4.0, 6.0, 7.0, 8.0, 10.0]))


@pytest.fixture(scope="function")
def loaded_registry(cbod_5cell, water_temp_5cell, depth_5cell, dox_5cell):
    """Stub registry pre-loaded with cbod / water_temperature / depth / DOX."""
    reg = _StubRegistry()
    reg.register("cbod", cbod_5cell)
    reg.register("water_temperature", water_temp_5cell)
    reg.register("depth", depth_5cell)
    reg.register("oxygen_dissolved", dox_5cell)
    return reg


@pytest.fixture(scope="function")
def time_zero() -> datetime:
    return datetime(2026, 1, 1)


def test_cbod_oxidation_matches_v1_CBOD_oxidation(
    loaded_registry, cbod_5cell, water_temp_5cell, dox_5cell, time_zero
):
    """v3 cached ``cbod_oxidation_rate`` matches v1
    ``CBOD_oxidation = (DOX / (KsOxbod + DOX)) * kbod_tc * CBOD``.
    """
    cbod = CBOD(
        parameters={
            "kbod_20": 0.12,
            "kbod_theta": 1.047,
            "ksbod_20": 0.0,
            "ksbod_theta": 1.047,
            "KsOxbod": 0.5,
        },
        time_step=timedelta(minutes=5),
    )
    cbod.use_DOX = True
    cbod.run(time_zero, loaded_registry)

    kbod_tc = v1.kbod_tc(water_temp_5cell, 0.12, 1.047)
    v1_oxidation = v1.CBOD_oxidation(
        DOX=dox_5cell,
        CBOD=cbod_5cell,
        kbod_tc=kbod_tc,
        KsOxbod=0.5,
        use_DOX=True,
    )

    np.testing.assert_allclose(
        np.asarray(cbod.cbod_oxidation_rate),
        np.asarray(v1_oxidation),
        rtol=1e-6,
    )


def test_cbod_oxidation_no_dox_matches_v1_first_order(
    loaded_registry, cbod_5cell, water_temp_5cell, time_zero
):
    """When ``use_DOX=False`` v3 ``cbod_oxidation_rate`` matches v1's
    first-order branch ``kbod_tc * CBOD`` (DOX attenuation off)."""
    cbod = CBOD(
        parameters={
            "kbod_20": 0.12,
            "kbod_theta": 1.047,
            "ksbod_20": 0.0,
            "ksbod_theta": 1.047,
            "KsOxbod": 0.5,
        },
        time_step=timedelta(minutes=5),
    )
    cbod.use_DOX = False
    cbod.run(time_zero, loaded_registry)

    kbod_tc = v1.kbod_tc(water_temp_5cell, 0.12, 1.047)
    # v1.CBOD_oxidation with use_DOX=False -> kbod_tc * CBOD.
    v1_oxidation = v1.CBOD_oxidation(
        DOX=xr.DataArray(np.array([4.0, 6.0, 7.0, 8.0, 10.0])),
        CBOD=cbod_5cell,
        kbod_tc=kbod_tc,
        KsOxbod=0.5,
        use_DOX=False,
    )

    np.testing.assert_allclose(
        np.asarray(cbod.cbod_oxidation_rate),
        np.asarray(v1_oxidation),
        rtol=1e-6,
    )


def test_cbod_settling_matches_v1_CBOD_sedimentation_per_depth(
    loaded_registry, cbod_5cell, water_temp_5cell, depth_5cell, time_zero
):
    """v3 cached ``cbod_settling_rate`` matches v1
    ``CBOD_sedimentation = CBOD * ksbod_tc`` divided by depth.

    See module-level deviation note: v1 treats ksbod_tc as 1/d while v3
    treats it as m/d and divides by depth. The parity assertion divides
    the v1 reference by depth to compare apples-to-apples.
    """
    # Use a non-zero ksbod_20 so the term actually exercises the
    # settling branch (the v3 default ksbod_20=0 would zero it out).
    ksbod_20 = 0.05
    ksbod_theta = 1.047
    cbod = CBOD(
        parameters={
            "kbod_20": 0.12,
            "kbod_theta": 1.047,
            "ksbod_20": ksbod_20,
            "ksbod_theta": ksbod_theta,
            "KsOxbod": 0.5,
        },
        time_step=timedelta(minutes=5),
    )
    cbod.use_DOX = True
    cbod.run(time_zero, loaded_registry)

    ksbod_tc = v1.ksbod_tc(water_temp_5cell, ksbod_20, ksbod_theta)
    v1_sedimentation = v1.CBOD_sedimentation(CBOD=cbod_5cell, ksbod_tc=ksbod_tc)
    v1_settling_per_depth = v1_sedimentation / depth_5cell

    np.testing.assert_allclose(
        np.asarray(cbod.cbod_settling_rate),
        np.asarray(v1_settling_per_depth),
        rtol=1e-6,
    )
