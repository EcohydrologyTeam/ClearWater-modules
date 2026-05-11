"""v3 CBOD kinetic regression against frozen v1 reference values.

Each test constructs a v3 ``CBOD`` instance, drives ``run`` against an
in-memory registry, and compares the cached step-scoped rate
attributes (``cbod_oxidation_rate``, ``cbod_settling_rate``) to v1
reference arrays captured from the legacy ``clearwater_modules.nsm1``
helpers (``kbod_tc``, ``ksbod_tc``, ``CBOD_oxidation``,
``CBOD_sedimentation``) using the same 5-cell fixture inputs.

Migration history (Phase 2 of the v1 retirement plan,
``design/v1_retirement_plan.md``):

- Originally ``tests/test_5_cbod_calculations_v2.py``: imported v1
  directly to compute the reference inline. That required keeping
  ``src/clearwater_modules/nsm1/`` source alive.
- This v3 file: the v1-side reference values are frozen as numpy
  array literals (``V1_*_REFERENCE``) captured by running v1 against
  the test fixtures once. Decouples the test from v1 source, so v1
  can be retired without breaking the regression.

If a future v3 change intentionally diverges from v1 (e.g. a corrected
formula in a Phase 9.x audit), update the reference literal here and
document why in the commit message. The intent is bit-exact match at
``rtol=1e-6`` until a deliberate divergence is approved.

v3 deviation note (settling): v1 ``CBOD_sedimentation = CBOD * ksbod_tc``
(units treat ``ksbod_tc`` as 1/d). v3 ``cbod_settling_rate =
ksbod_tc / depth * cbod`` (treats ``ksbod_tc`` as m/d and divides by
depth to get 1/d). The two forms differ by ``1/depth``. The
``V1_SETTLING_PER_DEPTH_REFERENCE`` reference below is v1's
sedimentation flux divided by the 5-cell depth fixture so the
assertion holds. The default ``ksbod_20=0.0`` makes this term
identically zero in production use anyway.
"""
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.cbod import CBOD


# ---------------------------------------------------------------------------
# Frozen v1 reference values (captured 2026-05-10 from
# clearwater_modules.nsm1.processes against the 5-cell fixtures below)
# ---------------------------------------------------------------------------

V1_OXIDATION_WITH_DOX_REFERENCE = np.array([
    0.1695607429727131,
    0.2526188682295775,
    0.336,
    0.4333249694117647,
    0.5751555921142887,
])

V1_OXIDATION_NO_DOX_REFERENCE = np.array([
    0.19075583584430222,
    0.2736704405820423,
    0.36,
    0.46040777999999993,
    0.6039133717200031,
])

V1_SETTLING_PER_DEPTH_REFERENCE = np.array([
    0.15896319653691854,
    0.11402935024251763,
    0.10000000000000002,
    0.09591828749999998,
    0.08387685718333376,
])


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


class _StubRegistry:
    """Minimal stand-in for VariableRegistry; supports get_at_time /
    set_at_time / __contains__."""

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
    reg = _StubRegistry()
    reg.register("cbod", cbod_5cell)
    reg.register("water_temperature", water_temp_5cell)
    reg.register("depth", depth_5cell)
    reg.register("oxygen_dissolved", dox_5cell)
    return reg


@pytest.fixture(scope="function")
def time_zero() -> datetime:
    return datetime(2026, 1, 1)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_cbod_oxidation_matches_v1_CBOD_oxidation(loaded_registry, time_zero):
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

    np.testing.assert_allclose(
        np.asarray(cbod.cbod_oxidation_rate),
        V1_OXIDATION_WITH_DOX_REFERENCE,
        rtol=1e-6,
    )


def test_cbod_oxidation_no_dox_matches_v1_first_order(loaded_registry, time_zero):
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

    np.testing.assert_allclose(
        np.asarray(cbod.cbod_oxidation_rate),
        V1_OXIDATION_NO_DOX_REFERENCE,
        rtol=1e-6,
    )


def test_cbod_settling_matches_v1_CBOD_sedimentation_per_depth(
    loaded_registry, time_zero
):
    """v3 cached ``cbod_settling_rate`` matches v1
    ``CBOD_sedimentation = CBOD * ksbod_tc`` divided by depth (see
    module-level deviation note)."""
    cbod = CBOD(
        parameters={
            "kbod_20": 0.12,
            "kbod_theta": 1.047,
            "ksbod_20": 0.05,
            "ksbod_theta": 1.047,
            "KsOxbod": 0.5,
        },
        time_step=timedelta(minutes=5),
    )
    cbod.use_DOX = True
    cbod.run(time_zero, loaded_registry)

    np.testing.assert_allclose(
        np.asarray(cbod.cbod_settling_rate),
        V1_SETTLING_PER_DEPTH_REFERENCE,
        rtol=1e-6,
    )
