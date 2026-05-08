import pytest
from clearwater_modules_v3.processes.temperature import Temperature

import xarray as xr


@pytest.fixture(scope="module")
def temperature_instance():
    # Wind constants are not used for Richardson tests
    return Temperature(
        1,
        1,
        1,
    )


def test_rnumber_zero(temperature_instance):
    # if density_air_sat == density_air then rnumber should be zero
    # and richardson function should be calculated to stable (1)
    rnumber, rfunction = temperature_instance.richardson_number(1, 1, 1)
    assert rnumber == 0.0
    assert rfunction == 1.0


def test_rnumber_zero_array(temperature_instance):
    # if density_air_sat == density_air then rnumber should be zero
    # and richardson function should be calculated to stable (1)
    data = xr.DataArray([1, 1, 1])
    rnumber, rfunction = temperature_instance.richardson_number(data, data, data)
    xr.testing.assert_equal(rnumber, xr.DataArray([0.0, 0.0, 0.0]))
    xr.testing.assert_equal(rfunction, xr.DataArray([1.0, 1.0, 1.0]))
