"""
=========================================================================================
Unit Tests - Water Temperature Simulation Module (TSM): Temperature Energy Budget Module
=========================================================================================

Developed by:
Dr. Todd E. Steissberg (ERDC-EL)

Date: April 11, 2021
"""
import pytest
from clearwater_python_modules.gsm import GeneralConstituentKinetics

@pytest.fixture(scope='module')
def tolerance() -> float:
    """Return the tolerance for comparing floating point numbers."""
    return 0.0001


@pytest.fixture(scope='module')
def params() -> dict:
    """Return a dictionary of parameters for the GeneralConstituentKinetics class.

    Key definitions:
    GC: General constituent concentration
    TwaterC: Water temperature in degrees Celsius.
    order: Compute 1st order kinetics.
    k_rc20: Reaction rate at 20 degrees Celsius, decay.
    k_theta: Arrhenius temperature correction factor, decay.
    rgc_rc20: Reaction rate at 20 degrees Celsius, settling.
    rgc_theta: Arrhenius temperature correction factor, settling.
    release: Turn suspension on.
    settling: Turn settling on.
    depth: Bed depth.
    settling_rate: Settling rate.
    """
    return {
        "GC": 10.0,         
        "TwaterC": 25.0,
        "order": 1,
        "k_rc20": 0.5,
        "k_theta": 1.047,
        "rgc_rc20": 0.5,
        "rgc_theta": 1.047,
        "release": True,
        "settling": True,
        "depth": 1.0,
        "settling_rate": 0.0002,
    }


@pytest.fixture(scope='module')
def make_gck_instance(params: dict) -> GeneralConstituentKinetics:
    """Create a GeneralConstituentKinetics instance with the given parameters."""
    return GeneralConstituentKinetics(
        params["GC"],
        params["TwaterC"],
        params["order"],
        params["k_rc20"],
        params["k_theta"],
        params["rgc_rc20"],
        params["rgc_theta"],
        release=params["release"],
        settling=params["settling"],
        depth=params["depth"],
        settling_rate=params["settling_rate"]
    )


def test_GC(make_gck_instance, params) -> None:
    """"Test the GeneralConstituentKinetics class with different GC values."""
    params["GC"] = 0.0
    assert pytest.approx(
        make_gck_instance(params), 
        tolerance,
    ) == 0.6291

    params["GC"] = 10.0
    assert pytest.approx(
        make_gck_instance(
            params,
        ), 
        tolerance,
    ) == -5.6637

    params["GC"] = 20.0
    assert pytest.approx(
        make_gck_instance(
            params,
        ), 
        tolerance,
    ) == -11.9565

    params["GC"] = 30.0
    assert pytest.approx(
        make_gck_instance(
            params,
        ), 
        tolerance,
    ) == -18.2492

def test_TwaterC(make_gck_instance, params):
    """Test the GeneralConstituentKinetics class with different TwaterC values."""
    params["TwaterC"] = 0.0

    assert pytest.approx(
        make_gck_instance(
            params,
        ), 
        tolerance,
    ) == -1.7979
    
    params["TwaterC"] = 10.0 
    assert pytest.approx(
        make_gck_instance(
            params,
        ), 
        tolerance,
    ) == -2.84484

    params["TwaterC"] = 20.0
    assert pytest.approx(
        make_gck_instance(
            params,
        ),
        tolerance,
    ) == -4.5020

    params["TwaterC"] = 30.0
    assert pytest.approx(
        make_gck_instance(
            params,
        ),
        tolerance,
    ) == -7.1253 

