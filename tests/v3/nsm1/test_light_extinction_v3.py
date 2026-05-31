"""FloatingAlgae computes the light-extinction coefficient lambda.

Guards `design/clearwater_modules_v3_light_extinction.md`: FloatingAlgae's
`limit_light` previously used a constant `light_attenuation_coefficient`
(1.0 /m). It now computes lambda each step via `utils.light.L`, a verified
port of the NSM1-I Fortran assembly
(`fortran/NSM1/02_global/nsmi_global_params.f90:421-427`):

    lambda = lambda0 + lambdas*Solid
                     + (use_POC)   lambdam*POC/fcom
                     + (use_Algae) lambda1*Ap + lambda2*Ap**(2/3)

with the NSM1-I default coefficients (lambda0=0.02, lambdas=0.052,
lambdam=0.174, lambda1=0.0088, lambda2=0.054, fcom=0.4). The constant
`light_attenuation_coefficient` is retained as a scalar override via
`use_computed_light_extinction=False` (tests / didactic runs).
"""
from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.examples import InMemoryRegistry
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.utils.light import L as compute_light_extinction


N = 4
START = datetime(2026, 1, 1, 0, 0, 0)

# NSM1-I optical defaults (must match _LIGHT_EXT_DEFAULTS / the Fortran).
L0, L1, L2, LS, LM, FCOM = 0.02, 0.0088, 0.054, 0.052, 0.174, 0.4


def _da(value) -> xr.DataArray:
    return xr.DataArray(np.full(N, value, dtype=float), dims=["cell"])


def _run(*, algae, solid=None, poc=None, use_computed=True, scalar=1.0):
    proc = FloatingAlgae(
        use_computed_light_extinction=use_computed,
        light_attenuation_coefficient=scalar,
    )
    reg = InMemoryRegistry()
    reg.register("algae_floating", _da(algae))
    reg.register("ammonium", _da(0.10))
    reg.register("nitrate", _da(0.20))
    reg.register("tip", _da(0.10))
    reg.register("depth", _da(1.0))
    reg.register("water_temperature", _da(20.0))
    reg.register("solar_radiation", _da(300.0))
    if solid is not None:
        reg.register("Solid", _da(solid))
    if poc is not None:
        reg.register("poc", _da(poc))
    proc.run(START, reg)
    return proc


def test_computed_is_the_default():
    """use_computed_light_extinction defaults to True."""
    assert FloatingAlgae().use_computed_light_extinction is True


def test_optical_defaults_match_nsm1():
    """The optical coefficients on the instance match the NSM1-I defaults."""
    p = FloatingAlgae()
    assert (p.lambda0, p.lambda1, p.lambda2, p.lambdas, p.lambdam, p.fcom) == (
        L0, L1, L2, LS, LM, FCOM
    )
    assert p.use_POC is True and p.use_Algae is True


def test_computed_lambda_matches_utils_light_L():
    """The per-step lambda equals utils.light.L for the same inputs
    (which is the Fortran-verified formula)."""
    algae, solid, poc = 40.0, 6.0, 4.0
    proc = _run(algae=algae, solid=solid, poc=poc)

    expected = compute_light_extinction(
        lambda0=L0, lambda1=L1, lambda2=L2, lambdas=LS, lambdam=LM,
        Solid=_da(solid), POC=_da(poc), fcom=FCOM, Ap=_da(algae),
        use_Algae=True, use_POC=True,
    )
    np.testing.assert_allclose(
        np.asarray(proc._light_extinction), np.asarray(expected.values),
        rtol=1e-12, atol=0.0,
    )
    # And the closed-form value, for documentation.
    hand = L0 + LS * solid + LM * poc / FCOM + L1 * algae + L2 * algae ** 0.66667
    np.testing.assert_allclose(np.asarray(proc._light_extinction), hand, rtol=1e-6)


def test_lambda_increases_with_algae_self_shading():
    """More algae -> larger lambda (self-shading term)."""
    lo = np.asarray(_run(algae=5.0, solid=0.0, poc=0.0)._light_extinction)
    hi = np.asarray(_run(algae=80.0, solid=0.0, poc=0.0)._light_extinction)
    assert np.all(hi > lo)


def test_poc_absent_drops_the_poc_term():
    """With no POC registered, the POC attenuation term is absent (POC=0)."""
    proc = _run(algae=20.0, solid=2.0, poc=None)  # poc not registered
    expected = L0 + LS * 2.0 + L1 * 20.0 + L2 * 20.0 ** 0.66667  # no POC term
    np.testing.assert_allclose(np.asarray(proc._light_extinction), expected, rtol=1e-6)


def test_solid_absent_uses_scalar_fallback_in_lambda():
    """With no Solid registered, the lambda uses the self.Solid fallback
    (1.0 mg/L) for the suspended-solids term."""
    proc = _run(algae=10.0, solid=None, poc=0.0)
    expected = L0 + LS * 1.0 + L1 * 10.0 + L2 * 10.0 ** 0.66667  # Solid fallback 1.0
    np.testing.assert_allclose(np.asarray(proc._light_extinction), expected, rtol=1e-6)


def test_scalar_override_when_computed_disabled():
    """use_computed_light_extinction=False -> _light_extinction is the
    constant scalar, regardless of the optical state."""
    proc = _run(algae=40.0, solid=6.0, poc=4.0, use_computed=False, scalar=0.7)
    assert np.all(np.asarray(proc._light_extinction) == 0.7)


def test_computed_vs_override_change_light_limitation():
    """The computed lambda (>> 1 for a bloom) yields a different light
    limitation than the constant-lambda override -> the wiring is live."""
    computed = _run(algae=40.0, solid=6.0, poc=4.0, use_computed=True)
    override = _run(algae=40.0, solid=6.0, poc=4.0, use_computed=False, scalar=1.0)
    a = np.asarray(computed._light_extinction)
    b = np.asarray(override._light_extinction)
    assert np.all(a > b), "computed bloom lambda should exceed the 1.0 /m constant"
