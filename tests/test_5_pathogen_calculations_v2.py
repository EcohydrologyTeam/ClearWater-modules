"""Parity tests: v3 Pathogen sub-rate methods vs v1 nsm1.processes helpers.

Each test constructs a v3 ``Pathogen`` instance, calls one of its rate
sub-methods directly, and compares to the equivalent v1 helper-function
output computed with the same inputs.

Scope: the v1-equivalent sub-rate computations on Pathogen
(``_rate_natural_decay``, ``_rate_light_decay``, ``_rate_settling``)
plus the lumped ``rate`` aggregator. The integrator branch in
``Pathogen.run`` (Forward Euler + clip-with-log + registry write) is
exercised in ``tests/v3/nsm1/test_pathogen_tier1.py``.

v1 reference: ``clearwater_modules.nsm1.processes`` ``kdx_tc``,
``PathogenDeath``, ``PathogenDecay``, ``PathogenSettling``, ``dPXdt``.

v3 light-decay note: Phase 3.1 originally substituted
``I0 = q_solar * Fr_PAR`` for the surface irradiance in
``_rate_light_decay``. Phase 9.F.B reverted that substitution because
pathogen inactivation is largely UVA/UVB-mediated (not PAR-mediated)
and the canonical Auer & Niehaus (1993) / Chapra (1997) calibration
operates on total broadband solar radiation. v3 now uses raw
``q_solar`` directly, matching v1 ``PathogenDecay`` exactly at the
kinetics level. The fixture below still pins ``Fr_PAR=1.0`` for
defensive consistency (it is now a no-op in ``_rate_light_decay``).

Synthetic mesh: 5-cell numpy arrays, single time step.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.pathogen import Pathogen
from clearwater_modules_v3.utils.light import L


@pytest.fixture(scope="module")
def px_5cell():
    return xr.DataArray(np.array([1.0e3, 5.0e3, 1.0e4, 5.0e4, 1.0e5]))


@pytest.fixture(scope="module")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="module")
def depth_5cell():
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))


@pytest.fixture(scope="module")
def q_solar_5cell():
    return xr.DataArray(np.array([200.0, 250.0, 300.0, 350.0, 400.0]))


@pytest.fixture(scope="module")
def solid_5cell():
    return xr.DataArray(np.array([10.0, 12.0, 15.0, 18.0, 20.0]))


@pytest.fixture(scope="module")
def poc_5cell():
    return xr.DataArray(np.array([1.0, 1.2, 1.4, 1.6, 1.8]))


@pytest.fixture(scope="module")
def ap_5cell():
    return xr.DataArray(np.array([5.0, 6.0, 7.0, 8.0, 10.0]))


@pytest.fixture(scope="function")
def pathogen_instance():
    """Pathogen instance with v1-aligned defaults.

    Pins ``apx=1.0`` and ``vx=1.0`` (v1 placeholder values) so the
    parity tests compare like-with-like against v1 helper outputs at
    those same placeholder values. The Phase 9.F.B canonical defaults
    (``apx=0.017``, ``vx=1.38``) are exercised in their own pinning
    tests at the bottom of this module.

    ``Fr_PAR=1.0`` is set defensively; Phase 9.F.B reverted the
    Phase 3.1 PAR substitution in ``_rate_light_decay``, so this is
    now a no-op.
    """
    return Pathogen(
        parameters={
            "kdx_20": 0.8,
            "kdx_theta": 1.07,
            "apx": 1.0,
            "vx": 1.0,
            "Fr_PAR": 1.0,
        },
        time_step=timedelta(minutes=5),
    )


def test_pathogen_natural_decay_matches_v1_PathogenDeath(
    pathogen_instance, px_5cell, water_temp_5cell
):
    """v3 ``_rate_natural_decay`` == v1 ``PathogenDeath(kdx_tc, PX)``."""
    v3_rate = pathogen_instance._rate_natural_decay(px_5cell, water_temp_5cell)

    kdx_tc = v1.kdx_tc(water_temp_5cell, 0.8, 1.07)
    v1_rate = v1.PathogenDeath(kdx_tc, px_5cell)

    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6)


def test_pathogen_light_decay_matches_v1_PathogenDecay(
    pathogen_instance,
    px_5cell,
    depth_5cell,
    q_solar_5cell,
    solid_5cell,
    poc_5cell,
    ap_5cell,
):
    """v3 ``_rate_light_decay`` == v1 ``PathogenDecay(apx, q_solar, L, depth, PX)``.

    With ``Fr_PAR=1.0`` (set in the fixture) the v3 effective surface
    irradiance equals ``q_solar``, matching v1 exactly at the kinetics
    level.
    """
    v3_rate = pathogen_instance._rate_light_decay(
        px=px_5cell,
        depth=depth_5cell,
        q_solar=q_solar_5cell,
        solid=solid_5cell,
        poc=poc_5cell,
        ap=ap_5cell,
    )

    # Compute v1 KEXT via the v3 utility (same Beer-Lambert formula as the
    # v1 shared L). Inputs match the Pathogen instance attributes.
    kext = L(
        lambda0=pathogen_instance.lambda0,
        lambda1=pathogen_instance.lambda1,
        lambda2=pathogen_instance.lambda2,
        lambdas=pathogen_instance.lambdas,
        lambdam=pathogen_instance.lambdam,
        Solid=solid_5cell,
        POC=poc_5cell,
        fcom=pathogen_instance.fcom,
        Ap=ap_5cell,
        use_Algae=pathogen_instance.use_Algae,
        use_POC=pathogen_instance.use_POC,
    )
    v1_rate = v1.PathogenDecay(
        apx=1.0,
        q_solar=q_solar_5cell,
        L=kext,
        depth=depth_5cell,
        PX=px_5cell,
    )

    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6)


def test_pathogen_settling_matches_v1_PathogenSettling(
    pathogen_instance, px_5cell, depth_5cell
):
    """v3 ``_rate_settling`` == v1 ``PathogenSettling(vx, depth, PX)``."""
    v3_rate = pathogen_instance._rate_settling(px_5cell, depth_5cell)
    v1_rate = v1.PathogenSettling(vx=1.0, depth=depth_5cell, PX=px_5cell)

    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6)


def test_pathogen_total_rate_matches_v1_dPXdt(
    pathogen_instance,
    px_5cell,
    water_temp_5cell,
    depth_5cell,
    q_solar_5cell,
    solid_5cell,
    poc_5cell,
    ap_5cell,
):
    """v3 ``rate`` == v1 ``dPXdt(PathogenDeath, PathogenDecay, PathogenSettling)``.

    v3 ``rate`` returns the *signed* rate of change ``dPX/dt`` (the
    sum of the three negative loss terms). v1 ``dPXdt`` is the same
    quantity computed by negating the three positive loss helpers.
    """
    v3_rate = pathogen_instance.rate(
        px=px_5cell,
        water_temperature=water_temp_5cell,
        depth=depth_5cell,
        q_solar=q_solar_5cell,
        solid=solid_5cell,
        poc=poc_5cell,
        ap=ap_5cell,
    )

    kdx_tc = v1.kdx_tc(water_temp_5cell, 0.8, 1.07)
    death = v1.PathogenDeath(kdx_tc, px_5cell)
    kext = L(
        lambda0=pathogen_instance.lambda0,
        lambda1=pathogen_instance.lambda1,
        lambda2=pathogen_instance.lambda2,
        lambdas=pathogen_instance.lambdas,
        lambdam=pathogen_instance.lambdam,
        Solid=solid_5cell,
        POC=poc_5cell,
        fcom=pathogen_instance.fcom,
        Ap=ap_5cell,
        use_Algae=pathogen_instance.use_Algae,
        use_POC=pathogen_instance.use_POC,
    )
    decay = v1.PathogenDecay(
        apx=1.0,
        q_solar=q_solar_5cell,
        L=kext,
        depth=depth_5cell,
        PX=px_5cell,
    )
    settling = v1.PathogenSettling(vx=1.0, depth=depth_5cell, PX=px_5cell)
    v1_rate = v1.dPXdt(
        PathogenDeath=death,
        PathogenDecay=decay,
        PathogenSettling=settling,
    )

    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6)


def test_phase9f_q_solar_units_are_w_per_m2():
    """Documentation-anchored regression: ``q_solar`` is W/m^2 in v3.

    Phase 9.F resolved a v1 docstring defect that mislabeled ``q_solar``
    as having units of ``1/d``. The v3 default value (500), the v3
    consumption pattern (Beer-Lambert PAR scaling via
    ``utils.light.PAR``), and the v3 documentation (inline comment in
    ``parameters/global_vars.py``, ``utils/light.py:PAR`` Args block,
    and ``processes/pathogen.py:_rate_light_decay`` units note) are all
    consistent on ``W/m^2``. This test pins that consistency by:

    1. Asserting the canonical W/m^2 string appears in each of the
       three v3 documentation sites that reference ``q_solar`` units.
    2. Asserting the legacy ``FIXME(phase1-audit)`` tag has been
       cleared from the ``q_solar`` line in ``global_vars.py`` (Phase
       9.F cleanup).
    3. Numerically exercising ``utils.light.PAR`` to confirm the
       expected ``q_solar * Fr_PAR`` scaling (the implementation that
       cements ``q_solar``'s W/m^2 semantics in the kinetics path).

    See ``parameter_defaults_corrections.md`` Section 2.7.
    """
    from pathlib import Path

    from clearwater_modules_v3.parameters.global_vars import DEFAULTS
    from clearwater_modules_v3.utils.light import PAR

    repo_root = Path(__file__).resolve().parent.parent
    src_root = repo_root / "src" / "clearwater_modules_v3"
    global_vars_text = (src_root / "parameters" / "global_vars.py").read_text()
    light_text = (src_root / "utils" / "light.py").read_text()
    pathogen_text = (src_root / "processes" / "pathogen.py").read_text()

    # 1. Documentation sites mention W/m^2.
    assert "W/m^2" in global_vars_text, (
        "global_vars.py inline comment should state q_solar units as W/m^2"
    )
    assert "W/m^2" in light_text, (
        "utils/light.py PAR Args block should state q_solar as W/m^2"
    )
    assert "W/m^2" in pathogen_text, (
        "processes/pathogen.py _rate_light_decay should mention q_solar W/m^2 units"
    )

    # 2. Phase 9.F: the FIXME on the q_solar line is cleared. The
    #    file may still have FIXME tags on other parameters (e.g.,
    #    vb), so this assertion narrows to the q_solar entry only.
    q_solar_line = next(
        line for line in global_vars_text.splitlines()
        if line.lstrip().startswith("'q_solar'")
    )
    assert "FIXME" not in q_solar_line, (
        f"q_solar line should no longer carry a FIXME tag (Phase 9.F): {q_solar_line!r}"
    )
    assert "W/m^2" in q_solar_line, (
        f"q_solar inline comment should state W/m^2 explicitly: {q_solar_line!r}"
    )

    # 3. Numerical: PAR(q_solar, Fr_PAR) = q_solar * Fr_PAR. The
    #    default q_solar (W/m^2) and Fr_PAR (dimensionless) yield a
    #    PAR irradiance of q_solar * Fr_PAR W/m^2.
    q_solar_default = DEFAULTS["q_solar"]
    fr_par_default = DEFAULTS["Fr_PAR"]
    par = PAR(
        q_solar=xr.DataArray(np.array([q_solar_default])),
        Fr_PAR=xr.DataArray(np.array([fr_par_default])),
    )
    expected = q_solar_default * fr_par_default
    np.testing.assert_allclose(np.asarray(par), np.array([expected]), rtol=1e-12)


def test_phase9f_lambdas_fixme_removed():
    """Documentation-anchored regression: ``lambdas`` FIXME tag cleared.

    Phase 9.C verified via three-way audit that v1
    ``shared/processes.py:232`` applies ``lambdas * Solid``
    unconditionally in the Beer-Lambert sum (matching Fortran
    ``modGlobalParam.f90:LightExtCoefficient`` and v3
    ``utils/light.py``). The earlier Phase 0 framing that the term was
    "commented out / defined but not used" was a documentation defect.

    Phase 9.F removed the now-obsolete ``FIXME(phase1-audit):`` inline
    comment from the ``lambdas`` line in
    ``parameters/global_vars.py``. This test pins that cleanup and
    asserts the corrections doc Section 2.8 is marked RESOLVED.

    See ``parameter_defaults_corrections.md`` Section 2.8.
    """
    from pathlib import Path

    repo_root = Path(__file__).resolve().parent.parent
    src_root = repo_root / "src" / "clearwater_modules_v3"
    global_vars_text = (src_root / "parameters" / "global_vars.py").read_text()
    corrections_text = (src_root / "parameter_defaults_corrections.md").read_text()

    # The lambdas line in global_vars.py no longer carries a FIXME tag.
    lambdas_line = next(
        line for line in global_vars_text.splitlines()
        if line.lstrip().startswith("'lambdas'")
    )
    assert "FIXME" not in lambdas_line, (
        f"lambdas line should no longer carry a FIXME tag (Phase 9.F): {lambdas_line!r}"
    )
    # And the inline comment now affirmatively records the Phase 9.C
    # verification ("active" / "applied unconditionally").
    assert "applied unconditionally" in lambdas_line or "active" in lambdas_line, (
        f"lambdas inline comment should record the Phase 9.C verification: {lambdas_line!r}"
    )

    # The corrections doc no longer carries the false "commented out"
    # framing as a current claim. The phrase may still appear in the
    # historical-context narrative inside Section 2.8 (the audit-history
    # record), but only as something that was *corrected*, not as a
    # current statement of fact. Pin the resolution markers instead.
    assert "Section 2.8" in corrections_text or "lambdas" in corrections_text
    section_2_8_start = corrections_text.find("### 2.8 `lambdas`")
    assert section_2_8_start != -1, "Section 2.8 should exist in corrections doc"
    section_2_8_end = corrections_text.find("\n### ", section_2_8_start + 1)
    section_2_8 = corrections_text[section_2_8_start:section_2_8_end]
    assert "RESOLVED" in section_2_8, (
        "Section 2.8 should be marked RESOLVED"
    )
    assert "Phase 9.F" in section_2_8 or "Phase 9.C" in section_2_8, (
        "Section 2.8 should reference the Phase 9.C verification "
        "and/or Phase 9.F FIXME cleanup"
    )


# ---------------------------------------------------------------------------
# Phase 9.F.B regression tests — apx, vx canonical values + Fr_PAR revert
# ---------------------------------------------------------------------------
# Phase 9.F.B replaced the v1 placeholder defaults apx=1.0 and vx=1.0
# with the canonical Auer & Niehaus (1993) / Chapra (1997) literature
# values, and reverted the Phase 3.1 substitution
# I0 = q_solar * Fr_PAR in _rate_light_decay so that the v3 kinetics
# tie directly to the canonical broadband-solar calibration.
#
#   Source                                      apx    vx (m/d)
#   ----                                        ---    --------
#   Auer & Niehaus 1993 (Onondaga Lake)         0.017  1.38
#     (alpha = 0.00824 cm^2/cal in cgs ->
#      0.017 (W/m^2)^-1 d^-1 in SI)
#   Chapra 1997, Surface Water-Quality          ~0.017 ~1
#     Modeling, McGraw-Hill, Ch. 33 (cites
#     Auer & Niehaus 1993)
#   QUAL2K v2.11b8 (Chapra et al. 2008)         "user" "user"
#     §5.5.20.1 (formulation only)
#   Bowie et al. 1985 (compilation range)       --     0.5-2.5
#   v1 / Fortran / pre-9.F.B v3 (placeholder)   1.0    1.0
#   Phase 9.F.B v3                              0.017  1.38

def test_phase9fb_apx_canonical_value_pinned():
    """Phase 9.F.B: apx = 0.017 (W/m^2)^-1 d^-1 (Auer & Niehaus 1993
    canonical via Chapra 1997 / QUAL2K). Pin so any future change
    requires explicit reconciliation against the literature."""
    from clearwater_modules_v3.parameters.pathogen import (
        DEFAULTS as PATHOGEN_DEFAULTS,
    )

    apx = PATHOGEN_DEFAULTS["apx"]
    np.testing.assert_allclose(
        apx, 0.017, rtol=1e-6,
        err_msg=(
            "Phase 9.F.B canonical: apx = 0.017 (W/m^2)^-1 d^-1 "
            "(Auer & Niehaus 1993; cited by Chapra 1997 Ch. 33)"
        ),
    )
    # And confirm the value is no longer the v1/Fortran placeholder.
    assert apx != 1.0


def test_phase9fb_vx_canonical_value_pinned():
    """Phase 9.F.B: vx = 1.38 m/d (Auer & Niehaus 1993 sediment-trap
    canonical via Chapra 1997). Pin so any future change requires
    explicit reconciliation against the literature."""
    from clearwater_modules_v3.parameters.pathogen import (
        DEFAULTS as PATHOGEN_DEFAULTS,
    )

    vx = PATHOGEN_DEFAULTS["vx"]
    np.testing.assert_allclose(
        vx, 1.38, rtol=1e-6,
        err_msg=(
            "Phase 9.F.B canonical: vx = 1.38 m/d "
            "(Auer & Niehaus 1993; cited by Chapra 1997 Ch. 33)"
        ),
    )
    # And confirm the value is no longer the v1/Fortran placeholder.
    assert vx != 1.0
    # Bowie et al. 1985 typical range 0.5 - 2.5 m/d.
    assert 0.5 <= vx <= 2.5


def test_phase9fb_light_decay_uses_raw_q_solar():
    """Phase 9.F.B: ``_rate_light_decay`` uses raw broadband ``q_solar``
    rather than ``q_solar * Fr_PAR``.

    Constructs two pathogen instances with the same ``apx`` and
    ``q_solar`` but different ``Fr_PAR`` values. Pre-9.F.B the
    light-decay rate scaled with Fr_PAR; Phase 9.F.B reverted that
    substitution, so the two instances should now produce identical
    rates.
    """
    pa = Pathogen(
        parameters={
            "apx": 0.017, "vx": 1.38, "kdx_20": 0.8, "kdx_theta": 1.07,
            "Fr_PAR": 0.47,  # default
        },
        time_step=timedelta(minutes=5),
    )
    pb = Pathogen(
        parameters={
            "apx": 0.017, "vx": 1.38, "kdx_20": 0.8, "kdx_theta": 1.07,
            "Fr_PAR": 1.0,  # pre-9.F.B test trick
        },
        time_step=timedelta(minutes=5),
    )

    px = xr.DataArray(np.array([1.0e3, 5.0e3, 1.0e4, 5.0e4, 1.0e5]))
    depth = xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))
    q_solar = xr.DataArray(np.array([200.0, 250.0, 300.0, 350.0, 400.0]))
    solid = xr.DataArray(np.array([10.0, 12.0, 15.0, 18.0, 20.0]))
    poc = xr.DataArray(np.array([1.0, 1.2, 1.4, 1.6, 1.8]))
    ap = xr.DataArray(np.array([5.0, 6.0, 7.0, 8.0, 10.0]))

    rate_a = pa._rate_light_decay(
        px=px, depth=depth, q_solar=q_solar,
        solid=solid, poc=poc, ap=ap,
    )
    rate_b = pb._rate_light_decay(
        px=px, depth=depth, q_solar=q_solar,
        solid=solid, poc=poc, ap=ap,
    )
    np.testing.assert_allclose(
        np.asarray(rate_a), np.asarray(rate_b), rtol=1e-12,
        err_msg=(
            "Phase 9.F.B: _rate_light_decay should be insensitive to "
            "Fr_PAR after the PAR substitution was reverted"
        ),
    )


def test_phase9fb_light_decay_matches_v1_with_canonical_apx():
    """Phase 9.F.B: with apx=0.017 (canonical) and raw q_solar, the v3
    light-decay rate should match v1 ``PathogenDecay`` exactly,
    confirming the Phase 9.F.B revert restored kinetics-level parity
    with the canonical broadband formulation."""
    p = Pathogen(
        parameters={
            "apx": 0.017, "vx": 1.38, "kdx_20": 0.8, "kdx_theta": 1.07,
        },
        time_step=timedelta(minutes=5),
    )

    px = xr.DataArray(np.array([1.0e3, 5.0e3, 1.0e4, 5.0e4, 1.0e5]))
    depth = xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))
    q_solar = xr.DataArray(np.array([200.0, 250.0, 300.0, 350.0, 400.0]))
    solid = xr.DataArray(np.array([10.0, 12.0, 15.0, 18.0, 20.0]))
    poc = xr.DataArray(np.array([1.0, 1.2, 1.4, 1.6, 1.8]))
    ap = xr.DataArray(np.array([5.0, 6.0, 7.0, 8.0, 10.0]))

    v3_rate = p._rate_light_decay(
        px=px, depth=depth, q_solar=q_solar,
        solid=solid, poc=poc, ap=ap,
    )

    kext = L(
        lambda0=p.lambda0, lambda1=p.lambda1, lambda2=p.lambda2,
        lambdas=p.lambdas, lambdam=p.lambdam,
        Solid=solid, POC=poc, fcom=p.fcom, Ap=ap,
        use_Algae=p.use_Algae, use_POC=p.use_POC,
    )
    v1_rate = v1.PathogenDecay(
        apx=0.017, q_solar=q_solar, L=kext, depth=depth, PX=px,
    )
    np.testing.assert_allclose(
        np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6,
    )
