"""v3 Pathogen kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_pathogen_calculations_v2.py``.
"""
from datetime import timedelta
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.pathogen import Pathogen
from clearwater_modules_v3.utils.light import L


V1_NATURAL_DECAY_REFERENCE = np.array([
    570.3889435869346,
    3493.7549130928464,
    8000.0,
    45796.0,
    112204.13845600003,
])

V1_LIGHT_DECAY_REFERENCE = np.array([
    151181.25993312732,
    673435.4824252209,
    1100709.6986337688,
    4420511.790555701,
    6102804.008810551,
])

V1_SETTLING_REFERENCE = np.array([
    2000.0,
    5000.0,
    6666.666666666666,
    25000.0,
    33333.33333333333,
])

V1_DPXDT_REFERENCE = np.array([
    -153751.64887671426,
    -681929.2373383137,
    -1115376.3653004356,
    -4491307.790555701,
    -6248341.480599884,
])

V1_LIGHT_DECAY_CANONICAL_REFERENCE = np.array([
    2570.081418863165,
    11448.403201228757,
    18712.064876774068,
    75148.70043944691,
    103747.66814977938,
])


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
    return Pathogen(
        parameters={
            "kdx_20": 0.8, "kdx_theta": 1.07,
            "apx": 1.0, "vx": 1.0, "Fr_PAR": 1.0,
        },
        time_step=timedelta(minutes=5),
    )


def test_pathogen_natural_decay_matches_v1_PathogenDeath(
    pathogen_instance, px_5cell, water_temp_5cell
):
    v3_rate = pathogen_instance._rate_natural_decay(px_5cell, water_temp_5cell)
    np.testing.assert_allclose(
        np.asarray(v3_rate), V1_NATURAL_DECAY_REFERENCE, rtol=1e-6
    )


def test_pathogen_light_decay_matches_v1_PathogenDecay(
    pathogen_instance, px_5cell, depth_5cell,
    q_solar_5cell, solid_5cell, poc_5cell, ap_5cell,
):
    v3_rate = pathogen_instance._rate_light_decay(
        px=px_5cell, depth=depth_5cell, q_solar=q_solar_5cell,
        solid=solid_5cell, poc=poc_5cell, ap=ap_5cell,
    )
    np.testing.assert_allclose(
        np.asarray(v3_rate), V1_LIGHT_DECAY_REFERENCE, rtol=1e-6
    )


def test_pathogen_settling_matches_v1_PathogenSettling(
    pathogen_instance, px_5cell, depth_5cell
):
    v3_rate = pathogen_instance._rate_settling(px_5cell, depth_5cell)
    np.testing.assert_allclose(
        np.asarray(v3_rate), V1_SETTLING_REFERENCE, rtol=1e-6
    )


def test_pathogen_total_rate_matches_v1_dPXdt(
    pathogen_instance, px_5cell, water_temp_5cell, depth_5cell,
    q_solar_5cell, solid_5cell, poc_5cell, ap_5cell,
):
    v3_rate = pathogen_instance.rate(
        px=px_5cell, water_temperature=water_temp_5cell, depth=depth_5cell,
        q_solar=q_solar_5cell, solid=solid_5cell, poc=poc_5cell, ap=ap_5cell,
    )
    np.testing.assert_allclose(
        np.asarray(v3_rate), V1_DPXDT_REFERENCE, rtol=1e-6
    )


def test_phase9f_q_solar_units_are_w_per_m2():
    """Documentation-anchored regression: ``q_solar`` is W/m^2 in v3."""
    from clearwater_modules_v3.parameters.global_vars import DEFAULTS
    from clearwater_modules_v3.utils.light import PAR

    repo_root = Path(__file__).resolve().parents[3]
    src_root = repo_root / "src" / "clearwater_modules_v3"
    global_vars_text = (src_root / "parameters" / "global_vars.py").read_text()
    light_text = (src_root / "utils" / "light.py").read_text()
    pathogen_text = (src_root / "processes" / "pathogen.py").read_text()

    assert "W/m^2" in global_vars_text
    assert "W/m^2" in light_text
    assert "W/m^2" in pathogen_text

    q_solar_line = next(
        line for line in global_vars_text.splitlines()
        if line.lstrip().startswith("'q_solar'")
    )
    assert "FIXME" not in q_solar_line
    assert "W/m^2" in q_solar_line

    q_solar_default = DEFAULTS["q_solar"]
    fr_par_default = DEFAULTS["Fr_PAR"]
    par = PAR(
        q_solar=xr.DataArray(np.array([q_solar_default])),
        Fr_PAR=xr.DataArray(np.array([fr_par_default])),
    )
    expected = q_solar_default * fr_par_default
    np.testing.assert_allclose(np.asarray(par), np.array([expected]), rtol=1e-12)


def test_phase9f_lambdas_fixme_removed():
    """Phase 9.F: ``lambdas`` FIXME tag cleared."""
    repo_root = Path(__file__).resolve().parents[3]
    src_root = repo_root / "src" / "clearwater_modules_v3"
    global_vars_text = (src_root / "parameters" / "global_vars.py").read_text()
    corrections_text = (src_root / "parameter_defaults_corrections.md").read_text()

    lambdas_line = next(
        line for line in global_vars_text.splitlines()
        if line.lstrip().startswith("'lambdas'")
    )
    assert "FIXME" not in lambdas_line
    assert "applied unconditionally" in lambdas_line or "active" in lambdas_line

    section_2_8_start = corrections_text.find("### 2.8 `lambdas`")
    assert section_2_8_start != -1
    section_2_8_end = corrections_text.find("\n### ", section_2_8_start + 1)
    section_2_8 = corrections_text[section_2_8_start:section_2_8_end]
    assert "RESOLVED" in section_2_8
    assert "Phase 9.F" in section_2_8 or "Phase 9.C" in section_2_8


def test_phase9fb_apx_canonical_value_pinned():
    """Phase 9.F.B: apx = 0.017 (Auer & Niehaus 1993 canonical)."""
    from clearwater_modules_v3.parameters.pathogen import (
        DEFAULTS as PATHOGEN_DEFAULTS,
    )
    apx = PATHOGEN_DEFAULTS["apx"]
    np.testing.assert_allclose(apx, 0.017, rtol=1e-6)
    assert apx != 1.0


def test_phase9fb_vx_canonical_value_pinned():
    """Phase 9.F.B: vx = 1.38 m/d (Auer & Niehaus 1993 canonical)."""
    from clearwater_modules_v3.parameters.pathogen import (
        DEFAULTS as PATHOGEN_DEFAULTS,
    )
    vx = PATHOGEN_DEFAULTS["vx"]
    np.testing.assert_allclose(vx, 1.38, rtol=1e-6)
    assert vx != 1.0
    assert 0.5 <= vx <= 2.5


def test_phase9fb_light_decay_uses_raw_q_solar():
    """Phase 9.F.B: light decay uses raw broadband q_solar, not q_solar*Fr_PAR."""
    pa = Pathogen(
        parameters={
            "apx": 0.017, "vx": 1.38, "kdx_20": 0.8, "kdx_theta": 1.07,
            "Fr_PAR": 0.47,
        },
        time_step=timedelta(minutes=5),
    )
    pb = Pathogen(
        parameters={
            "apx": 0.017, "vx": 1.38, "kdx_20": 0.8, "kdx_theta": 1.07,
            "Fr_PAR": 1.0,
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
        px=px, depth=depth, q_solar=q_solar, solid=solid, poc=poc, ap=ap
    )
    rate_b = pb._rate_light_decay(
        px=px, depth=depth, q_solar=q_solar, solid=solid, poc=poc, ap=ap
    )
    np.testing.assert_allclose(
        np.asarray(rate_a), np.asarray(rate_b), rtol=1e-12
    )


def test_phase9fb_light_decay_matches_v1_with_canonical_apx():
    """Phase 9.F.B: v3 light-decay rate with canonical apx matches frozen
    v1 ``PathogenDecay`` reference computed at apx=0.017."""
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
        px=px, depth=depth, q_solar=q_solar, solid=solid, poc=poc, ap=ap
    )
    np.testing.assert_allclose(
        np.asarray(v3_rate), V1_LIGHT_DECAY_CANONICAL_REFERENCE, rtol=1e-6
    )
