"""Phase 9.C parameter-correction regression tests.

Pins the canonical default values for parameters changed by the Phase 9.C
three-way audit (v1 / Fortran NSM1 / v3) documented in
``design/clearwater_modules_v3_nsm1_audit_utilities_params.md`` and
``src/clearwater_modules_v3/parameter_defaults_corrections.md`` Sections 1.8
and 1.9.

These tests guard against accidental regression of the corrections back to
v1's flawed values:

* ``vson_20``: corrected from 0.1 to 0.01 m/d (matches Fortran
  ``modGlobalParam.f90:92`` and v1 ``GlobalVars.vson``). Prior 0.1 was an
  internal v3 inconsistency: 10x the Fortran/v1 value and 10x v3's own
  ``global_vars.vson``.
* ``lambdam``: corrected from 0.0174 to 0.174 L/(mg*m) (matches Fortran
  ``modGlobalParam.f90:68`` and QUAL2K Table 6). Prior 0.0174 was inherited
  from a v1 typo; the legacy v1 NSM test suite already overrides with
  ``lambdam=0.174`` (e.g., ``test_7_nsm_algae_calculations.py:340``).

The tests also verify that the inline ``_LIGHT_DEFAULTS`` fallback in
``processes/pathogen.py`` agrees with the canonical ``global_vars`` value.
"""
from clearwater_modules_v3.parameters.global_vars import DEFAULTS as GLOBAL_VARS_DEFAULTS
from clearwater_modules_v3.parameters.nitrogen import DEFAULTS as NITROGEN_DEFAULTS


def test_vson_20_matches_fortran_and_v1():
    """``vson_20`` must equal 0.01 m/d (Fortran modGlobalParam.f90:92, v1
    GlobalVars.vson). Phase 9.C audit fix; see corrections doc Section 1.8.
    """
    assert NITROGEN_DEFAULTS["vson_20"] == 0.01


def test_vson_20_matches_global_vars_vson():
    """The migrated ``vson_20`` (nitrogen group) must agree with the legacy
    ``global_vars.vson`` to maintain internal v3 consistency. Phase 9.C
    audit found a 10x mismatch (0.1 vs 0.01) which has now been corrected.
    """
    assert NITROGEN_DEFAULTS["vson_20"] == GLOBAL_VARS_DEFAULTS["vson"]


def test_lambdam_matches_fortran_and_qual2k():
    """``lambdam`` must equal 0.174 L/(mg*m) (Fortran modGlobalParam.f90:68
    and QUAL2K Table 6). Phase 9.C audit fix; see corrections doc Section
    1.9. Prior 0.0174 was a v1 typo inherited by v3.
    """
    assert GLOBAL_VARS_DEFAULTS["lambdam"] == 0.174


def test_pathogen_inline_lambdam_matches_global_vars():
    """The inline ``_LIGHT_DEFAULTS`` fallback in ``processes/pathogen.py``
    must agree with the canonical ``global_vars.lambdam`` to prevent the
    Pathogen Process from using a stale value when constructed without
    user overrides. Phase 9.C audit confirmed the inline value was 0.0174
    while the canonical was being corrected to 0.174.
    """
    from clearwater_modules_v3.processes.pathogen import _LIGHT_DEFAULTS

    assert _LIGHT_DEFAULTS["lambdam"] == GLOBAL_VARS_DEFAULTS["lambdam"]


def test_lambdas_active_in_light_extinction():
    """The ``lambdas`` parameter is applied unconditionally in v3's
    light-extinction utility (matching v1 ``shared/processes.py:232`` and
    Fortran ``modGlobalParam.f90`` ``LightExtCoefficient``). Phase 9.C
    audit corrected an earlier corrections-doc claim that ``lambdas * Solid``
    was "commented out". Verify the parameter is non-zero (default 0.052).
    """
    assert GLOBAL_VARS_DEFAULTS["lambdas"] == 0.052


def test_lambdam_smoke_via_L_utility():
    """End-to-end smoke test: the ``lambdam`` value is consumed by the
    ``utils.light.L`` utility and the corrected 0.174 value produces a
    measurably-different extinction coefficient than the prior 0.0174.
    """
    import numpy as np
    import xarray as xr

    from clearwater_modules_v3.utils.light import L

    # 5-cell smoke registry. POC=2 mg/L makes the lambdam term meaningful.
    cells = 5
    lambda0 = xr.DataArray(np.full(cells, GLOBAL_VARS_DEFAULTS["lambda0"]))
    lambda1 = xr.DataArray(np.full(cells, GLOBAL_VARS_DEFAULTS["lambda1"]))
    lambda2 = xr.DataArray(np.full(cells, GLOBAL_VARS_DEFAULTS["lambda2"]))
    lambdas = xr.DataArray(np.full(cells, GLOBAL_VARS_DEFAULTS["lambdas"]))
    lambdam = xr.DataArray(np.full(cells, GLOBAL_VARS_DEFAULTS["lambdam"]))
    Solid = xr.DataArray(np.full(cells, 1.0))
    POC = xr.DataArray(np.full(cells, 2.0))
    fcom = xr.DataArray(np.full(cells, GLOBAL_VARS_DEFAULTS["fcom"]))
    Ap = xr.DataArray(np.full(cells, 0.0))
    use_Algae = xr.DataArray(np.full(cells, False))
    use_POC = xr.DataArray(np.full(cells, True))

    extinction = L(
        lambda0=lambda0,
        lambda1=lambda1,
        lambda2=lambda2,
        lambdas=lambdas,
        lambdam=lambdam,
        Solid=Solid,
        POC=POC,
        fcom=fcom,
        Ap=Ap,
        use_Algae=use_Algae,
        use_POC=use_POC,
    )

    # Expected: lambda0 + lambdas*Solid + lambdam * POC / fcom
    #         = 0.02 + 0.052*1.0 + 0.174 * 2.0 / 0.4
    #         = 0.02 + 0.052 + 0.87 = 0.942
    expected = (
        GLOBAL_VARS_DEFAULTS["lambda0"]
        + GLOBAL_VARS_DEFAULTS["lambdas"] * 1.0
        + GLOBAL_VARS_DEFAULTS["lambdam"] * 2.0 / GLOBAL_VARS_DEFAULTS["fcom"]
    )
    np_extinction = np.asarray(extinction)
    np.testing.assert_allclose(np_extinction, np.full(cells, expected), rtol=1e-12)

    # Sanity: with the corrected lambdam=0.174 the POC term dominates;
    # if the value reverted to 0.0174 the extinction would drop by ~0.78
    # (an order of magnitude).
    assert np_extinction[0] > 0.5
