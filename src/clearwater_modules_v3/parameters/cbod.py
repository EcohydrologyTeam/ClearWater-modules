"""v3 NSM1 CBOD parameter defaults.

Consumed by: ``CBOD`` Process (multi-group CBOD; defaults apply per group).
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``CBODStaticVariables`` /
``DEFAULT_CBOD``.

Corrections applied:

* ``ksbod_20`` FIXME cleared (Phase 9.F.C). Phase 9.F research
  confirmed that the zero default is the intentional, defensible
  modern-convention value:

  - QUAL2K v2.11b8 (Chapra, Pelletier & Tao 2008), QUAL2Kw, WASP7
    EUTRO, and CE-QUAL-W2 all treat CBOD as a **dissolved-only**
    state variable and provide **no CBOD settling parameter at all**;
    particulate organic matter is carried separately (detritus,
    LPOM/RPOM) with its own settling velocity.
  - QUAL2E (Brown & Barnwell 1987, EPA/600/3-87/007), the
    direct ancestor of NSM1, defaults ``K_3 = 0`` for the same reason.
  - EPA TMDL Technical Guidance Book II (Sample Calc B-3) explicitly
    assumes ``K_s = 0`` for treated effluent.
  - Yamuna River QUAL2E case (Parmar & Keshari, citing Kazmi &
    Agrawal 2005) calibrated ``K_3 = 0.9`` /d uniformly across 16
    reaches in a heavily polluted urban stretch where particulate-
    laden CBOD settling dominated removal — illustrating that
    nonzero values are site-, source-, and treatment-specific.

  Two related defects are flagged for follow-up (Section 4 / future
  audit, NOT addressed in Phase 9.F.C):

  - **Units form**: v3's ``processes/cbod.py`` divides ``ksbod_tc`` by
    depth (``ksbod_tc / depth * cbod``) implementing it as a
    settling-velocity (m/d), but Fortran NSM1 (``modCBOD.f90:114``,
    no depth division) and QUAL2E both treat it as a 1/d rate
    constant. With ``ksbod_20 = 0`` the form difference is silent;
    nonzero user values would diverge by a factor of 1/depth.
  - **Theta**: v3 ``ksbod_theta = 1.047`` differs from
    Fortran/QUAL2E ``1.024`` (the canonical settling-coefficient
    Arrhenius value per Bowie 1985 / QUAL2E).

  See ``parameter_defaults_corrections.md`` Section 2.3 and research
  record in ``docs/clearwater_modules_v3_nsm1_research_2_3_ksbod.md``.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'KsOxbod': 0.5,         # mg-O2/L; oxygen half-saturation for CBOD oxidation
    'kbod_20': 0.12,        # 1/d; CBOD oxidation rate at 20 C
    'ksbod_20': 0.0,        # CBOD settling rate at 20 C. v3/v1/Fortran/QUAL2E default; intentional zero per modern dissolved-CBOD convention (QUAL2K, WASP, CE-QUAL-W2 omit the parameter entirely). Nonzero values are site/source-specific. Phase 9.F.C documentation fix; see corrections doc Section 2.3.
    'kbod_theta': 1.047,    # unitless; Arrhenius coefficient for CBOD decay
    'ksbod_theta': 1.047,   # unitless; Arrhenius coefficient for CBOD settling (note: Fortran/QUAL2E use 1.024 for settling; v3 inherits v1's 1.047, flagged for follow-up - see corrections doc Section 2.3)
}
