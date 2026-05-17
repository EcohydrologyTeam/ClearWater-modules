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

  Two related defects — **RESOLVED 2026-05-16 (NSM1-SCI-CB1,
  gold-standard spec C2; research doc
  ``docs/clearwater_modules_v3_nsm1_research_2_3_ksbod.md``):**

  - **Units form**: v3's ``processes/cbod.py`` divided ``ksbod_tc`` by
    depth (``ksbod_tc / depth * cbod``), implementing it as a
    settling-velocity (m/d), but Fortran NSM1 (``modCBOD.f90:114``,
    no depth division) and QUAL2E both treat it as a 1/d rate
    constant. **Fixed** to ``ksbod_tc * cbod`` (no depth division;
    1/d at 20 °C). With ``ksbod_20 = 0`` the form difference was
    silent; nonzero user values would have diverged by 1/depth.
  - **Theta**: v3 ``ksbod_theta`` was ``1.047`` (the oxidation
    coefficient); the canonical **settling** value is ``1.024``
    (Bowie 1985 / QUAL2E). **Fixed** to ``1.024``.

  See ``parameter_defaults_corrections.md`` Section 2.3 and research
  record in ``docs/clearwater_modules_v3_nsm1_research_2_3_ksbod.md``.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'KsOxbod': 0.5,         # mg-O2/L; oxygen half-saturation for CBOD oxidation
    'kbod_20': 0.12,        # 1/d; CBOD oxidation rate at 20 C
    'ksbod_20': 0.0,        # 1/d at 20 C; CBOD settling first-order RATE constant (NOT a velocity). v3/v1/Fortran/QUAL2E default; intentional zero per modern dissolved-CBOD convention (QUAL2K, WASP, CE-QUAL-W2 omit the parameter entirely). Nonzero values are site/source-specific. NSM1-SCI-CB1 (spec C2): processes/cbod.py applies ksbod_tc*cbod with NO depth division (Fortran modCBOD.f90:114 / QUAL2E). See corrections doc Section 2.3.
    'kbod_theta': 1.047,    # unitless; Arrhenius coefficient for CBOD decay
    'ksbod_theta': 1.024,   # unitless; Arrhenius coefficient for CBOD SETTLING. NSM1-SCI-CB1 (spec C2): 1.047 -> 1.024, the canonical settling-coefficient value (Bowie 1985 / QUAL2E); the prior 1.047 was the oxidation coefficient mis-applied to settling. See corrections doc Section 2.3.
}
