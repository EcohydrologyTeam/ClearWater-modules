"""v3 NSM1 particulate organic matter (POM) parameter defaults.

Consumed by: ``POM`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``POMStaticVariables`` /
``DEFAULT_POM``.

Corrections applied:

* ``h2`` FIXME cleared (Phase 9.F.C). Phase 0 flagged ``h2`` with
  "unclear physical role" but the v1 ``static_variables.py:921``
  declaration and Fortran ``modGlobalParam.f90:38`` both define
  ``h2`` unambiguously as the **active sediment layer thickness (m)**.
  ``h2 = 0.1`` matches the Di Toro (2001) / QUAL2K convention for the
  lower anaerobic sediment layer thickness ``H_2`` (10 cm; QUAL2K
  v2.11 §5.6 Eq. 214). NSM1's POM state variable represents bed-
  sediment POM (Fortran ``POM2`` — the "2" suffix denoting Di Toro's
  layer 2), and ``h2`` is the divisor that converts areal water-column
  fluxes (m * mg/L/d) into volumetric concentration changes (mg/L/d)
  in that bed layer. The Phase 0 audit comment reflected a
  documentation gap, not a substantive issue. See
  ``parameter_defaults_corrections.md`` Section 2.5.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'kpom_20': 0.1,         # 1/d; POM dissolution rate at 20 C
    'h2': 0.1,              # m; active sediment layer thickness (Di Toro 2001 / QUAL2K H_2 anaerobic-layer convention, ~10 cm). Divides areal water-column fluxes (m * mg/L/d) into bed volumetric concentrations (mg/L/d). Phase 9.F.C documentation fix; see corrections doc Section 2.5.
    'kpom_theta': 1.047,    # unitless; Arrhenius coefficient for POM dissolution
}
