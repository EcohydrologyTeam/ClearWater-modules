"""v3 NSM1 pathogen parameter defaults.

Consumed by: ``Pathogen`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``PathogenStaticVariables`` /
``DEFAULT_PATHOGEN``.

Corrections applied:

* ``apx`` corrected from placeholder ``1.0`` to ``0.017`` (W/m^2)^-1 d^-1
  (Phase 9.F.B). Canonical value of Auer & Niehaus (1993, *Wat. Res.*
  27(4):693-701) for fecal coliform sunlight inactivation in Onondaga
  Lake, equivalent to alpha = 0.00824 cm^2/cal in cgs units (Chapra
  1997, *Surface Water-Quality Modeling*, McGraw-Hill, Ch. 33; QUAL2K
  v2.11b8 §5.5.20.1). The v1 docstring claim that ``apx`` is
  "dimensionless" was dimensionally incorrect: the rate-balance
  ``[1/d] = apx * q_solar * (dimensionless optical factor)`` requires
  ``apx`` to carry units ``(W/m^2)^-1 d^-1`` because ``q_solar`` is
  W/m^2 in v3. See ``parameter_defaults_corrections.md`` Section 1.15.

* ``vx`` corrected from placeholder ``1.0`` to ``1.38`` m/d (Phase
  9.F.B). Auer & Niehaus (1993) sediment-trap measurement of
  particle-associated fecal coliform in Onondaga Lake; the canonical
  value cited in Chapra (1997, Ch. 33), QUAL2K v2.11b8, and subsequent
  modeling studies. Bowie et al. (1985) compilation reports a 0.5-2.5
  m/d range across studies, bracketing the canonical 1.38 m/d.
  See ``parameter_defaults_corrections.md`` Section 1.16.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'kdx_20': 0.8,          # 1/d; pathogen decay (inactivation) rate at 20 C
    'kdx_theta': 1.07,      # unitless; Arrhenius coefficient for pathogen decay
    'apx': 0.017,           # (W/m^2)^-1 d^-1; pathogen sunlight-inactivation efficiency. Auer & Niehaus 1993 / Chapra 1997 canonical (alpha = 0.00824 cm^2/cal in cgs); was 1.0 placeholder in v1/v3 pre-9.F.B. See corrections doc Section 1.15.
    'vx': 1.38,             # m/d; pathogen net settling velocity. Auer & Niehaus 1993 / Chapra 1997 canonical (sediment-trap measurement, Onondaga Lake); was 1.0 placeholder in v1/v3 pre-9.F.B. See corrections doc Section 1.16.
}
