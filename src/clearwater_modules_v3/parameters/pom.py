"""v3 NSM1 particulate organic matter (POM) parameter defaults.

Consumed by: ``POM`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``POMStaticVariables`` /
``DEFAULT_POM``.
Corrections applied: none.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'kpom_20': 0.1,         # 1/d; POM dissolution rate at 20 C
    'h2': 0.1,              # FIXME(phase1-audit): m; sediment burial/sedimentation depth denominator, unclear physical role
    'kpom_theta': 1.047,    # unitless; Arrhenius coefficient for POM dissolution
}
