"""v3 NSM1 pathogen parameter defaults.

Consumed by: ``Pathogen`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``PathogenStaticVariables`` /
``DEFAULT_PATHOGEN``.
Corrections applied: none.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'kdx_20': 0.8,          # 1/d; pathogen decay (inactivation) rate at 20 C
    'kdx_theta': 1.07,      # unitless; Arrhenius coefficient for pathogen decay
    'apx': 1.0,             # FIXME(phase1-audit): unitless; pathogen sunlight inactivation/algal coupling, placeholder without literature basis
    'vx': 1.0,              # FIXME(phase1-audit): m/d; pathogen settling velocity, placeholder without literature basis
}
