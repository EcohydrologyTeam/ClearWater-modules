"""v3 NSM1 CBOD parameter defaults.

Consumed by: ``CBOD`` Process (multi-group CBOD; defaults apply per group).
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``CBODStaticVariables`` /
``DEFAULT_CBOD``.
Corrections applied: none.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'KsOxbod': 0.5,         # mg-O2/L; oxygen half-saturation for CBOD oxidation
    'kbod_20': 0.12,        # 1/d; CBOD oxidation rate at 20 C
    'ksbod_20': 0.0,        # FIXME(phase1-audit): CBOD never settles, confirm with LimnoTech. m/d at 20 C.
    'kbod_theta': 1.047,    # unitless; Arrhenius coefficient for CBOD decay
    'ksbod_theta': 1.047,   # unitless; Arrhenius coefficient for CBOD settling
}
