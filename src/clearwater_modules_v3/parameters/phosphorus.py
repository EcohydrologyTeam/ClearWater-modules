"""v3 NSM1 phosphorus parameter defaults.

Consumed by: ``Phosphorus`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``PhosphorusStaticVariables`` /
``DEFAULT_PHOSPHORUS``, plus ``vsop`` and ``vs`` migrated from v1 ``GlobalVars``
since they are phosphorus-specific settling velocities.
Corrections applied: ``vsop=0.1`` (was 999), ``vs=0.1`` (was 999).
See ``parameter_defaults_corrections.md`` Section 1.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'kop_20': 0.1,          # 1/d; OrgP decomposition rate at 20 C
    'rpo4_20': 0.0,         # FIXME(phase1-audit): sediment P release silently disabled; verify gated by use_SedFlux. 1/d at 20 C.
    'kop_theta': 1.047,     # unitless; Arrhenius coefficient for OrgP decay
    'rpo4_theta': 1.074,    # unitless; Arrhenius coefficient for sediment P release
    'kdpo4': 0.0,           # FIXME(phase1-audit): TIP partitioning feature disabled (NSM2 territory). L/kg.
    'vsop': 0.1,            # m/d; OrgP settling velocity (was 999 in v1, see corrections doc Section 1)
    'vs': 0.1,              # m/d; TIP settling velocity (was 999 in v1, see corrections doc Section 1)
}
