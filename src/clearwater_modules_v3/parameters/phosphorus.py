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
    'rpo4_20': 0.0,         # g-P/m^2/d at 20 C; sediment P release rate. Held at zero by design: v3 Phosphorus.run() gates dip_from_bed by use_TIP only, NOT use_SedFlux, so the zero default is the de facto gate. Phase 9.F.C added a defensive NotImplementedError guard in Phosphorus.__init__ that fires if use_SedFlux=True is passed (sediment-flux feature requires NSM2 path). See corrections doc Section 2.1.
    'kop_theta': 1.047,     # unitless; Arrhenius coefficient for OrgP decay
    'rpo4_theta': 1.074,    # unitless; Arrhenius coefficient for sediment P release
    'kdpo4': 0.0,           # FIXME(phase1-audit): TIP partitioning feature disabled (NSM2 territory). L/kg.
    'vsop': 0.1,            # m/d; OrgP settling velocity (was 999 in v1, see corrections doc Section 1)
    'vs': 0.1,              # m/d; TIP settling velocity (was 999 in v1, see corrections doc Section 1)
}
