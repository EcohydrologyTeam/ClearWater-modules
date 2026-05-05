"""v3 NSM1 nitrogen parameter defaults.

Consumed by: ``Nitrogen`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``NitrogenStaticVariables`` /
``DEFAULT_NITROGEN``.
Corrections applied: none.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'KNR': 0.6,             # mg-O2/L; oxygen half-saturation for nitrification inhibition (1 - exp(-KNR*DOX))
    'knit_20': 0.1,         # 1/d; nitrification rate at 20 C
    'kon_20': 0.1,          # 1/d; organic-N hydrolysis (ammonification) rate at 20 C
    'kdnit_20': 0.002,      # 1/d; denitrification rate at 20 C
    'rnh4_20': 0.0,         # FIXME(phase1-audit): sediment NH4 release silently disabled; verify gated by use_SedFlux. 1/d at 20 C.
    'vno3_20': 0.0,         # FIXME(phase1-audit): sediment NO3 denitrification silently disabled; verify gated by use_SedFlux. 1/d at 20 C.
    'vson_20': 0.1,         # m/d; OrgN settling velocity at 20 C (v1 default; design spec Section 3.4)
    'knit_theta': 1.083,    # unitless; Arrhenius coefficient for nitrification
    'kon_theta': 1.074,     # unitless; Arrhenius coefficient for OrgN hydrolysis
    'kdnit_theta': 1.08,    # unitless; Arrhenius coefficient for denitrification
    'rnh4_theta': 1.047,    # unitless; Arrhenius coefficient for sediment NH4 release
    'vno3_theta': 1.045,    # unitless; Arrhenius coefficient for sediment NO3 denitrification
    'vson_theta': 1.024,    # unitless; Arrhenius coefficient for OrgN settling
    'KsOxdn': 0.1,          # mg-O2/L; oxygen half-saturation for denitrification inhibition
    'PN': 0.5,              # unitless; algal preference fraction for NH4 over NO3
    'PNb': 0.5,             # unitless; benthic algal preference fraction for NH4 over NO3
    'use_OrgN': True,       # bool; enable organic-nitrogen state variable
}
