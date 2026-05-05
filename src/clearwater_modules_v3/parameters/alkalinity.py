"""v3 NSM1 alkalinity parameter defaults.

Consumed by: ``Alkalinity`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``AlkalinityStaticVariables`` /
``DEFAULT_ALKALINITY``.
Corrections applied: none.

All entries are stoichiometric ratios for alkalinity sources/sinks:
algal photosynthesis, algal N uptake, nitrification, denitrification, and
their benthic counterparts.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'r_alkaa': 14.0 / 106.0 / 12.0 / 1000.0,    # eq/mg-C; alk stoich for algal photosynthesis (Alk source)
    'r_alkan': 18.0 / 106.0 / 12.0 / 1000.0,    # eq/mg-C; alk stoich for algal N uptake (Alk sink)
    'r_alkn': 2.0 / 14.0 / 1000.0,              # eq/mg-N; alk stoich for nitrification (Alk sink)
    'r_alkden': 4.0 / 14.0 / 1000.0,            # eq/mg-N; alk stoich for denitrification (Alk source)
    'r_alkba': 14.0 / 106.0 / 12.0 / 1000.0,    # eq/mg-C; alk stoich for benthic algae photosynthesis
    'r_alkbn': 18.0 / 106.0 / 12.0 / 1000.0,    # eq/mg-C; alk stoich for benthic algae N uptake
}
