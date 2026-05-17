"""v3 NSM1 alkalinity parameter defaults.

Consumed by: ``Alkalinity`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``AlkalinityStaticVariables`` /
``DEFAULT_ALKALINITY``.

Corrections applied:

* ``r_alkden`` corrected ``4/14/1000`` -> ``1/14/1000`` eq/mg-N
  (NSM1-SCI-N1, gold-standard spec A2, 2026-05-16). Denitrification
  produces **1 eq of alkalinity per mole of NO3-N reduced**
  (NO3- -> 1/2 N2), per CE-QUAL-W2 ``water-quality.f90:3157`` and
  Stumm & Morgan, *Aquatic Chemistry* 3rd ed. (Ch. 4, alkalinity
  changes from N redox). The shipped ``4/14/1000`` is a **deliberate,
  reference-anchored divergence from the upstream NSM1 Fortran**
  (``modAlkalinity.f90:54``), which is 4x the stoichiometric value;
  the error is shared Fortran = v1 = v3 (invisible to v1<->v3 parity)
  and is reported upstream in
  ``design/clearwater_modules_v3_nsm1_upstream_fortran_defects.md``.
  Cross-check: ``r_alkn = 2/14/1000`` (nitrification consumes 2
  eq/mol-N) is the textbook value on the same basis, confirming the
  1:2 denitrification:nitrification alkalinity ratio.

All entries are stoichiometric ratios for alkalinity sources/sinks:
algal photosynthesis, algal N uptake, nitrification, denitrification, and
their benthic counterparts.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'r_alkaa': 14.0 / 106.0 / 12.0 / 1000.0,    # eq/mg-C; alk stoich for algal photosynthesis (Alk source)
    'r_alkan': 18.0 / 106.0 / 12.0 / 1000.0,    # eq/mg-C; alk stoich for algal N uptake (Alk sink)
    'r_alkn': 2.0 / 14.0 / 1000.0,              # eq/mg-N; alk stoich for nitrification (Alk sink)
    'r_alkden': 1.0 / 14.0 / 1000.0,            # eq/mg-N; alk stoich for denitrification (Alk source). 1 eq/mol-N (NO3- -> 1/2 N2); CE-QUAL-W2 water-quality.f90:3157 + Stumm & Morgan. NSM1-SCI-N1: deliberate divergence from upstream Fortran modAlkalinity.f90:54 (=4/14/1000, 4x too high; Fortran=v1=v3).
    'r_alkba': 14.0 / 106.0 / 12.0 / 1000.0,    # eq/mg-C; alk stoich for benthic algae photosynthesis
    'r_alkbn': 18.0 / 106.0 / 12.0 / 1000.0,    # eq/mg-C; alk stoich for benthic algae N uptake
}
