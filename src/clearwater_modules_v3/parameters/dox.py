"""v3 NSM1 dissolved oxygen (DOX) parameter defaults.

Consumed by: ``DOX`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``DOXStaticVariables`` /
``DEFAULT_DOX``, plus SOD and reaeration parameters (``SOD_20``, ``SOD_theta``,
``kaw_20_user``, ``kah_20_user``, ``kaw_theta``, ``kah_theta``,
``hydraulic_reaeration_option``, ``wind_reaeration_option``) migrated from v1
``GlobalVars`` since they are DOX-specific kinetics.
Corrections applied: ``SOD_20=1.0`` (was 999), ``SOD_theta=1.060`` (was 999),
``kaw_20_user=0.0`` (was 999), ``kah_20_user=0.0`` (was 999),
``hydraulic_reaeration_option=5`` (Phase 9.E; was 1 = user-supplied path
with silent constant default; new default matches QUAL2K's documented
"Internal/Covar 1976" depth-piecewise empirical formula).
See ``parameter_defaults_corrections.md`` Section 1.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'ron': 2.0 * 32.0 / 14.0,           # mg-O2/mg-N; stoichiometry of O2 consumed per N nitrified
    'KsSOD': 1.0,                       # mg-O2/L; oxygen half-saturation for SOD
    'SOD_20': 1.0,                      # g-O2/m^2/d; SOD at 20 C (was 999 in v1, see corrections doc Section 1)
    'SOD_theta': 1.060,                 # unitless; Arrhenius coefficient for SOD (was 999 in v1, see corrections doc Section 1)
    'kaw_20_user': 0.0,                 # m/d; user-override wind reaeration at 20 C, only consulted when wind_reaeration_option == 1 (was 999 in v1, see corrections doc Section 1)
    'kah_20_user': 0.0,                 # 1/d; user-override hydraulic reaeration at 20 C, only consulted when hydraulic_reaeration_option == 1 (was 999 in v1, see corrections doc Section 1)
    'kaw_theta': 1.024,                 # unitless; Arrhenius coefficient for wind reaeration
    'kah_theta': 1.024,                 # unitless; Arrhenius coefficient for hydraulic reaeration
    'min_reaeration_ka': 0.0,           # 1/d; OPT-IN minimum reaeration floor on ka_tc (NSM1-DOX-F2, spec C4). 0.0 = OFF (preserves v1/Fortran parity; DEFAULT). >0 follows the CE-QUAL-W2 MINKL precedent: ka_tc = max(ka_tc, min_reaeration_ka), preventing silent zero atmospheric reaeration when hydraulic_reaeration_option==1 with kah_20_user==0.
    'min_reaeration_depth': 0.0,        # m; OPT-IN floor on the cell mean depth fed to the reaeration formulas (kah_20 + ka_tc). 0.0 = OFF (DEFAULT; byte-identical). >0 clamps depth = max(depth, min_reaeration_depth) so the inverse-depth terms (depth**-1.85 etc.) cannot blow up at a sub-physical newly-wet (~0 m) cell in a coupled HEC-RAS-2D run. Shared by DOX/N2/Carbon, which all source this reaeration menu. Companion to min_reaeration_ka (floors the result); this floors the input depth, which the kah_20 hydraulic term needs since min_reaeration_ka cannot un-blow an already-exploded ka_tc cheaply.

    'hydraulic_reaeration_option': 5,   # selector: 1-9 hydraulic reaeration formula. Default 5 = Cover 1976 / Internal (depth-piecewise blend of Owens-Gibbs / O'Connor-Dobbins / Churchill); matches QUAL2K Chapra & Pelletier 2008 manual p56 documented default. Phase 9.E corrected from 1 (user-supplied path) to 5; see corrections doc Section 1.6.
    'wind_reaeration_option': 1,        # selector: 1-13 wind reaeration formula. Default 1 (user-supplied with kaw_20_user=0.0) corresponds to QUAL2K wind option 1 "omitted" (manual p57); appropriate for stream/river-focused NSM1 applications where hydraulic dominates. Lake/reservoir users should opt into option 4 (Banks-Herrera) or 5 (Wanninkhof).
}
