"""v3 NSM1 dissolved oxygen (DOX) parameter defaults.

Consumed by: ``DOX`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``DOXStaticVariables`` /
``DEFAULT_DOX``, plus SOD and reaeration parameters (``SOD_20``, ``SOD_theta``,
``kaw_20_user``, ``kah_20_user``, ``kaw_theta``, ``kah_theta``,
``hydraulic_reaeration_option``, ``wind_reaeration_option``) migrated from v1
``GlobalVars`` since they are DOX-specific kinetics.
Corrections applied: ``SOD_20=1.0`` (was 999), ``SOD_theta=1.060`` (was 999),
``kaw_20_user=0.0`` (was 999), ``kah_20_user=0.0`` (was 999).
See ``parameter_defaults_corrections.md`` Section 1.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'ron': 2.0 * 32.0 / 14.0,           # mg-O2/mg-N; stoichiometry of O2 consumed per N nitrified
    'KsSOD': 1.0,                       # mg-O2/L; oxygen half-saturation for SOD
    'SOD_20': 1.0,                      # g-O2/m^2/d; SOD at 20 C (was 999 in v1, see corrections doc Section 1)
    'SOD_theta': 1.060,                 # unitless; Arrhenius coefficient for SOD (was 999 in v1, see corrections doc Section 1)
    'kaw_20_user': 0.0,                 # m/d; user-override wind reaeration at 20 C, disabled unless opted in (was 999 in v1, see corrections doc Section 1)
    'kah_20_user': 0.0,                 # 1/d; user-override hydraulic reaeration at 20 C, disabled unless opted in (was 999 in v1, see corrections doc Section 1)
    'kaw_theta': 1.024,                 # unitless; Arrhenius coefficient for wind reaeration
    'kah_theta': 1.024,                 # unitless; Arrhenius coefficient for hydraulic reaeration
    'hydraulic_reaeration_option': 1,   # selector: 1-9 hydraulic reaeration formula
    'wind_reaeration_option': 1,        # selector: 1-13 wind reaeration formula
}
