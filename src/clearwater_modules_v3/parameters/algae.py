"""v3 NSM1 floating algae parameter defaults.

Consumed by: ``FloatingAlgae`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``AlgaeStaticVariables`` /
``DEFAULT_ALGAE``.
Corrections applied: none.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'AWd': 100.0,                   # mg-D/ug-Chla; algal dry weight ratio (stoichiometry)
    'AWc': 40.0,                    # mg-C/ug-Chla; algal carbon ratio
    'AWn': 7.2,                     # mg-N/ug-Chla; algal nitrogen ratio
    'AWp': 1.0,                     # mg-P/ug-Chla; algal phosphorus ratio
    'AWa': 1000.0,                  # ug-Chla/ug-Chla; chlorophyll per algal unit (rda = AWd/AWa)
    'KL': 10.0,                     # W/m^2; light limitation half-saturation
    'KsN': 0.04,                    # mg-N/L; nitrogen half-saturation for algal growth
    'KsP': 0.0012,                  # mg-P/L; phosphorus half-saturation for algal growth
    'mu_max_20': 1.0,               # 1/d; max algal growth rate at 20 C
    'kdp_20': 0.15,                 # 1/d; algal death rate at 20 C
    'krp_20': 0.2,                  # 1/d; algal respiration rate at 20 C
    'mu_max_theta': 1.047,          # unitless; Arrhenius coefficient for growth
    'kdp_theta': 1.047,             # unitless; Arrhenius coefficient for death
    'krp_theta': 1.047,             # unitless; Arrhenius coefficient for respiration
    'vsap': 0.15,                   # m/d; algal settling velocity
    'growth_rate_option': 1,        # selector: 1=Multiplicative, 2=Min, 3=Harmonic
    'light_limitation_option': 1,   # selector: 1=Half-Saturation, 2=Smith, 3=Steele
}
