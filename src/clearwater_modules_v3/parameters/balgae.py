"""v3 NSM1 benthic algae parameter defaults.

Consumed by: ``BenthicAlgae`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``BalgaeStaticVariables`` /
``DEFAULT_BALGAE``.
Corrections applied: none.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'BWd': 100.0,                       # mg-D/g-D; benthic algae dry weight ratio
    'BWc': 40.0,                        # mg-C/g-D; benthic algae carbon content
    'BWn': 7.2,                         # mg-N/g-D; benthic algae nitrogen content
    'BWp': 1.0,                         # mg-P/g-D; benthic algae phosphorus content
    'BWa': 3500.0,                      # g-D/m^2; benthic algae areal density reference
    'KLb': 10.0,                        # W/m^2; benthic light limitation half-saturation
    'KsNb': 0.25,                       # mg-N/L; nitrogen half-saturation for benthic algae
    'KsPb': 0.125,                      # mg-P/L; phosphorus half-saturation for benthic algae
    'Ksb': 10.0,                        # g-D/m^2; benthic biomass half-saturation (FSb space limitation)
    'mub_max_20': 0.4,                  # 1/d; max benthic algae growth rate at 20 C
    'krb_20': 0.2,                      # 1/d; benthic algae respiration rate at 20 C
    'kdb_20': 0.3,                      # 1/d; benthic algae death rate at 20 C
    'mub_max_theta': 1.047,             # unitless; Arrhenius coefficient for growth
    'krb_theta': 1.06,                  # unitless; Arrhenius coefficient for respiration (note: 1.06 differs from typical 1.047)
    'kdb_theta': 1.047,                 # unitless; Arrhenius coefficient for death
    'b_growth_rate_option': 1,          # selector: 1=Multiplicative, 2=Min, 3=Harmonic
    'b_light_limitation_option': 1,     # selector: 1=Half-Saturation, 2=Smith, 3=Steele
    'Fw': 0.9,                          # unitless; fraction of benthic photosynthesis released as CO2 (1-Fw routed to DOC)
    'Fb': 0.9,                          # unitless; fraction of benthic biomass released in death (1-Fb routed to DOC)
}
