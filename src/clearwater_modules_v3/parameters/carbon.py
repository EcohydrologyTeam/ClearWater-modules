"""v3 NSM1 carbon parameter defaults.

Consumed by: ``Carbon`` Process (POC, DOC, DIC kinetics).
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``CarbonStaticVariables`` /
``DEFAULT_CARBON``.
Corrections applied: none.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'f_pocp': 0.8,          # unitless; fraction of algal death routed to POC (1-f_pocp -> DOC). NSM1-SCI-A2 (spec C1): 0.9->0.8 to match the operative FloatingAlgae _FDP_DEFAULTS value (Carbon consumes the FloatingAlgae-cached mortality rate; this carbon-side copy is kept consistent, not an independent knob). CE-QUAL-W2 APOM ~0.8.
    'kdoc_20': 0.01,        # 1/d; DOC decomposition rate at 20 C
    'kdoc_theta': 1.047,    # unitless; Arrhenius coefficient for DOC decay
    'f_pocb': 0.8,          # unitless; fraction of benthic algae death routed to POC (1-f_pocb -> DOC). NSM1-SCI-A2 (spec C1): 0.9->0.8 to match the operative BenthicAlgae _BENTHIC_FDP_DEFAULTS value (Carbon consumes the BenthicAlgae-cached rate; kept consistent). CE-QUAL-W2 APOM ~0.8.
    'kpoc_20': 0.005,       # 1/d; POC decomposition rate at 20 C
    'kpoc_theta': 1.047,    # unitless; Arrhenius coefficient for POC decay
    'KsOxmc': 1.0,          # mg-O2/L; oxygen half-saturation for organic-matter mineralization
    'pCO2': 383.0,          # ppm; atmospheric CO2 partial pressure (~2024 value)
    'FCO2': 0.2,            # unitless; fraction of DIC as aqueous CO2 (placeholder; full carbonate speciation is NSM2 territory)
    'roc': 32.0 / 12.0,     # mg-O2/mg-C; respiration stoichiometry (O2 consumed per C oxidized)
}
