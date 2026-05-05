"""v3 NSM1 global feature-flag parameter defaults.

Consumed by: most NSM1 Process classes (each consults the relevant ``use_*``
flag to gate constituent-specific kinetics).
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``GlobalParameters`` /
``DEFAULT_GLOBALPARAMETERS``, plus ``pressure_mb`` migrated from v1
``GlobalVars`` since it is a model-level environmental scalar consumed by
DOX (O2sat) and N2 (N2sat) processes.
Corrections applied: ``pressure_mb=1013.25`` (was 2026.5).
See ``parameter_defaults_corrections.md`` Section 1.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'use_NH4': True,            # bool; enable NH4 (ammonium) state variable
    'use_NO3': True,            # bool; enable NO3 (nitrate) state variable
    'use_OrgN': True,           # bool; enable organic-N state variable
    'use_OrgP': True,           # bool; enable organic-P state variable
    'use_TIP': True,            # bool; enable total inorganic phosphorus state variable
    'use_SedFlux': False,       # bool; enable sediment flux features (rnh4_20, vno3_20, rpo4_20 are zero unless overridden)
    'use_POC': True,            # bool; enable particulate organic carbon
    'use_DOC': True,            # bool; enable dissolved organic carbon
    'use_DOX': True,            # bool; enable dissolved oxygen
    'use_DIC': True,            # bool; enable dissolved inorganic carbon
    'use_Algae': True,          # bool; enable floating algae
    'use_Balgae': True,         # bool; enable benthic algae
    'use_N2': True,             # bool; enable dissolved nitrogen gas
    'use_Pathogen': True,       # bool; enable pathogen tracer
    'use_Alk': True,            # bool; enable alkalinity
    'use_POM': True,            # bool; enable particulate organic matter
    'pressure_mb': 1013.25,     # hPa; standard sea-level pressure (was 2026.5 in v1, see corrections doc Section 1)
}
