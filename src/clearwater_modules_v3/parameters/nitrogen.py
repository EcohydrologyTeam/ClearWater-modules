"""v3 NSM1 nitrogen parameter defaults.

Consumed by: ``Nitrogen`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``NitrogenStaticVariables`` /
``DEFAULT_NITROGEN``.
Corrections applied:

* ``vson_20`` corrected from ``0.1`` to ``0.01`` m/d (Phase 9.C audit fix). The
  prior 0.1 was an internal v3 inconsistency: v1 GlobalVars, v3 ``global_vars``,
  and Fortran ``modGlobalParam.f90:92`` all use ``0.01`` m/d. The 0.1 in this
  module did not match any reference and was 10x the canonical value.
* ``vson_theta=1.024`` is a v3-only addition. v1 uses ``vson`` raw (no Arrhenius
  correction) and Fortran has no ``vson_theta``. Phase 2.B added the Arrhenius
  correction here for consistency with the other settling-velocity Arrhenius
  forms; the Process consumes ``arrhenius_correction(T, vson_20, vson_theta)``.

See ``parameter_defaults_corrections.md`` Sections 1.8 and 4.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'KNR': 0.6,             # mg-O2/L; oxygen half-saturation for nitrification inhibition (1 - exp(-KNR*DOX))
    'knit_20': 0.1,         # 1/d; nitrification rate at 20 C
    'kon_20': 0.1,          # 1/d; organic-N hydrolysis (ammonification) rate at 20 C
    'kdnit_20': 0.002,      # 1/d; denitrification rate at 20 C
    'rnh4_20': 0.0,         # FIXME(phase1-audit): sediment NH4 release silently disabled; verify gated by use_SedFlux. 1/d at 20 C.
    'vno3_20': 0.0,         # FIXME(phase1-audit): sediment NO3 denitrification silently disabled; verify gated by use_SedFlux. 1/d at 20 C.
    'vson_20': 0.01,        # m/d; OrgN settling velocity at 20 C (matches Fortran modGlobalParam.f90:92 vson=0.01 and v1 GlobalVars vson=0.01; corrected from 0.1 in Phase 9.C, see corrections doc Section 1.8)
    'knit_theta': 1.083,    # unitless; Arrhenius coefficient for nitrification (matches Fortran)
    'kon_theta': 1.074,     # FIXME(phase9c-audit): unitless; v3/v1 use 1.074 but Fortran modNitrogen.f90:89 uses kon%theta=1.047. Possible v1 typo; left as-is pending LimnoTech reconciliation.
    'kdnit_theta': 1.08,    # FIXME(phase9c-audit): unitless; v3/v1 use 1.08 but Fortran modNitrogen.f90:95 uses kdnit%theta=1.045. Possible v1 typo; left as-is pending LimnoTech reconciliation.
    'rnh4_theta': 1.047,    # FIXME(phase9c-audit): unitless; v3/v1 use 1.047 but Fortran modNitrogen.f90:82 uses rnh4%theta=1.074. Possible v1 swap; left as-is pending LimnoTech reconciliation.
    'vno3_theta': 1.045,    # FIXME(phase9c-audit): unitless; v3/v1 use 1.045 but Fortran modNitrogen.f90:100 uses vno3%theta=1.08. Possible v1 swap; left as-is pending LimnoTech reconciliation.
    'vson_theta': 1.024,    # unitless; Arrhenius coefficient for OrgN settling (v3 addition; v1 uses vson raw, no Fortran counterpart; see module docstring)
    'KsOxdn': 0.1,          # mg-O2/L; oxygen half-saturation for denitrification inhibition
    'PN': 0.5,              # unitless; algal preference fraction for NH4 over NO3
    'PNb': 0.5,             # unitless; benthic algal preference fraction for NH4 over NO3
    'use_OrgN': True,       # bool; enable organic-nitrogen state variable
}
