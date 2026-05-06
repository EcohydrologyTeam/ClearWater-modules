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
* **Nitrogen theta transposition fix (Phase 9.E)**. The four nitrogen Arrhenius
  theta values were transposed in pairs in v1's port from Fortran. Evidence:
  (1) Fortran ``modNitrogen.f90:82, 89, 95, 100`` initializes
  ``kon%theta=1.047``, ``rnh4%theta=1.074``, ``kdnit%theta=1.045``,
  ``vno3%theta=1.08`` (lines 82, 89, 95, 100); v1 has the values transposed
  within each pair (kon-rnh4 swap and kdnit-vno3 swap). (2) The phosphorus
  parallel-process check confirms the convention: Fortran/v1/v3 all agree on
  ``kop_theta=1.047`` (organic-P hydrolysis, parallel to ``kon``) and
  ``rpo4_theta=1.074`` (sediment-P release, parallel to ``rnh4``); the
  nitrogen pair should mirror this and does in Fortran. (3) Literature
  convention (Chapra 1997, QUAL2K manual): organic-matter hydrolysis uses
  ``theta=1.047`` (the universal NSM1 organic-matter Arrhenius default,
  matching ``mu_max_theta``, ``kdp_theta``, ``krp_theta``, ``kpoc_theta``,
  ``kdoc_theta``, ``kop_theta``, ``kpom_theta``, ``kbod_theta``);
  sediment-water exchange velocities use steeper temperature dependence
  (~1.074-1.08). v3 now adopts the Fortran-aligned values.

See ``parameter_defaults_corrections.md`` Sections 1.8, 1.10 (Phase 9.E
nitrogen theta), and the Phase 9.E commit message.
"""

DEFAULTS: dict[str, float | int | bool] = {
    'KNR': 0.6,             # mg-O2/L; oxygen half-saturation for nitrification inhibition (1 - exp(-KNR*DOX))
    'knit_20': 0.1,         # 1/d; nitrification rate at 20 C
    'kon_20': 0.1,          # 1/d; organic-N hydrolysis (ammonification) rate at 20 C
    'kdnit_20': 0.002,      # 1/d; denitrification rate at 20 C
    'rnh4_20': 0.0,         # FIXME(phase1-audit): sediment NH4 release silently disabled; verify gated by use_SedFlux. 1/d at 20 C.
    'vno3_20': 0.0,         # FIXME(phase1-audit): sediment NO3 denitrification silently disabled; verify gated by use_SedFlux. 1/d at 20 C.
    'vson_20': 0.01,        # m/d; OrgN settling velocity at 20 C (matches Fortran modGlobalParam.f90:92 vson=0.01 and v1 GlobalVars vson=0.01; corrected from 0.1 in Phase 9.C, see corrections doc Section 1.8)
    'knit_theta': 1.083,    # unitless; Arrhenius coefficient for nitrification (matches Fortran modNitrogen.f90)
    'kon_theta': 1.047,     # unitless; OrgN hydrolysis Arrhenius (Phase 9.E correction; was 1.074 in v1/v3, transposed with rnh4_theta during v1 port; matches Fortran kon%theta=1.047 modNitrogen.f90:89 and the universal NSM1 organic-matter convention)
    'kdnit_theta': 1.045,   # unitless; denitrification Arrhenius (Phase 9.E correction; was 1.08 in v1/v3, transposed with vno3_theta during v1 port; matches Fortran kdnit%theta=1.045 modNitrogen.f90:95 and Chapra 1997)
    'rnh4_theta': 1.074,    # unitless; sediment NH4 release Arrhenius (Phase 9.E correction; was 1.047 in v1/v3, transposed with kon_theta during v1 port; matches Fortran rnh4%theta=1.074 modNitrogen.f90:82 and the rpo4_theta=1.074 phosphorus parallel)
    'vno3_theta': 1.08,     # unitless; sediment denitrification Arrhenius (Phase 9.E correction; was 1.045 in v1/v3, transposed with kdnit_theta during v1 port; matches Fortran vno3%theta=1.08 modNitrogen.f90:100)
    'vson_theta': 1.024,    # unitless; Arrhenius coefficient for OrgN settling (v3 addition; v1 uses vson raw, no Fortran counterpart; see module docstring)
    'KsOxdn': 0.1,          # mg-O2/L; oxygen half-saturation for denitrification inhibition
    'PN': 0.5,              # unitless; algal preference fraction for NH4 over NO3
    'PNb': 0.5,             # unitless; benthic algal preference fraction for NH4 over NO3
    'use_OrgN': True,       # bool; enable organic-nitrogen state variable
}
