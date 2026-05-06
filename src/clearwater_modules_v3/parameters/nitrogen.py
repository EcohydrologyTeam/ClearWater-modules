"""v3 NSM1 nitrogen parameter defaults.

Consumed by: ``Nitrogen`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``NitrogenStaticVariables`` /
``DEFAULT_NITROGEN``.
Corrections applied:

* ``vson_20`` corrected from ``0.1`` to ``0.01`` m/d (Phase 9.C audit fix). The
  prior 0.1 was an internal v3 inconsistency: v1 GlobalVars, v3 ``global_vars``,
  and Fortran ``modGlobalParam.f90:92`` all use ``0.01`` m/d. The 0.1 in this
  module did not match any reference and was 10x the canonical value.
* ``vson_theta`` removed in Phase 9.E. The parameter was originally added in
  Phase 1.2 by analogy with rate-constant theta values. Phase 2.B's
  ``organic_nitrogen_settling`` consumed it via
  ``arrhenius_correction(T, vson_20, vson_theta)`` with a docstring claiming
  "parity with v1" -- but the parity claim was false. Both v1
  (``processes.py:1333`` ``OrgN_Settling = vson / depth * OrgN``) and Fortran
  (``modNitrogen.f90:233`` ``OrgN_Settling = vson(r) / depth * OrgN``) use raw
  ``vson`` without temperature correction. Fortran's deliberate type
  distinction (rate constants are ``TempCorrectionStruct`` with ``%rc20`` and
  ``%theta``; settling velocities are plain ``real(R8)``) reflects the
  physical convention that biochemical reaction rates scale with Arrhenius
  temperature dependence (Q10 approx 2-3) but settling velocities depend on
  water viscosity (theta approx 1.009, much smaller than the 1.024 v3 had
  applied). Phase 9.E removes the parameter; the Process now uses raw
  ``vson_20`` directly, matching v1 and Fortran exactly.
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
    'rnh4_20': 0.0,         # 1/d at 20 C; sediment NH4 release rate. Held at zero by design: v2/v3 Nitrogen does NOT gate ammonium_from_bed by use_SedFlux, so the zero default is the de facto gate. Phase 9.F.C added a defensive NotImplementedError guard in Nitrogen.__init__ that fires if use_SedFlux=True is passed (sediment-flux feature requires NSM2 path). See corrections doc Section 2.1.
    'vno3_20': 0.0,         # 1/d at 20 C; sediment NO3 denitrification rate. Same gating story as rnh4_20 (held at zero by design; ungated nitrate_bed_denitrification consumer; defensive guard added in Phase 9.F.C). See corrections doc Section 2.1.
    'vson_20': 0.01,        # m/d; OrgN settling velocity at 20 C (matches Fortran modGlobalParam.f90:92 vson=0.01 and v1 GlobalVars vson=0.01; corrected from 0.1 in Phase 9.C, see corrections doc Section 1.8)
    'knit_theta': 1.083,    # unitless; Arrhenius coefficient for nitrification (matches Fortran modNitrogen.f90)
    'kon_theta': 1.047,     # unitless; OrgN hydrolysis Arrhenius (Phase 9.E correction; was 1.074 in v1/v3, transposed with rnh4_theta during v1 port; matches Fortran kon%theta=1.047 modNitrogen.f90:89 and the universal NSM1 organic-matter convention)
    'kdnit_theta': 1.045,   # unitless; denitrification Arrhenius (Phase 9.E correction; was 1.08 in v1/v3, transposed with vno3_theta during v1 port; matches Fortran kdnit%theta=1.045 modNitrogen.f90:95 and Chapra 1997)
    'rnh4_theta': 1.074,    # unitless; sediment NH4 release Arrhenius (Phase 9.E correction; was 1.047 in v1/v3, transposed with kon_theta during v1 port; matches Fortran rnh4%theta=1.074 modNitrogen.f90:82 and the rpo4_theta=1.074 phosphorus parallel)
    'vno3_theta': 1.08,     # unitless; sediment denitrification Arrhenius (Phase 9.E correction; was 1.045 in v1/v3, transposed with kdnit_theta during v1 port; matches Fortran vno3%theta=1.08 modNitrogen.f90:100)
    # Phase 9.E removed `vson_theta` (was 1.024). v1 and Fortran both use
    # raw `vson` for OrgN settling without Arrhenius correction; the
    # parameter was an unjustified v3 addition. See module docstring.
    'KsOxdn': 0.1,          # mg-O2/L; oxygen half-saturation for denitrification inhibition
    'PN': 0.5,              # unitless; algal preference fraction for NH4 over NO3
    'PNb': 0.5,             # unitless; benthic algal preference fraction for NH4 over NO3
    'use_OrgN': True,       # bool; enable organic-nitrogen state variable
}
