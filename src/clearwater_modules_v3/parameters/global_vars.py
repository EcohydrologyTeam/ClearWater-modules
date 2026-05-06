"""v3 NSM1 global runtime/environmental defaults and IC placeholders.

Consumed by: most NSM1 Process classes that need shared environmental scalars
(temperature, depth, velocity, light extinction parameters, settling velocities
for organic matter not routed to a specific constituent, etc.).
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``GlobalVars`` /
``DEFAULT_GLOBALVARS``, with the following entries migrated to other groups:

* ``vsop``, ``vs`` -> ``phosphorus``
* ``SOD_20``, ``SOD_theta``, ``kaw_20_user``, ``kah_20_user``, ``kaw_theta``,
  ``kah_theta``, ``hydraulic_reaeration_option``, ``wind_reaeration_option``
  -> ``dox``
* ``pressure_mb`` -> ``global_parameters``

Corrections applied:

* ``lambdam`` corrected from ``0.0174`` to ``0.174`` L/(mg*m) (Phase 9.C audit
  fix). Fortran ``modGlobalParam.f90:68`` initializes ``lambdam = 0.174``; v1
  GlobalVars used ``0.0174`` (likely typo, 10x lower than canonical) and v3
  inherited the v1 typo. The 0.174 value matches QUAL2K Table 6 and is used
  throughout the legacy v1 NSM test suite (e.g., ``test_7_nsm_algae_calculations``,
  ``test_10_nsm_carbon_calculations``, ``test_17_nsm_N2_calculations``).
  See ``parameter_defaults_corrections.md`` Section 1.9.

The seven other critical sentinel-999 corrections live in the modules the
affected parameters were migrated to (``dox``, ``phosphorus``,
``global_parameters``).

Many runtime scalars (``dt``, ``depth``, ``TwaterC``, ``velocity``, ``flow``,
``topwidth``, ``slope``, ``shear_velocity``, ``wind_speed``, ``Solid``) are
toy placeholder values; production simulations override them per cell/time via
the model's coupling to the hydraulic driver.
"""

DEFAULTS: dict[str, float | int | bool] = {
    # Settling velocities for organic-matter constituents not specific to a single nutrient
    'vson': 0.01,                       # m/d; organic-N settling velocity
    'vsoc': 0.01,                       # m/d; POC settling velocity
    # Generic thermal Arrhenius coefficient
    'theta': 1.047,                     # unitless; default Arrhenius coefficient for processes lacking an explicit theta
    # Burial / sediment composition
    'vb': 6.85e-6,                      # m/d; sediment burial velocity (= 0.0025 m/yr = 0.25 cm/yr). Phase 9.F.A correction; was 0.01 m/d in v1/v3 pre-9.F, a 1460x v1 unit-conversion bug (v1 dropped Fortran's runtime /365 conversion without rescaling the numerical default from m/yr to m/d). Canonical value matches WASP7/WASP8 Appendix A, Fortran modGlobalParam.f90:138 (0.0025 m/yr / 365), and Di Toro 2001 sediment-flux model. See corrections doc Section 1.14.
    'fcom': 0.4,                        # unitless; fraction of sediment as combustible organic matter
    # Simulation control (toy placeholders; overridden at runtime)
    'dt': 1.0,                          # d; default timestep (toy value; overridden by model)
    'depth': 1.5,                       # m; default water depth (toy; overridden per cell)
    'TwaterC': 20.0,                    # C; reference water temperature (typically overridden per cell/time)
    'velocity': 1.0,                    # m/s; toy velocity (overridden per cell)
    'flow': 2.0,                        # m^3/s; toy flow (overridden per cell)
    'topwidth': 1.0,                    # m; toy top width (overridden per cell)
    'slope': 2.0,                       # unitless; toy slope (overridden per cell)
    'shear_velocity': 4.0,              # m/s; toy shear velocity (overridden per cell)
    'wind_speed': 4.0,                  # m/s; toy wind speed (overridden per cell/time)
    'q_solar': 500.0,                   # W/m^2; total incident solar radiation at the water surface. Note: v1's docstring incorrectly labeled this parameter as 1/d, but the value (500) and the consumption pattern (Beer-Lambert PAR scaling via utils/light.py:PAR and processes/pathogen.py:_rate_light_decay) are unambiguously W/m^2. Resolved in Phase 9.F; see corrections doc Section 2.7.
    'Solid': 1,                         # mg/L; suspended solids concentration (toy; overridden per cell)
    # Light attenuation (Beer-Lambert composite extinction coefficient)
    'lambda0': 0.02,                    # 1/m; background light extinction (clear water)
    'lambda1': 0.0088,                  # (1/m)/(ug-Chla/L); linear self-shading by chlorophyll
    'lambda2': 0.054,                   # unitless; non-linear chlorophyll extinction coefficient
    'lambdas': 0.052,                   # L/(mg*m); ISS (suspended-solids) extinction parameter; active per Phase 9.C three-way audit verification: applied unconditionally in utils/light.py (matches v1 shared/processes.py:232 and Fortran modGlobalParam.f90 LightExtCoefficient). The earlier Phase 0 "commented out / defined but not used" claim was a documentation defect, corrected in Phase 9.C and Phase 9.F. Multi-solid-class generalization out of scope for 1.0.0; see corrections doc Section 2.8.
    'lambdam': 0.174,                   # L/(mg*m); POM extinction coefficient (matches Fortran modGlobalParam.f90:68 and QUAL2K Table 6; corrected from v1 typo 0.0174 in Phase 9.C, see corrections doc Section 1.9)
    'Fr_PAR': 0.47,                     # unitless; PAR fraction of total solar radiation
}
