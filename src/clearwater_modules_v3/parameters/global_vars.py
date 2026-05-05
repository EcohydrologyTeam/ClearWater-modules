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

Corrections applied: none (the seven critical corrections live in the modules
the affected parameters were migrated to).

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
    'vb': 0.01,                         # FIXME(phase1-audit): m/d; burial velocity, magnitude not validated
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
    'q_solar': 500.0,                   # FIXME(phase1-audit): W/m^2 (v1 docstring says 1/d but value/usage is W/m^2)
    'Solid': 1,                         # mg/L; suspended solids concentration (toy; overridden per cell)
    # Light attenuation (Beer-Lambert composite extinction coefficient)
    'lambda0': 0.02,                    # 1/m; background light extinction (clear water)
    'lambda1': 0.0088,                  # (1/m)/(ug-Chla/L); linear self-shading by chlorophyll
    'lambda2': 0.054,                   # unitless; non-linear chlorophyll extinction coefficient
    'lambdas': 0.052,                   # FIXME(phase1-audit): L/(mg*m); ISS extinction parameter currently disabled in code path
    'lambdam': 0.0174,                  # L/(mg*m); POM extinction coefficient
    'Fr_PAR': 0.47,                     # unitless; PAR fraction of total solar radiation
}
