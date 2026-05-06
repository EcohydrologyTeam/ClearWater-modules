# TSM v3: per-component flux diagnostics and equilibrium temperature

The v3 ``Temperature`` process exposes seven per-component heat-flux
pathway outputs and a Newton-Raphson equilibrium-temperature
diagnostic. Both are **opt-in via registry pre-registration** — they
cost nothing when not requested and, when requested, are written to
the registry alongside ``water_temperature`` each substep.

This guide shows how to enable them, what the values mean, and the
conventions for interpreting their signs.

## Quick reference: the seven flux components

| Registry key | Description | Sign | Composition into `q_net` |
|---|---|---|---|
| `q_sensible` | Sensible heat flux | signed by `T_air − T_water` | `+ q_sensible` |
| `q_latent` | Latent (evaporative) heat flux **magnitude** | signed by `e_sat − e_air` | `− q_latent` |
| `q_longwave_up` | Upwelling longwave (Stefan-Boltzmann) **magnitude** | always ≥ 0 | `− q_longwave_up` |
| `q_longwave_down` | Atmospheric (downwelling) longwave **magnitude** | always ≥ 0 | `+ q_longwave_down` |
| `q_solar` | Solar input (passthrough of `solar_radiation` forcing) | always ≥ 0 | `+ q_solar` |
| `q_sediment` | Sediment heat flux | signed by `T_sed − T_water` | `+ q_sediment` |
| `q_net` | Net surface heat flux at the water column | composition | — |

All values are in **W/m²**. The convention is that **positive
contributions to `q_net` heat the water column**; `q_latent` and
`q_longwave_up` are stored as positive magnitudes and subtracted at
composition. This matches v1 ``tsm/processes.py:q_net`` and Fortran-A
``modTemperature.f90:257``.

`q_sensible` and `q_sediment` carry sign through their respective
temperature-gradient arguments — a positive value means the
corresponding source is heating the water column.

## Enabling the diagnostics

The pattern is identical to N2's optional ``total_dissolved_gas``
output: register the variable upfront with the same shape as
``water_temperature``, and the ``Temperature`` process will write to
it each substep.

```python
from clearwater_data import Variable
from clearwater_modules_v3 import Model

# ... build your model and registry as usual ...

# Pre-register the diagnostic outputs you care about.
# Each variable must have the same spatial shape as water_temperature.
shape_template = registry.get_variable("water_temperature").get()
zeros = xr.zeros_like(shape_template)

for diagnostic_name in (
    "q_net",
    "q_latent",
    "q_sensible",
    "q_longwave_up",
    "q_longwave_down",
    "q_solar",
    "q_sediment",
    "equilibrium_temperature",
):
    registry.register(diagnostic_name, Variable(name=diagnostic_name, data=zeros))

# Then run the model normally; the diagnostics will be filled in
# at every substep.
model.run()
```

If a diagnostic is **not** pre-registered, the ``Temperature`` process
silently skips writing it (no error, no overhead). This keeps the cost
of an unused diagnostic at zero — useful when you want
``q_latent`` only for a specific calibration run but the rest of the
production runs do not need it.

## Reading from sibling processes

The seven diagnostics are also cached as instance attributes on the
``Temperature`` process after each ``run`` call:

```python
temperature_process = model.get_process("Temperature")
# After model.run() executes the Temperature substep:
print(temperature_process.q_net)         # most recent net flux per cell
print(temperature_process.q_latent)      # most recent evaporative magnitude
print(temperature_process.q_sensible)    # most recent sensible flux
# ... etc
```

This is the same caching pattern that DOX uses to read
``cbod_oxidation_rate`` from CBOD. You do not need to register the
diagnostic in the registry to read the cached value from another
process.

## Equilibrium temperature

``Temperature.equilibrium_temperature(...)`` solves for the water
temperature ``T_eq`` at which the surface heat budget vanishes
under the current met conditions:

```text
q_net(T_eq | atmospheric forcing, sediment forcing) = 0
```

The implementation uses a vectorized Newton-Raphson loop with up to
10 iterations (configurable) and a 0.01 K per-iteration tolerance
(configurable), matching Fortran-A
``modTemperature.f90:209-263``. Analytic derivatives are computed for
upwelling longwave (`4εσT³`), latent (via the Brutsaert
`de_sat / dT_K` polynomial), sensible, and sediment fluxes.

**Cost:** computed only when ``equilibrium_temperature`` is
pre-registered in the registry. Off the hot path otherwise. Each
iteration is one full flux evaluation plus four analytic derivatives;
typical convergence is in 3–6 iterations.

```python
# To get T_eq alongside the standard model state:
registry.register(
    "equilibrium_temperature",
    Variable(name="equilibrium_temperature", data=xr.zeros_like(shape_template)),
)
model.run()

# After each substep:
teq_array = registry.get_at_time("equilibrium_temperature", current_time)
```

For a single-shot equilibrium-temperature evaluation outside the model
substep loop, call the method directly:

```python
teq = temperature_process.equilibrium_temperature(
    cloudiness=0.3,
    air_temperature=20.0,
    solar_flux=400.0,
    wind_speed=3.0,
    atmospheric_pressure=1013.0,
    atmospheric_vapor_pressure=15.0,
    sediment_temperature=20.0,
    sediment_thickness=0.1,
    max_iterations=20,    # tighter convergence
    tolerance_kelvin=1e-4,
)
```

## Sign-convention reference

The audit on 2026-05-05 (see
``design/clearwater_modules_v3_tsm_audit_2026-05-05.md`` finding
F-sign-convention) refactored the v3 flux methods to the
"magnitudes-only, signs-at-composition" convention. **Do not
pre-negate inside any flux method**:

```python
# Correct (v3 1.0+):
q_net = (
    sensible
    + solar
    + sediment
    + atmospheric_lw          # downwelling LW, magnitude
    - upwelling_lw            # upwelling LW, magnitude
    - latent                  # evaporative, magnitude
)
```

The seven per-component diagnostics are written in this convention
(`q_latent`, `q_longwave_up` are magnitudes; the rest are signed by
their physical gradients). When summing for analysis, apply the same
signs as the composition above to recover `q_net`.

## See also

- `tests/v3/test_tsm_diagnostics_v3.py` — 9 regression tests covering
  the diagnostic contract.
- `tests/v3/test_tsm_sign_convention_v3.py` — 7 tests pinning the
  magnitudes-only sign convention.
- `design/clearwater_modules_v3_tsm_audit_2026-05-05.md` — full audit
  documenting the design decisions behind these features (open
  questions 2 and 3, finding F-sign-convention).
