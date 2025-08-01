Temperature Simulation Module (TSM)
===================================

Overview
--------

The Temperature Simulation Module (TSM), formerly known as TEMP, simulates water temperature dynamics in aquatic systems. It accounts for various heat exchange processes including solar radiation, atmospheric heat exchange, and bed heat conduction.

Key Features
------------

- Comprehensive heat budget calculations
- Surface heat exchange modeling
- Bed heat conduction
- Ice formation and melting processes
- Integration with hydrodynamic models

Module Components
-----------------

Constants
~~~~~~~~~

The TSM uses various physical constants defined in ``tsm.constants``:

- Stefan-Boltzmann constant
- Specific heat of water
- Latent heat of vaporization
- Thermal conductivity values

State Variables
~~~~~~~~~~~~~~~

Primary state variables tracked by TSM:

- **TwaterC**: Water temperature (°C)
- **surface_area**: Water surface area (m²)
- **volume**: Water volume (m³)

Processes
~~~~~~~~~

The module implements several key processes:

1. **Solar Radiation**
   
   - Calculation of incident solar radiation
   - Absorption through the water column
   - Reflection at the water surface

2. **Atmospheric Heat Exchange**
   
   - Long-wave radiation
   - Sensible heat transfer
   - Latent heat transfer (evaporation/condensation)

3. **Bed Heat Conduction**
   
   - Heat exchange with sediments
   - Temperature gradient calculations

Usage Example
-------------

Basic usage of the TSM module:

.. code-block:: python

   from clearwater_modules.tsm import model
   from clearwater_modules.tsm import state_variables, static_variables
   
   # Initialize state variables
   state_vars = state_variables.StateVariables(
       water_temp_c=20.0,
       surface_area=1000.0,
       volume=10000.0
   )
   
   # Set static parameters
   static_vars = static_variables.StaticVariables(
       latitude=40.0,
       longitude=-74.0,
       elevation=100.0
   )
   
   # Run the model
   model = model.TSM()
   results = model.run(
       state_vars=state_vars,
       static_vars=static_vars,
       timestep=3600  # 1 hour
   )

Input Requirements
------------------

Meteorological Data
~~~~~~~~~~~~~~~~~~~

- Air temperature
- Wind speed
- Relative humidity
- Cloud cover
- Solar radiation (optional)

Physical Parameters
~~~~~~~~~~~~~~~~~~~

- Water depth
- Surface area
- Volume
- Sediment properties

Output Variables
----------------

The TSM provides the following outputs:

- Water temperature time series
- Heat flux components
- Energy budget terms
- Ice formation indicators

References
----------

The TSM algorithms are based on established heat transfer principles and have been validated against field data from various aquatic systems. Key references include:

- Heat exchange formulations from CE-QUAL-W2
- Solar radiation algorithms from EPA models
- Bed heat conduction methods from scientific literature