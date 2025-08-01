General Constituent Simulation Module (GSM)
===========================================

Overview
--------

The General Constituent Simulation Module (GSM) provides a flexible framework for simulating water quality constituents that don't require the specialized processes of the nutrient or contaminant modules. It's designed to handle conservative and non-conservative substances with first-order decay.

Key Features
------------

- Generic constituent tracking
- First-order decay kinetics
- Temperature-dependent reactions
- Settling and resuspension
- Simple sorption processes
- Flexible configuration for multiple constituents

Supported Constituents
----------------------

The GSM can simulate various general constituents including:

- Total Suspended Solids (TSS)
- Salinity/Total Dissolved Solids (TDS)
- Generic tracers
- Simple organic compounds
- Conservative substances
- User-defined constituents

Module Components
-----------------

Core Processes
~~~~~~~~~~~~~~

1. **Transport**
   
   - Advection with flow
   - Dispersion/diffusion
   - Settling (for particulates)

2. **Transformation**
   
   - First-order decay
   - Temperature correction
   - Photolysis (optional)

3. **Interactions**
   
   - Sorption to particles
   - Resuspension from bed
   - Exchange with sediments

Configuration
~~~~~~~~~~~~~

The GSM uses a configuration-based approach:

.. code-block:: python

   constituent_config = {
       'name': 'TSS',
       'units': 'mg/L',
       'decay_rate': 0.0,  # Conservative
       'settling_velocity': 0.1,  # m/day
       'temperature_coefficient': 1.0,
       'initial_concentration': 50.0
   }

Usage Example
-------------

Basic usage for simulating multiple constituents:

.. code-block:: python

   from clearwater_modules.gsm import general_constituents
   
   # Define constituents
   constituents = [
       {
           'name': 'TSS',
           'decay_rate': 0.0,
           'settling_velocity': 0.5,
           'initial_concentration': 100.0
       },
       {
           'name': 'Tracer',
           'decay_rate': 0.0,
           'settling_velocity': 0.0,
           'initial_concentration': 1.0
       },
       {
           'name': 'Generic_Organic',
           'decay_rate': 0.1,  # 1/day
           'temperature_coefficient': 1.047,
           'initial_concentration': 10.0
       }
   ]
   
   # Initialize GSM
   gsm = general_constituents.GSM(constituents=constituents)
   
   # Set environmental conditions
   gsm.set_temperature(20.0)  # °C
   gsm.set_flow_velocity(0.5)  # m/s
   
   # Run simulation
   results = gsm.simulate(
       timestep=3600,  # 1 hour
       duration=86400 * 7  # 7 days
   )

Advanced Features
-----------------

Sorption Modeling
~~~~~~~~~~~~~~~~~

For constituents that sorb to particles:

.. code-block:: python

   constituent_config = {
       'name': 'Sorbing_Chemical',
       'partition_coefficient': 100.0,  # L/kg
       'fraction_dissolved': 0.7,
       'sorption_kinetics': 'equilibrium'
   }

Temperature Effects
~~~~~~~~~~~~~~~~~~~

Temperature correction for reaction rates:

.. code-block:: python

   # Arrhenius temperature correction
   rate_T = rate_20 * theta ** (T - 20)
   
   # Where theta is the temperature coefficient

Settling and Resuspension
~~~~~~~~~~~~~~~~~~~~~~~~~

For particulate constituents:

.. code-block:: python

   settling_config = {
       'settling_velocity': 1.0,  # m/day
       'critical_shear_resuspension': 0.1,  # N/m²
       'resuspension_rate': 0.01  # kg/m²/day
   }

Input Requirements
------------------

Physical Parameters
~~~~~~~~~~~~~~~~~~~

- Water depth
- Flow velocity
- Temperature
- Turbulence/mixing coefficients

Constituent Properties
~~~~~~~~~~~~~~~~~~~~~~

- Initial concentrations
- Decay rates
- Settling velocities
- Partition coefficients (if applicable)

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

- Upstream concentrations
- Point source loads
- Distributed loads

Output Variables
----------------

The GSM provides:

- Concentration time series for each constituent
- Mass balance information
- Settling fluxes
- Transformation rates
- Spatial distributions (if coupled with transport model)

Best Practices
--------------

1. **Constituent Selection**
   
   - Use GSM for simple constituents
   - Consider NSM for nutrients
   - Use CSM for complex contaminants

2. **Parameter Estimation**
   
   - Literature values for common constituents
   - Field calibration for site-specific applications
   - Sensitivity analysis for uncertain parameters

3. **Model Validation**
   
   - Compare with field measurements
   - Check mass balance
   - Verify reaction kinetics

Integration with Transport Models
---------------------------------

The GSM is designed to couple with various transport models:

.. code-block:: python

   # Example coupling with hypothetical transport model
   transport_model = TransportModel()
   gsm = general_constituents.GSM(constituents)
   
   for timestep in simulation_period:
       # Get transport terms from hydrodynamic model
       velocity = transport_model.get_velocity()
       dispersion = transport_model.get_dispersion()
       
       # Update GSM with transport information
       gsm.update_transport(velocity, dispersion)
       
       # Calculate constituent changes
       gsm.step(timestep)
       
       # Update transport model with new concentrations
       transport_model.update_concentrations(gsm.get_concentrations())

Limitations
-----------

- First-order kinetics only (no complex reactions)
- No biological interactions
- Simple sorption models
- No pH or redox effects

For more complex processes, consider using the specialized modules (NSM, CSM, MSM).