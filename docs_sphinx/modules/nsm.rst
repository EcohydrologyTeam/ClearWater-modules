Nutrient Simulation Modules (NSM)
=================================

Overview
--------

The Nutrient Simulation Modules (NSM) model the cycling of nutrients in aquatic systems. Two versions are available:

- **NSM-I**: A simplified nutrient model focusing on essential processes
- **NSM-II**: An advanced model with comprehensive nutrient cycling capabilities

Both modules simulate nitrogen, phosphorus, and carbon dynamics along with associated biological and chemical processes.

NSM-I: Simplified Nutrient Model
--------------------------------

Key Features
~~~~~~~~~~~~

- Basic nitrogen cycle (ammonia, nitrate, organic nitrogen)
- Phosphorus cycling (orthophosphate, organic phosphorus)
- Phytoplankton growth and decay
- Dissolved oxygen dynamics
- CBOD (Carbonaceous Biochemical Oxygen Demand)

Primary Processes
~~~~~~~~~~~~~~~~~

1. **Nitrogen Cycle**
   
   - Nitrification (NH₄ → NO₃)
   - Denitrification
   - Organic nitrogen mineralization
   - Algal uptake

2. **Phosphorus Cycle**
   
   - Organic phosphorus mineralization
   - Algal uptake
   - Sorption/desorption

3. **Phytoplankton Dynamics**
   
   - Growth limited by nutrients, light, and temperature
   - Respiration and mortality
   - Settling

NSM-II: Advanced Nutrient Model
--------------------------------

Enhanced Features
~~~~~~~~~~~~~~~~~

NSM-II includes all NSM-I capabilities plus:

- Multiple algae groups
- Benthic algae
- Sediment-water interactions
- pH and alkalinity dynamics
- Pathogen modeling
- Detailed organic matter pools

Additional Processes
~~~~~~~~~~~~~~~~~~~~

1. **Benthic Processes**
   
   - Sediment oxygen demand
   - Nutrient flux from sediments
   - Benthic algae growth

2. **Carbon Dynamics**
   
   - Dissolved inorganic carbon
   - pH calculations
   - Alkalinity

3. **Advanced Biological Processes**
   
   - Multiple phytoplankton groups
   - Zooplankton grazing
   - Macrophyte interactions

Module Structure
----------------

Both NSM versions follow a similar structure:

Constants
~~~~~~~~~

Physical and biological rate constants:

- Growth rates
- Decay rates
- Half-saturation constants
- Temperature coefficients

State Variables
~~~~~~~~~~~~~~~

NSM-I core variables:

- Ammonia nitrogen (NH₄-N)
- Nitrate nitrogen (NO₃-N)
- Organic nitrogen
- Orthophosphate (PO₄-P)
- Organic phosphorus
- Phytoplankton biomass
- Dissolved oxygen
- CBOD

NSM-II additional variables:

- Multiple algae groups
- Benthic algae
- Alkalinity
- pH
- Pathogens
- Particulate organic matter

Usage Example
-------------

NSM-I Example:

.. code-block:: python

   from clearwater_modules.nsm1 import model
   import numpy as np
   
   # Initialize the model
   nsm = model.NSM1()
   
   # Set initial conditions
   initial_state = {
       'NH4': 0.5,  # mg/L
       'NO3': 2.0,  # mg/L
       'OrgN': 1.0,  # mg/L
       'PO4': 0.1,  # mg/L
       'OrgP': 0.05,  # mg/L
       'Algae': 10.0,  # mg/L
       'DO': 8.0,  # mg/L
       'CBOD': 5.0  # mg/L
   }
   
   # Run simulation
   results = nsm.run(
       initial_state=initial_state,
       timestep=3600,  # 1 hour
       duration=86400 * 30  # 30 days
   )

NSM-II Example:

.. code-block:: python

   from clearwater_modules.nsm2 import nsm2
   
   # Initialize with multiple algae groups
   model = nsm2.NSM2(
       algae_groups=['diatoms', 'green_algae', 'cyanobacteria']
   )
   
   # Configure benthic processes
   model.enable_benthic_processes = True
   model.enable_pH_calculations = True
   
   # Run simulation with advanced features
   results = model.simulate(
       initial_conditions=initial_conditions,
       boundary_conditions=boundary_conditions,
       meteorological_data=met_data
   )

Input Requirements
------------------

Environmental Conditions
~~~~~~~~~~~~~~~~~~~~~~~~

- Water temperature
- Solar radiation
- Flow velocity
- Depth

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

- Nutrient loadings
- Upstream concentrations
- Point source inputs

Kinetic Parameters
~~~~~~~~~~~~~~~~~~

- Growth rates
- Decay rates
- Settling velocities
- Temperature coefficients

Output Variables
----------------

Standard Outputs
~~~~~~~~~~~~~~~~

- Nutrient concentrations over time
- Phytoplankton biomass
- Dissolved oxygen
- Organic matter

Advanced Outputs (NSM-II)
~~~~~~~~~~~~~~~~~~~~~~~~~

- pH and alkalinity
- Sediment fluxes
- Multiple algae group dynamics
- Pathogen concentrations

Model Selection Guide
---------------------

**Use NSM-I when:**

- Basic nutrient cycling is sufficient
- Limited data availability
- Faster computation is needed
- Teaching or demonstration purposes

**Use NSM-II when:**

- Detailed ecosystem modeling is required
- pH/alkalinity is important
- Benthic processes are significant
- Multiple algae groups need to be tracked
- Sediment-water interactions are critical

References
----------

The NSM modules are based on:

- QUAL2K model algorithms
- CE-QUAL-W2 kinetic formulations
- EPA water quality models
- Recent scientific literature on nutrient cycling

Key technical reports:

- Zhang, Z. and Johnson, B.E. (2016). Aquatic nutrient simulation modules (NSMs) developed for hydrologic and hydraulic models. ERDC/EL TR-16-1.