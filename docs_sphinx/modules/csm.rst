Contaminant Simulation Module (CSM)
===================================

Overview
--------

The Contaminant Simulation Module (CSM) simulates the fate and transport of chemical contaminants in aquatic systems. It includes advanced processes for organic chemicals, metals, and emerging contaminants, accounting for complex environmental interactions.

Key Features
------------

- Multi-phase partitioning (water, particles, DOC, sediment)
- Transformation processes (hydrolysis, photolysis, biodegradation)
- Bioaccumulation modeling
- Volatilization
- Advanced sorption kinetics
- Metal speciation (planned enhancement)

Supported Contaminants
----------------------

The CSM can simulate:

- **Organic Chemicals**
  
  - Pesticides
  - PCBs
  - PAHs
  - Pharmaceuticals
  - Industrial chemicals

- **Metals** (basic implementation)
  
  - Heavy metals
  - Trace elements

- **Emerging Contaminants**
  
  - PFAS compounds
  - Microplastics (simplified)

Module Architecture
-------------------

Core Components
~~~~~~~~~~~~~~~

The CSM consists of:

- ``csm.py``: Main module interface
- ``csm_functions.py``: Process calculations
- ``csm_globals.py``: Constants and parameters

Process Framework
~~~~~~~~~~~~~~~~~

1. **Partitioning**
   
   - Equilibrium partitioning
   - Kinetic sorption/desorption
   - DOC binding
   - Multi-phase distribution

2. **Transformation**
   
   - Abiotic degradation
   - Biodegradation
   - Photochemical reactions
   - Hydrolysis

3. **Transport**
   
   - Advection/dispersion
   - Settling with particles
   - Volatilization
   - Sediment-water exchange

Usage Example
-------------

Basic contaminant simulation:

.. code-block:: python

   from clearwater_modules.csm import csm
   
   # Define contaminant properties
   contaminant = {
       'name': 'Atrazine',
       'molecular_weight': 215.68,  # g/mol
       'log_kow': 2.7,
       'henrys_constant': 2.48e-9,  # atm-m³/mol
       'koc': 100,  # L/kg
       'biodegradation_rate': 0.01,  # 1/day
       'photolysis_rate': 0.005,  # 1/day
       'hydrolysis_rate': 0.001  # 1/day
   }
   
   # Initialize CSM
   model = csm.CSM(contaminant=contaminant)
   
   # Set environmental conditions
   model.set_conditions(
       temperature=20,  # °C
       pH=7.5,
       tss=50,  # mg/L
       doc=5,  # mg/L
       solar_radiation=300  # W/m²
   )
   
   # Run simulation
   results = model.simulate(
       initial_concentration=1.0,  # μg/L
       timestep=3600,  # 1 hour
       duration=86400 * 30  # 30 days
   )

Advanced Features
-----------------

Multi-Phase Partitioning
~~~~~~~~~~~~~~~~~~~~~~~~

The CSM calculates contaminant distribution across phases:

.. code-block:: python

   # Phase distribution calculation
   phases = model.calculate_phase_distribution(
       total_concentration=10.0,  # μg/L
       tss=100,  # mg/L
       doc=5,  # mg/L
       poc_fraction=0.1  # fraction of TSS that is POC
   )
   
   print(f"Dissolved: {phases['dissolved']} μg/L")
   print(f"Particle-bound: {phases['particulate']} μg/L")
   print(f"DOC-bound: {phases['doc_bound']} μg/L")

Transformation Processes
~~~~~~~~~~~~~~~~~~~~~~~~

Configure multiple degradation pathways:

.. code-block:: python

   transformation_config = {
       'biodegradation': {
           'rate': 0.05,  # 1/day at 20°C
           'temperature_coefficient': 1.047,
           'half_saturation': 0.5  # mg/L DOC
       },
       'photolysis': {
           'quantum_yield': 0.001,
           'absorption_coefficient': 0.01,
           'light_screening': True
       },
       'hydrolysis': {
           'neutral_rate': 0.001,
           'acid_catalyzed': 0.01,
           'base_catalyzed': 0.1,
           'activation_energy': 50  # kJ/mol
       }
   }

Bioaccumulation
~~~~~~~~~~~~~~~

Model contaminant uptake by organisms:

.. code-block:: python

   # Bioaccumulation factor calculation
   baf = model.calculate_bioaccumulation(
       organism_lipid_fraction=0.05,
       metabolism_rate=0.01,
       growth_dilution=0.005
   )

Input Parameters
----------------

Contaminant Properties
~~~~~~~~~~~~~~~~~~~~~~

Required:

- Molecular weight
- Log Kₒw or Kₒc
- Henry's constant
- Degradation rates

Optional:

- pKa (for ionizable compounds)
- Diffusion coefficients
- Quantum yield
- Absorption spectra

Environmental Conditions
~~~~~~~~~~~~~~~~~~~~~~~~

- Temperature
- pH
- Dissolved oxygen
- Total suspended solids
- Dissolved organic carbon
- Solar radiation
- Water depth
- Wind speed

Output Variables
----------------

The CSM provides comprehensive output:

Concentrations
~~~~~~~~~~~~~~

- Total concentration
- Dissolved phase
- Particulate phase
- DOC-associated
- Sediment concentration

Process Rates
~~~~~~~~~~~~~

- Degradation rates
- Volatilization flux
- Settling flux
- Sediment-water exchange

Mass Balance
~~~~~~~~~~~~

- Input loads
- Transformation losses
- Transport fluxes
- Accumulation

Advanced Modeling Options
-------------------------

Stochastic Modeling
~~~~~~~~~~~~~~~~~~~

Account for parameter uncertainty:

.. code-block:: python

   # Monte Carlo simulation
   results = model.monte_carlo_simulation(
       parameter_distributions={
           'log_kow': ('normal', 2.7, 0.2),
           'biodeg_rate': ('lognormal', 0.01, 0.5)
       },
       n_iterations=1000
   )

Spatial Modeling
~~~~~~~~~~~~~~~~

Couple with transport models:

.. code-block:: python

   # 2D spatial simulation
   spatial_model = csm.SpatialCSM(
       contaminant=contaminant,
       grid=spatial_grid,
       transport_model=transport
   )

Time-Variable Inputs
~~~~~~~~~~~~~~~~~~~~

Handle dynamic conditions:

.. code-block:: python

   # Time-series inputs
   model.set_time_series(
       solar_radiation=hourly_radiation,
       temperature=hourly_temperature,
       flow_rate=daily_flow
   )

Best Practices
--------------

1. **Parameter Selection**
   
   - Use measured values when available
   - Apply QSAR estimates carefully
   - Consider site-specific calibration

2. **Model Validation**
   
   - Compare with field measurements
   - Check mass balance closure
   - Sensitivity analysis for key parameters

3. **Uncertainty Analysis**
   
   - Identify sensitive parameters
   - Propagate uncertainty through model
   - Report confidence intervals

Integration Example
-------------------

Coupling CSM with hydrodynamic model:

.. code-block:: python

   # Example integration
   from clearwater_modules.csm import csm
   from hypothetical import HydrodynamicModel
   
   # Initialize models
   hydro = HydrodynamicModel()
   contaminant_model = csm.CSM(contaminant_properties)
   
   # Simulation loop
   for time in simulation_period:
       # Get flow field
       velocity = hydro.get_velocity(time)
       mixing = hydro.get_dispersion(time)
       
       # Update CSM transport
       contaminant_model.set_transport(velocity, mixing)
       
       # Calculate contaminant fate
       contaminant_model.step(timestep)
       
       # Feedback to hydro (if needed)
       density_change = contaminant_model.get_density_effect()
       hydro.update_density(density_change)

References
----------

The CSM algorithms are based on:

- WASP contaminant modules
- EFDC toxics routines
- EPA's EXAMS model
- Current literature on contaminant fate

Key technical report:

- Zhang, Z. and Johnson, B.E. (2016). Aquatic contaminant and mercury simulation modules developed for hydrologic and hydraulic models. ERDC/EL TR-16-8.