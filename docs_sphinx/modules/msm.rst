Mercury Simulation Module (MSM)
===============================

Overview
--------

The Mercury Simulation Module (MSM) is a specialized module for simulating mercury cycling in aquatic environments. It addresses the unique chemistry of mercury, including speciation, methylation, and bioaccumulation processes that make mercury a particularly concerning environmental contaminant.

Key Features
------------

- Mercury speciation (Hg⁰, Hg²⁺, MeHg)
- Methylation and demethylation processes
- Multi-phase partitioning
- Bioaccumulation and biomagnification
- Atmospheric exchange
- Sediment-water interactions
- Redox transformations

Mercury Species Simulated
-------------------------

The MSM tracks three primary mercury species:

1. **Elemental Mercury (Hg⁰)**
   
   - Volatile form
   - Low solubility
   - Subject to oxidation/reduction

2. **Inorganic Mercury (Hg²⁺)**
   
   - Most common dissolved form
   - Strong binding to particles and DOC
   - Precursor to methylmercury

3. **Methylmercury (MeHg)**
   
   - Highly toxic form
   - Bioaccumulates in food chains
   - Produced by microbial methylation

Module Components
-----------------

Process Framework
~~~~~~~~~~~~~~~~~

1. **Speciation and Transformation**
   
   - Oxidation: Hg⁰ → Hg²⁺
   - Reduction: Hg²⁺ → Hg⁰
   - Methylation: Hg²⁺ → MeHg
   - Demethylation: MeHg → Hg²⁺

2. **Partitioning**
   
   - Dissolved phase
   - Particle-bound
   - DOC complexation
   - Sediment sorption

3. **Transport Processes**
   
   - Volatilization of Hg⁰
   - Settling with particles
   - Sediment-water diffusion
   - Resuspension

4. **Bioaccumulation**
   
   - Uptake by phytoplankton
   - Transfer through food web
   - Biomagnification factors

Usage Example
-------------

Basic mercury simulation:

.. code-block:: python

   from clearwater_modules.msm import msm
   
   # Initialize MSM
   mercury_model = msm.MSM()
   
   # Set initial conditions
   initial_conditions = {
       'Hg0': 0.1,      # ng/L
       'Hg2': 5.0,      # ng/L
       'MeHg': 0.5,     # ng/L
       'Hg_particulate': 2.0,  # ng/L
       'MeHg_particulate': 0.2  # ng/L
   }
   
   # Environmental parameters
   environment = {
       'temperature': 20,    # °C
       'pH': 7.0,
       'dissolved_oxygen': 8.0,  # mg/L
       'tss': 10,           # mg/L
       'doc': 5,            # mg/L
       'sulfide': 0.001,    # mg/L
       'light_intensity': 300  # W/m²
   }
   
   # Run simulation
   results = mercury_model.simulate(
       initial_conditions=initial_conditions,
       environment=environment,
       timestep=3600,      # 1 hour
       duration=86400 * 365  # 1 year
   )

Advanced Features
-----------------

Methylation Modeling
~~~~~~~~~~~~~~~~~~~~

The MSM includes detailed methylation kinetics:

.. code-block:: python

   # Configure methylation parameters
   methylation_params = {
       'base_rate': 0.01,          # 1/day
       'sulfate_half_sat': 10,     # mg/L
       'optimal_sulfide': 0.01,    # mg/L
       'temperature_coeff': 1.08,
       'pH_optimal': 6.5,
       'microbial_activity': 'high'
   }
   
   mercury_model.set_methylation_parameters(methylation_params)

Bioaccumulation Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Model mercury accumulation through food webs:

.. code-block:: python

   # Define food web structure
   food_web = {
       'phytoplankton': {
           'trophic_level': 1,
           'uptake_rate': 0.1,
           'elimination_rate': 0.05
       },
       'zooplankton': {
           'trophic_level': 2,
           'diet': {'phytoplankton': 1.0},
           'assimilation_efficiency': 0.7
       },
       'small_fish': {
           'trophic_level': 3,
           'diet': {'zooplankton': 0.8, 'phytoplankton': 0.2},
           'growth_rate': 0.01
       },
       'predator_fish': {
           'trophic_level': 4,
           'diet': {'small_fish': 1.0},
           'elimination_rate': 0.001
       }
   }
   
   # Calculate bioaccumulation factors
   baf = mercury_model.calculate_food_web_accumulation(food_web)

Sediment Interactions
~~~~~~~~~~~~~~~~~~~~~

Detailed sediment-water exchange:

.. code-block:: python

   # Sediment mercury parameters
   sediment_config = {
       'bulk_density': 1.5,        # g/cm³
       'porosity': 0.8,
       'organic_content': 0.05,    # fraction
       'mixing_depth': 0.1,        # m
       'bioturbation_rate': 0.01,  # 1/day
       'methylation_depth': 0.05   # m
   }
   
   mercury_model.configure_sediment(sediment_config)

Input Requirements
------------------

Mercury Sources
~~~~~~~~~~~~~~~

- Atmospheric deposition (wet and dry)
- Point source discharges
- Watershed runoff
- Sediment release

Environmental Conditions
~~~~~~~~~~~~~~~~~~~~~~~~

Critical parameters affecting mercury cycling:

- Temperature
- pH
- Dissolved oxygen
- Sulfate/sulfide
- Organic carbon (DOC and POC)
- Redox potential
- Microbial activity

Biological Data
~~~~~~~~~~~~~~~

For bioaccumulation modeling:

- Food web structure
- Growth rates
- Feeding rates
- Lipid content

Output Variables
----------------

Concentrations
~~~~~~~~~~~~~~

- Mercury species in water column
- Particulate mercury
- Sediment mercury
- Biota mercury concentrations

Process Rates
~~~~~~~~~~~~~

- Methylation/demethylation rates
- Volatilization flux
- Settling rates
- Bioaccumulation factors

Mass Balance
~~~~~~~~~~~~

- Total mercury budget
- Species transformations
- Atmospheric exchange
- Sediment accumulation

Special Considerations
----------------------

Sulfur Cycling
~~~~~~~~~~~~~~

Mercury methylation is strongly linked to sulfur cycling:

.. code-block:: python

   # Link to sulfur cycle
   mercury_model.couple_sulfur_cycle(
       sulfate_reduction_rate=0.1,
       sulfide_oxidation_rate=0.5
   )

Seasonal Variations
~~~~~~~~~~~~~~~~~~~

Account for seasonal changes:

.. code-block:: python

   # Seasonal methylation patterns
   seasonal_factors = {
       'spring': 1.5,
       'summer': 2.0,
       'fall': 1.2,
       'winter': 0.5
   }
   
   mercury_model.set_seasonal_methylation(seasonal_factors)

Model Validation
----------------

The MSM should be validated against:

- Water column mercury measurements
- Sediment mercury profiles
- Fish tissue concentrations
- Methylation rate measurements

Best Practices
--------------

1. **Initial Conditions**
   
   - Use site-specific data when available
   - Consider historical contamination
   - Account for background levels

2. **Parameter Selection**
   
   - Methylation rates vary widely
   - Site-specific calibration essential
   - Consider seasonal variations

3. **Uncertainty Analysis**
   
   - Mercury cycling highly uncertain
   - Perform sensitivity analysis
   - Use ensemble runs

Integration with Other Modules
------------------------------

The MSM can be coupled with other ClearWater modules:

.. code-block:: python

   # Example coupling with NSM for nutrient effects
   from clearwater_modules.nsm2 import nsm2
   from clearwater_modules.msm import msm
   
   # Initialize modules
   nutrient_model = nsm2.NSM2()
   mercury_model = msm.MSM()
   
   # Link primary productivity to methylation
   for time in simulation:
       # Get algae biomass from NSM
       algae = nutrient_model.get_algae_biomass()
       
       # Adjust methylation based on productivity
       mercury_model.set_productivity_factor(algae)
       
       # Run mercury calculations
       mercury_model.step(timestep)

Limitations
-----------

- Simplified food web representation
- Methylation processes are empirical
- Limited metal-ligand complexation
- No explicit microbial community dynamics

Future Enhancements
-------------------

Planned improvements include:

- Explicit sulfate-reducing bacteria modeling
- Enhanced metal speciation
- Stable isotope tracking
- Coupled contaminant interactions

References
----------

The MSM is based on:

- EPA's MCM (Mercury Cycling Model)
- WASP mercury module
- SERAFM mercury components
- Recent mercury cycling research

Key reference:

- Zhang, Z. and Johnson, B.E. (2016). Aquatic contaminant and mercury simulation modules developed for hydrologic and hydraulic models. ERDC/EL TR-16-8.