Examples
========

This section provides practical examples of using ClearWater modules for various water quality modeling scenarios.

Getting Started
---------------

Basic Module Usage
~~~~~~~~~~~~~~~~~~

Here's a simple example of running the Temperature Simulation Module:

.. code-block:: python

   from clearwater_modules.tsm import model, state_variables, static_variables
   import numpy as np
   
   # Initialize state variables
   state = state_variables.StateVariables()
   state.water_temp_c = 20.0  # Initial water temperature
   state.surface_area = 10000.0  # m²
   state.volume = 50000.0  # m³
   
   # Set static parameters
   static = static_variables.StaticVariables()
   static.latitude = 40.0
   static.longitude = -74.0
   
   # Create and run model
   tsm = model.TSM()
   results = tsm.run(state, static, timestep=3600, duration=86400)

Complete Workflows
------------------

Multi-Module Integration
~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to run TSM and NSM together:

.. code-block:: python

   from clearwater_modules.tsm import model as tsm_model
   from clearwater_modules.nsm1 import model as nsm_model
   import pandas as pd
   
   # Initialize both modules
   tsm = tsm_model.TSM()
   nsm = nsm_model.NSM1()
   
   # Set up shared environmental conditions
   timestep = 3600  # 1 hour
   duration = 86400 * 30  # 30 days
   
   # Time loop
   results = []
   for t in range(0, duration, timestep):
       # Run temperature module
       temp_results = tsm.step(timestep)
       
       # Use temperature in nutrient module
       nsm.set_temperature(temp_results['water_temp'])
       
       # Run nutrient module
       nutrient_results = nsm.step(timestep)
       
       # Store results
       results.append({
           'time': t,
           'temperature': temp_results['water_temp'],
           'ammonia': nutrient_results['NH4'],
           'nitrate': nutrient_results['NO3'],
           'phosphate': nutrient_results['PO4']
       })
   
   # Convert to DataFrame for analysis
   df = pd.DataFrame(results)

Spatial Modeling
~~~~~~~~~~~~~~~~

Example of applying modules to multiple grid cells:

.. code-block:: python

   import numpy as np
   from clearwater_modules.nsm2 import nsm2
   
   # Define spatial grid
   nx, ny = 10, 10  # 10x10 grid
   
   # Initialize NSM for each cell
   models = [[nsm2.NSM2() for j in range(ny)] for i in range(nx)]
   
   # Set spatially varying parameters
   for i in range(nx):
       for j in range(ny):
           # Depth decreases from left to right
           depth = 10.0 - i * 0.5
           
           # Temperature varies with location
           temp = 20.0 + i * 0.2 + j * 0.1
           
           models[i][j].set_depth(depth)
           models[i][j].set_temperature(temp)
   
   # Run simulation with transport
   for t in range(simulation_steps):
       # Calculate transport between cells
       # (simplified - actual transport would use hydrodynamic model)
       
       # Update each cell
       for i in range(nx):
           for j in range(ny):
               models[i][j].step(timestep)

Real-World Applications
-----------------------

River Water Quality
~~~~~~~~~~~~~~~~~~~

Modeling a river reach with point source inputs:

.. code-block:: python

   from clearwater_modules.nsm2 import nsm2
   from clearwater_modules.tsm import model as tsm_model
   
   # River segments
   n_segments = 20
   segment_length = 1000  # meters
   
   # Initialize modules for each segment
   segments = []
   for i in range(n_segments):
       seg = {
           'nsm': nsm2.NSM2(),
           'tsm': tsm_model.TSM(),
           'volume': segment_length * 100 * 3,  # length * width * depth
           'flow': 10.0  # m³/s
       }
       segments.append(seg)
   
   # Add point source at segment 5
   point_source = {
       'segment': 5,
       'flow': 0.5,  # m³/s
       'NH4': 10.0,  # mg/L
       'PO4': 2.0,   # mg/L
       'temperature': 25.0  # °C
   }
   
   # Simulation loop
   for t in range(simulation_steps):
       # Transport between segments
       for i in range(n_segments-1):
           # Simple advective transport
           travel_time = segment_length / segments[i]['flow']
           
           # Pass downstream
           segments[i+1]['nsm'].set_upstream(
               segments[i]['nsm'].get_concentrations()
           )
       
       # Add point source
       ps = point_source
       segments[ps['segment']]['nsm'].add_load(
           flow=ps['flow'],
           concentrations={'NH4': ps['NH4'], 'PO4': ps['PO4']}
       )
       
       # Run each segment
       for seg in segments:
           seg['tsm'].step(timestep)
           seg['nsm'].step(timestep)

Lake Eutrophication
~~~~~~~~~~~~~~~~~~~

Seasonal lake simulation with stratification:

.. code-block:: python

   from clearwater_modules.nsm2 import nsm2
   from clearwater_modules.tsm import model as tsm_model
   import numpy as np
   
   # Lake layers (epilimnion, metalimnion, hypolimnion)
   layers = {
       'epilimnion': {'depth': 5, 'volume': 1e6},
       'metalimnion': {'depth': 3, 'volume': 0.5e6},
       'hypolimnion': {'depth': 12, 'volume': 2e6}
   }
   
   # Initialize modules for each layer
   for name, layer in layers.items():
       layer['tsm'] = tsm_model.TSM()
       layer['nsm'] = nsm2.NSM2()
       layer['nsm'].enable_benthic_processes = (name == 'hypolimnion')
   
   # Seasonal simulation
   for day in range(365):
       # Calculate stratification strength
       if 120 < day < 270:  # Summer stratification
           mixing_rate = 0.01
       else:  # Mixed conditions
           mixing_rate = 1.0
       
       # Exchange between layers
       for i, (name1, layer1) in enumerate(layers.items()):
           for j, (name2, layer2) in enumerate(layers.items()):
               if i < j:  # Only mix adjacent layers
                   exchange_nutrients(layer1, layer2, mixing_rate)
       
       # Run models
       for layer in layers.values():
           layer['tsm'].step(86400)  # Daily timestep
           layer['nsm'].step(86400)

Advanced Features
-----------------

Uncertainty Analysis
~~~~~~~~~~~~~~~~~~~~

Monte Carlo simulation for parameter uncertainty:

.. code-block:: python

   import numpy as np
   from clearwater_modules.nsm1 import model
   import matplotlib.pyplot as plt
   
   # Define parameter distributions
   param_distributions = {
       'algae_growth_rate': ('normal', 2.0, 0.3),
       'algae_death_rate': ('lognormal', 0.1, 0.5),
       'nitrification_rate': ('uniform', 0.05, 0.15)
   }
   
   # Monte Carlo runs
   n_runs = 100
   results = []
   
   for i in range(n_runs):
       # Sample parameters
       params = {}
       for param, (dist, *args) in param_distributions.items():
           if dist == 'normal':
               params[param] = np.random.normal(*args)
           elif dist == 'lognormal':
               params[param] = np.random.lognormal(np.log(args[0]), args[1])
           elif dist == 'uniform':
               params[param] = np.random.uniform(*args)
       
       # Run model
       nsm = model.NSM1()
       nsm.set_parameters(params)
       output = nsm.run(duration=86400*30)
       
       results.append(output)
   
   # Analyze results
   algae_biomass = [r['algae'][-1] for r in results]
   plt.hist(algae_biomass, bins=20)
   plt.xlabel('Final Algae Biomass (mg/L)')
   plt.ylabel('Frequency')
   plt.title('Uncertainty in Algae Predictions')

Coupling with Hydrodynamic Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Example interface for HEC-RAS coupling:

.. code-block:: python

   from clearwater_modules.nsm2 import nsm2
   from clearwater_modules.tsm import model as tsm_model
   
   class HECRASCoupler:
       def __init__(self, ras_model, wq_modules):
           self.ras_model = ras_model
           self.wq_modules = wq_modules
           self.cells = self._initialize_cells()
       
       def _initialize_cells(self):
           cells = {}
           for cell_id in self.ras_model.get_cell_ids():
               cells[cell_id] = {
                   'tsm': tsm_model.TSM(),
                   'nsm': nsm2.NSM2()
               }
           return cells
       
       def step(self, timestep):
           # Get hydraulic conditions from RAS
           for cell_id, modules in self.cells.items():
               hydraulics = self.ras_model.get_cell_hydraulics(cell_id)
               
               # Update water quality modules
               modules['tsm'].set_depth(hydraulics['depth'])
               modules['tsm'].set_velocity(hydraulics['velocity'])
               modules['nsm'].set_depth(hydraulics['depth'])
               
               # Get transport fluxes
               for neighbor_id in self.ras_model.get_neighbors(cell_id):
                   flux = self.ras_model.get_flux(cell_id, neighbor_id)
                   self._transport_constituents(cell_id, neighbor_id, flux)
               
               # Run water quality calculations
               modules['tsm'].step(timestep)
               modules['nsm'].step(timestep)
       
       def _transport_constituents(self, from_cell, to_cell, flux):
           # Transport logic here
           pass

Best Practices
--------------

1. **Start Simple**: Begin with single modules before attempting integration
2. **Validate Often**: Compare results with analytical solutions or field data
3. **Monitor Mass Balance**: Ensure conservation of mass in your simulations
4. **Document Parameters**: Keep track of parameter sources and assumptions
5. **Use Version Control**: Track changes to your model configurations

Additional Resources
--------------------

- Jupyter notebooks in the `examples/` directory
- Test cases in the `tests/` directory
- Module-specific documentation in the API reference