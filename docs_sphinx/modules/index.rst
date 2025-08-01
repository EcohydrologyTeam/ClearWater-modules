Water Quality Modules
=====================

This section provides detailed documentation for each water quality module in the ClearWater system.

.. toctree::
   :maxdepth: 2
   :caption: Available Modules:

   tsm
   nsm
   gsm
   csm
   msm

Module Overview
---------------

Each water quality module in ClearWater is designed to simulate specific aspects of aquatic systems:

**Temperature Simulation Module (TSM)**
   Simulates water temperature dynamics including heat exchange processes.

**Nutrient Simulation Modules (NSM)**
   Models nutrient cycling (nitrogen, phosphorus, carbon) in aquatic systems.

**General Constituent Simulation Module (GSM)**
   Handles general water quality constituents.

**Contaminant Simulation Module (CSM)**
   Simulates contaminant fate and transport.

**Mercury Simulation Module (MSM)**
   Specialized module for mercury cycling and bioaccumulation.

Module Structure
----------------

Each module follows a consistent structure:

- **Constants**: Physical and chemical constants used in calculations
- **State Variables**: Variables that change over time
- **Static Variables**: Parameters that remain constant during simulation
- **Dynamic Variables**: Intermediate variables calculated during simulation
- **Processes**: Functions that implement the scientific algorithms
- **Model**: Main interface for running the module

Integration
-----------

The modules are designed to work independently or in combination, allowing users to select only the components needed for their specific application.