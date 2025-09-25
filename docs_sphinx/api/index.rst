API Reference
=============

This section contains the complete API reference for all ClearWater modules.

.. toctree::
   :maxdepth: 2
   :caption: Module APIs:

   clearwater_modules
   tsm
   nsm1
   nsm2
   gsm
   csm
   msm
   shared
   base

Module Overview
---------------

The ClearWater API is organized into the following main components:

**Core Modules**
   - :mod:`clearwater_modules.tsm` - Temperature Simulation Module
   - :mod:`clearwater_modules.nsm1` - Nutrient Simulation Module I
   - :mod:`clearwater_modules.nsm2` - Nutrient Simulation Module II
   - :mod:`clearwater_modules.gsm` - General Constituent Simulation Module
   - :mod:`clearwater_modules.csm` - Contaminant Simulation Module
   - :mod:`clearwater_modules.msm` - Mercury Simulation Module

**Shared Components**
   - :mod:`clearwater_modules.shared` - Shared utilities and types
   - :mod:`clearwater_modules.base` - Base classes for all modules
   - :mod:`clearwater_modules.utils` - General utilities

**Supporting Components**
   - :mod:`clearwater_modules.sorter` - Topological sorting for equations