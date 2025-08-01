Overview
========

Introduction
------------

The ClearWater-modules package is a comprehensive suite of water quality and vegetation process simulation modules designed for flexible integration with various water transport models. Developed by the U.S. Army Engineer Research and Development Center (ERDC), Environmental Laboratory (EL), these modules provide state-of-the-art capabilities for water quality modeling.

Purpose
-------

The primary goals of the ClearWater system are:

1. **Coupling water quality simulation capabilities** to state-of-the-art hydrologic and hydraulic modeling tools such as HEC-RAS-2D, CE-Qual-W2, and GSSHA
2. **Enabling users to leverage existing models** for river, reservoir, and watershed water quality studies
3. **Providing modern Python tools** that build on community-developed scientific workflows, standards, and libraries

Available Modules
-----------------

Temperature Simulation Module (TSM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The TSM (formerly TEMP) simulates water temperature dynamics in aquatic systems, accounting for various heat exchange processes including solar radiation, atmospheric heat exchange, and bed heat conduction.

Nutrient Simulation Modules (NSM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NSM consists of two versions:

- **NSM-I**: A simplified nutrient simulation module
- **NSM-II**: An advanced nutrient simulation module with enhanced capabilities

Both modules simulate the cycling of nutrients including nitrogen, phosphorus, and carbon in aquatic systems.

General Constituent Simulation Module (GSM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GSM provides capabilities for simulating general water quality constituents that don't require the complexity of the nutrient or contaminant modules.

Contaminant Simulation Module (CSM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CSM simulates the fate and transport of contaminants in aquatic systems, including processes such as sorption, decay, and transformation.

Mercury Simulation Module (MSM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MSM specifically addresses mercury cycling in aquatic environments, including methylation and bioaccumulation processes.

Architecture
------------

The ClearWater modules are designed with a modular architecture that allows:

- **Flexible coupling** with different hydrodynamic models
- **Independent module execution** for specific water quality constituents
- **Scalable computation** from small streams to large watersheds
- **Modern Python implementation** for ease of use and integration

Historical Background
---------------------

This Python library is a port and modernization of algorithms originally written in Fortran 95, released as version 1.0 in 2021. The original work is documented in several ERDC technical reports that provide detailed descriptions of the underlying science and algorithms.

Key Features
------------

- **Open-source implementation** in modern Python
- **Extensive documentation** and examples
- **Unit testing** for quality assurance
- **Integration with scientific Python ecosystem** (NumPy, pandas, xarray)
- **Jupyter notebook support** for interactive modeling