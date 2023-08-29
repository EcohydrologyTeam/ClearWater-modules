ClearWater Modules in Python
=====

The [ClearWater-modules-python](https://github.com/EcohydrologyTeam/ClearWater-modules-python) package is a collection of water quality and vegetation process simulation modules developed by the [U.S. Army Engineer Research and Development Center (ERDC)](https://www.erdc.usace.army.mil), [Environmental Laboratory (EL)](https://www.erdc.usace.army.mil/Locations/EL/).

- TSM: Temperature Simulation Module (formerly TEMP)
- NSM: Nutrient Simulation Modules (NSM-I and NSM-II)
- GSM: General Constituent Simulation Module
- CSM: Contaminant Simulation Module
- MSM: Mercury Simulation Module (Fortran only)
- SSM: Solids Simulation Module (Fortran only)
- RVSM: Riparian Vegetation Simulation Module (Fortran only)

These water quality modules form the central capabilities of the ClearWater (Corps Library for Environmental Analysis and Restoration of Watersheds) software system. The overall goal of the ClearWater system is to couple these water quality simulation capabilites to state-of-the art hydrologic and hydraulic modeling tools, such as HEC-RAS-2D, CE-Qual-W2, and GSSHA, allowing users to leverage existing river, reservoir, and waterhed models for water quality studies. A secondary goal is to develop a suite of easy-to-use modern Python tools that build on community-developed scientific workflows, standards, and libraries to automate model setup, prepare input datasets, store output data, and visualize results using Python-based user interfaces such as Jupyter Notebooks.

This Python library is a port and modernization of the algorithms and structures originally written in Fortran 95,  released as version 1.0 in 2021, and described in:

- Zhang, Z. and Johnson, B.E., 2016. Aquatic nutrient simulation modules (NSMs) developed for hydrologic and hydraulic models. Ecosystem Management and Restoration Research Program (EMRRP). ERDC/EL Technical Report 16-1. https://hdl.handle.net/11681/10112
- Zhang, Z. and Johnson, B.E., 2016. Aquatic contaminant and mercury simulation modules developed for hydrologic and hydraulic models. Environmental Quality Technology Research Program (EQTRP). ERDC/EL Technical Report 16-8. https://hdl.handle.net/11681/20249

# Getting Started

## Installation

## Examples


# Contributing


# Acknowlgements

The vision for modernizing this library, including the initial port to Python from Fortran, was develped by:

- Dr. Todd E. Steissberg (ERDC-EL)

The algorithms and structure of this program were adapted from the Fortran 95 version 1.0 of these modules, originally developed by:

- Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
- Dr. Zhonglong Zhang (Portland State University, LimnoTech)
- Mr. Mark Jensen (USACE HEC)
