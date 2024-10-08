[![Tests Status](https://github.com/EcohydrologyTeam/ClearWater-modules/actions/workflows/tests.yml/badge.svg)](https://github.com/EcohydrologyTeam/ClearWater-modules/actions/workflows/tests.yml)
[![Coverage](https://codecov.io/gh/EcohydrologyTeam/ClearWater-modules/graph/badge.svg)](https://codecov.io/gh/EcohydrologyTeam/ClearWater-modules)

# ClearWater Modules in Python

The [ClearWater-modules](https://github.com/EcohydrologyTeam/ClearWater-modules) package is a collection of water quality and vegetation process simulation modules written in modern Python and designed to flexibily couple with a variety of water transport models, such as RAS-2D, GSSHA, CE-Qual-W2, [AdH](https://www.erdc.usace.army.mil/Locations/CHL/AdH/), and others. These modules have been developed by the [U.S. Army Engineer Research and Development Center (ERDC)](https://www.erdc.usace.army.mil), [Environmental Laboratory (EL)](https://www.erdc.usace.army.mil/Locations/EL/).

- [TSM: Temperature Simulation Module](src/clearwater_modules/tsm) (formerly TEMP)
- [NSM: Nutrient Simulation Modules](src/clearwater_modules/nsm1) ([NSM-I](src/clearwater_modules/nsm1) and [NSM-II](src/clearwater_modules/nsm2))
- [GSM: General Constituent Simulation Module](src/clearwater_modules/gsm)
- [CSM: Contaminant Simulation Module](src/clearwater_modules/csm)
- [MSM: Mercury Simulation Module](src/clearwater_modules/msm)
- SSM: Solids Simulation Module (Fortran only)
- RVSM: Riparian Vegetation Simulation Module (Fortran only)

These water quality modules form the central capabilities of the [ClearWater (Corps Library for Environmental Analysis and Restoration of Watersheds)](https://ui.adsabs.harvard.edu/abs/2023EGUGA..2512470S/abstract) software system. The overall goal of the ClearWater system is to couple these water quality simulation capabilites to state-of-the art hydrologic and hydraulic modeling tools, such as HEC-RAS-2D, CE-Qual-W2, and GSSHA, allowing users to leverage existing river, reservoir, and waterhed models for water quality studies. At present, TSM and NSM have been successfully coupled to HEC-RAS-2D models via the [ClearWater-riverine]([url](https://github.com/EcohydrologyTeam/ClearWater-riverine)) package.

A secondary goal is to develop a suite of easy-to-use modern Python tools that build on community-developed scientific workflows, standards, and libraries to automate model setup, prepare input datasets, store output data, and visualize results using Python-based user interfaces such as Jupyter Notebooks.

This Python library is a port and modernization of the algorithms and structures originally written in Fortran 95, released as version 1.0 in 2021, and described in:

- Zhang, Zhonglong and Billy E. Johnson. 2016. Aquatic nutrient simulation modules (NSMs) developed for hydrologic and hydraulic models. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Ecosystem Management and Restoration Research Program (EMRRP). ERDC/EL Technical Report 16-1. https://hdl.handle.net/11681/10112
- Zhang, Zhonglong and Billy E. Johnson. 2016. Aquatic contaminant and mercury simulation modules developed for hydrologic and hydraulic models. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Environmental Quality Technology Research Program (EQTRP). ERDC/EL Technical Report 16-8. https://hdl.handle.net/11681/20249
- Johnson, Billy E. and Zhonglong Zhang. 2016. Testing and Validation Studies of the NSMII-Benthic Sediment Diagenesis Module. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Ecosystem Management and Restoration Research Program (EMRRP). ERDC/EL Technical Report 16-11. https://hdl.handle.net/11681/20343

## Repository Directories


## Getting Started

### Installation

Clearwater-modules was developed with **Python 3.11**. 

Follow these steps to install.

#### 1. Install Miniconda or Anaconda Distribution

We recommend installing the light-weight [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) that includes Python, the [conda](https://conda.io/docs/) environment and package management system, and their dependencies.

NOTE: Follow conda defaults to install in your local user director. DO NOT install for all users, to avoid substantial headaches with permissions.

If you have already installed the [**Anaconda Distribution**](https://www.anaconda.com/download), you can use it to complete the next steps, but you may need to [update to the latest version](https://docs.anaconda.com/free/anaconda/install/update-version/).

#### 2. Clone or Download this `ClearWater-modules` repository

From this Github site, click on the green "Code" dropdown button near the upper right. Select to either Open in GitHub Desktop (i.e. git clone) or "Download ZIP". We recommend using GitHub Desktop, to most easily receive updates.

Place your copy of this repo folder in any convenient location on your computer.

#### 3. Create a Conda Environment for this Repository (optional) 

We recommend creating a custom virtual environment with the [Conda](https://conda.io/docs/) package, dependency, and environment management for any language (i.e. easily install C++ packages such as GDAL).

We provide an [`environment.yml`](environment.yml) file that specifies for [Conda](https://conda.io/docs/) how to create a virtual environment that contains the same software dependencies that we've used in development and testing.

Create a `ClearWater-modules` environment using this [conda](https://conda.io/docs/) command in your terminal or Anaconda Prompt console. If necessary, replace `environment.yml` with the full file pathway to the `environment.yml` file in the local cloned repository.

```shell
conda env create --file environment.yml
```

Alternatively, use the faster [`libmamba` solver](https://conda.github.io/conda-libmamba-solver/getting-started/) with:

```shell
conda env create -f environment.yml --solver=libmamba
```

Activate the environment using the instructions printed by conda after the environment is created successfully.

To update your environment to the latest versions of dependencies and/or add additional dependencies to your environment (by first editting [`environment.yml`](environment.yml)), run the following command:

```shell
conda env update -f environment.yml --solver=libmamba --prune
```

or to recreate from scratch:

```shell
conda env create -f environment.yml --solver=libmamba --force
```

For additional information on managing conda environments, see [Conda's User Guide on Managing Environments](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html).

#### 4. Add your `ClearWater-modules` Path to Miniconda/Anaconda sites-packages

To have access to the `clearwater_riverine` module in your Python environments, it is necessary to have a path to your copy of Clearwater Riverine in Anaconda's `sites-packages` directory (i.e. something like `$HOME/path/to/anaconda/lib/pythonX.X/site-packages` or `$HOME/path/to/anaconda/lib/site-packages` similar).

The easiest way to do this is to use the [conda develop](https://docs.conda.io/projects/conda-build/en/latest/resources/commands/conda-develop.html) command in the console or terminal like this, replacing `/path/to/module/` with the full file pathway to the local cloned Clearwater-riverine repository:

```console
conda develop /path/to/module/
```

You should now be able to run the examples and create your own Jupyter Notebooks!


### Examples


## Contributing


## Acknowlgements

The vision for modernizing this library, including the initial port to Python from Fortran, was develped by:

- Dr. Todd E. Steissberg (ERDC-EL)

The algorithms and structure of this program were adapted from the Fortran 95 version 1.0 of these modules, originally developed by:

- Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
- Dr. Zhonglong Zhang (Portland State University, LimnoTech)
- Mr. Mark Jensen (USACE HEC)
