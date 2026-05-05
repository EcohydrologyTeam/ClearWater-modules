[![Tests Status](https://github.com/EcohydrologyTeam/ClearWater-modules/actions/workflows/tests.yml/badge.svg)](https://github.com/EcohydrologyTeam/ClearWater-modules/actions/workflows/tests.yml)
[![Coverage](https://codecov.io/gh/EcohydrologyTeam/ClearWater-modules/graph/badge.svg)](https://codecov.io/gh/EcohydrologyTeam/ClearWater-modules)

# ClearWater Modules

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

**[`src`](src)** contains the source code to create and run the `clearwater_modules`.

**[`examples`](examples)** contains tutorials and useful Juptyer Notebooks.

**[`docs`](docs)** contains relevant reference documentation.

**[`tests`](tests)** will contain clearwater_riverine tests once they are developed. 

## Getting Started

### Installation

Follow these steps to install the ClearWater Modeling System and its dependancies in a custom Python environment.

We recommend using [pixi](https://pixi.prefix.dev/latest/), the next-generation reproducible package management tool built on [conda](https://docs.conda.io/projects/conda/en/stable/) tooling.

If you are new to pixi but familiar with conda, this [Switching from Conda](https://pixi.prefix.dev/latest/switching_from/conda/) documentation succinctly compares similarities and differences.

Alternately, use a conda environment with the same dependencies.

#### 1. Install Pixi

Follow [Pixi Installation](https://pixi.prefix.dev/latest/installation/) instructions for your platform.

#### 2. Clone or Download the ClearWater family of repositories

There are three repositories which house the dependencies for this project. Navigate to each of the repositories listed below and follow the instructions to clone them to your local machine.

- [ClearWater-riverine](https://github.com/EcohydrologyTeam/ClearWater-riverine) transport process simulator
- [ClearWater-data](https://github.com/EcohydrologyTeam/ClearWater-data) perforamant data access and storage protocols
- [ClearWater-modules](https://github.com/EcohydrologyTeam/ClearWater-modules) (optional) water quality reaction process simulator

From these Github sites, click on the green "Code" dropdown button near the upper right. Select to either "Open in GitHub Desktop" (i.e. git clone) or "Download ZIP". 

We recommend using GitHub Desktop, to most easily manage git workflows by providing excellent visuals for stagging commits, exploring commit histories, comparing branches, and resolving merge conflicts in tight integration with Visual Studio Code.

Place your copy of these repos in any convenient location on your computer. Make sure that all are stored in the same directory OR you will need to update the `pyproject.toml` so that the `[pixi.tool.feature.dev.pypi-dependencies]` point to the correct location of your local clones.

### 3. Create Clearwater Workspace and Python Environments 

Create a project-specific Pixi workspace and Python enviornment(s) from the `pyproject.toml` manifest file.

#### Production 

Pardon our mess! Production instructions coming soon. For now you can install using the [Developer](#developer) instructions. 

#### Developers

Developer instructions install all ClearWater repos in "editable" mode in a second `dev` environment.

From your terminal or console, navigate to the directory of your `Clearwater-modules` clone and execute the following command to create a `dev` environment:

```shell
pixi install -e dev
```

To activate this environment in your shell, run the following:

```shell
pixi shell -e dev
```

You should now be able to run the examples and create your own Jupyter Notebooks!

### Examples

Try running our [03_Example_Coupled_TSM_and _Riverine.ipynb](examples/03_Example_Coupled_TSM_and _Riverine.ipynb) Jupyter Notebook.


## Contributing

We welcome your pull request.

## Acknowledgements

This library is developed by ERDC-EL through funding from the ECOMOD project.

Dr. Todd E. Steissberg (ERDC-EL) developed the vision for modernizing TSM, NSM, and related modules into this library, including coding the initial port form Fortran to Python.

The algorithms and structure of this program were adapted from the Fortran 95 version 1.0 of these modules, originally developed by:

- Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
- Dr. Zhonglong Zhang (Portland State University, LimnoTech)
- Mr. Mark Jensen (USACE HEC)
