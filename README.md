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

## Repository Directories


# Getting Started

## Installation

Clearwater-modules-python was developed with **Python 3.11**. 

Follow these steps to install.

#### 1. Install Miniconda or Anaconda Distribution

We recommend installing the light-weight [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) that includes Python, the [conda](https://conda.io/docs/) environment and package management system, and their dependencies.

If you have already installed the [**Anaconda Distribution**](https://www.anaconda.com/download), you can use it to complete the next steps, but you may need to [update to the latest version](https://docs.anaconda.com/free/anaconda/install/update-version/).



#### 2. Clone or Download this `ClearWater-modules-python` repository

From this Github site, click on the green "Code" dropdown button near the upper right. Select to either Open in GitHub Desktop (i.e. git clone) or "Download ZIP". We recommend using GitHub Desktop, to most easily receive updates.

Place your copy of this repo folder in any convenient location on your computer.

#### 3. Create a Conda Environment for this Repository (optional) 

We recommend creating a custom virtual environment with the same software dependencies that we've used in development and testing, as listed in the [`environment.yml`](environment.yml) file. 

Create a `ClearWater-modules-python` environment using this [conda](https://conda.io/docs/) command in your terminal or Anaconda Prompt console. If necessary, replace `environment.yml` with the full file pathway to the `environment.yml` file in the local cloned repository.


```shell
conda env create --file environment.yml
```

Alternatively, use the faster [`libmamba` solver](https://conda.github.io/conda-libmamba-solver/getting-started/) with:
```shell
conda env create -f environment.yml --solver=libmamba
```

Activate the environment using the instructions printed by conda after the environment is created successfully.

For additional information on managing conda environments, see [Conda's User Guide on Managing Environments](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html).


#### 4. Add your `ClearWater-modules-python` Path to Miniconda/Anaconda sites-packages

To have access to the `ClearWater-modules-python` module in your Python environment, it is necessary to have add it's path to enviornment's PATH variable.

The easiest way to do this is to use the [conda develop](https://docs.conda.io/projects/conda-build/en/latest/resources/commands/conda-develop.html) command in the console or terminal like this, replacing `/path/to/module/` with the full file pathway to the local cloned Clearwater-riverine repository:

```console
conda-develop /path/to/module/
```

You should now be able to run the examples and create your own Jupyter Notebooks!


## Examples


# Contributing


# Acknowlgements

The vision for modernizing this library, including the initial port to Python from Fortran, was develped by:

- Dr. Todd E. Steissberg (ERDC-EL)

The algorithms and structure of this program were adapted from the Fortran 95 version 1.0 of these modules, originally developed by:

- Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
- Dr. Zhonglong Zhang (Portland State University, LimnoTech)
- Mr. Mark Jensen (USACE HEC)
=======
This repository contains the ClearWater (Corps Library for Environmental Analysis and Restoration of Watersheds) modules. The following modules have been developed in Fortran. These are being rewritten in Python, starting with TSM and NSM.

- TSM: Temperature Simulation Module
- NSM: Nutrient Simulation Module (NSM-I and NSM-II)
- GSM: General Constituent Simulation Module
- CSM: Contaminant Simulation Module
- MSM: Mercury Simulation Module
- SSM: Solids Simulation Module

## Module Descriptions

### TSM: Temperature Simulation Module
The Temperature Simulation Module (TSM) (Zhang and Johnson, 2016) is an essential component of ClearWater. This module plays a crucial role in simulating and predicting water temperature within aquatic ecosystems. TSM utilizes a comprehensive energy balance approach to account for various factors contributing to heat inputs and outputs in the water environment. It considers both external forcing functions and heat exchanges occurring at the water surface and the sediment-water interface. The primary contributors to heat exchange at the water surface include shortwave solar radiation, longwave atmospheric radiation, heat conduction from the atmosphere to the water, and direct heat inputs. Conversely, the primary factors that remove heat from the system are longwave radiation emitted by the water, evaporation, and heat conduction from the water to the atmosphere.

The core principle behind TSM is the application of the laws of conservation of energy to compute water temperature. This means that the change in heat content of the water is directly related to changes in temperature, which, in turn, are influenced by various heat flux components. The specific heat of water is employed to establish this relationship. Each term of the heat flux equation can be calculated based on the input provided by the user, allowing for flexibility in modeling different environmental conditions.

For a more detailed understanding of the equations and underlying mechanisms used in TSM, users can refer to relevant literature sources such as Water Resources Engineers Inc. (1967), Brown and Barnwell (1987), and Deas and Lowney (2000). Much of the content in this module has been derived from HEC (2023) to ensure consistency with the original water temperature model and to facilitate accurate kinetic implementations. TSM is a valuable tool for environmental scientists and water resource managers to simulate and analyze temperature dynamics in aquatic systems, aiding in the preservation and restoration of watersheds.

### NSM: Nutrient Simulation Module (NSM-I and NSM-II)

### GSM: General Constituent Simulation Module

### CSM: Contaminant Simulation Module

### MSM: Mercury Simulation Module

### SSM: Solids Simulation Module

## References

- Andrews, R. 1980. Wärmeaustausch zwischen Wasser und Wattboden (Heat exchange between water and tidal flats). Deutsche Gewässerkundliche Mitteilungen. 24. 57-65.
- Bejan, A. 1993. Heat Transfer. Wiley, New York, NY. 
- Brown, L. C., and T. O. Barnwell. 1987. The enhanced stream water quality models QUAL2E and QUAL2E-UNCAS. EPA/600/3-87-007. Athens, GA: U.S. Environmental Protection Agency.
- Carslaw, H.S. and Jaeger, J.C. 1959. Conduction of Heat in Solids, Oxford Press, Oxford, UK, 510 pp.
- Cengel, Y.A. 1998 Heat Transfer: A Practical Approach. New York, McGraw-Hill.
- Chow, V.T., Maidment, D.R., and Mays, L.W. 1988. Applied Hydrology. New York, McGraw- Hill, 592 pp.
- Deas, M. L. and C. L. Lowney. 2000. Water temperature modeling review. Central Valley, CA.
- Geiger, R. 1965. The climate near the ground. Harvard University Press. Cambridge, MA.
- Grigull, U. and Sandner, H. 1984. Heat Conduction. Springer-Verlag, New York, NY.
- Hutchinson, G.E. 1957. A Treatise on Limnology, Vol. 1, Physics and Chemistry. Wiley, New York, NY.
- Hydrologic Engineering Center (HEC). 2023. HEC-RAS: River Analysis System User’s Reference Manual, Version 6.1. Davis, CA: Hydrologic Engineering Center. U.S. Corps of Engineers.
- Jobson, H.E. 1977. Bed Conduction Computation for Thermal Models. J. Hydraul. Div. ASCE. 103(10):1213-1217.
- Johnson, Billy E. and Zhonglong Zhang. 2016. Testing and Validation Studies of the NSMII-Benthic Sediment Diagenesis Module. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). http://www.dtic.mil/docs/citations/AD1012495.
- Kreith, F. and Bohn, M.S. 1986. Principles of Heat Transfer, 4th Ed. Harper and Row, New York, NY. 
- Likens, G. E., and Johnson, N. M. (1969). Measurements and analysis of the annual heat budget for sediments of two Wisconsin lakes. Limnol. Oceanogr., 14(1):115-135.
- Mills, A.F. 1992. Heat Transfer. Irwin, Homewood, IL. 
- Nakshabandi, G.A. and H. Kohnke. 1965. Thermal conductivity and diffusivity of soils as related to moisture tension and other physical properties. Agr. Met. Vol 2.
- Water Resources Engineers Inc. 1967. Prediction of thermal distribution in streams and reservoirs, report to California Department of Fish and Game, Walnut Creek, CA.
- Zhang, Zhonglong and Billy E. Johnson. 2016. Aquatic Contaminant and Mercury Simulation Modules Developed for Hydrologic and Hydraulic Models. Technical Report. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). DOI: 10.21236/AD1013220. http://www.dtic.mil/docs/citations/AD1013220.
- Zhang, Zhonglong and Billy E Johnson. 2016. Aquatic Nutrient Simulation Modules (NSMs) Developed for Hydrologic and Hydraulic Models. Technical Report. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC)
