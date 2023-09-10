# ClearWater-modules-python

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
