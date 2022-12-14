# Temperature Simulation Module (TSM)

## CLEARWATER (Corps Library for Environmental Analysis and Restoration of Watersheds) Version
## Version 1.0
## April 11, 2021

This is version 1.0 of a new Temperature Simulation Module (TSM) built from the energy budget algorithms in the original temperature simulation module. 

Developed by:
* Dr. Todd E. Steissberg (ERDC-EL)
* Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)
* Mr. Isaac Mudge (MVN)

This module computes the water quality of a single computational cell. The algorithms 
and structure of this program were adapted from the Fortran 95 version of this module, 
developed by:
* Dr. Billy E. Johnson (ERDC-EL)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

This module contains one submodule:
* Temperature Energy Budget Module
  * This module computes the complete energy budget.


The submodule compute the values for a single cell only. The original version contained two separate modules that computed the values for an array of cells for each 
of multiple regions. The objective of the new module is to compute values for a single cell only, so each module can be called in sequence by an external framework. Furthermore, the original module contained extensive code to handle XML strings for interoperability with a .Net GUI. The external framework will handle the array processing, sub-regions, etc. This allows the module to be more portable while simplifying the essential logic to focus on the water quality calculations.

Note: All constants and initial values are specified using a dictionary. (The C++ version uses an unordered_map to simulate a Python dictionary.) The dictionary enables extention of this code to replace the values from external sources, such as a UI, or import from an external file (CSV, JSON, XML, etc.)

To compile the sub-modules into executable programs:
g++ -std=c++14 temperature_equilibrium.cpp -o temperature_equilibrium
g++ -std=c++14 temperature_energy_budget.cpp -o temperature_energy_budget
