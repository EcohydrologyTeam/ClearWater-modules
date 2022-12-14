# General Constituent Simulation Module (GSM)

## CLEARWATER (Corps Library for Environmental Analysis and Restoration of Watersheds) Version
## Version 1.0
## April 11, 2021

Developed by:
* Dr. Todd E. Steissberg (ERDC-EL)
* Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

This module computes the water quality of a single computational cell. The algorithms 
and structure of this program were adapted from the Fortran 95 version of this module, 
developed by:
* Dr. Billy E. Johnson (ERDC-EL)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

This module computes the change of concentration of a general constituent for a single cell only. The original version contained two separate modules that computed the values for an array of cells for each of multiple regions. The objective of the new module is to compute values for a single cell only, so each module can be called in sequence by an external framework. Furthermore, the original module contained extensive code to handle XML strings for interoperability with a .Net GUI. The external framework will handle the array processing, sub-regions, etc. This allows the module to be more portable while simplifying the essential logic to focus on the water quality calculations.
