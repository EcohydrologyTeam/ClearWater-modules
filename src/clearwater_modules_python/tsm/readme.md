# Temperature Simulation Module (TSM)

## CLEARWATER (Corps Library for Environmental Analysis and Restoration of Watersheds)

This is a Python version of a new Temperature Simulation Module (TSM) built from the energy budget algorithms in the original temperature simulation module. 
This module is compatible with `xarray` based inputs, and can run rapidly across a grid of input values since each cell is computed independent from it's neighbors.

This module is a component of CLEARWATER (Corps Library for Environmental Analysis and Restoration of Watersheds).


# Credits

Developed by:
* Dr. Todd E. Steissberg (ERDC-EL)
* Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)
* Mr. Isaac Mudge (MVN)
* Mr. Xavier Nogueira (LimnoTech)

The algorithms and structure of this program were adapted from the Fortran 95 version of this module, 
developed by:
* Dr. Billy E. Johnson (ERDC-EL)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)
