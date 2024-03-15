"""
=======================================================================================
Nutrient Simulation Module 1 (NSM1)
=======================================================================================

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

Version 1.0

Initial Version: June 5, 2021
"""

from clearwater_modules.nsm1 import algae
from clearwater_modules.nsm1 import alkalinity
from clearwater_modules.nsm1 import balgae
from clearwater_modules.nsm1 import carbon
from clearwater_modules.nsm1 import CBOD
from clearwater_modules.nsm1 import DOX
from clearwater_modules.nsm1 import nitrogen
from clearwater_modules.nsm1 import POM


class NSM1:
    def __init__(self):
        pass
