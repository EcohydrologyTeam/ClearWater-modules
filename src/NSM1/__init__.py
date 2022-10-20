'''
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
'''

from ._algae_split import Algae
from ._benthic_algae import BenthicAlgae
from ._carbon import Carbon
from ._cbod import CBOD
from ._dox import DOX
from ._n2 import N2
from ._nitrogen import Nitrogen
from ._pathogen import Pathogen
from ._phosphorus import Phosphorus
from ._pom import POM
from ._sed_flux import SedFlux
from ._globals import Globals
from ._temp_correction import TempCorrection

from ._kelsey_global import Kelsey_Global
from ._kelsey_algae import Algae

class NSM1:
    def __init__(self):
        pass
