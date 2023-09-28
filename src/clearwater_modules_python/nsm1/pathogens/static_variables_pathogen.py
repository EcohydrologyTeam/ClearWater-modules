"""
File includes static variables only used in Pathogen module
"""

import clearwater_modules_python.base as base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.pathogens.pathogen_processes as pathogen_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

Variable(
    name='kdx',
    long_name='Pathogen death rate',
    units='1/d',
    description='Pathogen death rate',
    use='static',
)

Variable(
    name='apx',
    long_name='Light efficiency factor for pathogen decay',
    units='unitless',
    description='Light efficiency factor for pathogen decay',
    use='static',
)

Variable(
    name='vx',
    long_name='Pathogen net settling velocity',
    units='unitless',
    description='Pathogen net settling velocity',
    use='static',
)
