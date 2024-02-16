"""
File includes static variables only used in Pathogen module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.pathogens.processes as processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

Variable(
    name='kdx_20',
    long_name='Pathogen death rate at 20C',
    units='1/d',
    description='Pathogen death rate at 20C',
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
