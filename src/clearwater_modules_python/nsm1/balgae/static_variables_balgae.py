"""
File includes static variables only used in Algae module
"""

import clearwater_modules_python.base as base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.algae.algae_processes as algae_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

#Only Balgae Variables 

Variable(
    name='Fw',
    long_name='Fraction of benthic algae mortality into water column',
    units='unitless',
    description='Fraction of benthic algae mortality into water column',
    use='static',
)

Variable(
    name='Fb',
    long_name='Fraction of bottom area available for benthic algae',
    units='unitless',
    description='Fraction of bottom area available for benthic algae',
    use='static',
)
