"""
File includes static variables only used in Algae module
"""

import clearwater_modules_python.base as base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.phosphorus.phosphorus_processes as phosphorus_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

#Only phosphorus variables

Variable(
    name='kop_20',
    long_name='Decay rate of organic P to DIP',
    units='1/d',
    description='Decay rate of organic P to DIP',
    use='static',
)

Variable(
    name='rpo4_20',
    long_name='Benthic sediment release rate of DIP',
    units='g-P/m^2/d',
    description='Benthic sediment release rate of DIP',
    use='static',
)
