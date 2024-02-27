"""
File contains static variables related to the Phosphorus module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget

@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...


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


Variable(
    name='kdop4',
    long_name='solid partitioning coeff. of PO4',
    units='L/kg',
    description='solid partitioning coeff. of PO4',
    use='static',
)