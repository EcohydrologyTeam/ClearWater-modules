"""
File contains static variables related to the DOX module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...


Variable(
    name='ron',
    long_name='O2:N ratio for nitrification',
    units='mg-O2/mg-N',
    description='2*32/14',
    use='static'
)


Variable(
    name='KsSOD',
    long_name='half saturation oxygen attenuation constant for SOD',
    units='mg/L',
    description='half saturation oxygen attenuation constant for SOD',
    use='static'
)
