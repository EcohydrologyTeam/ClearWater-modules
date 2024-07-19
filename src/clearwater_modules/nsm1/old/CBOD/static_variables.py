"""
File contains static variables related to the CBOD module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...


Variable(
    name='kbod_20',
    long_name='CBOD oxidation rate at 20C',
    units='1/d',
    description='CBOD oxidation rate at 20C',
    use='static'
)

Variable(
    name='ksbod_20',
    long_name='CBOD sedimentation rate at 20C',
    units='m/d',
    description='CBOD sedimentation rate at 20C',
    use='static'
)

Variable(
    name='KsOxbod',
    long_name='Half saturation oxygen attenuation constant for CBOD oxidation',
    units='mg-O2/L',
    description='Half saturation oxygen attenuation constant for CBOD oxidation',
    use='static'
)

