# TODO: figure out what 'model' to import
import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


# CBOD variables for each CBOD group - array
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
