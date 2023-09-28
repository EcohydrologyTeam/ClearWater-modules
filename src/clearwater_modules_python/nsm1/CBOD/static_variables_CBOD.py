# TODO: figure out what 'model' to import
import clearwater_modules_python.base as base
from clearwater_modules_python.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...

#CBOD variables for each CBOD group - array
Variable(
    name='kbod_i_20',
    long_name='CBOD oxidation rate for each CBOD group at 20C',
    units='/d',
    description='CBOD oxidation rate for each CBOD group at 20C, array',
    use='static'
)

Variable(
    name='ksbod_i_20',
    long_name='CBOD sedimentation rate for each CBOD group at 20C',
    units='m/d',
    description='CBOD sedimentation rate for each CBOD group at 20C, array',
    use='static'
)

Variable(
    name='ksOxbod_i',
    long_name='Half saturation oxygen attenuation constant for CBOD oxidation for each CBOD group',
    units='mg-O2/L',
    description='Half saturation oxygen attenuation constant for CBOD oxidation for each CBOD group, array',
    use='static'
)