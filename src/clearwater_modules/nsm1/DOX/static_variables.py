# TODO: figure out what 'model' to import
import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...

Variable(
    name='patm',
    long_name='Atmospheric pressure',
    units='atm',
    description='Atmospheric pressure',
    use='static'
)

Variable(
    name='KsSOD',
    long_name='half saturation oxygen attenuation constant for SOD',
    units='mg/L',
    description='half saturation oxygen attenuation constant for SOD',
    use='static'
)
