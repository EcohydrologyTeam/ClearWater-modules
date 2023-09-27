# TODO: figure out what 'model' to import
import clearwater_modules_python.base as base
from clearwater_modules_python.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...

Variable(
    name='kaw_20',
    long_name='Wind oxygen reaeration velocity at 20C',
    units='m/d',
    description='Wind oxygen reaeration velocity at 20C',
    use='static'
)

Variable(
    name='kah_20',
    long_name='Hydraulic oxygen reaeration rate at 20C',
    units='/d',
    description='Hydraulic oxygen reaeration rate at 20C',
    use='static'
)

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