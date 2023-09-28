# TODO: figure out what 'model' to import
import clearwater_modules_python.base as base
from clearwater_modules_python.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...

Variable(
    name='kaw_20_user',
    long_name='Wind oxygen reaeration velocity at 20C',
    units='m/d',
    description='Wind oxygen reaeration velocity at 20C',
    use='static'
)

Variable(
    name='kah_20_user',
    long_name='Hydraulic oxygen reaeration rate at 20C',
    units='/d',
    description='Hydraulic oxygen reaeration rate at 20C',
    use='static'
)

Variable(
    name='hydraulic_reaeration_option',
    long_name='Option for chosing the method by which O2 reaeration rate is calculated',
    units='unitless',
    description='Selects method for computing O2 reaeration rate',
    use='static'
)

Variable(
    name='wind_reaeration_option',
    long_name='Option for chosing the method by which wind reaeration is calculated',
    units='unitless',
    description='Selects method for computing O2 reaeration due to wind',
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