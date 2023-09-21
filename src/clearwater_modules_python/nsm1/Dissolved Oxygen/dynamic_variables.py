# TODO: figure out imports

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.tsm.model import EnergyBudget
from clearwater_modules_python.nsm1 import processes


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...

Variable(
    name='DOX_sat',
    long_name='DO saturation concentration',
    units='mg/L',
    description='DO saturation concentration in water as a function of water temperature (K)',
    use='dynamic',
    process=processes.DOX_sat
)

Variable(
    name='pwv',
    long_name='Partial pressure of water vapor',
    units='atm',
    description='Partial pressure of water vapor',
    use='dynamic',
    process=processes.pwv
)

Variable(
    name='DOs_atm_alpha',
    long_name='DO saturation atmospheric correction coefficient',
    units='unitless',
    description='DO saturation atmospheric correction coefficient',
    use='dynamic',
    process=processes.DOs_atm_alpha
)

Variable(
    name = 'kah_20_r',
    long_name='Hydraulic oxygen reaeration rate adjusted for hydraulics',
    units='/d',
    description='Hydraulic oxygen reaeration rate adjusted for hydraulic parameters according to XX lit',
    use='dynamic',
    process=processes.kah_20_r
)

Variable(
    name='kah_T_r',
    long_name='Hydraulic oxygen reaeration rate adjusted for temperature',
    units='/d',
    description='Hydraulic oxygen reaeration rate adjusted for temperature',
    use='dynamic',
    process=processes.kah_T_r
)

Variable(
    name = 'kaw_20_r',
    long_name='Wind oxygen reaeration velocity adjusted for hydraulics',
    units='m/d',
    description='Wind oxygen reaeration velocity adjusted for hydraulic parameters according to XX lit',
    use='dynamic',
    process=processes.kaw_20_r
)

Variable(
    name='kaw_T_r',
    long_name='Wind oxygen reaeration velocity adjusted for temperature',
    units='m/d',
    description='Wind oxygen reaeration velocity adjusted for temperature',
    use='dynamic',
    process=processes.kaw_T_r
)