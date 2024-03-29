import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget

import clearwater_modules.tsm as tsm
import clearwater_modules.tsm.processes as tsm_processes

import clearwater_modules.shared.processes as shared_processes

@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...


Variable(
    name='depth',
    long_name='Average water depth in cell',
    units='m',
    description='Average water depth in cell computed by dividing volume by surface area',
    use='dynamic',
    process=shared_processes.compute_depth
)

Variable(
    name='TwaterC',
    long_name='Water Temperature',
    units='C',
    description='Water Temperature Degree Celsius',
    use='dynamic',
    process=tsm_processes.t_water_c
)

Variable(
    name='TwaterK',
    long_name='Water Temperature K',
    units='K',
    description='Water temperature degree kelvin',
    use='dynamic',
    process=shared_processes.TwaterK
)

Variable(
    name='SOD_tc',
    long_name='Sediment Oxygen Demand at water temperature tc',
    units='mg/L',
    description='Sediment Oxygen Demand at water temperature tc',
    use='dynamic',    
    process=shared_processes.SOD_tc
)

Variable(
    name='kah_20',
    long_name='Hydraulic oxygen reaeration rate adjusted for hydraulics',
    units='1/d',
    description='Hydraulic oxygen reaeration rate adjusted for hydraulic parameters according to XX lit',
    use='dynamic',
    process=shared_processes.kah_20
)

Variable(
    name='kah_tc',
    long_name='Hydraulic oxygen reaeration rate adjusted for temperature',
    units='1/d',
    description='Hydraulic oxygen reaeration rate adjusted for temperature',
    use='dynamic',
    process=shared_processes.kah_tc
)

Variable(
    name='kaw_20',
    long_name='Wind oxygen reaeration velocity adjusted for hydraulics',
    units='m/d',
    description='Wind oxygen reaeration velocity adjusted for hydraulic parameters according to XX lit',
    use='dynamic',
    process=shared_processes.kaw_20
)

Variable(
    name='kaw_tc',
    long_name='Wind oxygen reaeration velocity adjusted for temperature',
    units='m/d',
    description='Wind oxygen reaeration velocity adjusted for temperature',
    use='dynamic',
    process=shared_processes.kaw_tc
)

Variable(
    name='ka_tc',
    long_name='Oxygen reaeration rate',
    units='1/d',
    description='Oxygen reaeration rate',
    use='dynamic',
    process=shared_processes.ka_tc
)


Variable(
    name='L',
    long_name='Light attenuation coefficient',
    units='unitless',
    description='Light attenuation coefficient',
    use='dynamic',
    process=shared_processes.L
)

Variable(
    name='PAR',
    long_name='surface light intensity',
    units='W/m2',
    description='surface light intensity',
    use='dynamic',
    process=shared_processes.PAR
)

Variable(
    name='fdp',
    long_name='Fraction phosphorus dissolved',
    units='Unitless',
    description='Fraction phosphorus dissolved',
    use='dynamic',
    process=shared_processes.fdp
)
