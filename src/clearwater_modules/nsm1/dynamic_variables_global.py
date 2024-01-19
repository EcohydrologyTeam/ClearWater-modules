import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget
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
    name='SOD_tc',
    long_name='Sediment oxygen demand adjusted for temperature',
    units='mg-O2/L/d',
    description='Sediment oxygen demand adjusted for temperature',
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
    name='kah_T',
    long_name='Hydraulic oxygen reaeration rate adjusted for temperature',
    units='1/d',
    description='Hydraulic oxygen reaeration rate adjusted for temperature',
    use='dynamic',
    process=shared_processes.kah_T
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
    name='kaw_T',
    long_name='Wind oxygen reaeration velocity adjusted for temperature',
    units='m/d',
    description='Wind oxygen reaeration velocity adjusted for temperature',
    use='dynamic',
    process=shared_processes.kaw_T
)

Variable(
    name='ka_T',
    long_name='Oxygen reaeration rate',
    units='1/d',
    description='Oxygen reaeration rate',
    use='dynamic',
    process=shared_processes.ka_T
)