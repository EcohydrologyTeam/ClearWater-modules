import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.carbon.model import CarbonBudget
from clearwater_modules.nsm1.alkalinity import processes


@base.register_variable(models=CarbonBudget)
class Variable(base.Variable):
    ...


Variable(
    name='kdp_T',
    long_name='',
    units='',
    description='',
    use='dynamic',
    process=
)

Variable(
    name='kdb_T',
    long_name='',
    units='',
    description='',
    use='dynamic',
    process=
)

Variable(
    name='kpom_T',
    long_name='',
    units='',
    description='',
    use='dynamic',
    process=
)

Variable(
    name='depth',
    long_name='',
    units='',
    description='',
    use='dynamic',
    process=shared_processes.compute_depth
)

Variable(
    name='',
    long_name='',
    units='',
    description='',
    use='dynamic',
    process=
)

Variable(
    name='',
    long_name='',
    units='',
    description='',
    use='dynamic',
    process=
)