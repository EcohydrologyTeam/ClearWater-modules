import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


Variable(
    name='kdp_20',
    long_name='',
    units='',
    description='',
    use='static'
)

Variable(
    name='kdb_20',
    long_name='',
    units='',
    description='',
    use='static'
)

Variable(
    name='kpom_20',
    long_name='',
    units='',
    description='',
    use='static'
)