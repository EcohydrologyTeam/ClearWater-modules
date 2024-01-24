import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


Variable(
    name='kpom_20',
    long_name='POM dissolution rate at 20C',
    units='1/d',
    description='POM dissolution rate at 20C',
    use='static'
)

