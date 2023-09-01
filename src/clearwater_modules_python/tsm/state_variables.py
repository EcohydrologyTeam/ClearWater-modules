from clearwater_modules_python import base
from clearwater_modules_python.tsm.model import EnergyBudget

@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    """TSM state variables."""
    ...

Variable(
    name='t_water_c',
    long_name='Water temperature',
    units='degC',
    description='TSM state variable for water temperature',
    use='state',
    process=None,
)
Variable(
    name='surface_area',
    long_name='Surface area',
    units='m^2',
    description='Surface area',
    use='state',
    process=None,
)
Variable(
    name='volume',
    long_name='Volume',
    units='m^3',
    description='Volume',
    use='state',
    process=None,
)

