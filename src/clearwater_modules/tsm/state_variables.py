from clearwater_modules import base
from clearwater_modules.tsm.model import EnergyBudget
from clearwater_modules.tsm import processes


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    """TSM state variables."""
    ...


Variable(
    name='water_temp_c',
    long_name='Water temperature',
    units='degC',
    description='TSM state variable for water temperature',
    use='state',
    process=processes.t_water_c,
)

# TODO: remove mock_equation


def mock_equation(water_temp_c: float) -> float:
    return water_temp_c ** 2


Variable(
    name='surface_area',
    long_name='Surface area',
    units='m^2',
    description='Surface area',
    use='state',
    process=mock_equation,
)
Variable(
    name='volume',
    long_name='Volume',
    units='m^3',
    description='Volume',
    use='state',
    process=mock_equation,
)
