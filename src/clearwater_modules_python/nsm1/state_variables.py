from clearwater_modules_python import base
from clearwater_modules_python.nsm1.model import NutrientBudget
from clearwater_modules_python.nsm1 import algae_processes
from clearwater_modules_python.tsm import processes

@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    """NSM1 state variables."""
    ...

# TODO: remove mock_equation
def mock_equation(water_temp_c: float) -> float:
    return water_temp_c ** 2

Variable(
    name='Ap',
    long_name='Algae Concentration',
    units='ug-Chla/L',
    description='Algal Concentration',
    use='state',
    process=mock_equation
)

Variable(
    name='NH4',
    long_name='Ammonium Concentration',
    units='mg-N/L',
    description='Ammonium Concentration',
    use='state',
    process=mock_equation
)

Variable(
    name='NO3',
    long_name='Nitrate Concentration',
    units='mg-N/L',
    description='Nitrate Concentration',
    use='state',
    process=mock_equation
)

Variable(
    name='TIP',
    long_name='Total Inorganic Phosphrous',
    units='mg-P/L',
    description='Total Inorganic Phosphrous Concentration',
    use='state',
    process=mock_equation
)

#TODO not sure the order of calling with tsm
Variable(
    name='TwaterC',
    long_name='Water Temperature',
    units='C',
    description='Water Temperature Degree Celcius',
    use='state',
    process=processes.t_water_c
)

