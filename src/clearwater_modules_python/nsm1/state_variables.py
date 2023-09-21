from clearwater_modules_python import base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.algae.algae_processes as algae_processes
import clearwater_modules_python.nsm1.nitrogen.nitrogen_processes as nitrogen_processes
import clearwater_modules_python.shared.processes as shared_processes
import clearwater_modules_python.tsm.processes as tsm_processes

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
    process=algae_processes.Ap 
)

Variable(
    name='NH4',
    long_name='Ammonium Concentration',
    units='mg-N/L',
    description='Ammonium Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with nitrogen module
)

Variable(
    name='NO3',
    long_name='Nitrate Concentration',
    units='mg-N/L',
    description='Nitrate Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with nitrogen module
)

Variable(
    name='TIP',
    long_name='Total Inorganic Phosphrous',
    units='mg-P/L',
    description='Total Inorganic Phosphrous Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with phosphorous module
)

Variable(
    name='DOX',
    long_name='Dissolved Oxygen',
    units='mg-O2/L',
    description='Dissolved Oxygen',
    use='state',
    process=mock_equation #TODO this variable only changes with DOX module
)

Variable(
    name='OrgN',
    long_name='Organic Nitrogen Concentration',
    units='mg-N/L',
    description='Organic Nitrogen Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with DOX module
)

#TODO not sure the order of calling with tsm
Variable(
    name='TwaterC',
    long_name='Water Temperature',
    units='C',
    description='Water Temperature Degree Celcius',
    use='state',
    process=tsm_processes.t_water_c 
)

Variable(
    name='depth',
    long_name='Water Depth',
    units='m',
    description='Water depth from surface',
    use='state',
    process=shared_processes.depth_calc
)