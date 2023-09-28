from clearwater_modules_python import base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.algae.algae_processes as algae_processes
import clearwater_modules_python.nsm1.nitrogen.nitrogen_processes as nitrogen_processes
import clearwater_modules_python.nsm1.n2.n2_processes as n2_processes
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
    name='Ab',
    long_name='Benthic Algae Concentration',
    units='g-D/m^2',
    description='Benthic Algae Concentration',
    use='state',
    process=mock_equation #TODO depends on benthic algae module 
)

Variable(
    name='NH4',
    long_name='Ammonium Concentration',
    units='mg-N/L',
    description='Ammonium Concentration',
    use='state',
    process=nitrogen_processes.NH4_new
)

Variable(
    name='NO3',
    long_name='Nitrate Concentration',
    units='mg-N/L',
    description='Nitrate Concentration',
    use='state',
    process=nitrogen_processes.NO3_new
)

Variable(
    name='OrgN',
    long_name='Organic Nitrogen Concentration',
    units='mg-N/L',
    description='Organic Nitrogen Concentration',
    use='state',
    process=nitrogen_processes.OrgN_new
)

Variable(
    name='N2',
    long_name='Nitrogen concentration air',
    units='mg-N/L',
    description='Nitrogen concentration air',
    use='state',
    process=n2_processes.N2_new
)

Variable(
    name='TIP',
    long_name='Total Inorganic Phosphorus',
    units='mg-P/L',
    description='Total Inorganic Phosphorus Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with phosphorous module
)

Variable(
    name='OrgP',
    long_name='Total Organic Phosphorus',
    units='mg-P/L',
    description='Total Organic Phosphorus Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with phosphorous module
)

Variable(
    name='POC',
    long_name='Particulate Organic Carbon',
    units='mg-C/L',
    description='Particulate Organic Carbon Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with carbon module
)

Variable(
    name='DOC',
    long_name='Dissolved Organic Carbon',
    units='mg-C/L',
    description='Dissolved Organic Carbon Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with carbon module
)

Variable(
    name='DIC',
    long_name='Dissolved Inorganic Carbon',
    units='mg-C/L',
    description='Dissolved Inorganic Carbon Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with carbon module
)

Variable(
    name='POM',
    long_name='Particulate Organic Matter',
    units='mg-D/L',
    description='Particulate Organic Matter Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with pom module
)

Variable(
    name='POM2',
    long_name='Sediment Particulate Organic Matter',
    units='mg-D/L',
    description='Sediment Particulate Organic Matter Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with pom module
)

Variable(
    name='CBOD',
    long_name='Carbonaceous Biochemical Oxygen Demand',
    units='mg-O2/L',
    description='Carbonaceous Biochemical Oxygen Demand Concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with cbod module
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
    name='PX',
    long_name='Pathogen',
    units='cfu/100mL',
    description='Pathogen concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with pathogen module
)

Variable(
    name='Alk',
    long_name='Alkalinity',
    units='mg-CaCO3/L',
    description='Alkalinity concentration',
    use='state',
    process=mock_equation #TODO this variable only changes with alkalinity module
)

#TODO not sure the order of calling with tsm
Variable(
    name='TwaterC',
    long_name='Water Temperature',
    units='C',
    description='Water Temperature Degree Celsius',
    use='state',
    process=tsm_processes.t_water_c 
)

Variable(
    name='TwaterK',
    long_name='Water Temperature',
    units='K',
    description='Water Temperature Degree Kelvin',
    use='state',
    process=mock_equation #TODO not sure where to call this from shows up in TSM and N2_processes
)

Variable(
    name='depth',
    long_name='Water Depth',
    units='m',
    description='Water depth from surface',
    use='state',
    process=shared_processes.depth_calc
)