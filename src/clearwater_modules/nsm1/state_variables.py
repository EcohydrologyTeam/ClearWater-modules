from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.algae.processes as algae_processes
import clearwater_modules.nsm1.alkalinity.processes as alkalinity_processes
import clearwater_modules.nsm1.balgae.processes as balgae_processes
import clearwater_modules.nsm1.carbon.processes as carbon_processes
import clearwater_modules.nsm1.CBOD.processes as CBOD_processes
import clearwater_modules.nsm1.DOX.processes as DOX_processes
import clearwater_modules.nsm1.n2.processes as n2_processes
import clearwater_modules.nsm1.nitrogen.processes as nitrogen_processes
import clearwater_modules.nsm1.pathogens.processes as pathogens_processes
import clearwater_modules.nsm1.phosphorus.processes as phosphorus_processes
import clearwater_modules.nsm1.POM.processes as POM_processes
import clearwater_modules.shared.processes as shared_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    """NSM1 state variables."""
    ...

# TODO: remove mock_equation
def mock_equation(water_temp_c: float) -> float:
    return water_temp_c ** 2

# TODO: import state variables from CWR such as surface_area volume, and timestep, as well as kah inputs

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
    process=balgae_processes.Ab
)

Variable(
    name='NH4',
    long_name='Ammonium Concentration',
    units='mg-N/L',
    description='Ammonium Concentration',
    use='state',
    process=nitrogen_processes.NH4
)

Variable(
    name='NO3',
    long_name='Nitrate Concentration',
    units='mg-N/L',
    description='Nitrate Concentration',
    use='state',
    process=nitrogen_processes.NO3
)

Variable(
    name='OrgN',
    long_name='Organic Nitrogen Concentration',
    units='mg-N/L',
    description='Organic Nitrogen Concentration',
    use='state',
    process=nitrogen_processes.OrgN
)

Variable(
    name='N2',
    long_name='Nitrogen concentration air',
    units='mg-N/L',
    description='Nitrogen concentration air',
    use='state',
    process=n2_processes.N2
)

Variable(
    name='TIP',
    long_name='Total Inorganic Phosphorus',
    units='mg-P/L',
    description='Total Inorganic Phosphorus Concentration',
    use='state',
    process=phosphorus_processes.TIP
)

Variable(
    name='OrgP',
    long_name='Total Organic Phosphorus',
    units='mg-P/L',
    description='Total Organic Phosphorus Concentration',
    use='state',
    process=phosphorus_processes.OrgP
)

Variable(
    name='POC',
    long_name='Particulate Organic Carbon',
    units='mg-C/L',
    description='Particulate Organic Carbon Concentration',
    use='state',
    process=carbon_processes.POC
)

Variable(
    name='DOC',
    long_name='Dissolved Organic Carbon',
    units='mg-C/L',
    description='Dissolved Organic Carbon Concentration',
    use='state',
    process=carbon_processes.DOC
)

Variable(
    name='DIC',
    long_name='Dissolved Inorganic Carbon',
    units='mg-C/L',
    description='Dissolved Inorganic Carbon Concentration',
    use='state',
    process=carbon_processes.DIC
)

Variable(
    name='POM',
    long_name='Particulate Organic Matter',
    units='mg-D/L',
    description='Particulate Organic Matter Concentration',
    use='state',
    process=POM_processes.POM
)

Variable(
    name='POM2',
    long_name='Sediment Particulate Organic Matter',
    units='mg-D/L',
    description='Sediment Particulate Organic Matter Concentration',
    use='state',
    process=mock_equation#TODO might be Sedflux
)

Variable(
    name='CBOD',
    long_name='Carbonaceous Biochemical Oxygen Demand',
    units='mg-O2/L',
    description='Carbonaceous Biochemical Oxygen Demand Concentration',
    use='state',
    process=CBOD_processes.CBOD
)

Variable(
    name='DOX',
    long_name='Dissolved Oxygen',
    units='mg-O2/L',
    description='Dissolved Oxygen',
    use='state',
    process=DOX_processes.DOX
)

Variable(
    name='PX',
    long_name='Pathogen',
    units='cfu/100mL',
    description='Pathogen concentration',
    use='state',
    process=pathogens_processes.PX
)

Variable(
    name='Alk',
    long_name='Alkalinity',
    units='mg-CaCO3/L',
    description='Alkalinity concentration',
    use='state',
    process=alkalinity_processes.Alk 
)