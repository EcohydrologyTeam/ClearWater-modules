from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.processes as processes


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
    process=processes.Ap 
)

Variable(
    name='Ab',
    long_name='Benthic Algae Concentration',
    units='g-D/m^2',
    description='Benthic Algae Concentration',
    use='state',
    process=processes.Ab
)

Variable(
    name='NH4',
    long_name='Ammonium Concentration',
    units='mg-N/L',
    description='Ammonium Concentration',
    use='state',
    process=processes.NH4
)

Variable(
    name='NO3',
    long_name='Nitrate Concentration',
    units='mg-N/L',
    description='Nitrate Concentration',
    use='state',
    process=processes.NO3
)

Variable(
    name='OrgN',
    long_name='Organic Nitrogen Concentration',
    units='mg-N/L',
    description='Organic Nitrogen Concentration',
    use='state',
    process=processes.OrgN
)

Variable(
    name='N2',
    long_name='Nitrogen concentration air',
    units='mg-N/L',
    description='Nitrogen concentration air',
    use='state',
    process=processes.N2
)

Variable(
    name='TIP',
    long_name='Total Inorganic Phosphorus',
    units='mg-P/L',
    description='Total Inorganic Phosphorus Concentration',
    use='state',
    process=processes.TIP
)

Variable(
    name='OrgP',
    long_name='Total Organic Phosphorus',
    units='mg-P/L',
    description='Total Organic Phosphorus Concentration',
    use='state',
    process=processes.OrgP
)

Variable(
    name='POC',
    long_name='Particulate Organic Carbon',
    units='mg-C/L',
    description='Particulate Organic Carbon Concentration',
    use='state',
    process=processes.POC
)

Variable(
    name='DOC',
    long_name='Dissolved Organic Carbon',
    units='mg-C/L',
    description='Dissolved Organic Carbon Concentration',
    use='state',
    process=processes.DOC
)

Variable(
    name='DIC',
    long_name='Dissolved Inorganic Carbon',
    units='mg-C/L',
    description='Dissolved Inorganic Carbon Concentration',
    use='state',
    process=processes.DIC
)

Variable(
    name='POM',
    long_name='Particulate Organic Matter',
    units='mg-D/L',
    description='Particulate Organic Matter Concentration',
    use='state',
    process=processes.POM
)

Variable(
    name='CBOD',
    long_name='Carbonaceous Biochemical Oxygen Demand',
    units='mg-O2/L',
    description='Carbonaceous Biochemical Oxygen Demand Concentration',
    use='state',
    process=processes.CBOD
)

Variable(
    name='DOX',
    long_name='Dissolved Oxygen',
    units='mg-O2/L',
    description='Dissolved Oxygen',
    use='state',
    process=processes.DOX
)

Variable(
    name='PX',
    long_name='Pathogen',
    units='cfu/100mL',
    description='Pathogen concentration',
    use='state',
    process=processes.PX
)

Variable(
    name='Alk',
    long_name='Alkalinity',
    units='mg-CaCO3/L',
    description='Alkalinity concentration',
    use='state',
    process=processes.Alk 
)