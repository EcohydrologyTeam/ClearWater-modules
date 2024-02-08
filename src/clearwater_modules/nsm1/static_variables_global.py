import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.algae.algae_processes as algae_processes
import clearwater_modules.nsm1.nitrogen.nitrogen_processes as nitrogen_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

# Global parameters



Variable(
    name='fdp',
    long_name='Fraction P dissolved',
    units='unitless',
    description='Fraction P dissolved',
    use='static',
)

Variable(
    name='PAR',
    long_name='Surface light intensity',
    units='W/m^2',
    description='Surface light intensity',
    use='static',
)

Variable(
    name='vson',
    long_name='Organic N settling velocity',
    units='m/d',
    description='Organic N settling velocity',
    use='static',
)

Variable(
    name='vsoc',
    long_name='POC settling velocity',
    units='m/d',
    description='POC settling velocity',
    use='static'
)

Variable(
    name='SOD_20',
    long_name='Sediment oxygen demand at 20 degrees C',
    units='g-O2/m/d',
    description='Sediment oxygen demand at 20 degrees C',
    use='static'
)

Variable(
    name='SOD_theta',
    long_name='Arrhenius coefficient for sediment oxygen demand',
    units='unitless',
    description='Arrhenius coefficient for sediment oxygen demand',
    use='static'
)

# Global module options

Variable(
    name='use_NH4',
    long_name='Use ammonium module',
    units='unitless',
    description='True/False use ammonium module',
    use='static',
)

Variable(
    name='use_NO3',
    long_name='Use nitrate module',
    units='unitless',
    description='True/False use nitrate module',
    use='static',
)

Variable(
    name='use_OrgN',
    long_name='Use organic nitrogen module',
    units='unitless',
    description='True/False use organic nitrogen module',
    use='static',
)

Variable(
    name='use_SedFlux',
    long_name='Use sediment flux module',
    units='unitless',
    description='True/False use sediment flux module',
    use='static',
)

Variable(
    name='use_DOX',
    long_name='Use dissolved oxygen module',
    units='unitless',
    description='True/False use dissolved oxygen module',
    use='static',
)

Variable(
    name='use_Algae',
    long_name='Use algae module',
    units='unitless',
    description='True/False use algae module',
    use='static',
)

Variable(
    name='use_Balgae',
    long_name='Use benthic algae module',
    units='unitless',
    description='True/False use benthic algae module',
    use='static',
)

Variable(
    name='use_TIP',
    long_name='Use total inorganic phosphorus module',
    units='unitless',
    description='True/False use total inorganic phosphorus module',
    use='static',
)

Variable(
    name='use_OrgP',
    long_name='Use total organic phosphorus module',
    units='unitless',
    description='True/False use total organic phosphorus module',
    use='static',
)

Variable(
    name='use_POC',
    long_name='Use particulate organic carbon module',
    units='unitless',
    description='True/False use particulate organic carbon module',
    use='static',
)

Variable(
    name='use_DOC',
    long_name='Use dissolved organic carbon module',
    units='unitless',
    description='True/False use dissolved organic carbon module',
    use='static',
)

Variable(
    name='use_DIC',
    long_name='Use dissolved inorganic carbon module',
    units='unitless',
    description='True/False use dissolved inorganic carbon module',
    use='static',
)

Variable(
    name='use_N2',
    long_name='Use dissolved N2 module',
    units='unitless',
    description='True/False use N2 module',
    use='static',
)

Variable(
    name='use_Pathogen',
    long_name='Use pathogen module',
    units='unitless',
    description='True/False use pathogen module',
    use='static',
)

Variable(
    name='use_Alk',
    long_name='Use alkalinity module',
    units='unitless',
    description='True/False use alkalinity module',
    use='static',
)

Variable(
    name='use_POM2',
    long_name='Use particulate organic matter module',
    units='unitless',
    description='True/False use particulate organic matter module',
    use='static',
)
