import clearwater_modules_python.base as base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.algae.algae_processes as algae_processes
import clearwater_modules_python.nsm1.nitrogen.nitrogen_processes as nitrogen_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

#Global parameters

Variable(
    name='lambda',
    long_name='Light attenuation coefficient',
    units='unitless',
    description='Light attenuation coefficient',
    use='static',
)

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
#Global module options

Variable(
    name='use_NH4',
    long_name='Use ammonium module',
    units='unitless',
    description='True/Fasle use ammonium module',
    use='static',
)

Variable(
    name='use_NO3',
    long_name='Use nitrate module',
    units='unitless',
    description='True/Fasle use nitrate module',
    use='static',
)

Variable(
    name='use_OrgN',
    long_name='Use organic nitrogen module',
    units='unitless',
    description='True/Fasle use organic nitrogen module',
    use='static',
)

Variable(
    name='use_SedFlux',
    long_name='Use sediment flux module',
    units='unitless',
    description='True/Fasle use sediment flux module',
    use='static',
)

Variable(
    name='use_DOX',
    long_name='Use dissolved oxygen module',
    units='unitless',
    description='True/Fasle use dissovled oxygen module',
    use='static',
)

Variable(
    name='use_Algae',
    long_name='Use algae module',
    units='unitless',
    description='True/Fasle use algae module',
    use='static',
)

Variable(
    name='use_Balgae',
    long_name='Use benthic algae module',
    units='unitless',
    description='True/Fasle use benthic algae module',
    use='static',
)

Variable(
    name='use_TIP',
    long_name='Use total organic phosphrous module',
    units='unitless',
    description='True/Fasle use total organic phosphrous module',
    use='static',
)

