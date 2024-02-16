import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.algae.processes as processes
import clearwater_modules.nsm1.nitrogen.processes as processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

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

Variable(
    name='timestep',
    long_name='timestep',
    units='d',
    description='calculation timestep',
    use='static',
)

#Assume these are static within the one cell calculation but is pulled from the flow model

Variable(
    name='velocity',
    long_name='velocity',
    units='m/s',
    description='Average water velocity in cell',
    use='static',
)

Variable(
    name='flow',
    long_name='flow',
    units='m3/s',
    description='Average flow rate in cell',
    use='static',
)

Variable(
    name='topwidth',
    long_name='topwidth',
    units='m',
    description='Average topwidth of cell',
    use='static',
)

#TODO find units for slope
Variable(
    name='slope',
    long_name='slope',
    units='TODO',
    description='Average slope of bottom surface',
    use='static',
)

Variable(
    name='shear_velocity',
    long_name='shear_velocity',
    units='TODO',
    description='Average shear velocity on bottom surface',
    use='static',
)

Variable(
    name='pressure_atm',
    long_name='pressure_atm',
    units='TODO',
    description='atmospheric pressure in atm',
    use='static',
)

Variable(
    name='wind_speed',
    long_name='Wind speed at 10 meters above the water surface',
    units='m/s',
    description='Wind speed at 10 meters above the water surface',
    use='static',
)

Variable(
    name='wind_reaeration_option',
    long_name='wind reaeration option',
    units='unitless',
    description='Integer value which selects method for computing wind oxygen reaeration velocity',
    use='static',
)

#TODO figure out what Solid is
Variable(
    name='Solid',
    long_name='Solid reaeration option',
    units='Unknown',
    description='Solid',
    use='static',
)
