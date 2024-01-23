"""
File includes static variables only used in Nitrogen module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.nitrogen.nitrogen_processes as nitrogen_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

# Only Nitrogen Variables


Variable(
    name='KNR',
    long_name='Oxygen inhabitation factor for nitrification',
    units='mg-O2/L',
    description='Oxygen inhabitation factor for nitrification',
    use='static',
)

Variable(
    name='knit_20',
    long_name='Nitrification Rate Ammonia decay at 20C',
    units='1/d',
    description='Nitrification Rate Ammonia NH4 -> NO3 decay at 20C',
    use='static',
)

Variable(
    name='kon_20',
    long_name='Decay Rate of OrgN to NH4 at 20C',
    units='1/d',
    description='Decay Rate of OrgN to NH4 at 20C',
    use='static',
)

Variable(
    name='kdnit_20',
    long_name='Denitrification rate at 20C',
    units='1/d',
    description='Denitrification rate at 20C',
    use='static',
)

Variable(
    name='rnh4_20',
    long_name='Sediment release rate of NH4 at 20C',
    units='g-N/m^2/d',
    description='Sediment release rate of NH4 at 20C',
    use='static'
)

Variable(
    name='vno3_20',
    long_name='Sediment denitrification velocity at 20C',
    units='m/d',
    description='Sediment denitrification velocity at 20C',
    use='static',
)

Variable(
    name='KsOxdn',
    long_name='Half-saturation oxygen inhibition constant for denitrification',
    units='mg-O2/L',
    description='Half-saturation oxygen inhibition constant for denitrification',
    use='static',
)


Variable(
    name='PN',
    long_name='NH4 preference factor algae',
    units='unitless',
    description='NH4 preference factor algae (1=full NH4, 0=full NO3)',
    use='static',
)

Variable(
    name='PNb',
    long_name='NH4 preference factor benthic algae',
    units='unitless',
    description='NH4 preference factor benthic algae (1=full NH4, 0=full NO3)',
    use='static',
)
