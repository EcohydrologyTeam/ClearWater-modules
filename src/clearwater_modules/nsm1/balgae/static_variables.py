"""
File includes static variables only used in Algae module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.algae.processes as processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

#Only Balgae Variables 

Variable(
    name='Fw',
    long_name='Fraction of benthic algae mortality into water column',
    units='unitless',
    description='Fraction of benthic algae mortality into water column',
    use='static',
)

Variable(
    name='Fb',
    long_name='Fraction of bottom area available for benthic algae',
    units='unitless',
    description='Fraction of bottom area available for benthic algae',
    use='static',
)

Variable(
    name='BWd',
    long_name='Benthic algae dry weight',
    units='unitless',
    description='Benthic algae dry weight',
    use='static',
)

Variable(
    name='BWc',
    long_name='Benthic algae carbon',
    units='unitless',
    description='Benthic algae carbon',
    use='static',
)

Variable(
    name='BWn',
    long_name='Benthic algae nitrogen',
    units='unitless',
    description='Benthic algae nitrogen',
    use='static',
)

Variable(
    name='BWp',
    long_name='Benthic algae phosphorus',
    units='unitless',
    description='Benthic algae phosphorus',
    use='static',
)

Variable(
    name='BWa',
    long_name='Benthic algae Chla',
    units='unitless',
    description='Benthic algae Chla',
    use='static',
)

Variable(
    name='KLb',
    long_name='Light limiting constant for benthic algae growth',
    units='W/m^2',
    description='Light limiting constant for benthic algae growth',
    use='static',
)

Variable(
    name='KsNb',
    long_name='Half-Saturation N limiting constant for Benthic algae',
    units='mg-N/L',
    description='Half-Saturation N limiting constant for Benthic algae',
    use='static',
)

Variable(
    name='KsPb',
    long_name='Half-Saturation P limiting constant for Benthic algae',
    units='mg-P/L',
    description='Half-Saturation P limiting constant for Benthic algae',
    use='static',
)

Variable(
    name='Ksb',
    long_name='Half-Saturation density constant for benthic algae growth',
    units='g-D/m^2',
    description='Half-Saturation density constant for benthic algae growth',
    use='static',
)

Variable(
    name='mub_max_20',
    long_name='Maximum benthic algal growth rate',
    units='1/d',
    description='maximum benthic algal growth rate',
    use='static',
)

Variable(
    name='krb_20',
    long_name='Benthic algal respiration rate',
    units='1/d',
    description='Benthic algal respiration rate',
    use='static',
)

Variable(
    name='kdb_20',
    long_name='Benthic algal mortality rate',
    units='1/d',
    description='Benthic algal mortality rate',
    use='static',
)

Variable(
    name='b_growth_rate_option',
    long_name='Benthic Algal growth rate options',
    units='unitless',
    description='Benthic Algal growth rate with two options: 1) Multiplicative, 2) Limiting Nutritent',
    use='static',
)

Variable(
    name='b_light_limitation_option',
    long_name='Benthic Algal light limitation rate options',
    units='unitless',
    description='Benthic Algal light limitation rate with three options: 1) Half-saturation formulation, 2) Smiths Model, 3) Steeles Model'
    use='static',
)

Variable(
    name='Fb',
    long_name='Fraction of bottom area available for benthic algae growth',
    units='unitless',
    description='Fraction of bottom area available for benthic algae growth',
    use='static'
)

Variable(
    name='Fw',
    long_name='Fraction of benthic algae mortality into water column',
    units='unitless',
    description='Fraction of benthic algae mortality into water column',
    use='static'
)
