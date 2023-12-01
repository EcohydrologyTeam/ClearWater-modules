"""
File includes static variables only used in Algae module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.algae.algae_processes as algae_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

# Only Algae Variables


Variable(
    name='Awd',
    long_name='Algal Dry Weight',
    units='mg',
    description='Algal Dry Weight',
    use='static',
)

Variable(
    name='AWc',
    long_name='Carbon Weight',
    units='mg',
    description='Carbon Weight',
    use='static',
)

Variable(
    name='AWn',
    long_name='Nitrogen Weight',
    units='mg',
    description='Nitrogen Weight',
    use='static',
)

Variable(
    name='AWp',
    long_name='Phosphorus Weight',
    units='mg',
    description='Phosphorus Weight',
    use='static',
)

Variable(
    name='AWa',
    long_name='Algal Chlorophyll',
    units='ug Chla',
    description='Algal Chlorophyll',
    use='static',
)


Variable(
    name='KL',
    long_name='Light Limiting Constant for Algal Growth',
    units='W/m^2',
    description='Light Limiting Constant for Algal Growth',
    use='static',
)

Variable(
    name='KsN',
    long_name='Half-Saturation N Limiting Constant for Algal Growth',
    units='mg-N/L',
    description='Half-Saturation N Limiting Constant for Algal Growth',
    use='static',
)

Variable(
    name='KsP',
    long_name='Half-Saturation P Limiting Constant for Algal Growth',
    units='mg-P/L',
    description='Half-Saturation P Limiting Constant for Algal Growth',
    use='static',
)

Variable(
    name='mu_max_20',
    long_name='Max Algae Growth',
    units='1/d',
    description='Max Algae Growth at 20C',
    use='static',
)

Variable(
    name='kdp_20',
    long_name='Algal Mortality Rate',
    units='1/d',
    description='Algal Mortality Rate at 20C',
    use='static',
)

Variable(
    name='krp_20',
    long_name='Algal Respiration Rate',
    units='1/d',
    description='Algal Respiration Rate at 20C',
    use='static',
)

Variable(
    name='vsap',
    long_name='Algal Setting Velocity',
    units='m/d',
    description='Algal Setting Velocity',
    use='static',
)

Variable(
    name='growth_rate_option',
    long_name='Growth Rate Option',
    units='1/d',
    description='Algal growth rate option 1) multiplicative, 2) Limiting Nutrient, 3) Harmonic Mean Option',
    use='static',
)

Variable(
    name='growth_rate_option',
    long_name='Growth Rate Option',
    units='1/d',
    description='Algal growth rate option 1) multiplicative, 2) Limiting Nutrient, 3) Harmonic Mean Option',
    use='static',
)

Variable(
    name='light_limitation_option',
    long_name='Light Limitation Option',
    units='1/d',
    description='Algal light limitation 1) half-saturation, 2) Smith model, 3) Steele model',
    use='static',
)
