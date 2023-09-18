import clearwater_modules_python.base as base
from clearwater_modules_python.nsm1.model import NutrientBudget
from clearwater_modules_python.nsm1 import algae_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

#Global Variables that are used in Algae
Variable(
    name='depth',
    long_name='Water Depth',
    units='m',
    description='Water depth from surface',
    use='static',
)

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
    name='use_TIP',
    long_name='Use total organic phosphrous module',
    units='unitless',
    description='True/Fasle use total organic phosphrous module',
    use='static',
)

#Only Algae Variables 

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
    name='rna',
    long_name='Algal N:Chla Ratio',
    units='mg-N/ug Chla',
    description='Algal N:Chla Ratio',
    use='static',
    process=algae_processes.rna
)

Variable(
    name='rpa',
    long_name='Algal P:Chla Ratio',
    units='mg-P/ug Chla',
    description='Algal P:Chla Ratio',
    use='static',
    process=algae_processes.rpa
)

Variable(
    name='rca',
    long_name='Algal C:Chla Ratio',
    units='mg-C/ug Chla',
    description='Algal C:Chla Ratio',
    use='static',
    process=algae_processes.rca
)

Variable(
    name='rda',
    long_name='Algal D:Chla Ratio',
    units='mg-D/ug Chla',
    description='Algal D:Chla Ratio',
    use='static',
    process=algae_processes.rda
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
    name='mu_max',
    long_name='Max Algae Growth',
    units='1/d',
    description='Max Algae Growth',
    use='static',
)

Variable(
    name='kdp',
    long_name='Algal Mortality Rate',
    units='1/d',
    description='Algal Mortality Rate',
    use='static',
)

Variable(
    name='krp',
    long_name='Algal Respiration Rate',
    units='1/d',
    description='Algal Respiration Rate',
    use='static',
)

Variable(
    name='vsap',
    long_name='Algal Setting Velocity',
    units='1/d',
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