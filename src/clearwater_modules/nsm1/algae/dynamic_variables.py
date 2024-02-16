"""
File contains dynamic variables related to the Algae module
"""

import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.algae.processes as processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

Variable(
    name='L',
    long_name='Light attenuation coefficient',
    units='unitless',
    description='Light attenuation coefficient',
    use='dynamic',
    process=processes.L
)

Variable(
    name='PAR',
    long_name='surface light intensity',
    units='W/m2',
    description='surface light intensity',
    use='dynamic',
    process=processes.PAR
)

Variable(
    name='rna',
    long_name='Algal N:Chla Ratio',
    units='mg-N/ug Chla',
    description='Algal N:Chla Ratio',
    use='dynamic',
    process=processes.rna
)

Variable(
    name='rpa',
    long_name='Algal P:Chla Ratio',
    units='mg-P/ug Chla',
    description='Algal P:Chla Ratio',
    use='dynamic',
    process=processes.rpa
)

Variable(
    name='rca',
    long_name='Algal C:Chla Ratio',
    units='mg-C/ug Chla',
    description='Algal C:Chla Ratio',
    use='dynamic',
    process=processes.rca
)

Variable(
    name='rda',
    long_name='Algal D:Chla Ratio',
    units='mg-D/ug Chla',
    description='Algal D:Chla Ratio',
    use='dynamic',
    process=processes.rda
)

Variable(
    name='mu_max_tc',
    long_name='Max Algae Growth with Temperature Correction',
    units='1/d',
    description='Max Algae Growth with Temperature Correction',
    use='dynamic',
    process=processes.mu_max_tc,
)

Variable(
    name='krp_tc',
    long_name='Algal Respiration Rate with Temperature Correction',
    units='1/d',
    description='Algal Respiration Rate with Temperature Correction',
    use='dynamic',
    process=processes.krp_tc,
)

Variable(
    name='kdp_tc',
    long_name='Algal Mortality Rate with Temperature Correction',
    units='1/d',
    description='Algal Mortality Rate with Temperature Correction',
    use='dynamic',
    process=processes.kdp_tc,
)

Variable(
    name='FL',
    long_name='Algal Light Limitation',
    units='unitless',
    description='Algal Light Limitation',
    use='dynamic',
    process=processes.FL,
)

Variable(
    name='FN',
    long_name='Algal Nitrogen Limitation',
    units='unitless',
    description='Algal Nitrogen Limitation',
    use='dynamic',
    process=processes.FN,
)

Variable(
    name='FP',
    long_name='Algal Phosphorus Limitation',
    units='unitless',
    description='Algal Phosphorus Limitation',
    use='dynamic',
    process=processes.FP,
)

Variable(
    name='mu',
    long_name='Algal Growth Rate',
    units='1/d',
    description='Algal Growth Rate',
    use='dynamic',
    process=processes.mu,
)

Variable(
    name='ApGrowth',
    long_name='Algal Growth',
    units='ug-Chala/L/d',
    description='Algal Growth',
    use='dynamic',
    process=processes.ApGrowth,
)

Variable(
    name='ApRespiration',
    long_name='Algal Respiration',
    units='ug-Chala/L/d',
    description='Algal Respiration',
    use='dynamic',
    process=processes.ApRespiration,
)

Variable(
    name='ApDeath',
    long_name='Algal Death',
    units='ug-Chala/L/d',
    description='Algal Death',
    use='dynamic',
    process=processes.ApDeath,
)

Variable(
    name='ApSettling',
    long_name='Algal Settling',
    units='ug-Chala/L/d',
    description='Algal Settling',
    use='dynamic',
    process=processes.ApSettling,
)

Variable(
    name='dApdt',
    long_name='Algal Biomass Concentration Change',
    units='ug-Chala/L/d',
    description='Algal Biomass Concentration Change',
    use='dynamic',
    process=processes.dApdt,
)

Variable(
    name='Ap',
    long_name='New algal biomass concentration',
    units='ug-Chala/L/d',
    description='New algal biomass concentration',
    use='dynamic',
    process=processes.Ap,
)
