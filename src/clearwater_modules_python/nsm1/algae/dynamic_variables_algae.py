"""
File includes dynamic variables computed in Algae module. Dynamic variables may be accessed by other modules.
"""

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.algae.algae_processes as algae_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

Variable(
    name='rna',
    long_name='Algal N:Chla Ratio',
    units='mg-N/ug Chla',
    description='Algal N:Chla Ratio',
    use='dynamic',
    process=algae_processes.rna
)

Variable(
    name='rpa',
    long_name='Algal P:Chla Ratio',
    units='mg-P/ug Chla',
    description='Algal P:Chla Ratio',
    use='dynamic',
    process=algae_processes.rpa
)

Variable(
    name='rca',
    long_name='Algal C:Chla Ratio',
    units='mg-C/ug Chla',
    description='Algal C:Chla Ratio',
    use='dynamic',
    process=algae_processes.rca
)

Variable(
    name='rda',
    long_name='Algal D:Chla Ratio',
    units='mg-D/ug Chla',
    description='Algal D:Chla Ratio',
    use='dynamic',
    process=algae_processes.rda
)

Variable(
    name='mu_max_tc', 
    long_name='Max Algae Growth with Temperature Correction',
    units='1/d',
    description='Max Algae Growth with Temperature Correction',
    use='dynamic',
    process=algae_processes.mu_max_tc,
)

Variable(
    name='krp_tc', 
    long_name='Algal Respiration Rate with Temperature Correction',
    units='1/d',
    description='Algal Respiration Rate with Temperature Correction',
    use='dynamic',
    process=algae_processes.krp_tc,
)

Variable(
    name='kdp_tc', 
    long_name='Algal Mortality Rate with Temperature Correction',
    units='1/d',
    description='Algal Mortality Rate with Temperature Correction',
    use='dynamic',
    process=algae_processes.kdp_tc,
)

Variable(
    name='FL', 
    long_name='Algal Light Limitation',
    units='unitless',
    description='Algal Light Limitation',
    use='dynamic',
    process=algae_processes.FL,
)

Variable(
    name='FN', 
    long_name='Algal Nitrogen Limitation',
    units='unitless',
    description='Algal Nitrogen Limitation',
    use='dynamic',
    process=algae_processes.FN,
)

Variable(
    name='FP', 
    long_name='Algal Phosphorus Limitation',
    units='unitless',
    description='Algal Phosphorus Limitation',
    use='dynamic',
    process=algae_processes.FP,
)

Variable(
    name='mu', 
    long_name='Algal Growth Rate',
    units='1/d',
    description='Algal Growth Rate',
    use='dynamic',
    process=algae_processes.mu,
)

Variable(
    name='ApGrowth', 
    long_name='Algal Growth',
    units='ug-Chala/L/d',
    description='Algal Growth',
    use='dynamic',
    process=algae_processes.ApGrowth,
)

Variable(
    name='ApRespiration', 
    long_name='Algal Respiration',
    units='ug-Chala/L/d',
    description='Algal Respiration',
    use='dynamic',
    process=algae_processes.ApRespiration,
)

Variable(
    name='ApDeath', 
    long_name='Algal Death',
    units='ug-Chala/L/d',
    description='Algal Death',
    use='dynamic',
    process=algae_processes.ApDeath,
)

Variable(
    name='ApSettling', 
    long_name='Algal Settling',
    units='ug-Chala/L/d',
    description='Algal Settling',
    use='dynamic',
    process=algae_processes.ApSettling,
)

Variable(
    name='dApdt',
    long_name='Algal Biomass Concentration Change',
    units='ug-Chala/L/d',
    description='Algal Biomass Concentration Change',
    use='dynamic',
    process=algae_processes.dApdt,
)

Variable(
    name='Ap_new',
    long_name='New algal biomass concentration',
    units='ug-Chala/L/d',
    description='New algal biomass concentration',
    use='dynamic',
    process=algae_processes.Ap_new,
)