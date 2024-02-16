"""
File contains dynamic variables related to the alkalinity module
"""

import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
from clearwater_modules.nsm1.alkalinity import processes

@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

Variable(
    name='Alk_denitrification',
    long_name='Alkalinity change due to denitrification',
    units='mg/L/d',
    description='Alkalinity change due to denitrification',
    use='dynamic',
    process=processes.Alk_denitrification
)

Variable(
    name='Alk_nitrification',
    long_name='Alkalinity change due to nitrification',
    units='mg/L/d',
    description='Alkalinity change due to nitrification',
    use='dynamic',
    process=processes.Alk_nitrification
)

Variable(
    name='Alk_algal_growth',
    long_name='Alkalinity change due to algal growth',
    units='mg/L/d',
    description='Alkalinity change due to algal growth',
    use='dynamic',
    process=processes.Alk_algal_growth
)

Variable(
    name='Alk_algal_respiration',
    long_name='Alkalinity change due to algal respiration',
    units='mg/L/d',
    description='Alkalinity change due to algal respiration',
    use='dynamic',
    process=processes.Alk_algal_respiration
)

Variable(
    name='Alk_benthic_algae_growth',
    long_name='Alkalinity change due to benthic algae growth',
    units='mg/L/d',
    description='Alkalinity change due to benthic algae growth',
    use='dynamic',
    process=processes.Alk_benthic_algae_growth
)

Variable(
    name='Alk_benthic_algae_respiration',
    long_name='Alkalinity change due to benthic algae growth',
    units='mg/L/d',
    description='Alkalinity change due to benthic algae growth',
    use='dynamic',
    process=processes.Alk_benthic_algae_respiration
)

Variable(
    name='dAlkdt',
    long_name='Alkalinity concentration change per timestep',
    units='mg/L/d',
    description='Alkalinity concentration change per timestep',
    use='dynamic',
    process=processes.dAlkdt
)



