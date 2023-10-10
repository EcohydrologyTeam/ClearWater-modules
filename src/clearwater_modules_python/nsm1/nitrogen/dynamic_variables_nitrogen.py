"""
File includes dynamic variables computed in Algae module. Dynamic variables may be accessed by other modules.
"""

import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.nitrogen.nitrogen_processes as nitrogen_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...


Variable(
    name='knit_tc',
    long_name='Nitrification rate ammonia decay',
    units='1/d',
    description='Nitrification rate ammonia decay temperature correction',
    use='dynamic',
    process=nitrogen_processes.knit_tc
)

Variable(
    name='rnh4_tc',
    long_name='Sediment release rate of NH4',
    units='1/d',
    description=' Sediment release rate of NH4 temperature correction',
    use='dynamic',
    process=nitrogen_processes.rnh4_tc
)

Variable(
    name='vno3_tc',
    long_name='Sediment denitrification velocity',
    units='m/d',
    description='Sediment denitrification velocity temperature correction',
    use='dynamic',
    process=nitrogen_processes.vno3_tc
)

Variable(
    name='kon_tc',
    long_name='Decay rate of OrgN to NH4',
    units='1/d',
    description='Decay rate of OrgN to NH4 temperature correction',
    use='dynamic',
    process=nitrogen_processes.kon_tc
)

Variable(
    name='kdnit_tc',
    long_name='Denitrification rate',
    units='1/d',
    description='Denitrification rate temperature correction',
    use='dynamic',
    process=nitrogen_processes.kdnit_tc
)

Variable(
    name='ApUptakeFr_NH4',
    long_name='Fraction of actual floating algal uptake from ammonia pool',
    units='unitless',
    description='Fraction of actual floating algal uptake from ammonia pool',
    use='dynamic',
    process=nitrogen_processes.ApUptakeFr_NH4
)

Variable(
    name='ApUptakeFr_NO3',
    long_name='Fraction of actual floating algal uptake from nitrate pool',
    units='unitless',
    description='Fraction of actual floating algal uptake from nitrate pool',
    use='dynamic',
    process=nitrogen_processes.ApUptakeFr_NO3
)

Variable(
    name='AbUptakeFr_NH4',
    long_name='Fraction of actual benthic algal uptake from ammonia pool',
    units='unitless',
    description='Fraction of actual benthic algal uptake from ammonia pool',
    use='dynamic',
    process=nitrogen_processes.AbUptakeFr_NH4
)

Variable(
    name='AbUptakeFr_NO3',
    long_name='Fraction of actual benthic algal uptake from nitrate pool',
    units='unitless',
    description='Fraction of actual benthic algal uptake from nitrate pool',
    use='dynamic',
    process=nitrogen_processes.AbUptakeFr_NO3
)

Variable(
    name='dOrgNdt',
    long_name='Change in organic nitrogen',
    units='mg-N/L',
    description='Change in organic nitrogen',
    use='dynamic',
    process=nitrogen_processes.dOrgNdt
)

Variable(
    name='dNH4dt',
    long_name='Change in ammonium concentration',
    units='mg-N/L',
    description='Change in ammonium concentration',
    use='dynamic',
    process=nitrogen_processes.dNH4dt
)

Variable(
    name='dNO3dt',
    long_name='Change in nitrate concentration',
    units='mg-N/L',
    description='Change in nitrate concentration',
    use='dynamic',
    process=nitrogen_processes.dNO3dt
)

Variable(
    name='DIN',
    long_name='Dissolve inorganic nitrogen',
    units='mg-N/L',
    description='Dissolve inorganic nitrogen',
    use='dynamic',
    process=nitrogen_processes.DIN
)

Variable(
    name='TON',
    long_name='Total organic nitrogen',
    units='mg-N/L',
    description='Total organic nitrogen',
    use='dynamic',
    process=nitrogen_processes.TON
)

Variable(
    name='TKN',
    long_name='Total kjeldhl nitrogen',
    units='mg-N/L',
    description='Total kjeldhl nitrogen',
    use='dynamic',
    process=nitrogen_processes.TKN
)

Variable(
    name='TN',
    long_name='Total nitrogen',
    units='mg-N/L',
    description='Total nitrogen',
    use='dynamic',
    process=nitrogen_processes.TN
)
