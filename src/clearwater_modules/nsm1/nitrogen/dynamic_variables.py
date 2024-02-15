"""
File includes dynamic variables computed in Algae module. Dynamic variables may be accessed by other modules.
"""

import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.nitrogen.processes as processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...


Variable(
    name='knit_tc',
    long_name='Nitrification rate ammonia decay',
    units='1/d',
    description='Nitrification rate ammonia decay temperature correction',
    use='dynamic',
    process=processes.knit_tc
)

Variable(
    name='rnh4_tc',
    long_name='Sediment release rate of NH4',
    units='1/d',
    description=' Sediment release rate of NH4 temperature correction',
    use='dynamic',
    process=processes.rnh4_tc
)

Variable(
    name='vno3_tc',
    long_name='Sediment denitrification velocity',
    units='m/d',
    description='Sediment denitrification velocity temperature correction',
    use='dynamic',
    process=processes.vno3_tc
)

Variable(
    name='kon_tc',
    long_name='Decay rate of OrgN to NH4',
    units='1/d',
    description='Decay rate of OrgN to NH4 temperature correction',
    use='dynamic',
    process=processes.kon_tc
)

Variable(
    name='kdnit_tc',
    long_name='Denitrification rate',
    units='1/d',
    description='Denitrification rate temperature correction',
    use='dynamic',
    process=processes.kdnit_tc
)

Variable(
    name='ApUptakeFr_NH4',
    long_name='Fraction of actual floating algal uptake from ammonia pool',
    units='unitless',
    description='Fraction of actual floating algal uptake from ammonia pool',
    use='dynamic',
    process=processes.ApUptakeFr_NH4
)

Variable(
    name='ApUptakeFr_NO3',
    long_name='Fraction of actual floating algal uptake from nitrate pool',
    units='unitless',
    description='Fraction of actual floating algal uptake from nitrate pool',
    use='dynamic',
    process=processes.ApUptakeFr_NO3
)

Variable(
    name='AbUptakeFr_NH4',
    long_name='Fraction of actual benthic algal uptake from ammonia pool',
    units='unitless',
    description='Fraction of actual benthic algal uptake from ammonia pool',
    use='dynamic',
    process=processes.AbUptakeFr_NH4
)

Variable(
    name='AbUptakeFr_NO3',
    long_name='Fraction of actual benthic algal uptake from nitrate pool',
    units='unitless',
    description='Fraction of actual benthic algal uptake from nitrate pool',
    use='dynamic',
    process=processes.AbUptakeFr_NO3
)

Variable(
    name='ApDeath_OrgN',
    long_name='Algae -> OrgN',
    units='mg-N/L/d',
    description='Algae conversion to Organic nitrogen',
    use='dynamic',
)

Variable(
    name='AbDeath_OrgN',
    long_name='Benthic Algae -> OrgN',
    units='mg-N/L/d',
    description='Benthic algae conversion to Organic nitrogen',
    use='dynamic',
)

Variable(
    name='OrgN_NH4_Decay',
    long_name='OrgN -> NH4',
    units='mg-N/L/d',
    description='Organic nitrogen to ammonium decay',
    use='dynamic',
)

Variable(
    name='OrgN_Settling',
    long_name='OrgN -> bed',
    units='mg-N/L/d',
    description='Organic nitrogen to bed settling',
    use='dynamic',
)

Variable(
    name='dOrgNdt',
    long_name='Change in organic nitrogen',
    units='mg-N/L',
    description='Change in organic nitrogen',
    use='dynamic',
    process=processes.dOrgNdt
)

Variable(
    name='NH4_Nitrification',
    long_name='NH4 -> NO3  Nitrification',
    units='mg-N/L/d',
    description='NH4 Nitrification',
    use='dynamic',
)

Variable(
    name='NH4fromBed',
    long_name='bed ->  NH4 (diffusion)',
    units='mg-N/L/d',
    description='Sediment bed release of NH4',
    use='dynamic',
)

Variable(
    name='NH4_ApRespiration',
    long_name='Floating algae -> NH4',
    units='mg-N/L/d',
    description='Floating algae to NH4',
    use='dynamic',
)

Variable(
    name='NH4_ApGrowth',
    long_name='NH4 -> Floating algae',
    units='mg-N/L/d',
    description='NH4 uptake to algae',
    use='dynamic',
)

Variable(
    name='NH4_AbRespiration',
    long_name='Benthic algae -> NH4',
    units='mg-N/L/d',
    description='Benthic algae release of NH4',
    use='dynamic',
)

Variable(
    name='NH4_AbGrowth',
    long_name='NH4 -> Benthic Algae',
    units='mg-N/L/d',
    description='Benthic algae uptake of NH4',
    use='dynamic',
)

Variable(
    name='dNH4dt',
    long_name='Change in ammonium concentration',
    units='mg-N/L',
    description='Change in ammonium concentration',
    use='dynamic',
    process=processes.dNH4dt
)

Variable(
    name='NO3_Denit',
    long_name='NO3 -> Loss',
    units='mg-N/L/d',
    description='NO3 loss from denitrification',
    use='dynamic',
)

Variable(
    name='NO3_BedDenit',
    long_name='Sediment denitrification',
    units='mg-N/L/d',
    description='Sediment denitrification',
    use='dynamic',
)

Variable(
    name='NO3_ApGrowth',
    long_name='NO3 -> Floating algae',
    units='mg-N/L/d',
    description='NO3 uptake to floating algae',
    use='dynamic',
)

Variable(
    name='NO3_AbGrowth',
    long_name='NO3 -> Benthic algae',
    units='mg-N/L/d',
    description='NO3 uptake to benthic algae',
    use='dynamic',
)

Variable(
    name='dNO3dt',
    long_name='Change in nitrate concentration',
    units='mg-N/L',
    description='Change in nitrate concentration',
    use='dynamic',
    process=processes.dNO3dt
)

Variable(
    name='DIN',
    long_name='Dissolve inorganic nitrogen',
    units='mg-N/L',
    description='Dissolve inorganic nitrogen',
    use='dynamic',
    process=processes.DIN
)

Variable(
    name='TON',
    long_name='Total organic nitrogen',
    units='mg-N/L',
    description='Total organic nitrogen',
    use='dynamic',
    process=processes.TON
)

Variable(
    name='TKN',
    long_name='Total kjeldhl nitrogen',
    units='mg-N/L',
    description='Total kjeldhl nitrogen',
    use='dynamic',
    process=processes.TKN
)

Variable(
    name='TN',
    long_name='Total nitrogen',
    units='mg-N/L',
    description='Total nitrogen',
    use='dynamic',
    process=processes.TN
)
