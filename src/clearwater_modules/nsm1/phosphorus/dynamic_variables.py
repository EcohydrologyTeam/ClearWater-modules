"""
File includes dynamic variables computed in phosphorus module. Dynamic variables may be accessed by other modules.
"""

import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.phosphorus.processes as processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

Variable(
    name='kop_tc',
    long_name='Decay rate of organic P to DIP',
    units='1/d',
    description='Decay rate of organic P to DIP temperature correction',
    use='dynamic',
    process=processes.kop_tc
)

Variable(
    name='rpo4_tc',
    long_name='Benthic sediment release rate of DIP',
    units='g-P/m2/d',
    description='Benthic sediment release rate of DIP temperature correction',
    use='dynamic',
    process=processes.rpo4_tc
)

Variable(
    name='OrgP_DIP_decay',
    long_name='Organic phosphorus decay to dissolve inorganic phosphorus',
    units='mg-P/L/d',
    description='Organic phosphorus decay to dissolve inorganic phosphorus',
    use='dynamic',
    process=processes.OrgP_DIP_decay
)

Variable(
    name='OrgP_Settling',
    long_name='Organic phosphorus settling to sediment',
    units='mg-P/L/d',
    description='Organic phosphorus settling to sediment',
    use='dynamic',
    process=processes.OrgP_Settling
)

Variable(
    name='ApDeath_OrgP',
    long_name='Algal death turning into organic phosphorus',
    units='mg-P/L/d',
    description='Algal death turning into organic phosphorus',
    use='dynamic',
    process=processes.ApDeath_OrgP
)

Variable(
    name='AbDeath_OrgP',
    long_name='Benthic algal death turning into organic phosphorus',
    units='mg-P/L/d',
    description='Benthic algal death turning into organic phosphorus',
    use='dynamic',
    process=processes.AbDeath_OrgP
)

Variable(
    name='dOrgPdt',
    long_name='Change in organic phosphorus concentration',
    units='mg-P/L/d',
    description='Change in organic phosphorus concentration',
    use='dynamic',
    process=processes.dOrgPdt
)

Variable(
    name='DIPfromBed_SedFlux',
    long_name='Dissolved Organic Phosphorus coming from Bed calculated using SedFlux modules',
    units='mg-P/L/d',
    description='Dissolved Organic Phosphorus coming from Bed calculated using SedFlux modules',
    use='dynamic',
    process=processes.DIPfromBed_SedFlux
)

Variable(
    name='DIPfromBed_NoSedFlux',
    long_name='Dissolved Organic Phosphorus coming from Bed calculated without SedFlux modules',
    units='mg-P/L/d',
    description='Dissolved Organic Phosphorus coming from Bed calculated without SedFlux modules',
    use='dynamic',
    process=processes.DIPfromBed_NoSedFlux
)

Variable(
    name='TIP_Settling',
    long_name='Total inorganic phosphorus settling from water to bed',
    units='mg-P/L/d',
    description='Total inorganic phosphorus settling from water to bed',
    use='dynamic',
    process=processes.TIP_Settling
)

Variable(
    name='OrgP_DIP_decay',
    long_name='Total organic phosphorus decaying to dissolved inorganic phosphrous',
    units='mg-P/L/d',
    description='Total organic phosphorus decaying to dissolved inorganic phosphrous',
    use='dynamic',
    process=processes.OrgP_DIP_decay
)

Variable(
    name='DIP_ApRespiration',
    long_name='Dissolved inorganic phosphorus released from algal respiration',
    units='mg-P/L/d',
    description='Dissolved inorganic phosphorus released from algal respiration',
    use='dynamic',
    process=processes.DIP_ApRespiration
)

Variable(
    name='DIP_ApGrowth',
    long_name='Dissolved inorganic phosphorus consumed for algal growth',
    units='mg-P/L/d',
    description='Dissolved inorganic phosphorus consumed for algal growth',
    use='dynamic',
    process=processes.DIP_ApGrowth
)

Variable(
    name='DIP_AbRespiration',
    long_name='Dissolved inorganic phosphorus released for benthic algal respiration',
    units='mg-P/L/d',
    description='Dissolved inorganic phosphorus released for benthic algal respiration',
    use='dynamic',
    process=processes.DIP_AbRespiration
)

Variable(
    name='DIP_AbGrowth',
    long_name='Dissolved inorganic phosphorus consumed for benthic algal growth',
    units='mg-P/L/d',
    description='Dissolved inorganic phosphorus consumed for benthic algal growth',
    use='dynamic',
    process=processes.DIP_AbGrowth
)

Variable(
    name='dTIPdt',
    long_name='Change in dissolved inorganic phosphorus water concentration',
    units='mg-P/L/d',
    description='Change in dissolved inorganic phosphorus water concentration',
    use='dynamic',
    process=processes.dTIPdt
)


Variable(
    name='TOP',
    long_name='Total organic phosphorus',
    units='mg-P/L',
    description='Total organic phosphorus',
    use='dynamic',
    process=processes.TOP
)

Variable(
    name='TP',
    long_name='Total phosphorus',
    units='mg-P/L',
    description='Total phosphorus',
    use='dynamic',
    process=processes.TP
)

Variable(
    name='DIP',
    long_name='Dissolve inorganich phosphorus',
    units='mg-P/L',
    description='Dissolve inorganich phosphorus',
    use='dynamic',
    process=processes.DIP
)