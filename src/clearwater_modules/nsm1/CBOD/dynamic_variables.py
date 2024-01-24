# TODO: figure out imports...

import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.CBOD.processes as processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...


Variable(
    name='kbod_T',
    long_name='Temperature adjusted oxidation rate',
    units='1/d',
    description='Temperature adjusted oxidation rate',
    use='dynamic',
    process=processes.kbod_T
)

Variable(
    name='ksbod_T',
    long_name='Temperature adjusted sedimentation rate',
    units='m/d',
    description='Temperature adjusted sedimentation rate',
    use='dynamic',
    process=processes.ksbod_T
)

Variable(
    name='CBOD_oxidation',
    long_name='CBOD oxidation',
    units='mg/L/d',
    description='CBOD oxidation',
    use='dynamic',
    process=processes.CBOD_oxidation
)

Variable(
    name='CBOD_sedimentation',
    long_name='CBOD sedimentation',
    units='mg/L/d',
    description='CBOD sedimentation',
    use='dynamic',
    process=processes.CBOD_sedimentation
)

Variable(
    name='dCBODdt',
    long_name='Change in CBOD concentration for the given timestep',
    units='mg/L/d',
    description='Change in CBOD concentration for the given timestep',
    use='dynamic',
    process=processes.dCBODdt
)
