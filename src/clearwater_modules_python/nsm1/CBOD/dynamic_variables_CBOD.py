# TODO: figure out imports

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.tsm.model import EnergyBudget
from clearwater_modules_python.nsm1 import processes


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...

Variable(
    name='kbod_i_T',
    long_name='Temperature adjusted oxidation rate for each CBOD group',
    units='/d',
    description='Temperature adjusted oxidation rate for each CBOD group, array',
    use='dynamic',
    process=processes.kbod_i_T
)

Variable(
    name='ksbod_i_T',
    long_name='Temperature adjusted sedimentation rate for each CBOD group',
    units='m/d',
    description='Temperature adjusted sedimentation rate for each CBOD group, array',
    use='dynamic',
    process=processes.ksbod_i_T
)

Variable(
    name='CBOD_oxidation',
    long_name='CBOD oxidation for all CBOD groups, array',
    units='mg/L/t',
    description='CBOD oxidation for all CBOD groups, array',
    use='dynamic',
    process=processes.CBOD_oxidation
)

Variable(
    name='CBOD_sedimentation',
    long_name='CBOD sedimentation for all CBOD groups, array',
    units='mg/L/t',
    description='CBOD sedimentation for all CBOD groups, array',
    use='dynamic',
    process=processes.CBOD_sedimentation
)

Variable(
    name='dCBOD_idt',
    long_name='Change in each CBOD group concentration for the given timestep, array',
    units='mg/L/t',
    description='Change in each CBOD group concentration for the given timestep, arrayChange in each CBOD group concentration for the given timestep, array',
    use='dynamic',
    process=processes.CBOD_change
)