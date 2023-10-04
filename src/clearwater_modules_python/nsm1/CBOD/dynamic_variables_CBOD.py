# TODO: figure out imports

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.tsm.model import EnergyBudget
from clearwater_modules_python.nsm1.CBOD import CBOD_processes


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...

Variable(
    name='kbod_i_T',
    long_name='Temperature adjusted oxidation rate for each CBOD group',
    units='1/d',
    description='Temperature adjusted oxidation rate for each CBOD group, array',
    use='dynamic',
    process=CBOD_processes.kbod_i_T
)

Variable(
    name='ksbod_i_T',
    long_name='Temperature adjusted sedimentation rate for each CBOD group',
    units='m/d',
    description='Temperature adjusted sedimentation rate for each CBOD group, array',
    use='dynamic',
    process=CBOD_processes.ksbod_i_T
)

Variable(
    name='CBOD_oxidation',
    long_name='CBOD oxidation for all CBOD groups, array',
    units='mg/L/d',
    description='CBOD oxidation for all CBOD groups, array',
    use='dynamic',
    process=CBOD_processes.CBOD_oxidation
)

Variable(
    name='CBOD_sedimentation',
    long_name='CBOD sedimentation for all CBOD groups, array',
    units='mg/L/d',
    description='CBOD sedimentation for all CBOD groups, array',
    use='dynamic',
    process=CBOD_processes.CBOD_sedimentation
)

Variable(
    name='dCBODdt',
    long_name='Change in each CBOD group concentration for the given timestep, array',
    units='mg/L/d',
    description='Change in each CBOD group concentration for the given timestep, arrayChange in each CBOD group concentration for the given timestep, array',
    use='dynamic',
    process=CBOD_processes.dCBODdt
)