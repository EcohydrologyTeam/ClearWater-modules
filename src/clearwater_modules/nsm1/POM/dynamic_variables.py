import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.carbon.model import CarbonBudget
from clearwater_modules.nsm1.POM import processes


@base.register_variable(models=CarbonBudget)
class Variable(base.Variable):
    ...


Variable(
    name='kdb_T',
    long_name='Benthic algal mortality rate adjusted for temperature',
    units='1/d',
    description='Benthic algal mortality rate adjusted for temperature',
    use='dynamic',
    process=processes.kdb_T
)

Variable(
    name='kpom_T',
    long_name='POM dissolution rate adjusted for temperature',
    units='1/d',
    description='POM dissolution rate adjusted for temperature',
    use='dynamic',
    process=processes.kpom_T
)

Variable(
    name='depth',
    long_name='Water depth in computation cell',
    units='m',
    description='Water depth in computation cell',
    use='dynamic',
    process=shared_processes.compute_depth
)

Variable(
    name='POM_algal_settling',
    long_name='POM concentration change due to algal settling',
    units='mg/L/d',
    description='POM concentration change due to algal settling',
    use='dynamic',
    process=processes.POM_algal_settling
)

Variable(
    name='POM_dissolution',
    long_name='POM concentration change due to dissolution',
    units='mg/L/d',
    description='POM concentration change due to dissolution',
    use='dynamic',
    process=processes.POM_dissolution
)

Variable(
    name='POM_POC_settling',
    long_name='POM concentration change due to POC settling',
    units='mg/L/d',
    description='POM concentration change due to POC settling',
    use='dynamic',
    process=processes.POM_POC_settling
)

Variable(
    name='POM_benthic_algae_mortality',
    long_name='POM concentration change due to algae mortality',
    units='mg/L/d',
    description='POM concentration change due to algae mortality',
    use='dynamic',
    process=processes.POM_benthic_algae_mortality
)

Variable(
    name='POM_burial',
    long_name='POM concentration change due to burial',
    units='mg/L/d',
    description='POM concentration change due to burial',
    use='dynamic',
    process=processes.POM_burial
)

Variable(
    name='dPOMdt',
    long_name='Change in POM concentration for one timestep',
    units='mg/L/d',
    description='Change in POM concentration for one timestep',
    use='dynamic',
    process=processes.dPOMdt
)
