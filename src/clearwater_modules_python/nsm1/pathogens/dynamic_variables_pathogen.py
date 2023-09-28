"""
File includes dynamic variables computed in pathogen module. Dynamic variables may be accessed by other modules.
"""

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.pathogens.pathogen_processes as pathogen_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

Variable(
    name='kdx_tc',
    long_name='Pathogen death rate',
    units='1/d',
    description='Pathogen death rate with temperature correction',
    use='dynamic',
    process=pathogen_processes.kdx_tc
)

Variable(
    name='PathogenDeath',
    long_name='Pathogen natural death',
    units='cfu/100mL/d',
    description='Pathogen natural death',
    use='dynamic',
    process=pathogen_processes.PathogenDeath
)

Variable(
    name='PathogenDecay',
    long_name='Pathogen death due to light',
    units='cfu/100mL/d',
    description='Pathogen death due to light',
    use='dynamic',
    process=pathogen_processes.PathogenDecay
)

Variable(
    name='PathogenSettling',
    long_name='Pathogen settling',
    units='cfu/100mL/d',
    description='Pathogen settling',
    use='dynamic',
    process=pathogen_processes.PathogenSettling
)

Variable(
    name='dPXdt',
    long_name='Change in pathogen concentration',
    units='cfu/100mL/d',
    description='Change in pathogen concentration',
    use='dynamic',
    process=pathogen_processes.dPXdt
)

Variable(
    name='PX_new',
    long_name='New pathogen concentration',
    units='cfu/100mL',
    description='New pathogen concentration',
    use='dynamic',
    process=pathogen_processes.PX_new
)