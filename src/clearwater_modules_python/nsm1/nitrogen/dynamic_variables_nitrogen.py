"""
File includes dynamic variables computed in Algae module. Dynamic variables may be accessed by other modules.
"""

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.nitrogen.nitrogen_processes as nitrogen_processes 


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

