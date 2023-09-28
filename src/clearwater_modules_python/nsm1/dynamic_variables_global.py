import clearwater_modules_python.base as base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.shared.processes as shared_processes

@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

Variable(
    name='depth',
    long_name='Average water depth in cell',
    units='m',
    description='Average water depth in cell computed by dividing volume by surface area',
    use='dynamic',
    process=shared_processes.compute_depth
)

Variable(
    name='SOD_tc',
    long_name='Sediment oxygen demand adjusted for temperature',
    units='mg-O2/L/d',
    description='Sediment oxygen demand adjusted for temperature',
    use='dynamic',
    process=shared_processes.arrhenius_correction
)