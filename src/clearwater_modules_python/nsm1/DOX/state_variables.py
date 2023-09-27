# TODO: figure out imports, determine if state variables can be vectorized

from clearwater_modules_python import base
from clearwater_modules_python.nsmI.model import CarbonBalance
from clearwater_modules_python.nsmI import processes


@base.register_variable(models=CarbonBalance)
class Variable(base.Variable):
    """Dissolved Oxygen state variables."""
    ...

Variable(
    name='DOX',
    long_name='Dissolved oxygen',
    units='mg/L',
    description='Dissolved oxygen',
    use='state',
    process=processes.update_DOX
)