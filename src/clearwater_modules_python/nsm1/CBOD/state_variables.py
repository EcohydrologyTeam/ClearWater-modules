# TODO: figure out imports, determine if state variables can be vectorized

from clearwater_modules_python import base
from clearwater_modules_python.nsmI.model import CarbonBalance
from clearwater_modules_python.nsmI import processes


@base.register_variable(models=CarbonBalance)
class Variable(base.Variable):
    """CBOD state variables."""

Variable(
    name='CBOD_i',
    long_name='Carbonaceous biochemical oxygen demand',
    units='mg/L',
    description='Array of CBOD concentrations representing the demand from different microbe groups',
    use='state',
    process=processes.update_CBOD
)