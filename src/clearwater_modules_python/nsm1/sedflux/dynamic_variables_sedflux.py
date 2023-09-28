"""
File includes dynamic variables computed in Algae module. Dynamic variables may be accessed by other modules.
"""

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.nsm1.model import NutrientBudget


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

def mock_equation(water_temp_c: float) -> float:
    return water_temp_c ** 2

Variable(
    name='JNH4',
    long_name='Sediment water flux of ammonia',
    units='g-N/m^2/d',
    description='Sediment water flux of ammonia',
    use='dynamic',
    process=mock_equation #TODO this depends on sedflux module
)

Variable(
    name='JNO3',
    long_name='Sediment water flux of nitrate',
    units='g-N/m^2/d',
    description='Sediment water flux of nitrate',
    use='dynamic',
    process=mock_equation #TODO this depends on sedflux module
)


Variable(
    name='JDIP',
    long_name='Sediment water flux of phosphate',
    units='g-P/m^2/d',
    description='Sediment water flux of phosphate',
    use='dynamic',
    process=mock_equation #TODO this depends on sedflux module
)