"""
File contains dynamic variables related to the SedFlux module
"""

import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.nsm1.model import NutrientBudget
from clearwater_modules.nsm1.sedflux import processes


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
    process=mock_equation  # TODO this depends on sedflux module
)

Variable(
    name='JNO3',
    long_name='Sediment water flux of nitrate',
    units='g-N/m^2/d',
    description='Sediment water flux of nitrate',
    use='dynamic',
    process=mock_equation  # TODO this depends on sedflux module
)


Variable(
    name='JDIP',
    long_name='Sediment water flux of phosphate',
    units='g-P/m^2/d',
    description='Sediment water flux of phosphate',
    use='dynamic',
    process=mock_equation #TODO this depends on sedflux module
)