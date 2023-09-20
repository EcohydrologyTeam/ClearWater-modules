from clearwater_modules_python import base
from clearwater_modules_python.nsmI.model import CarbonBalance
from clearwater_modules_python.nsmI import processes


@base.register_variable(models=CarbonBalance)
class Variable(base.Variable):
    """Carbon state variables."""
    ...


Variable(
    name='DOC',
    long_name='Dissolved organic carbon',
    units='mg/L',
    description='NSMI state variable for dissolved organic carbon',
    use='state',
    process=processes.update_DOC,
)

Variable(
    name='POC',
    long_name='Particulate organic carbon',
    units='mg/L',
    description='NSMI state variable for particulate organic carbon',
    use='state',
    process=processes.update_POC,
)

Variable(
    name='DIC',
    long_name='Dissolved inorganic carbon',
    units='mg/L',
    description='NSMI state variable for dissolved inorganic carbon',
    use='state',
    process=processes.update_DIC,
)


### Calculated outside of module but pulled in... from previous timestep.. so are these static? 

Variable(
    name='Algae',
    long_name='Algae concentration',
    units='mg/L',
    description='NSMI state variable for algae concentration',
    use='state',
    process=*ANOTHER MODULE?.processes.algae,
)

Variable(
    name='POM',
    long_name='Particulate organic matter',
    units='mg/L',
    description='NSMI state variable for particulate organic matter',
    use='state',
    process=*ANOTHER MODULE?.processes.POM,
)

Variable(
    name='TSS',
    long_name='Total suspended solids',
    units='mg/L',
    description='NSMI state variable for total suspended solids',
    use='state',
    process=*ANOTHER MODULE?.processes.TSS,
)

###

# TODO: remove mock_equation


def mock_equation(water_temp_c: float) -> float:
    return water_temp_c ** 2


Variable(
    name='surface_area',
    long_name='Surface area',
    units='m^2',
    description='Surface area',
    use='state',
    process=mock_equation,
)
Variable(
    name='volume',
    long_name='Volume',
    units='m^3',
    description='Volume',
    use='state',
    process=mock_equation,
)