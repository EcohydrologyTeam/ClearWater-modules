import xarray as xr
from clearwater_modules import base
from clearwater_modules.tsm.model import EnergyBudget
from clearwater_modules.tsm import processes


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    """TSM state variables."""
    ...

Variable(
    name='water_temp_c',
    long_name='Water temperature',
    units='degC',
    description='TSM state variable for water temperature',
    use='state',
    process=processes.t_water_c,
)

# TODO: remove mock_equation

def mock_surface_area(surface_area: xr.DataArray) -> xr.DataArray:
    return surface_area

def mock_volume(volume: xr.DataArray) -> xr.DataArray:
    return volume


Variable(
    name='surface_area',
    long_name='Surface area',
    units='m^2',
    description='Surface area',
    use='state',
    process=mock_surface_area,
)
Variable(
    name='volume',
    long_name='Volume',
    units='m^3',
    description='Volume',
    use='state',
    process=mock_volume,
)
