from clearwater_modules_python.base import (
    Variable,
)
from clearwater_modules_python.tsm import processes

mixing_ratio_air = Variable(
    name='mixing_ratio_air',
    long_name='Mixing ratio of air',
    units='unitless',
    description='Mixing ratio of air',
    process=processes.mixing_ratio_air,
)

density_air = Variable(
    name='density_air',
    long_name='Density of air',
    units='kg/m^3',
    description='Density of air',
    process=processes.density_air,
)

emissivity_air = Variable(
    name='emissivity_air',
    long_name='Emissivity of air',
    units='unitless',
    description='Emissivity of air',
    process=processes.emissivity_air,
)

wind_function = Variable(
    name='wind_function',
    long_name='Wind function',
    units='unitless',
    description='Wind function',
    process=processes.wind_function,
)

q_latent = Variable(
    name='q_latent',
    long_name='Latent heat flux',
    units='W/m^2',
    description='Latent heat flux',
    process=processes.q_latent,
)

q_sensible = Variable(
    name='q_sensible',
    long_name='Sensible heat flux',
    units='W/m^2',
    description='Sensible heat flux',
    process=processes.q_sensible,
)

q_sediment = Variable(
    name='q_sediment',
    long_name='Sediment heat flux',
    units='W/m^2',
    description='Sediment heat flux',
    process=processes.q_sediment,
)

dTdt_sediment_c = Variable(
    name='dTdt_sediment_c',
    long_name='Sediment temperature change',
    units='degC',
    description='Sediment temperature change',
    process=processes.dTdt_sediment_c,
)

q_net = Variable(
    name='q_net',
    long_name='Net heat flux',
    units='W/m^2',
    description='Net heat flux',
    process=processes.q_net,
)

dTdt_water_c = Variable(
    name='dTdt_water_c',
    long_name='Water temperature change',
    units='degC',
    description='Water temperature change',
    process=processes.dTdt_water_c,
)
