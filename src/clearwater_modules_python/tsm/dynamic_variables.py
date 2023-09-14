from clearwater_modules_python import base
from clearwater_modules_python.tsm.model import EnergyBudget
from clearwater_modules_python.tsm import processes

@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


Variable(
    name='mixing_ratio_air',
    long_name='Mixing ratio of air',
    units='unitless',
    description='Mixing ratio of air',
    use='dynamic',
    process=processes.mixing_ratio_air,
)
Variable(
    name='density_air',
    long_name='Density of air',
    units='kg/m^3',
    description='Density of air',
    use='dynamic',
    process=processes.density_air,
)
Variable(
    name='emissivity_air',
    long_name='Emissivity of air',
    units='unitless',
    description='Emissivity of air',
    use='dynamic',
    process=processes.emissivity_air,
)
Variable(
    name='wind_function',
    long_name='Wind function',
    units='unitless',
    description='Wind function',
    use='dynamic',
    process=processes.wind_function,
)
Variable(
    name='q_latent',
    long_name='Latent heat flux',
    units='W/m^2',
    description='Latent heat flux',
    use='dynamic',
    process=processes.q_latent,
)
Variable(
    name='q_sensible',
    long_name='Sensible heat flux',
    units='W/m^2',
    description='Sensible heat flux',
    use='dynamic',
    process=processes.q_sensible,
)
Variable(
    name='q_sediment',
    long_name='Sediment heat flux',
    units='W/m^2',
    description='Sediment heat flux',
    use='dynamic',
    process=processes.q_sediment,
)
Variable(
    name='dTdt_sediment_c',
    long_name='Sediment temperature change',
    units='degC',
    description='Sediment temperature change',
    use='dynamic',
    process=processes.dTdt_sediment_c,
)
Variable(
    name='q_net',
    long_name='Net heat flux',
    units='W/m^2',
    description='Net heat flux',
    use='dynamic',
    process=processes.q_net,
)
Variable(
    name='dTdt_water_c',
    long_name='Water temperature change',
    units='degC',
    description='Water temperature change',
    use='dynamic',
    process=processes.dTdt_water_c,
)
