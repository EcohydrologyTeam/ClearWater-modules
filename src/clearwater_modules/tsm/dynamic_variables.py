import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.tsm.model import EnergyBudget
from clearwater_modules.tsm import processes


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


Variable(
    name='air_temp_k',
    long_name='Air temperature',
    units='K',
    description='Air temperature',
    use='dynamic',
    process=processes.air_temp_k,
)
Variable(
    name='water_temp_k',
    long_name='Water temperature',
    units='K',
    description='Water temperature',
    use='dynamic',
    process=processes.water_temp_k,
)
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
    name='density_water',
    long_name='Density of water',
    units='kg/m^3',
    description='Density of water',
    use='dynamic',
    process=processes.mf_density_water,
)
Variable(
    name='esat_mb',
    long_name='Saturation vapor pressure',
    units='mb',
    description='Saturation vapor pressure',
    use='dynamic',
    process=processes.mf_esat_mb,
)
Variable(
    name='density_air_sat',
    long_name='Density of air at saturation',
    units='kg/m^3',
    description='Density of air at saturation',
    use='dynamic',
    process=processes.mf_density_air_sat,
)
Variable(
    name='ri_number',
    long_name='Richardson number',
    units='unitless',
    description='Richardson number',
    use='dynamic',
    process=processes.ri_number,
)
Variable(
    name='ri_function',
    long_name='Richardson function',
    units='unitless',
    description='Richardson function',
    use='dynamic',
    process=processes.ri_function,
)
Variable(
    name='lv',
    long_name='Latent heat of vaporization',
    units='J/kg',
    description='Latent heat of vaporization',
    use='dynamic',
    process=processes.mf_latent_heat_vaporization,
)
Variable(
    name='cp_water',
    long_name='Specific heat of water',
    units='J/kg/K',
    description='Specific heat of water',
    use='dynamic',
    process=processes.mf_cp_water,
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
    name='q_net',
    long_name='Net heat flux',
    units='W/m^2',
    description='Net heat flux',
    use='dynamic',
    process=processes.q_net,
)
Variable(
    name='q_longwave_down',
    long_name='Downwelling longwave radiation',
    units='W/m2',
    description='Downwelling longwave radiation',
    use='dynamic',
    process=processes.mf_q_longwave_down,
)

Variable(
    name='q_longwave_up',
    long_name='Upwelling longwave radiation',
    units='W/m2',
    description='Upwelling longwave radiation',
    use='dynamic',
    process=processes.mf_q_longwave_up,
)

Variable(
    name='dTdt_water_c',
    long_name='Water temperature change',
    units='degC',
    description='Water temperature change',
    use='dynamic',
    process=processes.dTdt_water_c,
)
