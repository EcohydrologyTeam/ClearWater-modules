"""Variable data classes used by one or more modules."""
from clearwater_modules.base import Variable
import clearwater_modules.shared.processes as processes

arrhenius_correction = Variable(
    name='arrhenius_correction',
    long_name='Arrhenius correction',
    units='unitless',
    description='Arrhenius correction',
    process=processes.arrhenius_correction,
)

mf_d_esat_dT = Variable(
    name='mf_d_esat_dT',
    long_name='Derivative of saturation vapor pressure',
    units='K',
    description='Derivative of saturation vapor pressure',
    process=processes.mf_d_esat_dT,
)


mf_esat_mb = Variable(
    name='esat_mb',
    long_name='Saturation vapor pressure',
    units='mb',
    description='Saturation vapor pressure',
    process=processes.mf_esat_mb,
)

ri_number = Variable(
    name='ri_number',
    long_name='Richardson number',
    units='unitless',
    description='Richardson number',
    process=processes.ri_number,
)

ri_function = Variable(
    name='ri_function',
    long_name='Richardson function',
    units='unitless',
    description='Richardson function',
    process=processes.ri_function,
)

mf_latent_heat_vaporization = Variable(
    name='mf_latent_heat_vaporization',
    long_name='Latent heat of vaporization',
    units='W/m2',
    description='Latent heat of vaporization',
    process=processes.mf_latent_heat_vaporization,
)

mf_density_water = Variable(
    name='mf_density_water',
    long_name='Density of water',
    units='kg/m3',
    description='Density of water',
    process=processes.mf_density_water,
)

mf_density_air_sat = Variable(
    name='mf_density_air_sat',
    long_name='Density of saturated air at water surface temperature',
    units='kg/m3',
    description='Density of air',
    process=processes.mf_density_air,
)

mf_cp_water = Variable(
    name='mf_cp_water',
    long_name='Specific heat of water as a function of water temperature (degC)',
    units='J/kg/K',
    description='Specific heat of water',
    process=processes.mf_cp_water,
)
