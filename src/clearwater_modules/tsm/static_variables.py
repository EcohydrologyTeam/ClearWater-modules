import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values


Variable(
    name='use_sed_temp',
    long_name='Use Sediment Temperature?',
    units='boolean',
    description='Controls whether to use/calculate sediment temperature or not.',
    use='static',
)

Variable(
    name='stefan_boltzmann',
    long_name='Stefan-Boltzmann Constant',
    units='W m-2 K-4',
    description='The Stefan-Boltzmann constant.',
    use='static',
)

Variable(
    name='cp_air',
    long_name='Specific Heat Capacity of Air',
    units='J kg-1 K-1',
    description='The specific heat capacity of air.',
    use='static',
)

Variable(
    name='emissivity_water',
    long_name='Emissivity of Water',
    units='1',
    description='The emissivity of water.',
    use='static',
)

Variable(
    name='gravity',
    long_name='Gravity',
    units='m s-2',
    description='The acceleration due to gravity.',
    use='static',
)

Variable(
    name='a0',
    long_name='Albedo of Water',
    units='unitless',
    description='The albedo of water.',
    use='static',
)

Variable(
    name='a1',
    long_name='Albedo of Water',
    units='unitless',
    description='The albedo of water.',
    use='static',
)
Variable(
    name='a2',
    long_name='Albedo of Water',
    units='unitless',
    description='The albedo of water.',
    use='static',
)
Variable(
    name='a3',
    long_name='Albedo of Water',
    units='unitless',
    description='The albedo of water.',
    use='static',
)
Variable(
    name='a4',
    long_name='Albedo of Water',
    units='unitless',
    description='The albedo of water.',
    use='static',
)
Variable(
    name='a5',
    long_name='Albedo of Water',
    units='unitless',
    description='The albedo of water.',
    use='static',
)
Variable(
    name='a6',
    long_name='Albedo of Water',
    units='unitless',
    description='The albedo of water.',
    use='static',
)
Variable(
    name='pb',
    long_name='Bulk Density of Sediment',
    units='kg m-3',
    description='The bulk density of sediment.',
    use='static',
)
Variable(
    name='cps',
    long_name='CPS',
    units='J kg-1 K-1',
    description='The CPS.',
    use='static',
)
Variable(
    name='h2',
    long_name='H2',
    units='m',
    description='The H2.',
    use='static',
)
Variable(
    name='alphas',
    long_name='Alphas',
    units='m-1',
    description='The alphas.',
    use='static',
)
Variable(
    name='richardson_option',
    long_name='Richardson Option',
    units='unitless',
    description='The Richardson option.',
    use='static',
)
Variable(
    name='air_temp_c',
    long_name='Air Temperature',
    units='C',
    description='The air temperature.',
    use='static',
)
Variable(
    name='q_solar',
    long_name='Solar Heat Flux',
    units='W m-2',
    description='The solar heat flux.',
    use='static',
)
Variable(
    name='sed_temp_c',
    long_name='Sediment Temperature',
    units='C',
    description='The sediment temperature.',
    use='static',
)
Variable(
    name='eair_mb',
    long_name='Vapor pressure',
    units='mb',
    description='Vapor pressure',
    use='static',
)
Variable(
    name='pressure_mb',
    long_name='Pressure',
    units='mb',
    description='The pressure.',
    use='static',
)
Variable(
    name='cloudiness',
    long_name='Cloudiness',
    units='unitless',
    description='The cloudiness.',
    use='static',
)
Variable(
    name='wind_speed',
    long_name='Wind Speed',
    units='m s-1',
    description='The wind speed.',
    use='static',
)
Variable(
    name='wind_a',
    long_name='Wind A',
    units='m s-1',
    description='The wind A.',
    use='static',
)
Variable(
    name='wind_b',
    long_name='Wind B',
    units='m s-1',
    description='The wind B.',
    use='static',
)
Variable(
    name='wind_c',
    long_name='Wind C',
    units='m s-1',
    description='The wind C.',
    use='static',
)
Variable(
    name='wind_kh_kw',
    long_name='Wind KH KW',
    units='m s-1',
    description='The wind KH KW.',
    use='static',
)
