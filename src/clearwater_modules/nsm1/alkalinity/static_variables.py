import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


Variable(
    name='kdnit_20',
    long_name='Denitrification rate at 20C',
    units='1/d',
    description='Denitrification rate at 20C',
    use='static',
)

Variable(
    name='knit_20',
    long_name='Nitrification Rate Ammonia decay at 20C',
    units='1/d',
    description='Nitrification Rate Ammonia NH4 -> NO3 decay at 20C',
    use='static',
)

Variable(
    name='KNR',
    long_name='Oxygen inhabitation factor for nitrification',
    units='mg-O2/L',
    description='Oxygen inhabitation factor for nitrification',
    use='static',
)

Variable(
    name='r_alkaa',
    long_name='Ratio translating algal growth into Alk if NH4 is the N source',
    units='eq/ug-Chla',
    description='Ratio translating algal growth into Alk if NH4 is the N source',
    use='static'
)

Variable(
    name='r_alkan',
    long_name='Ratio translating algal growth into Alk if NO3 is the N source',
    units='eq/ug-Chla',
    description='Ratio translating algal growth into Alk if NO3 is the N source',
    use='static'
)

Variable(
    name='r_alkn',
    long_name='Ratio translating NH4 nitrification into Alk',
    units='eq/mg-N',
    description='Ratio translating NH4 nitrification into Alk',
    use='static'
)

Variable(
    name='r_alkden',
    long_name='Ratio translating NO3 denitrification into Alk',
    units='eq/mg-N',
    description='Ratio translating NO3 denitrification into Alk',
    use='static'
)

Variable(
    name='r_alkba',
    long_name='Ratio translating benthic algae growth into Alk if NH4 is the N source',
    units='eq/mg-D',
    description='Ratio translating benthic algae growth into Alk if NH4 is the N source',
    use='static'
)

Variable(
    name='r_alkbn',
    long_name='Ratio translating benthic algae growth into Alk if NO3 is the N source',
    units='eq/mg-D',
    description='Ratio translating benthic algae growth into Alk if NO3 is the N source',
    use='static'
)
