import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


Variable(
    name='KsOxdn',
    long_name='',
    units='',
    description='',
    use='static'
)

Variable(
    name='kdnit_20',
    long_name='',
    units='',
    description='',
    use='static'
)

Variable(
    name='knit_20',
    long_name='',
    units='',
    description='',
    use='static'
)

Variable(
    name='KNR',
    long_name='',
    units='',
    description='',
    use='static'
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

Variable(
    name='F1',
    long_name='Preference fraction of algal N uptake from NH4',
    units='unitless',
    description='Preference fraction of algal N uptake from NH4',
    use='static'
)

Variable(
    name='F2',
    long_name='Preference fraction of benthic algae N uptake from NH4',
    units='unitless',
    description='Preference fraction of benthic algae N uptake from NH4',
    use='static'
)

Variable(
    name='Fb',
    long_name='Fraction of bottom area available for benthic algae growth',
    units='unitless',
    description='Fraction of bottom area available for benthic algae growth',
    use='static'
)

