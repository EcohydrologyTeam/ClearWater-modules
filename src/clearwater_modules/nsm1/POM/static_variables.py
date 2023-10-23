import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


Variable(
    name='kdb_20',
    long_name='Benthic algal mortality rate at 20C',
    units='1/d',
    description='Benthic algal mortality rate at 20C',
    use='static',
)

Variable(
    name='kpom_20',
    long_name='POM dissolution rate at 20C',
    units='1/d',
    description='POM dissolution rate at 20C',
    use='static'
)

Variable(
    name='rda',
    long_name='Algal D:Chla Ratio',
    units='mg-D/ug Chla',
    description='Algal D:Chla Ratio',
    use='static'
)

Variable(
    name='vsoc',
    long_name='POC settling velocity',
    units='m/d',
    description='POC settling velocity',
    use='static'
)

Variable(
    name='vsap',
    long_name='Algal Setting Velocity',
    units='m/d',
    description='Algal Setting Velocity',
    use='static',
)

Variable(
    name='Fw',
    long_name='Fraction of benthic algae mortality into water column',
    units='unitless',
    description='Fraction of benthic algae mortality into water column',
    use='static',
)

Variable(
    name='Fb',
    long_name='Fraction of bottom area available for benthic algae',
    units='unitless',
    description='Fraction of bottom area available for benthic algae',
    use='static',
)

Variable(
    name='fcom',
    long_name='Fraction of carbon in organic matter',
    units='mg-C/mg-D',
    description='Fraction of carbon in organic matter',
    use='static'
)

Variable(
    name='vb',
    long_name='Burial velocity',
    units='m/d',
    description='Rate at which constituents are buried on the bottom',
    use='static'
)
