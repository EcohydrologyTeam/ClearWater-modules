# TODO: figure out what 'model' to import
import clearwater_modules_python.base as base
from clearwater_modules_python.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

###
# To get from algae/benthic: rca?, rcb?, kdp_20 (static), kdp_T (dyn), Ap (state), depth (state), Ab (state), Fw, Fb
Variable(
    name='F_pocp',
    long_name='Fraction of algal mortality into POC',
    units='unitless',
    description='Fraction of dead algae that converts to particulate organic carbon',
    use='static'
)

Variable( ##regional?
    name='rca',
    long_name='algal C : Chla ratio',
    units='mg-C/ug-Chla',
    description='Ratio of algal carbon to chlorophyll-a',
    use='static'
)

Variable(
    name='rcb',
    long_name='benthic algae C : D ratio',
    units='mg-C/mg-D',
    description='Ratio of benthic algal carbon to algal biomass',
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
    name='F_pocb',
    long_name='fraction of benthic algal mortality into POC',
    units='unitless',
    description='fraction of benthic algal mortality into POC',
    use='static'
)

Variable(
    name='k_poc_20',
    long_name='POC hydrolysis rate at 20 degrees Celsius',
    units='/d',
    description='POC hydrolysis rate at 20 degrees Celsius',
    use='static'
)

#if we are using POM
Variable(
    name='fcom',
    long_name='fraction of carbon in organic matter',
    units='mg-C/mg-D',
    description='fraction of carbon in organic matter',
    use='static'
)

Variable(
    name='K_sOxmc',
    long_name='half saturation oxygen attenuation constant for DOC oxidation rate',
    units='mg-O2/L',
    description='half saturation oxygen attenuation constant for DOC oxidation rate',
    use='static'
)

Variable(
    name='kac_20',
    long_name='CO2 reaeration rate',
    units='/d',
    description='CO2 reaeration rate',
    use='static'
)

Variable( ## Regional
    name='pCO2',
    long_name='partial atmospheric CO2 pressure',
    units='ppm',
    description='partial pressure of CO2 in the atmosphere',
    use='static'
)

Variable( ## Regional
    name='FCO2',
    long_name='CO2 reaeration rate',
    units='/d',
    description='CO2 reaeration rate',
    use='static'
)








