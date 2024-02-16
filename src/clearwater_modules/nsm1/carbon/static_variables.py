# TODO: figure out what 'model' to import
import clearwater_modules.base as base
from clearwater_modules.nsm1.carbon.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...
# TODO: verify all these values

###
# To get from algae/benthic: rca?, rcb?, kdp_20 (static), kdp_tc (dyn), Ap (state), depth (state), Ab (state), Fw, Fb


Variable(
    name='f_pocp',
    long_name='Fraction of algal mortality into POC',
    units='unitless',
    description='Fraction of dead algae that converts to particulate organic carbon',
    use='static'
)


Variable(
    name='kdoc_20',
    long_name='Dissolved organic carbon oxidation rate',
    units='1/d',
    description='Dissolved organic carbon oxidation rate',
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
    name='f_pocb',
    long_name='fraction of benthic algal mortality into POC',
    units='unitless',
    description='fraction of benthic algal mortality into POC',
    use='static'
)

Variable(
    name='kpoc_20',
    long_name='POC hydrolysis rate at 20 degrees Celsius',
    units='1/d',
    description='POC hydrolysis rate at 20 degrees Celsius',
    use='static'
)

Variable(
    name='KsOxmc',
    long_name='half saturation oxygen attenuation constant for DOC oxidation rate',
    units='mg-O2/L',
    description='half saturation oxygen attenuation constant for DOC oxidation rate',
    use='static'
)


Variable(
    name='kac_20',
    long_name='CO2 reaeration rate',
    units='1/d',
    description='CO2 reaeration rate',
    use='static'
)

Variable(
    name='pCO2',
    long_name='partial atmospheric CO2 pressure',
    units='ppm',
    description='partial pressure of CO2 in the atmosphere',
    use='static'
)

Variable(
    name='FCO2',
    long_name='CO2 reaeration rate',
    units='1/d',
    description='CO2 reaeration rate',
    use='static'
)

Variable(
    name='kah_20',
    long_name='re-aeration rate (diffusion from atomsphere)',
    units='1/d',
    description='re-aeration rate (diffusion from atomsphere)',
    use='static'
)
Variable(
    name='kaw_tc',
    long_name='wind derived re-aeration transfer velocity',
    units='m/d',
    description='wind derived re-aeration transfer velocity',
    use='static'
)

#TODO define roc long name and description
Variable(
    name='roc',
    long_name='O2:C ratio for carbon oxidation',
    units='mg-O2/mg-C',
    description='32/12',
    use='static'
)