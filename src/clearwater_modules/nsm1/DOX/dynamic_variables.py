# TODO: figure out imports

import clearwater_modules.shared.processes as shared_processes
from clearwater_modules import base
from clearwater_modules.tsm.model import EnergyBudget
from clearwater_modules.nsm1.DOX import processes


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


Variable(
    name='DOX_sat',
    long_name='DO saturation concentration',
    units='mg/L',
    description='DO saturation concentration in water as a function of water temperature (K)',
    use='dynamic',
    process=processes.DOX_sat
)

Variable(
    name='pwv',
    long_name='Partial pressure of water vapor',
    units='atm',
    description='Partial pressure of water vapor',
    use='dynamic',
    process=processes.pwv
)

Variable(
    name='DOs_atm_alpha',
    long_name='DO saturation atmospheric correction coefficient',
    units='unitless',
    description='DO saturation atmospheric correction coefficient',
    use='dynamic',
    process=processes.DOs_atm_alpha
)

Variable(
    name='Atm_O2_reaeration',
    long_name='Atmospheric oxygen reaeration',
    units='mg/L/d',
    description='Atmospheric oxygen reaeration, can fluctuate both in and out of waterbody',
    use='dynamic',
    process=processes.Atm_O2_reaeration
)

# TODO: UPDATE BASED ON FORTRAN
Variable(
    name='DOX_ApGrowth',
    long_name='Dissolved oxygen flux due to algal photosynthesis',
    units='mg/L/d',
    description='Dissolved oxygen flux due to algal photosynthesis',
    use='dynamics',
    process=processes.DOX_ApGrowth
)

# TODO: UPDATE BASED ON FORTRAN
Variable(
    name='DOX_algal_respiration',
    long_name='Dissolved oxygen flux due to algal respiration',
    units='mg/L/d',
    description='Dissolved oxygen flux due to algal respiration',
    use='dynamic',
    process=processes.DOX_ApRespiration
)

Variable(
    name='DOX_Nitrification',
    long_name='Dissolved oxygen flux due to nitrification',
    units='mg/L/d',
    description='Dissolved oxygen flux due to nitrification',
    use='dynamic',
    process=processes.DOX_Nitrification
)

Variable(
    name='DOX_DOC_Oxidation',
    long_name='Dissolved oxygen flux due to DOC oxidation',
    units='mg/L/d',
    description='Dissolved oxygen flux due to DOC oxidation',
    use='dynamic',
    process=processes.DOX_DOC_Oxidation
)

Variable(
    name='DOX_CBOD_Oxidation',
    long_name='Dissolved oxygen flux due to CBOD oxidation',
    units='mg/L/d',
    description='Dissolved oxygen flux due to CBOD oxidation',
    use='dynamic',
    process=processes.DOX_CBOD_Oxidation
)

Variable(
    name='DOX_AbGrowth',
    long_name='Dissolved oxygen flux due to benthic algae photosynthesis',
    units='mg/L/d',
    description='Dissolved oxygen flux due to benthic algae photosynthesis',
    use='dynamics',
    process=processes.DOX_AbGrowth
)

Variable(
    name='DOX_AbRespiration',
    long_name='Dissolved oxygen flux due to benthic algae respiration',
    units='mg/L/d',
    description='Dissolved oxygen flux due to benthic algae respiration',
    use='dynamic',
    process=processes.DOX_AbRespiration
)


Variable(
    name='DOX_SOD',
    long_name='Dissolved oxygen flux due to sediment oxygen demand',
    units='mg/L/d',
    description='Dissolved oxygen flux due to sediment oxygen demand',
    use='dynamic',
    process=processes.DOX_SOD
)

Variable(
    name='dDOXdt',
    long_name='Change in dissolved oxygen concentration for one timestep',
    units='mg/L/d',
    description='Change in dissolved oxygen concentration for one timestep',
    use='dynamic',
    process=processes.dDOXdt
)
