# TODO: figure out imports

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.nsm1.carbon.model import CarbonBudget
from clearwater_modules_python.nsm1.carbon import carbon_processes


@base.register_variable(models=CarbonBudget)
class Variable(base.Variable):
    ...


Variable( #Shared dynamic variable
    name='depth',
    long_name='Average water depth in cell',
    units='m',
    description='Average water depth in cell computed by dividing volume by surface area',
    use='dynamic',
    process=shared_processes.compute_depth
)

Variable(
    name='kpoc_T',
    long_name='Temperature adjusted POC hydrolysis rate',
    units='d-1',
    description='Temperature adjusted POC hydrolysis rate',
    use='dynamic',
    process=carbon_processes.kpoc_T,
)

Variable(
    name='POC_settling',
    long_name='POC concentration removed from cell due to settling',
    units='mg/L/s',
    description='POC concentration removed from cell due to settling',
    use='dynamic',
    process=carbon_processes.POC_settling,
)

Variable(
    name='POC_hydrolysis',
    long_name='POC concentration removed from cell due to hydrolysis',
    units='mg/L/s',
    description='POC concentration removed from cell due to hydrolysis',
    use='dynamic',
    process=carbon_processes.POC_hydrolysis,
)

Variable(
    name='POC_algal_mortality',
    long_name='POC concentration added to cell due to algal mortality',
    units='mg/L/s',
    description='POC concentration added to cell due to algal mortality',
    use='dynamic',
    process=carbon_processes.POC_algal_mortality,
)

Variable(
    name='POC_benthic_algae_mortality',
    long_name='POC concentration added to cell due to benthic algae mortality',
    units='mg/L/s',
    description='POC concentration added to cell due to benthic algae mortality',
    use='dynamic',
    process=carbon_processes.POC_benthic_algae_mortality,
)

Variable(
    name='dPOCdt',
    long_name='POC concentration change per timestep',
    units='mg/L/s',
    description='POC concentration change per timestep',
    use='dynamic',
    process=carbon_processes.POC_change
)

##DOC

Variable(
    name='kdoc_T',
    long_name='Dissolved organic carbon oxidation rate adjusted for temperature',
    units='/d',
    description='Dissolved organic carbon oxidation rate adjusted for temperature',
    use='dynamic',
    process=carbon_processes.kdoc_T
)

Variable(
    name='kpom_T',
    long_name='POM dissolution rate adjusted for temperature',
    units='/d',
    description='POM dissolution rate adjusted for temperature',
    use='dynamic',
    process=carbon_processes.kpom_T
)

Variable(
    name='DOC_POM_dissolution',
    long_name='DOC concentration change due to POM dissolution',
    units='mg/L/d',
    description='DOC concentration change due to POM dissolution',
    use='dynamic',
    process=carbon_processes.DOC_POM_dissolution
)

Variable(
    name='DOC_denitrification',
    long_name='DOC concentration change due to denitrification',
    units='mg/L/d',
    description='DOC concentration change due to denitrification',
    use='dynamic',
    process=carbon_processes.DOC_denitrification
)

Variable(
    name='DOC_algal_mortality',
    long_name='DOC concentration added to cell due to algal mortality',
    units='mg/L/s',
    description='DOC concentration added to cell due to algal mortality',
    use='dynamic',
    process=carbon_processes.DOC_algal_mortality
)

Variable(
    name='DOC_benthic_algae_mortality',
    long_name='DOC concentration added to cell due to benthic algae mortality',
    units='mg/L/s',
    description='DOC concentration added to cell due to benthic algae mortality',
    use='dynamic',
    process=carbon_processes.DOC_benthic_algae_mortality,
)

Variable(
    name='DOC_oxidation',
    long_name='DOC concentration lost to cell due to oxidation',
    units='mg/L/s',
    description='DOC concentration lost to cell due to oxidation',
    use='dynamic',
    process=carbon_processes.DOC_oxidation    
)

Variable(
    name='dDOCdt',
    long_name='DOC concentration change per timestep',
    units='mg/L/s',
    description='DOC concentration change per timestep',
    use='dynamic',
    process=carbon_processes.DOC_change
)

#DIC

Variable(
    name='K_H',
    long_name='Henrys coefficient',
    units='mol/L-atm',
    description='Henrys coefficient controlling the relative proportion of gaseous and aqueous phase CO2',
    use='dynamic',
    process=carbon_processes.Henrys_k
)

Variable(
    name='kac_T',
    long_name='temperature dependent CO2 reaeration rate',
    units='/d',
    description='temperature dependent CO2 reaeration rate',
    use='dynamic',
    process=carbon_processes.kac_T
)

Variable(
    name='Atm_CO2_reaeration',
    long_name='Atmospheric CO2 reaeration',
    units='mg/L/t',
    description='Amount of DIC flux due to atmospheric exchange',
    use='dynamic',
    process=carbon_processes.Atmospheric_CO2_reaeration
)

Variable(
    name='DIC_algal_respiration',
    long_name='DIC generated by algal respiration',
    units='mg/L/t',
    description='DIC generated by algal respiration',
    use='dynamic',
    process=carbon_processes.DIC_algal_respiration
)

Variable(
    name='DIC_algal_photosynthesis',
    long_name='DIC consumed by algal photosynthesis',
    units='mg/L/t',
    description='DIC consumed by algal photosynthesis',
    use='dynamic',
    process=carbon_processes.DIC_algal_photosynthesis
)

Variable(
    name='DIC_benthic_algae_respiration',
    long_name='DIC generated by benthic algae respiration',
    units='mg/L/t',
    description='DIC generated by benthic algae respiration',
    use='dynamic',
    process=carbon_processes.DIC_benthic_algae_respiration
)

Variable(
    name='DIC_benthic_algae_photosynthesis',
    long_name='DIC consumed by benthic algae photosynthesis',
    units='mg/L/t',
    description='DIC consumed by benthic algae photosynthesis',
    use='dynamic',
    process=carbon_processes.DIC_benthic_algae_photosynthesis
)

Variable(
    name='DIC_CBOD_oxidation',
    long_name='DIC flux due to CBOD oxidation',
    units='mg/L/t',
    description='DIC flux due to CBOD oxidation',
    use='dynamic',
    process=carbon_processes.DIC_CBOD_oxidation
)

Variable(
    name='DIC_sed_release',
    long_name='DIC flux due to sediment release',
    units='mg/L/t',
    description='DIC flux due to sediment release',
    use='dynamic',
    process=carbon_processes.DIC_sed_release
)

Variable(
    name='dDICdt',
    long_name='DIC flux per timestep',
    units='mg/L/t',
    description='DIC flux per timestep',
    use='dynamic',
    process=carbon_processes.DIC_change
)

