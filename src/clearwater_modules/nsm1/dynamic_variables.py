import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget

import clearwater_modules.nsm1.processes as processes

@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

############################################ From dynamic_variables_global
Variable(
    name='depth',
    long_name='Average water depth in cell',
    units='m',
    description='Average water depth in cell computed by dividing volume by surface area',
    use='dynamic',
    process=processes.compute_depth
)

Variable(
    name='TwaterK',
    long_name='Water Temperature K',
    units='K',
    description='Water temperature degree kelvin',
    use='dynamic',
    process=processes.TwaterK
)

Variable(
    name='SOD_tc',
    long_name='Sediment Oxygen Demand at water temperature tc',
    units='mg/L',
    description='Sediment Oxygen Demand at water temperature tc',
    use='dynamic',    
    process=processes.SOD_tc
)

Variable(
    name='kah_20',
    long_name='Hydraulic oxygen reaeration rate adjusted for hydraulics',
    units='1/d',
    description='Hydraulic oxygen reaeration rate adjusted for hydraulic parameters according to XX lit',
    use='dynamic',
    process=processes.kah_20
)

Variable(
    name='kah_tc',
    long_name='Hydraulic oxygen reaeration rate adjusted for temperature',
    units='1/d',
    description='Hydraulic oxygen reaeration rate adjusted for temperature',
    use='dynamic',
    process=processes.kah_tc
)

Variable(
    name='kaw_20',
    long_name='Wind oxygen reaeration velocity adjusted for hydraulics',
    units='m/d',
    description='Wind oxygen reaeration velocity adjusted for hydraulic parameters according to XX lit',
    use='dynamic',
    process=processes.kaw_20
)

Variable(
    name='kaw_tc',
    long_name='Wind oxygen reaeration velocity adjusted for temperature',
    units='m/d',
    description='Wind oxygen reaeration velocity adjusted for temperature',
    use='dynamic',
    process=processes.kaw_tc
)

Variable(
    name='ka_tc',
    long_name='Oxygen reaeration rate',
    units='1/d',
    description='Oxygen reaeration rate',
    use='dynamic',
    process=processes.ka_tc
)

Variable(
    name='L',
    long_name='Light attenuation coefficient',
    units='unitless',
    description='Light attenuation coefficient',
    use='dynamic',
    process=processes.L
)

Variable(
    name='PAR',
    long_name='surface light intensity',
    units='W/m2',
    description='surface light intensity',
    use='dynamic',
    process=processes.PAR
)

Variable(
    name='fdp',
    long_name='Fraction phosphorus dissolved',
    units='Unitless',
    description='Fraction phosphorus dissolved',
    use='dynamic',
    process=processes.fdp
)

############################################ From algae
Variable(
    name='rna',
    long_name='Algal N:Chla Ratio',
    units='mg-N/ug Chla',
    description='Algal N:Chla Ratio',
    use='dynamic',
    process=processes.rna
)

Variable(
    name='rpa',
    long_name='Algal P:Chla Ratio',
    units='mg-P/ug Chla',
    description='Algal P:Chla Ratio',
    use='dynamic',
    process=processes.rpa
)

Variable(
    name='rca',
    long_name='Algal C:Chla Ratio',
    units='mg-C/ug Chla',
    description='Algal C:Chla Ratio',
    use='dynamic',
    process=processes.rca
)

Variable(
    name='rda',
    long_name='Algal D:Chla Ratio',
    units='mg-D/ug Chla',
    description='Algal D:Chla Ratio',
    use='dynamic',
    process=processes.rda
)

Variable(
    name='mu_max_tc',
    long_name='Max Algae Growth with Temperature Correction',
    units='1/d',
    description='Max Algae Growth with Temperature Correction',
    use='dynamic',
    process=processes.mu_max_tc,
)

Variable(
    name='krp_tc',
    long_name='Algal Respiration Rate with Temperature Correction',
    units='1/d',
    description='Algal Respiration Rate with Temperature Correction',
    use='dynamic',
    process=processes.krp_tc,
)

Variable(
    name='kdp_tc',
    long_name='Algal Mortality Rate with Temperature Correction',
    units='1/d',
    description='Algal Mortality Rate with Temperature Correction',
    use='dynamic',
    process=processes.kdp_tc,
)

Variable(
    name='FL',
    long_name='Algal Light Limitation',
    units='unitless',
    description='Algal Light Limitation',
    use='dynamic',
    process=processes.FL,
)

Variable(
    name='FN',
    long_name='Algal Nitrogen Limitation',
    units='unitless',
    description='Algal Nitrogen Limitation',
    use='dynamic',
    process=processes.FN,
)

Variable(
    name='FP',
    long_name='Algal Phosphorus Limitation',
    units='unitless',
    description='Algal Phosphorus Limitation',
    use='dynamic',
    process=processes.FP,
)

Variable(
    name='mu',
    long_name='Algal Growth Rate',
    units='1/d',
    description='Algal Growth Rate',
    use='dynamic',
    process=processes.mu,
)

Variable(
    name='ApGrowth',
    long_name='Algal Growth',
    units='ug-Chala/L/d',
    description='Algal Growth',
    use='dynamic',
    process=processes.ApGrowth,
)

Variable(
    name='ApRespiration',
    long_name='Algal Respiration',
    units='ug-Chala/L/d',
    description='Algal Respiration',
    use='dynamic',
    process=processes.ApRespiration,
)

Variable(
    name='ApDeath',
    long_name='Algal Death',
    units='ug-Chala/L/d',
    description='Algal Death',
    use='dynamic',
    process=processes.ApDeath,
)

Variable(
    name='ApSettling',
    long_name='Algal Settling',
    units='ug-Chala/L/d',
    description='Algal Settling',
    use='dynamic',
    process=processes.ApSettling,
)

Variable(
    name='dApdt',
    long_name='Algal Biomass Concentration Change',
    units='ug-Chala/L/d',
    description='Algal Biomass Concentration Change',
    use='dynamic',
    process=processes.dApdt,
)

############################################ From benthic algae
Variable(
    name='mub_max_tc',
    long_name='Maximum benthic algal growth rate',
    units='1/d',
    description='Maximum benthic algal growth rate with temperature correction',
    use='dynamic',
    process=processes.mub_max_tc
)

Variable(
    name='krb_tc',
    long_name='Benthic algae respiration rate',
    units='1/d',
    description='Benthic algae respiration rate with temperature correction',
    use='dynamic',
    process=processes.krb_tc
)

Variable(
    name='kdb_tc',
    long_name='Benthic algae mortality rate',
    units='1/d',
    description='Benthic algae mortality rate with temperature correction',
    use='dynamic',
    process=processes.kdb_tc
)

Variable(
    name='rnb',
    long_name='Ratio nitrogen to dry weight',
    units='mg-N/mg-D',
    description='Ratio benthic algae nitrogen to dry weight',
    use='dynamic',
    process=processes.rnb
)

Variable(
    name='rpb',
    long_name='Ratio benthic algae phosphorus to dry weight',
    units='mg-P/mg-D',
    description='Ratio benthic algae phosphorus to dry weight',
    use='dynamic',
    process=processes.rpb
)

Variable(
    name='rcb',
    long_name='Ratio benthic algae carbon to dry weight',
    units='mg-C/mg-D',
    description='Ratio benthic algae carbon to dry weight',
    use='dynamic',
    process=processes.rcb
)

Variable(
    name='rab',
    long_name='Ratio benthic algae chlorophyll-a to dry weight',
    units='ug-Chala-a/mg-D',
    description='Ratio benthic algae chlorophyll-a to dry weight',
    use='dynamic',
    process=processes.rab
)

Variable(
    name='FLb',
    long_name='Benthic algal light limitation factor',
    units='unitless',
    description='Benthic algal light limitation factor',
    use='dynamic',
    process=processes.FLb
)

Variable(
    name='FNb',
    long_name='Benthic algal nitrogen limitation factor',
    units='unitless',
    description='Benthic algal nitrogen limitation factor',
    use='dynamic',
    process=processes.FNb
)

Variable(
    name='FPb',
    long_name='Benthic algal phosphorous limitation factor',
    units='unitless',
    description='Benthic algal phosphorous limitation factor',
    use='dynamic',
    process=processes.FPb
)

Variable(
    name='FSb',
    long_name='Benthic algal density attenuation',
    units='unitless',
    description='Benthic algal density attenuation',
    use='dynamic',
    process=processes.FSb
)

Variable(
    name='mub',
    long_name='Benthic algae specific growth rate',
    units='1/d',
    description='Benthic algae specific growth rate',
    use='dynamic',
    process=processes.mub
)

Variable(
    name='AbGrowth',
    long_name='Benthic algae growth rate',
    units='g/m^2/d',
    description='Benthic algae growth rate',
    use='dynamic',
    process=processes.AbGrowth
)

Variable(
    name='AbRespiration',
    long_name='Benthic algae respiration rate',
    units='g/m^2/d',
    description='Benthic algae respiration rate',
    use='dynamic',
    process=processes.AbRespiration
)

Variable(
    name='AbDeath',
    long_name='Benthic algae death rate',
    units='g/m^2/d',
    description='Benthic algae death rate',
    use='dynamic',
    process=processes.AbDeath
)

Variable(
    name='dAbdt',
    long_name='Change in benthic algae concentration',
    units='g/m^2/d',
    description='Change in benthic algae concentration',
    use='dynamic',
    process=processes.dAbdt
)

Variable(#TODO: figure out what this is... 
    name='Chlb',
    long_name='Chlorophyll-a concentration',
    units='mg-Chla/m^2',
    description='Chlorophyll-a concentration',
    use='dynamic',
    process=processes.Chlb
)

############################################ From nitrogen
Variable(
    name='knit_tc',
    long_name='Nitrification rate ammonia decay',
    units='1/d',
    description='Nitrification rate ammonia decay temperature correction',
    use='dynamic',
    process=processes.knit_tc
)

Variable(
    name='rnh4_tc',
    long_name='Sediment release rate of NH4',
    units='1/d',
    description=' Sediment release rate of NH4 temperature correction',
    use='dynamic',
    process=processes.rnh4_tc
)

Variable(
    name='vno3_tc',
    long_name='Sediment denitrification velocity',
    units='m/d',
    description='Sediment denitrification velocity temperature correction',
    use='dynamic',
    process=processes.vno3_tc
)

Variable(
    name='kon_tc',
    long_name='Decay rate of OrgN to NH4',
    units='1/d',
    description='Decay rate of OrgN to NH4 temperature correction',
    use='dynamic',
    process=processes.kon_tc
)

Variable(
    name='kdnit_tc',
    long_name='Denitrification rate',
    units='1/d',
    description='Denitrification rate temperature correction',
    use='dynamic',
    process=processes.kdnit_tc
)

Variable(
    name='ApUptakeFr_NH4',
    long_name='Fraction of actual floating algal uptake from ammonia pool',
    units='unitless',
    description='Fraction of actual floating algal uptake from ammonia pool',
    use='dynamic',
    process=processes.ApUptakeFr_NH4
)

Variable(
    name='ApUptakeFr_NO3',
    long_name='Fraction of actual floating algal uptake from nitrate pool',
    units='unitless',
    description='Fraction of actual floating algal uptake from nitrate pool',
    use='dynamic',
    process=processes.ApUptakeFr_NO3
)

Variable(
    name='AbUptakeFr_NH4',
    long_name='Fraction of actual benthic algal uptake from ammonia pool',
    units='unitless',
    description='Fraction of actual benthic algal uptake from ammonia pool',
    use='dynamic',
    process=processes.AbUptakeFr_NH4
)

Variable(
    name='AbUptakeFr_NO3',
    long_name='Fraction of actual benthic algal uptake from nitrate pool',
    units='unitless',
    description='Fraction of actual benthic algal uptake from nitrate pool',
    use='dynamic',
    process=processes.AbUptakeFr_NO3
)

Variable(
    name='ApDeath_OrgN',
    long_name='Algae -> OrgN',
    units='mg-N/L/d',
    description='Algae conversion to Organic nitrogen',
    use='dynamic',
    process=processes.ApDeath_OrgN
)

Variable(
    name='AbDeath_OrgN',
    long_name='Benthic Algae -> OrgN',
    units='mg-N/L/d',
    description='Benthic algae conversion to Organic nitrogen',
    use='dynamic',
    process=processes.AbDeath_OrgN
)

Variable(
    name='OrgN_NH4_Decay',
    long_name='OrgN -> NH4',
    units='mg-N/L/d',
    description='Organic nitrogen to ammonium decay',
    use='dynamic',
    process=processes.OrgN_NH4_Decay
)

Variable(
    name='OrgN_Settling',
    long_name='OrgN -> bed',
    units='mg-N/L/d',
    description='Organic nitrogen to bed settling',
    use='dynamic',
    process=processes.OrgN_Settling
)

Variable(
    name='dOrgNdt',
    long_name='Change in organic nitrogen',
    units='mg-N/L',
    description='Change in organic nitrogen',
    use='dynamic',
    process=processes.dOrgNdt
)

Variable(
    name='NH4_Nitrification',
    long_name='NH4 -> NO3  Nitrification',
    units='mg-N/L/d',
    description='NH4 Nitrification',
    use='dynamic',
    process=processes.NH4_Nitrification
)

Variable(
    name='NH4fromBed',
    long_name='bed ->  NH4 (diffusion)',
    units='mg-N/L/d',
    description='Sediment bed release of NH4',
    use='dynamic',
    process=processes.NH4fromBed
)

Variable(
    name='NH4_ApRespiration',
    long_name='Floating algae -> NH4',
    units='mg-N/L/d',
    description='Floating algae to NH4',
    use='dynamic',
    process=processes.NH4_ApRespiration
)

Variable(
    name='NH4_ApGrowth',
    long_name='NH4 -> Floating algae',
    units='mg-N/L/d',
    description='NH4 uptake to algae',
    use='dynamic',
    process=processes.NH4_ApGrowth
)

Variable(
    name='NH4_AbRespiration',
    long_name='Benthic algae -> NH4',
    units='mg-N/L/d',
    description='Benthic algae release of NH4',
    use='dynamic',
    process=processes.NH4_AbRespiration
)

Variable(
    name='NH4_AbGrowth',
    long_name='NH4 -> Benthic Algae',
    units='mg-N/L/d',
    description='Benthic algae uptake of NH4',
    use='dynamic',
    process=processes.NH4_AbGrowth
)

Variable(
    name='dNH4dt',
    long_name='Change in ammonium concentration',
    units='mg-N/L',
    description='Change in ammonium concentration',
    use='dynamic',
    process=processes.dNH4dt
)

Variable(
    name='NO3_Denit',
    long_name='NO3 -> Loss',
    units='mg-N/L/d',
    description='NO3 loss from denitrification',
    use='dynamic',
    process=processes.NO3_Denit
)

Variable(
    name='NO3_BedDenit',
    long_name='Sediment denitrification',
    units='mg-N/L/d',
    description='Sediment denitrification',
    use='dynamic',
    process=processes.NO3_BedDenit
)

Variable(
    name='NO3_ApGrowth',
    long_name='NO3 -> Floating algae',
    units='mg-N/L/d',
    description='NO3 uptake to floating algae',
    use='dynamic',
    process=processes.NO3_ApGrowth
)

Variable(
    name='NO3_AbGrowth',
    long_name='NO3 -> Benthic algae',
    units='mg-N/L/d',
    description='NO3 uptake to benthic algae',
    use='dynamic',
    process=processes.NO3_AbGrowth
)

Variable(
    name='dNO3dt',
    long_name='Change in nitrate concentration',
    units='mg-N/L',
    description='Change in nitrate concentration',
    use='dynamic',
    process=processes.dNO3dt
)

Variable(
    name='DIN',
    long_name='Dissolve inorganic nitrogen',
    units='mg-N/L',
    description='Dissolve inorganic nitrogen',
    use='dynamic',
    process=processes.DIN
)

Variable(
    name='TON',
    long_name='Total organic nitrogen',
    units='mg-N/L',
    description='Total organic nitrogen',
    use='dynamic',
    process=processes.TON
)

Variable(
    name='TKN',
    long_name='Total kjeldhl nitrogen',
    units='mg-N/L',
    description='Total kjeldhl nitrogen',
    use='dynamic',
    process=processes.TKN
)

Variable(
    name='TN',
    long_name='Total nitrogen',
    units='mg-N/L',
    description='Total nitrogen',
    use='dynamic',
    process=processes.TN
)

Variable(
    name='NitrificationInhibition',
    long_name='Nitrification Inhibitation (limits nitrification under low DO conditions)',
    units='unitless',
    description='Nitrification Inhibitation (limits nitrification under low DO conditions)',
    use='dynamic',
    process=processes.NitrificationInhibition
)

############################################ From phosphorus
Variable(
    name='kop_tc',
    long_name='Decay rate of organic P to DIP',
    units='1/d',
    description='Decay rate of organic P to DIP temperature correction',
    use='dynamic',
    process=processes.kop_tc
)

Variable(
    name='rpo4_tc',
    long_name='Benthic sediment release rate of DIP',
    units='g-P/m2/d',
    description='Benthic sediment release rate of DIP temperature correction',
    use='dynamic',
    process=processes.rpo4_tc
)

Variable(
    name='OrgP_DIP_decay',
    long_name='Organic phosphorus decay to dissolve inorganic phosphorus',
    units='mg-P/L/d',
    description='Organic phosphorus decay to dissolve inorganic phosphorus',
    use='dynamic',
    process=processes.OrgP_DIP_decay
)

Variable(
    name='OrgP_Settling',
    long_name='Organic phosphorus settling to sediment',
    units='mg-P/L/d',
    description='Organic phosphorus settling to sediment',
    use='dynamic',
    process=processes.OrgP_Settling
)

Variable(
    name='ApDeath_OrgP',
    long_name='Algal death turning into organic phosphorus',
    units='mg-P/L/d',
    description='Algal death turning into organic phosphorus',
    use='dynamic',
    process=processes.ApDeath_OrgP
)

Variable(
    name='AbDeath_OrgP',
    long_name='Benthic algal death turning into organic phosphorus',
    units='mg-P/L/d',
    description='Benthic algal death turning into organic phosphorus',
    use='dynamic',
    process=processes.AbDeath_OrgP
)

Variable(
    name='dOrgPdt',
    long_name='Change in organic phosphorus concentration',
    units='mg-P/L/d',
    description='Change in organic phosphorus concentration',
    use='dynamic',
    process=processes.dOrgPdt
)

Variable(
    name='DIPfromBed_SedFlux',
    long_name='Dissolved Organic Phosphorus coming from Bed calculated using SedFlux modules',
    units='mg-P/L/d',
    description='Dissolved Organic Phosphorus coming from Bed calculated using SedFlux modules',
    use='dynamic',
    process=processes.DIPfromBed_SedFlux
)

Variable(#TODO: find correct process
    name='DIPfromBed',
    long_name='Dissolved Organic Phosphorus coming from Bed calculated without SedFlux modules',
    units='mg-P/L/d',
    description='Dissolved Organic Phosphorus coming from Bed calculated without SedFlux modules',
    use='dynamic',
    process=processes.DIPfromBed
)

Variable(
    name='TIP_Settling',
    long_name='Total inorganic phosphorus settling from water to bed',
    units='mg-P/L/d',
    description='Total inorganic phosphorus settling from water to bed',
    use='dynamic',
    process=processes.TIP_Settling
)

Variable(
    name='OrgP_DIP_decay',
    long_name='Total organic phosphorus decaying to dissolved inorganic phosphrous',
    units='mg-P/L/d',
    description='Total organic phosphorus decaying to dissolved inorganic phosphrous',
    use='dynamic',
    process=processes.OrgP_DIP_decay
)

Variable(
    name='DIP_ApRespiration',
    long_name='Dissolved inorganic phosphorus released from algal respiration',
    units='mg-P/L/d',
    description='Dissolved inorganic phosphorus released from algal respiration',
    use='dynamic',
    process=processes.DIP_ApRespiration
)

Variable(
    name='DIP_ApGrowth',
    long_name='Dissolved inorganic phosphorus consumed for algal growth',
    units='mg-P/L/d',
    description='Dissolved inorganic phosphorus consumed for algal growth',
    use='dynamic',
    process=processes.DIP_ApGrowth
)

Variable(
    name='DIP_AbRespiration',
    long_name='Dissolved inorganic phosphorus released for benthic algal respiration',
    units='mg-P/L/d',
    description='Dissolved inorganic phosphorus released for benthic algal respiration',
    use='dynamic',
    process=processes.DIP_AbRespiration
)

Variable(
    name='DIP_AbGrowth',
    long_name='Dissolved inorganic phosphorus consumed for benthic algal growth',
    units='mg-P/L/d',
    description='Dissolved inorganic phosphorus consumed for benthic algal growth',
    use='dynamic',
    process=processes.DIP_AbGrowth
)

Variable(
    name='dTIPdt',
    long_name='Change in dissolved inorganic phosphorus water concentration',
    units='mg-P/L/d',
    description='Change in dissolved inorganic phosphorus water concentration',
    use='dynamic',
    process=processes.dTIPdt
)

Variable(
    name='TOP',
    long_name='Total organic phosphorus',
    units='mg-P/L',
    description='Total organic phosphorus',
    use='dynamic',
    process=processes.TOP
)

Variable(
    name='TP',
    long_name='Total phosphorus',
    units='mg-P/L',
    description='Total phosphorus',
    use='dynamic',
    process=processes.TP
)

Variable(
    name='DIP',
    long_name='Dissolve inorganic phosphorus',
    units='mg-P/L',
    description='Dissolve inorganic phosphorus',
    use='dynamic',
    process=processes.DIP
)

############################################ From POM
Variable(
    name='kpom_tc',
    long_name='POM dissolution rate adjusted for temperature',
    units='1/d',
    description='POM dissolution rate adjusted for temperature',
    use='dynamic',
    process=processes.kpom_tc
)

Variable(
    name='POM_algal_settling',
    long_name='POM concentration change due to algal settling',
    units='mg/L/d',
    description='POM concentration change due to algal settling',
    use='dynamic',
    process=processes.POM_algal_settling
)

Variable(
    name='POM_dissolution',
    long_name='POM concentration change due to dissolution',
    units='mg/L/d',
    description='POM concentration change due to dissolution',
    use='dynamic',
    process=processes.POM_dissolution
)

Variable(
    name='POM_POC_settling',
    long_name='POM concentration change due to POC settling',
    units='mg/L/d',
    description='POM concentration change due to POC settling',
    use='dynamic',
    process=processes.POM_POC_settling
)

Variable(
    name='POM_benthic_algae_mortality',
    long_name='POM concentration change due to algae mortality',
    units='mg/L/d',
    description='POM concentration change due to algae mortality',
    use='dynamic',
    process=processes.POM_benthic_algae_mortality
)

Variable(
    name='POM_burial',
    long_name='POM concentration change due to burial',
    units='mg/L/d',
    description='POM concentration change due to burial',
    use='dynamic',
    process=processes.POM_burial
)

Variable(
    name='dPOMdt',
    long_name='Change in POM concentration for one timestep',
    units='mg/L/d',
    description='Change in POM concentration for one timestep',
    use='dynamic',
    process=processes.dPOMdt
)

############################################ From CBOD
Variable(
    name='kbod_tc',
    long_name='Temperature adjusted oxidation rate',
    units='1/d',
    description='Temperature adjusted oxidation rate',
    use='dynamic',
    process=processes.kbod_tc
)

Variable(
    name='ksbod_tc',
    long_name='Temperature adjusted sedimentation rate',
    units='m/d',
    description='Temperature adjusted sedimentation rate',
    use='dynamic',
    process=processes.ksbod_tc
)

Variable(
    name='CBOD_oxidation',
    long_name='CBOD oxidation',
    units='mg/L/d',
    description='CBOD oxidation',
    use='dynamic',
    process=processes.CBOD_oxidation
)

Variable(
    name='CBOD_sedimentation',
    long_name='CBOD sedimentation',
    units='mg/L/d',
    description='CBOD sedimentation',
    use='dynamic',
    process=processes.CBOD_sedimentation
)

Variable(
    name='dCBODdt',
    long_name='Change in CBOD concentration for the given timestep',
    units='mg/L/d',
    description='Change in CBOD concentration for the given timestep',
    use='dynamic',
    process=processes.dCBODdt
)

############################################ From carbon
Variable(
    name='kpoc_tc',
    long_name='Temperature adjusted POC hydrolysis rate',
    units='1/d',
    description='Temperature adjusted POC hydrolysis rate',
    use='dynamic',
    process=processes.kpoc_tc,
)

Variable(
    name='POC_settling',
    long_name='POC concentration removed from cell due to settling',
    units='mg/L/d',
    description='POC concentration removed from cell due to settling',
    use='dynamic',
    process=processes.POC_settling,
)

Variable(
    name='POC_hydrolysis',
    long_name='POC concentration removed from cell due to hydrolysis',
    units='mg/L/d',
    description='POC concentration removed from cell due to hydrolysis',
    use='dynamic',
    process=processes.POC_hydrolysis,
)

Variable(
    name='POC_algal_mortality',
    long_name='POC concentration added to cell due to algal mortality',
    units='mg/L/d',
    description='POC concentration added to cell due to algal mortality',
    use='dynamic',
    process=processes.POC_algal_mortality,
)

Variable(
    name='POC_benthic_algae_mortality',
    long_name='POC concentration added to cell due to benthic algae mortality',
    units='mg/L/d',
    description='POC concentration added to cell due to benthic algae mortality',
    use='dynamic',
    process=processes.POC_benthic_algae_mortality,
)

Variable(
    name='dPOCdt',
    long_name='POC concentration change per timestep',
    units='mg/L/d',
    description='POC concentration change per timestep',
    use='dynamic',
    process=processes.dPOCdt
)

Variable(
    name='kdoc_tc',
    long_name='Dissolved organic carbon oxidation rate adjusted for temperature',
    units='1/d',
    description='Dissolved organic carbon oxidation rate adjusted for temperature',
    use='dynamic',
    process=processes.kdoc_tc
)

Variable(
    name='DOC_algal_mortality',
    long_name='DOC concentration added to cell due to algal mortality',
    units='mg/L/d',
    description='DOC concentration added to cell due to algal mortality',
    use='dynamic',
    process=processes.DOC_algal_mortality
)

Variable(
    name='DOC_benthic_algae_mortality',
    long_name='DOC concentration added to cell due to benthic algae mortality',
    units='mg/L/d',
    description='DOC concentration added to cell due to benthic algae mortality',
    use='dynamic',
    process=processes.DOC_benthic_algae_mortality,
)

Variable(
    name='DOC_oxidation',
    long_name='DOC concentration lost to cell due to oxidation',
    units='mg/L/d',
    description='DOC concentration lost to cell due to oxidation',
    use='dynamic',
    process=processes.DOC_oxidation
)

Variable(
    name='dDOCdt',
    long_name='DOC concentration change per timestep',
    units='mg/L/d',
    description='DOC concentration change per timestep',
    use='dynamic',
    process=processes.dDOCdt
)

Variable(
    name='K_H',
    long_name='Henrys coefficient',
    units='mol/L-atm',
    description='Henrys coefficient controlling the relative proportion of gaseous and aqueous phase CO2',
    use='dynamic',
    process=processes.Henrys_k
)

Variable(
    name='Atm_CO2_reaeration',
    long_name='Atmospheric CO2 reaeration',
    units='mg/L/d',
    description='Amount of DIC concentration change due to atmospheric exchange',
    use='dynamic',
    process=processes.Atmospheric_CO2_reaeration
)

Variable(
    name='DIC_algal_respiration',
    long_name='DIC generated by algal respiration',
    units='mg/L/d',
    description='DIC generated by algal respiration',
    use='dynamic',
    process=processes.DIC_algal_respiration
)

Variable(
    name='DIC_algal_photosynthesis',
    long_name='DIC consumed by algal photosynthesis',
    units='mg/L/d',
    description='DIC consumed by algal photosynthesis',
    use='dynamic',
    process=processes.DIC_algal_photosynthesis
)

Variable(
    name='DIC_benthic_algae_respiration',
    long_name='DIC generated by benthic algae respiration',
    units='mg/L/d',
    description='DIC generated by benthic algae respiration',
    use='dynamic',
    process=processes.DIC_benthic_algae_respiration
)

Variable(
    name='DIC_benthic_algae_photosynthesis',
    long_name='DIC consumed by benthic algae photosynthesis',
    units='mg/L/d',
    description='DIC consumed by benthic algae photosynthesis',
    use='dynamic',
    process=processes.DIC_benthic_algae_photosynthesis
)

Variable(
    name='DIC_CBOD_oxidation',
    long_name='DIC concentration change due to CBOD oxidation',
    units='mg/L/d',
    description='DIC concentration change due to CBOD oxidation',
    use='dynamic',
    process=processes.DIC_CBOD_oxidation
)

Variable(
    name='DIC_sed_release',
    long_name='DIC concentration change due to sediment release',
    units='mg/L/d',
    description='DIC concentration change due to sediment release',
    use='dynamic',
    process=processes.DIC_sed_release
)

Variable(
    name='dDICdt',
    long_name='DIC concentration change per timestep',
    units='mg/L/d',
    description='DIC concentration change per timestep',
    use='dynamic',
    process=processes.dDICdt
)

############################################ From DOX
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

############################################ From pathogen
Variable(
    name='kdx_tc',
    long_name='Pathogen death rate',
    units='1/d',
    description='Pathogen death rate with temperature correction',
    use='dynamic',
    process=processes.kdx_tc
)

Variable(
    name='PathogenDeath',
    long_name='Pathogen natural death',
    units='cfu/100mL/d',
    description='Pathogen natural death',
    use='dynamic',
    process=processes.PathogenDeath
)

Variable(
    name='PathogenDecay',
    long_name='Pathogen death due to light',
    units='cfu/100mL/d',
    description='Pathogen death due to light',
    use='dynamic',
    process=processes.PathogenDecay
)

Variable(
    name='PathogenSettling',
    long_name='Pathogen settling',
    units='cfu/100mL/d',
    description='Pathogen settling',
    use='dynamic',
    process=processes.PathogenSettling
)

Variable(
    name='dPXdt',
    long_name='Change in pathogen concentration',
    units='cfu/100mL/d',
    description='Change in pathogen concentration',
    use='dynamic',
    process=processes.dPXdt
)

############################################ From alkalinity
Variable(
    name='Alk_denitrification',
    long_name='Alkalinity change due to denitrification',
    units='mg/L/d',
    description='Alkalinity change due to denitrification',
    use='dynamic',
    process=processes.Alk_denitrification
)

Variable(
    name='Alk_nitrification',
    long_name='Alkalinity change due to nitrification',
    units='mg/L/d',
    description='Alkalinity change due to nitrification',
    use='dynamic',
    process=processes.Alk_nitrification
)

Variable(
    name='Alk_algal_growth',
    long_name='Alkalinity change due to algal growth',
    units='mg/L/d',
    description='Alkalinity change due to algal growth',
    use='dynamic',
    process=processes.Alk_algal_growth
)

Variable(
    name='Alk_algal_respiration',
    long_name='Alkalinity change due to algal respiration',
    units='mg/L/d',
    description='Alkalinity change due to algal respiration',
    use='dynamic',
    process=processes.Alk_algal_respiration
)

Variable(
    name='Alk_benthic_algae_growth',
    long_name='Alkalinity change due to benthic algae growth',
    units='mg/L/d',
    description='Alkalinity change due to benthic algae growth',
    use='dynamic',
    process=processes.Alk_benthic_algae_growth
)

Variable(
    name='Alk_benthic_algae_respiration',
    long_name='Alkalinity change due to benthic algae growth',
    units='mg/L/d',
    description='Alkalinity change due to benthic algae growth',
    use='dynamic',
    process=processes.Alk_benthic_algae_respiration
)

Variable(
    name='dAlkdt',
    long_name='Alkalinity concentration change per timestep',
    units='mg/L/d',
    description='Alkalinity concentration change per timestep',
    use='dynamic',
    process=processes.dAlkdt
)

############################################ From N2
Variable(
    name='KHN2_tc',
    long_name='Henrys law constant',
    units='mol/L/atm',
    description='Henrys law constant temperature corrected',
    use='dynamic',
    process=processes.KHN2_tc
)

Variable(
    name='P_wv',
    long_name='Partial pressure water vapor',
    units='atm',
    description='Partial pressure water vapor',
    use='dynamic',
    process=processes.P_wv
)

Variable(
    name='N2sat',
    long_name='N2 at saturation',
    units='mg-N/L',
    description='N2 at saturation f(Twater and atm pressure)',
    use='dynamic',
    process=processes.N2sat
)

Variable(
    name='dN2dt',
    long_name='Change in N2 air concentration',
    units='mg-N/L/d',
    description='Change in N2 air concentration',
    use='dynamic',
    process=processes.dN2dt
)

Variable(
    name='TDG',
    long_name='Total dissolved gas',
    units='%',
    description='Total dissolved gas',
    use='dynamic',
    process=processes.TDG
)