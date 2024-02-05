import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
    

Variable(
    name='L',
    long_name='Light attenuation coefficient',
    units='unitless',
    description='Light attenuation coefficient',
    use='static',
)

Variable(
    name='fdp',
    long_name='Fraction P dissolved',
    units='unitless',
    description='Fraction P dissolved',
    use='static',
)

Variable(
    name='PAR',
    long_name='Surface light intensity',
    units='W/m^2',
    description='Surface light intensity',
    use='static',
)

Variable(
    name='vson',
    long_name='Organic N settling velocity',
    units='m/d',
    description='Organic N settling velocity',
    use='static',
)

Variable(
    name='vsoc',
    long_name='POC settling velocity',
    units='m/d',
    description='POC settling velocity',
    use='static'
)

Variable(
    name='vsop',
    long_name='Organic phosphorus settling velocity',
    units='m/d',
    description='Organic phosphorus settling velocity',
    use='static'
)

Variable(
    name='vs',
    long_name='Sediment settling velocity',
    units='m/d',
    description='Sediment settling velocity',
    use='static'
)

Variable(
    name='SOD_20',
    long_name='Sediment oxygen demand at 20 degrees C',
    units='g-O2/m/d',
    description='Sediment oxygen demand at 20 degrees C',
    use='static'
)

Variable(
    name='SOD_theta',
    long_name='Arrhenius coefficient for sediment oxygen demand',
    units='unitless',
    description='Arrhenius coefficient for sediment oxygen demand',
    use='static'
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


Variable(
    name='kaw_20_user',
    long_name='Wind oxygen reaeration velocity at 20C',
    units='m/d',
    description='Wind oxygen reaeration velocity at 20C',
    use='static'
)

Variable(
    name='kah_20_user',
    long_name='Hydraulic oxygen reaeration rate at 20C',
    units='1/d',
    description='Hydraulic oxygen reaeration rate at 20C',
    use='static'
)

Variable(
    name='hydraulic_reaeration_option',
    long_name='Option for chosing the method by which O2 reaeration rate is calculated',
    units='unitless',
    description='Selects method for computing O2 reaeration rate',
    use='static'
)

Variable(
    name='wind_reaeration_option',
    long_name='Option for chosing the method by which wind reaeration is calculated',
    units='unitless',
    description='Selects method for computing O2 reaeration due to wind',
    use='static'
)

# Global module options

Variable(
    name='use_NH4',
    long_name='Use ammonium module',
    units='unitless',
    description='True/False use ammonium module',
    use='static',
)

Variable(
    name='use_NO3',
    long_name='Use nitrate module',
    units='unitless',
    description='True/False use nitrate module',
    use='static',
)

Variable(
    name='use_OrgN',
    long_name='Use organic nitrogen module',
    units='unitless',
    description='True/False use organic nitrogen module',
    use='static',
)

Variable(
    name='use_SedFlux',
    long_name='Use sediment flux module',
    units='unitless',
    description='True/False use sediment flux module',
    use='static',
)

Variable(
    name='use_DOX',
    long_name='Use dissolved oxygen module',
    units='unitless',
    description='True/False use dissolved oxygen module',
    use='static',
)

Variable(
    name='use_Algae',
    long_name='Use algae module',
    units='unitless',
    description='True/False use algae module',
    use='static',
)

Variable(
    name='use_Balgae',
    long_name='Use benthic algae module',
    units='unitless',
    description='True/False use benthic algae module',
    use='static',
)

Variable(
    name='use_TIP',
    long_name='Use total inorganic phosphorus module',
    units='unitless',
    description='True/False use total inorganic phosphorus module',
    use='static',
)

Variable(
    name='use_OrgP',
    long_name='Use total organic phosphorus module',
    units='unitless',
    description='True/False use total organic phosphorus module',
    use='static',
)

Variable(
    name='use_POC',
    long_name='Use particulate organic carbon module',
    units='unitless',
    description='True/False use particulate organic carbon module',
    use='static',
)

Variable(
    name='use_DOC',
    long_name='Use dissolved organic carbon module',
    units='unitless',
    description='True/False use dissolved organic carbon module',
    use='static',
)

Variable(
    name='use_DIC',
    long_name='Use dissolved inorganic carbon module',
    units='unitless',
    description='True/False use dissolved inorganic carbon module',
    use='static',
)

Variable(
    name='use_N2',
    long_name='Use dissolved N2 module',
    units='unitless',
    description='True/False use N2 module',
    use='static',
)

Variable(
    name='use_Pathogen',
    long_name='Use pathogen module',
    units='unitless',
    description='True/False use pathogen module',
    use='static',
)

Variable(
    name='use_Alk',
    long_name='Use alkalinity module',
    units='unitless',
    description='True/False use alkalinity module',
    use='static',
)

Variable(
    name='use_POM',
    long_name='Use particulate organic matter module',
    units='unitless',
    description='True/False use particulate organic matter module',
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
    name='kop_20',
    long_name='Decay rate of organic P to DIP',
    units='1/d',
    description='Decay rate of organic P to DIP',
    use='static',
)

Variable(
    name='rpo4_20',
    long_name='Benthic sediment release rate of DIP',
    units='g-P/m^2/d',
    description='Benthic sediment release rate of DIP',
    use='static',
)

Variable(
    name='kdx',
    long_name='Pathogen death rate',
    units='1/d',
    description='Pathogen death rate',
    use='static',
)

Variable(
    name='apx',
    long_name='Light efficiency factor for pathogen decay',
    units='unitless',
    description='Light efficiency factor for pathogen decay',
    use='static',
)

Variable(
    name='vx',
    long_name='Pathogen net settling velocity',
    units='unitless',
    description='Pathogen net settling velocity',
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
    name='knit_20',
    long_name='Nitrification Rate Ammonia decay at 20C',
    units='1/d',
    description='Nitrification Rate Ammonia NH4 -> NO3 decay at 20C',
    use='static',
)

Variable(
    name='kon_20',
    long_name='Decay Rate of OrgN to NH4 at 20C',
    units='1/d',
    description='Decay Rate of OrgN to NH4 at 20C',
    use='static',
)

Variable(
    name='kdnit_20',
    long_name='Denitrification rate at 20C',
    units='1/d',
    description='Denitrification rate at 20C',
    use='static',
)

Variable(
    name='rnh4_20',
    long_name='Sediment release rate of NH4 at 20C',
    units='g-N/m^2/d',
    description='Sediment release rate of NH4 at 20C',
    use='static'
)

Variable(
    name='vno3_20',
    long_name='Sediment denitrification velocity at 20C',
    units='m/d',
    description='Sediment denitrification velocity at 20C',
    use='static',
)

Variable(
    name='KsOxdn',
    long_name='Half-saturation oxygen inhibition constant for denitrification',
    units='mg-O2/L',
    description='Half-saturation oxygen inhibition constant for denitrification',
    use='static',
)


Variable(
    name='PN',
    long_name='NH4 preference factor algae',
    units='unitless',
    description='NH4 preference factor algae (1=full NH4, 0=full NO3)',
    use='static',
)

Variable(
    name='PNb',
    long_name='NH4 preference factor benthic algae',
    units='unitless',
    description='NH4 preference factor benthic algae (1=full NH4, 0=full NO3)',
    use='static',
)

Variable(
    name='kbod_20',
    long_name='CBOD oxidation rate at 20C',
    units='1/d',
    description='CBOD oxidation rate at 20C',
    use='static'
)

Variable(
    name='ksbod_20',
    long_name='CBOD sedimentation rate at 20C',
    units='m/d',
    description='CBOD sedimentation rate at 20C',
    use='static'
)

Variable(
    name='ksOxbod',
    long_name='Half saturation oxygen attenuation constant for CBOD oxidation',
    units='mg-O2/L',
    description='Half saturation oxygen attenuation constant for CBOD oxidation',
    use='static'
)

Variable(
    name='F_pocp',
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
    name='F_pocb',
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
    name='K_sOxmc',
    long_name='half saturation oxygen attenuation constant for DOC oxidation rate',
    units='mg-O2/L',
    description='half saturation oxygen attenuation constant for DOC oxidation rate',
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
    name='BWd',
    long_name='Benthic algae dry weight',
    units='unitless',
    description='Benthic algae dry weight',
    use='static',
)

Variable(
    name='BWc',
    long_name='Benthic algae carbon',
    units='unitless',
    description='Benthic algae carbon',
    use='static',
)

Variable(
    name='BWn',
    long_name='Benthic algae nitrogen',
    units='unitless',
    description='Benthic algae nitrogen',
    use='static',
)

Variable(
    name='BWp',
    long_name='Benthic algae phosphorus',
    units='unitless',
    description='Benthic algae phosphorus',
    use='static',
)

Variable(
    name='BWa',
    long_name='Benthic algae Chla',
    units='unitless',
    description='Benthic algae Chla',
    use='static',
)

Variable(
    name='KLb',
    long_name='Light limiting constant for benthic algae growth',
    units='W/m^2',
    description='Light limiting constant for benthic algae growth',
    use='static',
)

Variable(
    name='KsNb',
    long_name='Half-Saturation N limiting constant for Benthic algae',
    units='mg-N/L',
    description='Half-Saturation N limiting constant for Benthic algae',
    use='static',
)

Variable(
    name='KsPb',
    long_name='Half-Saturation P limiting constant for Benthic algae',
    units='mg-P/L',
    description='Half-Saturation P limiting constant for Benthic algae',
    use='static',
)

Variable(
    name='Ksb',
    long_name='Half-Saturation density constant for benthic algae growth',
    units='g-D/m^2',
    description='Half-Saturation density constant for benthic algae growth',
    use='static',
)

Variable(
    name='mub_max_20',
    long_name='Maximum benthic algal growth rate',
    units='1/d',
    description='maximum benthic algal growth rate',
    use='static',
)

Variable(
    name='krb_20',
    long_name='Benthic algal respiration rate',
    units='1/d',
    description='Benthic algal respiration rate',
    use='static',
)

Variable(
    name='kdb_20',
    long_name='Benthic algal mortality rate',
    units='1/d',
    description='Benthic algal mortality rate',
    use='static',
)

Variable(
    name='b_growth_rate_option',
    long_name='Benthic Algal growth rate options',
    units='unitless',
    description='Benthic Algal growth rate with two options: 1) Multiplicative, 2) Limiting Nutritent',
    use='static',
)

Variable(
    name='b_light_limitation_option',
    long_name='Benthic Algal light limitation rate options',
    units='unitless',
    description='Benthic Algal light limitation rate with three options: 1) Half-saturation formulation, 2) Smiths Model, 3) Steeles Model',
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

Variable(
    name='Awd',
    long_name='Algal Dry Weight',
    units='mg',
    description='Algal Dry Weight',
    use='static',
)

Variable(
    name='AWc',
    long_name='Carbon Weight',
    units='mg',
    description='Carbon Weight',
    use='static',
)

Variable(
    name='AWn',
    long_name='Nitrogen Weight',
    units='mg',
    description='Nitrogen Weight',
    use='static',
)

Variable(
    name='AWp',
    long_name='Phosphorus Weight',
    units='mg',
    description='Phosphorus Weight',
    use='static',
)

Variable(
    name='AWa',
    long_name='Algal Chlorophyll',
    units='ug Chla',
    description='Algal Chlorophyll',
    use='static',
)


Variable(
    name='KL',
    long_name='Light Limiting Constant for Algal Growth',
    units='W/m^2',
    description='Light Limiting Constant for Algal Growth',
    use='static',
)

Variable(
    name='KsN',
    long_name='Half-Saturation N Limiting Constant for Algal Growth',
    units='mg-N/L',
    description='Half-Saturation N Limiting Constant for Algal Growth',
    use='static',
)

Variable(
    name='KsP',
    long_name='Half-Saturation P Limiting Constant for Algal Growth',
    units='mg-P/L',
    description='Half-Saturation P Limiting Constant for Algal Growth',
    use='static',
)

Variable(
    name='mu_max_20',
    long_name='Max Algae Growth',
    units='1/d',
    description='Max Algae Growth at 20C',
    use='static',
)

Variable(
    name='kdp_20',
    long_name='Algal Mortality Rate',
    units='1/d',
    description='Algal Mortality Rate at 20C',
    use='static',
)

Variable(
    name='krp_20',
    long_name='Algal Respiration Rate',
    units='1/d',
    description='Algal Respiration Rate at 20C',
    use='static',
)

Variable(
    name='vsap',
    long_name='Algal Setting Velocity',
    units='m/d',
    description='Algal Setting Velocity',
    use='static',
)

Variable(
    name='growth_rate_option',
    long_name='Growth Rate Option',
    units='1/d',
    description='Algal growth rate option 1) multiplicative, 2) Limiting Nutrient, 3) Harmonic Mean Option',
    use='static',
)

Variable(
    name='growth_rate_option',
    long_name='Growth Rate Option',
    units='1/d',
    description='Algal growth rate option 1) multiplicative, 2) Limiting Nutrient, 3) Harmonic Mean Option',
    use='static',
)

Variable(
    name='light_limitation_option',
    long_name='Light Limitation Option',
    units='1/d',
    description='Algal light limitation 1) half-saturation, 2) Smith model, 3) Steele model',
    use='static',
)











