"""
File includes dynamic variables computed in Balgae module. Dynamic variables may be accessed by other modules.
"""

import clearwater_modules_python.shared.processes as shared_processes
from clearwater_modules_python import base
from clearwater_modules_python.nsm1.model import NutrientBudget
import clearwater_modules_python.nsm1.balgae.balgae_processes as balgae_processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...

Variable(
    name='mub_max_tc',
    long_name='Maximum benthic algal growth rate',
    units='1/d',
    description='Maximum benthic algal growth rate with temperature correction',
    use='dynamic',
    process=balgae_processes.mub_max_tc
)

Variable(
    name='krb_tc',
    long_name='Benthic algae respiration rate',
    units='1/d',
    description='Benthic algae respiration rate with temperature correction',
    use='dynamic',
    process=balgae_processes.krb_tc
)

Variable(
    name='kdb_tc',
    long_name='Benthic algae mortality rate',
    units='1/d',
    description='Benthic algae mortality rate with temperature correction',
    use='dynamic',
    process=balgae_processes.kdb_tc
)

Variable(
    name='rnb',
    long_name='Ratio nitrogen to dry weight',
    units='mg-N/mg-D',
    description='Ratio benthic algae nitrogen to dry weight',
    use='dynamic',
    process=balgae_processes.rnb
)

Variable(
    name='rpb',
    long_name='Ratio benthic algae phosphorus to dry weight',
    units='mg-P/mg-D',
    description='Ratio benthic algae phosphorus to dry weight',
    use='dynamic',
    process=balgae_processes.rpb
)

Variable(
    name='rcb',
    long_name='Ratio benthic algae carbon to dry weight',
    units='mg-C/mg-D',
    description='Ratio benthic algae carbon to dry weight',
    use='dynamic',
    process=balgae_processes.rcb
)

Variable(
    name='rab',
    long_name='Ratio benthic algae chlorophyll-a to dry weight',
    units='ug-Chala-a/mg-D',
    description='Ratio benthic algae chlorophyll-a to dry weight',
    use='dynamic',
    process=balgae_processes.rab
)

Variable(
    name='FLb',
    long_name='Benthic algal light limitation factor',
    units='unitless',
    description='Benthic algal light limitation factor',
    use='dynamic',
    process=balgae_processes.FLb
)

Variable(
    name='FNb',
    long_name='Benthic algal nitrogen limitation factor',
    units='unitless',
    description='Benthic algal nitrogen limitation factor',
    use='dynamic',
    process=balgae_processes.FNb
)

Variable(
    name='FPb',
    long_name='Benthic algal phosphorous limitation factor',
    units='unitless',
    description='Benthic algal phosphorous limitation factor',
    use='dynamic',
    process=balgae_processes.FPb
)

Variable(
    name='FSb',
    long_name='Benthic algal density attenuation',
    units='unitless',
    description='Benthic algal density attenuation',
    use='dynamic',
    process=balgae_processes.FSb
)

Variable(
    name='mub',
    long_name='Benthic algae specific growth rate',
    units='1/d',
    description='Benthic algae specific growth rate',
    use='dynamic',
    process=balgae_processes.mub
)

Variable(
    name='AbGrowth',
    long_name='Benthic algae growth rate',
    units='g/m^2/d',
    description='Benthic algae growth rate',
    use='dynamic',
    process=balgae_processes.AbGrowth
)

Variable(
    name='AbRespiration',
    long_name='Benthic algae respiration rate',
    units='g/m^2/d',
    description='Benthic algae respiration rate',
    use='dynamic',
    process=balgae_processes.AbRespiration
)

Variable(
    name='AbDeath',
    long_name='Benthic algae death rate',
    units='g/m^2/d',
    description='Benthic algae death rate',
    use='dynamic',
    process=balgae_processes.AbDeath
)

Variable(
    name='dAbdt',
    long_name='Change in benthic algae concentration',
    units='g/m^2/d',
    description='Change in benthic algae concentration',
    use='dynamic',
    process=balgae_processes.dAbdt
)

Variable(
    name='Chlb',
    long_name='Chlorophyll-a concentration',
    units='mg-Chla/m^2',
    description='Chlorophyll-a concentration',
    use='dynamic',
    process=balgae_processes.Chlb
)