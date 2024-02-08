"""Variable data classes used by one or more modules."""
from clearwater_modules.base import Variable
import clearwater_modules.shared.processes as processes

arrhenius_correction = Variable(
    name='arrhenius_correction',
    long_name='Arrhenius correction',
    units='unitless',
    description='Arrhenius correction',
    process=processes.arrhenius_correction,
)

Variable(
    name='lambda0',
    long_name='lambda0',
    units='1/m',
    description='background portion',
    use='static',
)

Variable(
    name='lambda1',
    long_name='lambda1',
    units='1/m/(ug Chla/L)',
    description='linear self shading',
    use='static',
)

Variable(
    name='lambda2',
    long_name='lambda2',
    units='unitless',
    description='nonlinear',
    use='static',
)

Variable(
    name='lambdas',
    long_name='lambdas',
    units='L/mg/m',
    description='ISS portion',
    use='static',
)

Variable(
    name='lambdam',
    long_name='lambdam',
    units='L/mg/m',
    description='POM portion',
    use='static',
)
