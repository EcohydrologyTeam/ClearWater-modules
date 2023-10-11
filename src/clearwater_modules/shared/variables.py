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

