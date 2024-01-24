# TODO: figure out what 'model' to import
import clearwater_modules.base as base
from clearwater_modules.tsm.model import EnergyBudget


@base.register_variable(models=EnergyBudget)
class Variable(base.Variable):
    ...


