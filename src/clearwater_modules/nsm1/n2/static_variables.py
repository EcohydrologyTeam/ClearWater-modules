"""
File includes static variables only used in N2 module
"""

import clearwater_modules.base as base
from clearwater_modules.nsm1.model import NutrientBudget
import clearwater_modules.nsm1.n2.processes as processes


@base.register_variable(models=NutrientBudget)
class Variable(base.Variable):
    ...
