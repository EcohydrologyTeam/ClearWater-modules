from .riverine import Riverine
from .temperature import Temperature
from .benthic_algae import BenthicAlgae
from .floating_algae import FloatingAlgae

# Nitrogen has been moved to ``clearwater_modules_v3.processes.nitrogen``
# as part of the v2 retirement. Import from v3 directly.

__all__ = ["Riverine", "Temperature", "BenthicAlgae", "FloatingAlgae"]

RUN_ORDER = [
    Riverine,
    Temperature,
    BenthicAlgae,
    FloatingAlgae,
]
