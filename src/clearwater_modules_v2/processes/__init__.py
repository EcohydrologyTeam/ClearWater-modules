from .temperature import Temperature
from .benthic_algae import BenthicAlgae
from .floating_algae import FloatingAlgae

# Nitrogen and Riverine have been moved to ``clearwater_modules_v3.processes``
# as part of the v2 retirement. Import from v3 directly.

__all__ = ["Temperature", "BenthicAlgae", "FloatingAlgae"]

RUN_ORDER = [
    Temperature,
    BenthicAlgae,
    FloatingAlgae,
]
