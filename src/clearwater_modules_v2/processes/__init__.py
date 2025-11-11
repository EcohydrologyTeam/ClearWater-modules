from .riverine import Riverine
from .temperature import Temperature
from .benthic_algae import BenthicAlgae
from .floating_algae import FloatingAlgae
from .nitrogen import Nitrogen

__all__ = ["Riverine", "Temperature", "BenthicAlgae", "FloatingAlgae", "Nitrogen"]

RUN_ORDER = [
    Riverine,
    Temperature,
    BenthicAlgae,
    FloatingAlgae,
    Nitrogen,
]
