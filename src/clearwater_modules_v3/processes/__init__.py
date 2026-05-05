"""v3 processes package.

Phase 1: every process is an overlay re-export from v2. As later phases
replace each process with a v3-native implementation, the import in the
corresponding submodule is updated; this file does not change.

Run order matches v2's ``RUN_ORDER`` (Riverine first so transport state is
available before kinetics).
"""

from clearwater_modules_v3.processes.base import Process, ProcessFactory
from clearwater_modules_v3.processes.riverine import Riverine
from clearwater_modules_v3.processes.temperature import Temperature
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.nitrogen import Nitrogen

__all__ = [
    "Process",
    "ProcessFactory",
    "Riverine",
    "Temperature",
    "BenthicAlgae",
    "FloatingAlgae",
    "Nitrogen",
    "RUN_ORDER",
]

RUN_ORDER = [
    Riverine,
    Temperature,
    BenthicAlgae,
    FloatingAlgae,
    Nitrogen,
]
