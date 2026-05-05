"""v3 processes package.

Phase 1: every process is an overlay re-export from v2. As later phases
replace each process with a v3-native implementation, the import in the
corresponding submodule is updated; this file does not change.

m18 (review-findings 2026-05-04): the prior ``RUN_ORDER`` constant was
exported but never consulted by v3 orchestration code. Process firing
order is determined by the order of ``processes`` passed to ``Model``
(typically driven by the YAML ``processes:`` block), not by a
package-level constant. The constant has been removed to avoid the
false expectation that reordering it changes behavior.
"""

from clearwater_modules_v3.processes.base import Process, ProcessFactory
from clearwater_modules_v3.processes.riverine import Riverine
from clearwater_modules_v3.processes.temperature import Temperature
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.nitrogen import Nitrogen
from clearwater_modules_v3.processes.pathogen import Pathogen
from clearwater_modules_v3.processes.pom import POM
from clearwater_modules_v3.processes.cbod import CBOD
from clearwater_modules_v3.processes.n2 import N2
from clearwater_modules_v3.processes.phosphorus import Phosphorus
from clearwater_modules_v3.processes.carbon import Carbon
from clearwater_modules_v3.processes.dox import DOX
from clearwater_modules_v3.processes.alkalinity import Alkalinity

__all__ = [
    "Process",
    "ProcessFactory",
    "Riverine",
    "Temperature",
    "BenthicAlgae",
    "FloatingAlgae",
    "Nitrogen",
    "Pathogen",
    "POM",
    "CBOD",
    "N2",
    "Phosphorus",
    "Carbon",
    "DOX",
    "Alkalinity",
]
