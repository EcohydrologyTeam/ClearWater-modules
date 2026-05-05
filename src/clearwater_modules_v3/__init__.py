"""ClearWater Modules v3.

The v3 package converges the function-style v1 framework
(``clearwater_modules``) with the class-based v2 framework
(``clearwater_modules_v2``) into a single coherent codebase. v3 retains v2's
``Process`` abstraction, YAML configuration, per-process substepping, and
chunking execution path, and adds v1's kernel optimization, wet-mask gating,
hotstart, latent-heat unit fix, and thin-water stability guard.

Phase 1 (current): the package is a thin overlay; every symbol re-exports
its v2 counterpart unchanged. v3-native code lands in subsequent phases per
``design/clearwater_modules_v3_tsm_design_specification.md`` and the umbrella
``design/clearwater_modules_v3_architecture_specification.md``.
"""

from clearwater_modules_v3 import parameters
from clearwater_modules_v3 import utils
from clearwater_modules_v3.model import Model
from clearwater_modules_v3.config import init_from_file
from clearwater_modules_v3.processes import (
    Riverine,
    Temperature,
    BenthicAlgae,
    FloatingAlgae,
    Nitrogen,
)

__all__ = [
    "Model",
    "init_from_file",
    "Riverine",
    "Temperature",
    "BenthicAlgae",
    "FloatingAlgae",
    "Nitrogen",
]
