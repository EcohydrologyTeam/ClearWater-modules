"""ClearWater Modules v3.

The v3 package is a self-contained codebase that grew out of the
function-style v1 framework (``clearwater_modules``) and the class-based
v2 framework. v3 inherited v2's ``Process`` abstraction, YAML
configuration, per-process substepping, and chunking execution path,
and added v1's kernel optimization, wet-mask gating, hotstart,
latent-heat unit fix, and thin-water stability guard, plus
v3-native NSM1 process implementations with audit-fixed kinetics.

Originally implemented as a thin overlay re-exporting v2 symbols. The
v2 retirement work (PRs #1-#10 on the streaming branch) ported every
v2 process and helper in-tree and removed the ``clearwater_modules_v2``
package. Some kept-in-tree helpers note their v2 lineage in their
docstrings for archeological context.

See ``design/clearwater_modules_v3_tsm_design_specification.md`` and the
umbrella ``design/clearwater_modules_v3_architecture_specification.md``.
"""

from clearwater_modules_v3 import parameters
from clearwater_modules_v3 import utils
# Phase I-4 / F9 (2026-05-21): ``examples`` is no longer eagerly
# imported at package init. The ``examples`` subpackage
# transitively imports clearwater_riverine (via
# ``processes/riverine.py``), and pulling that into v3's import
# critical path made v3 Process unit tests fail any time
# clearwater_riverine had an import-time issue (cf. Phase G-1 / G-2).
# v3 Processes are intended to be usable standalone; users who want
# the examples should explicitly ``from clearwater_modules_v3
# import examples`` in their own code.
from clearwater_modules_v3.model import Model
from clearwater_modules_v3.config import init_from_file
from clearwater_modules_v3.processes import (
    Riverine,
    Temperature,
    BenthicAlgae,
    FloatingAlgae,
    Nitrogen,
    Pathogen,
    POM,
    CBOD,
    N2,
    Phosphorus,
    Carbon,
    DOX,
    Alkalinity,
)

__all__ = [
    "Model",
    "init_from_file",
    # Phase I-4 / F9 (2026-05-21): "build_nsm1_demo" and "examples"
    # were eagerly imported at package level; removed because they
    # transitively pulled clearwater_riverine into v3's import path
    # and broke standalone Process unit testing on import-time
    # defects in the consumer. Access them explicitly via
    # ``from clearwater_modules_v3 import examples`` or
    # ``from clearwater_modules_v3.examples import build_nsm1_demo``
    # in caller code.
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
