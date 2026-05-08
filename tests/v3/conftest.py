"""Test-environment shim for the v3 test suite.

The streaming repo's ``src/`` is the canonical home for
``clearwater_modules_v3`` (parameter library, ``utils/numerics``
clip-with-log, v3-native processes including the v2-retirement
ports). The conda env's editable ``clearwater_modules`` install,
however, points at a *vendor* copy under
``Publication-ClearWater-Riverine-01-Temperature/notebooks/vendor/`` that
does not carry the streaming-repo work and would shadow it on
``sys.path``.

Mitigation: prepend the streaming-repo ``src/`` to ``sys.path`` before
any tests in ``tests/v3/`` import ``clearwater_modules_v3``.
``sys.path.insert(0, ...)`` ensures the local source is preferred over
the vendor editable install.

Once the streaming repo's pyproject grows its own editable install
referencing ``./src``, this shim can be deleted.
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SRC = _REPO_ROOT / "src"

if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
