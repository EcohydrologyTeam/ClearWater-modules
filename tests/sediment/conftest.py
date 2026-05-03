"""Pytest configuration for SSM sediment tests.

The repo's `clearwater_modules_v2` package lives under `src/` in this
checkout, but the conda env's editable install of `clearwater_modules`
points at a vendor copy that does not contain the new
`processes.sediment` subpackage. To run these tests against the local
source tree without disturbing the rest of the repo, prepend the local
`src/` directory to `sys.path` here.

This file is scoped to `tests/sediment/` and does not affect other
test packages in the repo.
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SRC = _REPO_ROOT / "src"

if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
