"""Test-environment shim for the v3 test suite.

The streaming repo carries a stale local copy of ``clearwater_modules_v2``
(see Phase 0 gap analysis). The clean v2 lives in the modules repo and is
already on ``sys.path`` via the modules-repo's editable install. To make
v3 importable for tests without letting the stale local v2 shadow the
clean one, this conftest:

1. Imports v2 first while ``sys.path`` still favors the modules-repo
   path (so v2 modules get cached against the clean source).
2. Appends the streaming-repo ``src/`` directory to ``sys.path``
   *after* v2 is cached. v3's overlay imports of v2 will then hit the
   import cache and resolve to the modules-repo version.

This shim is **only** required for testing v3 from the streaming repo
against the modules-repo pixi env. Once the streaming repo grows its
own pixi env (architecture spec section 1) and the streaming-local v2
is brought up to upstream, this file can be deleted.
"""

import os
import sys

# Pre-cache clean v2 modules from the modules-repo's editable install
# before adding the streaming-repo path to sys.path.
import clearwater_modules_v2  # noqa: F401
import clearwater_modules_v2.config  # noqa: F401
import clearwater_modules_v2.config.init  # noqa: F401
import clearwater_modules_v2.config.read  # noqa: F401
import clearwater_modules_v2.model  # noqa: F401
import clearwater_modules_v2.processes  # noqa: F401
import clearwater_modules_v2.processes.base  # noqa: F401
import clearwater_modules_v2.processes.benthic_algae  # noqa: F401
import clearwater_modules_v2.processes.floating_algae  # noqa: F401
import clearwater_modules_v2.processes.nitrogen  # noqa: F401
import clearwater_modules_v2.processes.riverine  # noqa: F401
import clearwater_modules_v2.processes.temperature  # noqa: F401
import clearwater_modules_v2.utils  # noqa: F401
import clearwater_modules_v2.utils.constants  # noqa: F401
import clearwater_modules_v2.utils.conversions  # noqa: F401

# Now expose v3 by appending the streaming-repo src to sys.path.
_STREAMING_SRC = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "src",
)
if _STREAMING_SRC not in sys.path:
    sys.path.append(_STREAMING_SRC)
