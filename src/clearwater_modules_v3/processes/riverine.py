"""v3 ``Riverine`` process.

Phase 1: overlay re-export from v2.

The Riverine wrapper around ``clearwater_riverine`` is unmodified in v3
through at least the v3 1.0.0 release. NSM2 work may extend it later to
register additional state variables from the riverine mesh.
"""

from clearwater_modules_v2.processes.riverine import Riverine

__all__ = ["Riverine"]
