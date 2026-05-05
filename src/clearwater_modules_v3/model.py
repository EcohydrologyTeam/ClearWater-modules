"""v3 ``Model`` orchestration.

Phase 1: re-export v2's ``Model`` unchanged.

Phase 3 will replace this overlay with a v3-native ``Model`` that adds three
orchestration-level capabilities ported from v1:

- Cached compute plan and direct-array writes (kernel optimization).
- Wet-mask gating at the orchestration layer.
- Hotstart from ``xr.Dataset`` with optional per-process substep state.

It will also resolve the four chunking TODOs in v2's
``__process_loop_chunked`` by mirroring the chunking conventions used in
``clearwater_riverine``.
"""

from clearwater_modules_v2.model import Model

__all__ = ["Model"]
