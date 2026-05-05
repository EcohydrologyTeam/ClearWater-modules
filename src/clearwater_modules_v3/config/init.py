"""v3 ``init_from_file`` entry point.

Phase 1: re-export v2's ``init_from_file`` unchanged.

Phase 3 will replace this overlay with a v3-native function that accepts
two additional optional top-level YAML keys, in addition to everything v2
already accepts:

- ``hotstart``: ``{dataset_path: str, timestep: str|int}``
- ``wet_mask``: ``{variable: str, threshold: float}``

When neither is present, v3 behavior matches v2 exactly (backward
compatibility with all existing v2 configurations).
"""

from clearwater_modules_v2.config.init import init_from_file

__all__ = ["init_from_file"]
