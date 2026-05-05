"""v3 configuration package.

Phase 1: re-exports ``init_from_file`` and ``read_config`` from v2 unchanged.
Phase 3 will add support for two optional top-level YAML keys (``hotstart``
and ``wet_mask``) without breaking backward compatibility with v2 configs.
"""

from clearwater_modules_v3.config.init import init_from_file
from clearwater_modules_v2.config.read import read_config

__all__ = ["init_from_file", "read_config"]
