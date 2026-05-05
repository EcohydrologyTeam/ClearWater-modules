"""v3 utilities subpackage.

Phase 1: ``constants`` and ``conversions`` re-export from v2. Phase 1.1
adds the v3 NSM1 shared physics primitives (``reaeration``, ``sediment``,
``light``, ``partitioning``) and the new ``numerics`` module that
implements the resolved Q7 clip-with-log contract for the integrator.
"""

from clearwater_modules_v3.utils import (
    constants,
    conversions,
    light,
    numerics,
    partitioning,
    reaeration,
    sediment,
)

__all__ = [
    "constants",
    "conversions",
    "light",
    "numerics",
    "partitioning",
    "reaeration",
    "sediment",
]
