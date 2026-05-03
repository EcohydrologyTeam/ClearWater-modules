"""SSM — Sediment Simulation Module (EFDC SEDZLJ port).

Public API:
    from clearwater_modules_v2.processes.sediment import SSM, SedimentClass

The full design specification is at
``ClearWater-Riverine-streaming/design/ssm_design_spec.md``.
"""

from .classes import SedimentClass, SedimentClassRegistry
from .ssm import SSM

__all__ = ["SSM", "SedimentClass", "SedimentClassRegistry"]
