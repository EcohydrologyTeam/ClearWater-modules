"""v3 ``FloatingAlgae`` process.

Phase 1: overlay re-export from v2. NSM1 v3 work replaces this overlay
later (multiplicative-integrator and stray ``* 86400`` fixes, light
extinction, full nutrient coupling).
"""

from clearwater_modules_v2.processes.floating_algae import FloatingAlgae

__all__ = ["FloatingAlgae"]
