"""v3 ``Nitrogen`` process.

Phase 1: overlay re-export from v2.

The NSM1 v3 work (parallel session, ``v3-nsm1-merge`` branch) will replace
this overlay with a v3-native ``Nitrogen`` that fixes the v2 multiplicative
integrator bug and adds the legacy NSM1 nitrogen kinetics not currently in
v2 (e.g., OrgN hydrolysis and settling). See the NSM1 v3 design spec.
"""

from clearwater_modules_v2.processes.nitrogen import Nitrogen

__all__ = ["Nitrogen"]
