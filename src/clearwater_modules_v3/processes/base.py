"""v3 ``Process`` and ``ProcessFactory`` abstractions.

Re-export from v2 unchanged. The v2 ``Process`` contract is the
architectural baseline for v3 (architecture spec §3 non-goal: no new
framework). v3 establishes an *integrator-pattern contract* on top of this
class that future v3-native processes follow:

1. ``Process.run`` reads state from the registry at the current time.
2. Computes a net rate of change with units of ``[state] / second``.
3. Applies the rate via Forward Euler:
   ``state_new = state_old + rate * self.time_step.total_seconds()``.
4. Writes the updated state back via ``registry.set_at_time``.
5. Applies negative-state guards with ``xr.where`` where appropriate.

The contract is documented in the umbrella architecture spec §4.
"""

from clearwater_modules_v2.processes.base import Process, ProcessFactory

__all__ = ["Process", "ProcessFactory"]
