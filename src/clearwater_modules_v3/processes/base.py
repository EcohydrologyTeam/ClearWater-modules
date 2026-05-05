"""v3 ``Process`` and ``ProcessFactory`` abstractions.

Re-export from v2 unchanged. The v2 ``Process`` contract is the
architectural baseline for v3 (architecture spec section 3 non-goal: no
new framework). v3 establishes an *integrator-pattern contract* on top
of this class that future v3-native processes follow:

1. ``Process.run`` reads state from the registry at the current time.
2. Computes a net rate of change with units of ``[state] / second``.
3. Applies the rate via Forward Euler:
   ``state_new = state_old + rate * self.time_step.total_seconds()``.
4. Writes the updated state back via ``registry.set_at_time``.
5. Applies negative-state guards with ``xr.where`` where appropriate.

The contract is documented in the umbrella architecture spec section 4.

------------------------------------------------------------------
M5 ordering contract for ``init_process`` and ``from_hotstart``
------------------------------------------------------------------

The v3 ``Model.__init_model`` invokes per-process initialization in a
specific ordering that processes MUST respect when they implement
hotstart support:

1. Data sources (or first chunk) are loaded into the registry.
2. If a hotstart dataset is supplied, the registry is **seeded** from
   the saved dataset (``Model.__seed_from_hotstart``).
3. ``init_process`` is called on every process.
4. If a hotstart dataset is supplied, ``from_hotstart`` is called on
   every process that defines it.

The ordering means ``init_process`` runs **before** ``from_hotstart``,
and so:

- ``init_process`` MUST set the process's substep-internal state
  assuming a *fresh start*. It does not need to know whether a
  hotstart was supplied; it sets defaults appropriate to a fresh run.
- ``from_hotstart`` MUST override those fresh-start defaults with
  values restored from the saved dataset's ``attrs``. If a process
  adds new internal substep state in a future revision but forgets to
  also handle it in ``from_hotstart``, fresh-start runs and
  hotstart-resume runs will silently diverge starting at the first
  substep where the un-restored state matters.

Together, the two methods ensure that fresh-start and
hotstart-resume produce equivalent post-init state when the same
initial conditions are loaded. This invariant is a load-bearing
property of the v3 hotstart design (TSM design spec section 3.2).

The default ``Process.from_hotstart`` is implicit (``getattr`` -based
in the Model). Processes opt in by defining the method; processes
that don't have any internal substep state are free to omit it.
"""

from clearwater_modules_v2.processes.base import Process, ProcessFactory

__all__ = ["Process", "ProcessFactory"]
