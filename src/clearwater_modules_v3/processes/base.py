"""v3 ``Process`` and ``ProcessFactory`` abstractions.

Re-export from v2 unchanged. The v2 ``Process`` contract is the
architectural baseline for v3 (architecture spec section 3 non-goal: no
new framework).

------------------------------------------------------------------
Integration patterns (guideline, not enforced contract)
------------------------------------------------------------------

v3-native processes implement ``run(time, registry)`` to advance state
by one substep. Two integration patterns are common:

(a) Compute a per-second rate of change and apply Forward Euler
    ``state_new = state_old + rate * time_step_seconds``. This is
    suitable for processes with simple linear kinetics (e.g., NSM1
    reaction terms).
(b) Compute the per-substep ``delta_state`` directly and add it to the
    current state. This is suitable for processes whose update depends
    non-linearly on the substep length (e.g., the v3 ``Temperature``
    thin-water guards).

Both patterns are valid; the choice belongs to the process author. M16
(review-findings 2026-05-04) demoted an earlier numbered "5-step
contract" to this guideline because real v3 processes (notably
``Temperature``) follow pattern (b) and the contract wording did not
reflect actual practice.

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
