"""v3 NSM1 numerical safety-net utilities.

Implements the resolved Q7 clipping contract for the v3 NSM1 integrator
(design spec Section 5 step 4 and Section 14 Q7). Each ``Process.run``
calls ``clip_negative_state`` after the Forward Euler update; clips are
counted on a run-level ``Diagnostics`` object and structured log records
are appended with rate-limiting so that high-volume off-design events do
not flood logs. Tier 1 closed-system tests assert that ``clip_events`` is
empty under physically reasonable inputs.

The clip target is exactly 0, not a small epsilon. Monod ratios
``C / (C + K)`` are well-defined at ``C = 0``; any kinetic formula that
divides directly by a clipped state is reformulated rather than the clip
threshold being adjusted.

Defense-in-depth NaN/inf sanitization is provided by ``sanitize_rate``,
which catches both ``NaN`` and ``inf`` cells (e.g., ``x / depth`` at
``depth == 0`` produces ``inf``, not ``NaN``) and replaces them with 0
before they reach the Forward Euler update. The primary defense for dry
cells is the orchestration-layer wet-mask gating in ``Model``; the
sanitizer is the secondary defense that protects unconfigured-wet-mask
runs and also handles missing forcings that ``xr.DataArray`` operations
turn into ``NaN``.
"""

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import xarray as xr

from clearwater_data.custom_types import ArrayLike


_DEFAULT_DETAIL_LIMIT = 10


@dataclass
class Diagnostics:
    """Run-level diagnostics container for the v3 NSM1 integrator.

    Attributes:
        clip_events: per-state-variable count of cells clipped to zero,
            aggregated across all calls to ``clip_negative_state`` for the
            run. Tier 1 closed-system tests assert this dict is empty.
        clip_log: structured log records for the first
            ``detail_limit_per_call`` clips per call. Each record is a dict
            with keys ``name``, ``cell_index``, ``value_before``, and
            ``step``. Subsequent clips in the same call contribute only to
            the aggregate count, not to the detailed log.
        detail_limit_per_call: maximum number of structured detail records
            emitted per call to ``clip_negative_state``; defaults to 10 per
            the Q7 contract.
        current_step: integrator substep counter maintained by ``Model``;
            ``clip_negative_state`` reads this as the default value for
            the ``step`` field of every log record. Resolves
            pattern-alignment spec §10 Q1: the step index source is the
            ``Diagnostics`` field set by ``Model``'s substep loop, not
            a per-Process counter or a kwarg passed through every
            ``run`` signature. ``-1`` is the "no step set" sentinel; the
            log record then reads ``step=None`` for back-compat with the
            1.0.0 ``step=None`` convention.
    """

    clip_events: dict[str, int] = field(default_factory=dict)
    clip_log: list[dict[str, Any]] = field(default_factory=list)
    detail_limit_per_call: int = _DEFAULT_DETAIL_LIMIT
    current_step: int = -1


def clip_negative_state(
    state: ArrayLike,
    name: str,
    diagnostics: "Diagnostics | None" = None,
    step: int | None = None,
) -> ArrayLike:
    """Clip negative values of a state variable to exactly zero.

    Implements the v3 NSM1 clip-with-log contract. Negative cells are
    counted in ``diagnostics.clip_events[name]``; up to
    ``diagnostics.detail_limit_per_call`` per-cell records are appended to
    ``diagnostics.clip_log`` and any remaining clips contribute only to a
    single aggregate ``clipped_aggregate`` record. The clip target is
    exactly 0.

    Resolves pattern-alignment spec §10 Q2 (graceful no-op when
    ``diagnostics is None``) and Q1 (default step index comes from
    ``diagnostics.current_step``):

    * Callers may pass ``diagnostics=None`` to clip without
      counting/logging. The container-type return contract is preserved
      (``xr.DataArray`` in → ``xr.DataArray`` out; ``np.ndarray`` in →
      ``np.ndarray`` out; scalar in → scalar out). This lets every
      ``Process.run`` invoke ``clip_negative_state(state, name,
      self.diagnostics)`` without an ``isinstance`` / ``is not None``
      guard branch.
    * Callers may omit ``step``; when ``diagnostics`` is provided and
      ``step is None``, the default is taken from
      ``diagnostics.current_step`` (sentinel ``-1`` is normalised to
      ``None`` to match the 1.0.0 log convention). ``Model``'s substep
      loop sets ``current_step`` on every iteration.

    Args:
        state | ArrayLike | post-Forward-Euler state variable to clip.
            Accepts ``xr.DataArray``, ``np.ndarray``, or Python scalar.
        name | str | state-variable name used as the diagnostics key.
        diagnostics | Diagnostics or None | run-level diagnostics
            container, mutated in place. ``None`` skips counting/logging
            but still clips.
        step | int or None | optional integrator step index for log
            records. ``None`` defers to ``diagnostics.current_step``.

    Returns:
        Same container type as ``state``, with negative cells replaced
        by ``0``.
    """
    # Container-type-aware negative-mask + clip.
    if isinstance(state, xr.DataArray):
        values = state.values
        negative_mask = values < 0.0
    elif isinstance(state, np.ndarray):
        values = state
        negative_mask = state < 0.0
    else:
        # Python scalar path.
        if state < 0.0:
            if diagnostics is not None:
                diagnostics.clip_events[name] = (
                    diagnostics.clip_events.get(name, 0) + 1
                )
                effective_step = (
                    step
                    if step is not None
                    else (
                        diagnostics.current_step
                        if diagnostics.current_step >= 0
                        else None
                    )
                )
                diagnostics.clip_log.append(
                    {
                        "name": name,
                        "cell_index": (),
                        "value_before": float(state),
                        "step": effective_step,
                    }
                )
            return type(state)(0.0)
        return state

    n_clipped = int(negative_mask.sum())
    if n_clipped == 0:
        return state

    if diagnostics is not None:
        diagnostics.clip_events[name] = (
            diagnostics.clip_events.get(name, 0) + n_clipped
        )

        effective_step = (
            step
            if step is not None
            else (
                diagnostics.current_step
                if diagnostics.current_step >= 0
                else None
            )
        )

        detail_limit = diagnostics.detail_limit_per_call
        flat_indices = np.flatnonzero(negative_mask.ravel())
        detail_indices = flat_indices[:detail_limit]
        multi_index = np.unravel_index(detail_indices, values.shape)
        flat_values = values.ravel()
        for record_position, flat_index in enumerate(detail_indices):
            cell_index = tuple(int(axis[record_position]) for axis in multi_index)
            diagnostics.clip_log.append(
                {
                    "name": name,
                    "cell_index": cell_index,
                    "value_before": float(flat_values[flat_index]),
                    "step": effective_step,
                }
            )

        if n_clipped > detail_limit:
            diagnostics.clip_log.append(
                {
                    "name": name,
                    "cell_index": "clipped_aggregate",
                    "value_before": None,
                    "step": effective_step,
                    "n_suppressed": n_clipped - detail_limit,
                }
            )

    clipped_values = np.where(negative_mask, 0.0, values)
    if isinstance(state, xr.DataArray):
        return xr.DataArray(
            clipped_values,
            coords=state.coords,
            dims=state.dims,
            name=state.name,
            attrs=state.attrs,
        )
    # np.ndarray path.
    return clipped_values


def sanitize_rate(rate: ArrayLike) -> ArrayLike:
    """Replace NaN and inf cells with 0; preserves container type.

    Defense-in-depth guard for Process rate computations. Catches:

    * ``NaN`` from missing forcings or from undefined arithmetic
      (e.g., ``0 / 0`` propagated through ``xr.DataArray`` operations).
    * ``inf`` from division by zero (e.g., ``x / depth`` at
      ``depth == 0``).

    The primary defense for dry cells is the orchestration-layer
    wet-mask gating in ``Model.__apply_wet_mask`` (it overwrites
    Process outputs with ``NaN`` after the Process runs at any cell
    where the wet-mask threshold is not met). ``sanitize_rate`` is
    the secondary defense: it protects unconfigured-wet-mask runs
    and Tier 1 unit tests where wet-mask is not active. It also
    means that even at wet cells, a transient ``NaN`` or ``inf``
    in any sub-rate term cannot poison the Forward Euler update.

    Handles the three input container types Process rates may
    appear as in v3:

    * ``xr.DataArray`` -- canonical Process input/output container.
    * ``np.ndarray`` -- intermediate result of mixed numpy
      operations or a stand-alone test invocation.
    * native Python scalar (``float`` or ``int``) -- single-cell
      unit-test invocations or constant-input scenarios.

    Args:
        rate | ArrayLike | per-cell rate of change to sanitize.

    Returns:
        ArrayLike | the same container type as ``rate`` with any
        ``NaN`` or ``inf`` cells replaced by exactly 0.
    """
    if isinstance(rate, xr.DataArray):
        return xr.where(rate.isnull() | np.isinf(rate), 0.0, rate)
    if isinstance(rate, np.ndarray):
        return np.where(np.isnan(rate) | np.isinf(rate), 0.0, rate)
    # Native Python scalar (or any other duck-typed numeric value).
    if np.isnan(rate) or np.isinf(rate):
        return 0.0
    return rate
