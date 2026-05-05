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
"""

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import xarray as xr


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
    """

    clip_events: dict[str, int] = field(default_factory=dict)
    clip_log: list[dict[str, Any]] = field(default_factory=list)
    detail_limit_per_call: int = _DEFAULT_DETAIL_LIMIT


def clip_negative_state(
    state: xr.DataArray,
    name: str,
    diagnostics: Diagnostics,
    step: int | None = None,
) -> xr.DataArray:
    """Clip negative values of a state variable to exactly zero.

    Implements the v3 NSM1 clip-with-log contract. Negative cells are
    counted in ``diagnostics.clip_events[name]``; up to
    ``diagnostics.detail_limit_per_call`` per-cell records are appended to
    ``diagnostics.clip_log`` and any remaining clips contribute only to a
    single aggregate ``clipped_aggregate`` record. The clip target is
    exactly 0.

    Args:
        state | DataArray | post-Forward-Euler state variable to clip.
        name | str | state-variable name used as the diagnostics key.
        diagnostics | Diagnostics | run-level diagnostics container,
            mutated in place.
        step | int or None | optional integrator step index for log records.

    Returns:
        DataArray | clipped state with negative cells replaced by 0.
    """
    values = state.values
    negative_mask = values < 0.0
    n_clipped = int(negative_mask.sum())

    if n_clipped == 0:
        return state

    diagnostics.clip_events[name] = diagnostics.clip_events.get(name, 0) + n_clipped

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
                "step": step,
            }
        )

    if n_clipped > detail_limit:
        diagnostics.clip_log.append(
            {
                "name": name,
                "cell_index": "clipped_aggregate",
                "value_before": None,
                "step": step,
                "n_suppressed": n_clipped - detail_limit,
            }
        )

    clipped_values = np.where(negative_mask, 0.0, values)
    return xr.DataArray(
        clipped_values,
        coords=state.coords,
        dims=state.dims,
        name=state.name,
        attrs=state.attrs,
    )
