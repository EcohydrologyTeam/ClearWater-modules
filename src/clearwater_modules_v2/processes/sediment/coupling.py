"""ESM coupling helpers.

Reads optional vegetation feedback fields from the mesh dataset (when
ESM is in the run) and exposes them in a normalized form for the
erosion module. Producer-side logic (computing biostabilization,
exposing it on the mesh) lives in
``esm/processes/biostabilization.py``.

Reference: design spec §6.5, §10.4; ESM
``io/clearwater_interface.py:406–421`` (existing TODO note).
"""

from __future__ import annotations

import xarray as xr

from . import contracts


def _select_at_time(da: xr.DataArray, time) -> xr.DataArray:
    """Select ``da`` at ``time`` if it carries a time dimension; else return as-is.

    Vegetation feedback fields are spec'd as (time, nface) but we keep
    this tolerant of static (nface,) inputs, which simplifies tests and
    permits non-time-varying ESM outputs without special-casing them at
    the call sites.
    """
    if contracts.DIM_TIME in da.dims:
        return da.sel({contracts.DIM_TIME: time})
    return da


def read_vegetation_feedback(mesh: xr.Dataset, time):
    """Return (biostabilization, root_cohesion_pa, frontal_area) at ``time``.

    Each may be ``None`` if the corresponding ESM-side variable is not
    present on the mesh; the erosion module is responsible for
    treating ``None`` as "no feedback".
    """
    bio = (
        _select_at_time(mesh[contracts.VAR_VEGETATION_BIOSTABILIZATION], time)
        if contracts.VAR_VEGETATION_BIOSTABILIZATION in mesh.data_vars
        else None
    )
    root = (
        _select_at_time(mesh[contracts.VAR_VEGETATION_ROOT_COHESION], time)
        if contracts.VAR_VEGETATION_ROOT_COHESION in mesh.data_vars
        else None
    )
    frontal = (
        _select_at_time(mesh[contracts.VAR_VEGETATION_FRONTAL_AREA], time)
        if contracts.VAR_VEGETATION_FRONTAL_AREA in mesh.data_vars
        else None
    )
    return bio, root, frontal


def read_composite_manning_n(mesh: xr.Dataset, time, fallback: xr.DataArray) -> xr.DataArray:
    """Return ``composite_manning_n`` from ESM if present; else ``fallback``
    (typically the static :data:`contracts.VAR_MANNINGS_N` from RAS)."""
    if contracts.VAR_COMPOSITE_MANNINGS_N in mesh.data_vars:
        return _select_at_time(mesh[contracts.VAR_COMPOSITE_MANNINGS_N], time)
    return fallback
