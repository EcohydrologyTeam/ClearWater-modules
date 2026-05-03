"""Bed armoring: D50_avg, τ_crit interpolation, mass-fraction sorting.

Bed armoring emerges automatically from the per-class erosion gating
in :mod:`erosion` plus mass-conservative sorting in :mod:`bed`. The
helpers here compute the surface-layer D50 and look up τ_crit on the
SCND/TAUCRITE size-class interpolant table (s_sedzlj.f90:265).

Reference: SAND2008-5621; design spec §5.3, §5.6.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from . import contracts
from .classes import SedimentClassRegistry


def compute_d50_avg(
    class_fraction: xr.DataArray,    # (nface, ssm_class) PERSED of one layer
    d50_um_array: np.ndarray,        # (ssm_class,) μm
) -> xr.DataArray:
    """Mass-weighted mean D50 of a layer.

    .. math:: D_{50,{\\rm avg}} = \\sum_s f_s \\, D_{50,s}
    """
    d50_da = xr.DataArray(
        np.asarray(d50_um_array, dtype="float64"),
        dims=(contracts.DIM_CLASS,),
    )
    return (class_fraction * d50_da).sum(dim=contracts.DIM_CLASS)


def interpolate_taucrit_from_d50(
    d50_avg_um: xr.DataArray,         # (nface,) μm
    size_interpolants_um: np.ndarray,  # (NSICM,) SCND
    taucrit_per_size_pa: np.ndarray,   # (NSICM,) TAUCRITE
) -> xr.DataArray:
    """Linear interpolation of τ_crit on the SCND/TAUCRITE table.

    Mirrors s_sedzlj.f90:265. Out-of-range D50 values are clamped to
    the table endpoints; behaviour matches s_sedzlj.f90:242–263.
    """
    xp = np.asarray(size_interpolants_um, dtype="float64")
    fp = np.asarray(taucrit_per_size_pa, dtype="float64")

    if xp.ndim != 1 or fp.ndim != 1 or xp.shape != fp.shape:
        raise ValueError(
            "size_interpolants_um and taucrit_per_size_pa must be 1-D and same length"
        )
    if xp.size == 0:
        raise ValueError("size_interpolants_um must be non-empty")

    # np.interp requires monotonically increasing xp; sort defensively
    # because the SEDflume SCND table is conventionally ascending but we
    # do not want to rely on that invariant silently.
    if np.any(np.diff(xp) < 0):
        order = np.argsort(xp)
        xp = xp[order]
        fp = fp[order]

    # np.interp clamps to fp[0] and fp[-1] for out-of-range queries —
    # this matches the SEDZLJ s_sedzlj.f90:242–263 endpoint behaviour.
    return xr.apply_ufunc(
        lambda a: np.interp(a, xp, fp),
        d50_avg_um,
        dask="parallelized",
        output_dtypes=[float],
    )
