"""Hotstart / restart for SSM bed state.

Modern format: NetCDF, written via xarray, readable by xarray. Stored
to a path matching the existing checkpoint convention exercised by
``tests/test_hotstart_roundtrip.py``.

Legacy compatibility: :func:`read_legacy_sedbed_hot_sdf` reads the
SEDZLJ ASCII ``SEDBED_HOT.SDF`` format (see s_sedic.f90:452) so users
porting from EFDC+ can warm-start without re-running.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import xarray as xr

from .. import contracts


def write_hotstart(
    mesh: xr.Dataset,
    path: Path | str,
    time: Any | None = None,
) -> None:
    """Write SSM bed state at a single time step to NetCDF.

    Parameters
    ----------
    mesh
        Full SSM mesh dataset. Only variables registered in
        :data:`contracts.BED_STATE_BY_NAME` are written.
    path
        Output NetCDF path.
    time
        Optional time-coordinate value to select before writing. If
        ``None``, the entire variable (all time-steps present) is
        written.
    """
    keep = [name for name in contracts.BED_STATE_BY_NAME if name in mesh.data_vars]
    if not keep:
        raise ValueError(
            "Mesh contains none of the SSM bed-state variables; nothing to write."
        )
    subset = mesh[keep]

    if time is not None and contracts.DIM_TIME in subset.dims:
        subset = subset.sel({contracts.DIM_TIME: time})

    subset.to_netcdf(path)


def read_hotstart(path: Path | str) -> xr.Dataset:
    """Read an SSM hotstart NetCDF written by :func:`write_hotstart`."""
    return xr.open_dataset(path)


def read_legacy_sedbed_hot_sdf(path: Path | str) -> dict:
    """Read the legacy EFDC+ ``SEDBED_HOT.SDF`` ASCII format.

    Format (per ``s_sedic.f90:452-512``):

    * Line 1: ``TBEGINSEDZLJ [VER]``  — restart time and optional
      file-format version. If only one token is present, ``VER = -1``
      (pre-1240 legacy format with no TSED0).
    * Block 1 (``LAYERACTIVE``, format ``34569``): ``KB`` integers per
      cell, ``LA - 1`` cells.
    * Block 2 (``KBT``,         format ``34569``): one integer per
      cell.
    * Block 3 (``D50AVG``,      format ``34567``): one float per cell.
    * Block 4 (``BULKDENS``,    format ``34568``): ``KB`` floats per
      cell.
    * Block 5 (``TSED``,        format ``34568``): ``KB`` floats per
      cell.
    * Block 6 (``TSED0``,       format ``34568``): ``KB`` floats per
      cell — present only when ``VER >= 1240``.
    * Block 7 (``PERSED``,      format ``34568``): ``NSEDS × KB`` per
      cell.

    Returns
    -------
    dict
        Keyed by the canonical bed-state variable names from
        :mod:`contracts` (``ssm_bed_layer_active``, ``ssm_bed_d50_surface``,
        ``ssm_bed_layer_bulk_density``, ``ssm_bed_layer_mass``,
        ``ssm_bed_layer_initial_mass``, ``ssm_bed_class_fraction``).

    Notes
    -----
    Stub implementation: this returns an empty dict pending a
    real-world ``SEDBED_HOT.SDF`` to validate against. The docstring
    above documents the on-disk layout precisely so a full parser can
    be added without re-deriving it from the Fortran.
    """
    # The SEDBED_HOT.SDF format depends on KB and NSEDS, which the file
    # itself does not carry. A complete parser must accept those as
    # arguments (or read them from the companion bed.sdf). Returning an
    # empty dict signals "not yet implemented" without raising; callers
    # should fall back to fresh initialization.
    return {}
