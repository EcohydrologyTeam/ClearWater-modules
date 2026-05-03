"""Cell-to-core mapping for unstructured meshes.

The original SEDZLJ ``core_field.sdf`` is structured-grid (i, j)
addressed. For the ClearWater unstructured RAS mesh we need a
``core_id`` per ``nface``. Two paths supported:

1. CSV: ``core_field_unstructured.csv`` with columns
   ``Cell_Index, Core_ID``. Loaded by :func:`load_unstructured_core_map`.
2. Polygon overlay: GIS shapefile of core extents; cells assigned by
   point-in-polygon. **Phase 2 — not implemented in initial release.**
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def load_unstructured_core_map(path: Path | str, n_face: int) -> np.ndarray:
    """Read CSV mapping ``Cell_Index → Core_ID``.

    The CSV must have a header row whose first two columns are
    (case-insensitive) ``Cell_Index`` and ``Core_ID``; trailing columns
    are ignored. ``Cell_Index`` is interpreted as a 0-based ``nface``
    index.

    Returns
    -------
    np.ndarray  (n_face,) int
        Core ID for each cell. Cells absent from the CSV default to
        core 1 (matching s_sedic.f90:208 behaviour for VAR_BED=0).
    """
    out = np.ones(n_face, dtype=int)
    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh)
        rows = list(reader)
    if not rows:
        return out

    # Detect header. If the first row's first cell isn't an integer,
    # treat it as a header.
    start = 0
    try:
        int(rows[0][0])
    except (ValueError, IndexError):
        start = 1

    for row in rows[start:]:
        if len(row) < 2:
            continue
        try:
            cell = int(row[0])
            core = int(row[1])
        except ValueError:
            continue
        if 0 <= cell < n_face:
            out[cell] = core
    return out
