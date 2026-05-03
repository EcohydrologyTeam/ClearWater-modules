"""SEDflume input-file parsers (verbatim SEDZLJ format).

Three files per SAND2008-5621 §"SEDZLJ Sediment Input Files":

* ``bed.sdf``   — global parameters, per-class properties, per-interpolant
                  properties, ENRATE table or (A, n, max_rate) tuples.
                  See SAND2008-5621 Figure 3.
* ``erate.sdf`` — per-core SEDflume data: per-layer thicknesses, critical
                  shear stresses, bulk densities, water/sediment density,
                  particle size distribution, per-shear-level erosion
                  rates. See SAND2008-5621 Figure 1.
* ``core_field.sdf`` — integer matrix mapping (i, j) cell to core ID.
                  Both DSI standard and SNL standard formats are
                  supported per s_sedic.f90:151–180.

The bundled tests (``test_io_sedflume.py``) round-trip the SAND2008
example datasets to verify byte-for-byte fidelity isn't required, but
all numerical content survives load → save → reload.

Format notes / interpretive choices
-----------------------------------

The SAND2008-5621 manual and the production EFDC+ source
(``s_sedic.f90``) disagree on the column layout of the ``bed.sdf``
header card. We follow the **manual** (Figure 3) as the authoritative
reference because the user-facing input format is what we must remain
backward-compatible with; modern EFDC+ has reordered the same fields.
Specifically:

* ``bed.sdf`` line 2 (after the header comment) contains, per
  Figure 3's column-header comment:
  ``VAR_BED  Bedload  Nequil  KB  ISEDTIME  IMORPH  IFWAVE  MAXDEPLIMIT``.
  We map these into the bundle as
  ``var_bed, icalc_bl, _, n_layers, _, _, _, max_deposit_limit``.
* Comment lines start with ``#`` or ``*`` (after optional leading
  whitespace). Inline ``!`` and ``#`` comments are also stripped from
  data lines. Blank lines are skipped.
* All critical-shear values in the file are in dynes/cm² and converted
  to Pa (× 0.1) on read. D50 is left in μm. Densities are left in
  g/cm³.
* ``NSEDFLUME`` (table vs. power-law) is not present in the file
  itself — EFDC+ reads it from a separate ``efdc.inp`` card. We
  infer it from the trailing ENRATE-block layout: ``NSICM`` rows
  of ``ITBM`` columns ⇒ ``nsedflume = 1``; ``NSICM`` rows of 3
  columns ⇒ ``nsedflume = 2``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


# ---------------------------------------------------------------------------
# Unit conversion
# ---------------------------------------------------------------------------

#: Conversion factor: 1 dyne/cm² = 0.1 Pa.
DYNES_CM2_TO_PA: float = 0.1


# ---------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------


@dataclass
class SedflumeBundle:
    """All data parsed from the three SEDflume files, in SI / clean form.

    Critical-shear values are converted from dynes/cm² to Pa (× 0.1) on
    read; D50 is kept in μm; densities in g/cm³.
    """

    # From bed.sdf
    n_layers: int                                  # KB
    var_bed: int                                   # 0 or 1
    icalc_bl: int                                  # bedload flag
    nsedflume: int                                 # 1 (table) or 2 (power-law)
    zb_skin_um: float
    tau_const_pa: float
    bedload_cutoff_um: float
    max_deposit_limit: float
    d50_um: np.ndarray                             # (n_class,)
    tau_ce_pa: np.ndarray                          # (n_class,)
    tau_cs_pa: np.ndarray                          # (n_class,)
    settling_cm_s: np.ndarray                      # (n_class,) — -1 sentinel = compute
    size_interpolants_um: np.ndarray               # (NSICM,)   — SCND
    taucrit_per_size_pa: np.ndarray                # (NSICM,)   — TAUCRITE
    erate_active_table: np.ndarray | None          # (NSICM, ITBM) cm/s — for nsedflume=1
    actdep_a: np.ndarray | None                    # (NSICM,) — for nsedflume=2
    actdep_n: np.ndarray | None                    # (NSICM,)
    actdep_max: np.ndarray | None                  # (NSICM,)

    # From erate.sdf
    n_cores: int
    layer_thickness_cm: np.ndarray                 # (n_cores, n_layers)
    layer_taucrit_pa: np.ndarray                   # (n_cores, n_layers)
    bulk_density_g_cm3: np.ndarray                 # (n_cores, n_layers)
    water_density_g_cm3: float
    solid_density_g_cm3: float
    particle_size_distribution_pct: np.ndarray     # (n_cores, n_layers, n_class)
    tau_levels_pa: np.ndarray                      # (ITBM,)
    erate_per_core_cm_s: np.ndarray | None         # (n_cores, n_layers, ITBM) — nsedflume=1
    ea_per_core: np.ndarray | None                 # (n_cores, n_layers) — nsedflume=2
    en_per_core: np.ndarray | None                 # (n_cores, n_layers)
    max_rate_per_core_cm_s: np.ndarray | None      # (n_cores, n_layers)

    # From core_field.sdf (structured-grid form; unstructured mapping
    # produced separately by mesh_mapping.load_unstructured_core_map)
    core_field_ij: np.ndarray | None               # (jc, ic) int


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------


def _strip_inline_comment(line: str) -> str:
    """Remove inline ``!`` or ``#`` trailing comments while preserving
    leading content. Comments that span the whole line are detected
    upstream by :func:`_iter_data_lines`."""
    for marker in ("!", "#"):
        idx = line.find(marker)
        if idx >= 0:
            line = line[:idx]
    return line


def _iter_data_lines(path: Path | str):
    """Yield non-comment, non-blank lines from a SEDZLJ-format file.

    Comment lines start with ``#``, ``*``, ``!``, or ``C``/``c`` (the
    Fortran column-1 comment convention) after optional leading
    whitespace. Inline ``!`` and ``#`` comments are stripped from
    data lines.
    """
    with open(path, "r") as fh:
        for raw in fh:
            stripped = raw.strip()
            if not stripped:
                continue
            first = stripped[0]
            if first in ("#", "*", "!"):
                continue
            # Treat a leading "C " or "c " (Fortran fixed-form column-1
            # comment) as a comment, but NOT a leading "C123" token
            # which could be a data label.
            if first in ("C", "c") and (len(stripped) == 1 or not stripped[1].isdigit()):
                continue
            cleaned = _strip_inline_comment(stripped).strip()
            if not cleaned:
                continue
            yield cleaned


def _take_floats(line: str) -> list[float]:
    """Whitespace-tokenize, drop trailing comments, return floats."""
    tokens = _strip_inline_comment(line).split()
    return [float(t) for t in tokens]


def _take_ints(line: str) -> list[int]:
    """Whitespace-tokenize, drop trailing comments, return ints. Floats
    that happen to be written ``1.``  / ``1.0`` are accepted."""
    tokens = _strip_inline_comment(line).split()
    return [int(float(t)) for t in tokens]


# ---------------------------------------------------------------------------
# bed.sdf
# ---------------------------------------------------------------------------


def parse_bed_sdf(path: Path | str) -> dict:
    """Parse ``bed.sdf`` per SAND2008-5621 Figure 3 layout.

    Returns
    -------
    dict
        Mapping containing every :class:`SedflumeBundle` field that
        ``bed.sdf`` populates. Two top-level shapes are possible
        depending on the file's ``NSEDFLUME`` mode (inferred from the
        layout — see notes below).

    Notes
    -----
    The bed.sdf format does not carry an explicit ``NSEDFLUME`` flag.
    EFDC+ infers it from a separate ``efdc.inp`` card, which we don't
    consume. We use the following heuristic:

    * If the post-TAUCRITE block has ``KB`` lines whose token count
      matches ``ITBM``, NSEDFLUME=1.
    * Else if every line has exactly 3 tokens, NSEDFLUME=2.

    The number of sediment classes ``NSEDS`` is inferred from the
    length of the ``D50`` line.
    """
    lines = list(_iter_data_lines(path))
    cursor = 0

    # --- Header card ---------------------------------------------------
    # Per SAND2008-5621 Figure 3 column layout:
    #   VAR_BED  Bedload  Nequil  KB  ISEDTIME  IMORPH  IFWAVE  MAXDEPLIMIT
    header_tokens = _strip_inline_comment(lines[cursor]).split()
    cursor += 1
    var_bed = int(float(header_tokens[0]))
    icalc_bl = int(float(header_tokens[1]))
    # header_tokens[2] is NEQUIL — toxics flag, not stored.
    n_layers = int(float(header_tokens[3]))
    # MAXDEPLIMIT is the last number on the line per Figure 3.
    max_deposit_limit = float(header_tokens[-1])

    # --- ITBM, NSICM ---------------------------------------------------
    itbm_nsicm = _take_ints(lines[cursor])
    cursor += 1
    itbm = itbm_nsicm[0]
    nsicm = itbm_nsicm[1]

    # --- ZBSKIN, TAUCONST, ISSLOPE (optional), BEDLOAD_CUTOFF ----------
    zb_line = _take_floats(lines[cursor])
    cursor += 1
    zb_skin_um = zb_line[0]
    tau_const_dynes = zb_line[1] if len(zb_line) > 1 else 0.0
    # Layout in Figure 3: only ZBSKIN + TAUCONST. EFDC+ source reads
    # ZBSKIN, TAUCONST, ISSLOPE, BEDLOAD_CUTOFF. Default cutoff to 0
    # when not specified — caller can fall back to the 64 μm default.
    if len(zb_line) >= 4:
        bedload_cutoff_um = zb_line[3]
    else:
        bedload_cutoff_um = 0.0
    tau_const_pa = tau_const_dynes * DYNES_CM2_TO_PA

    # --- D50 line -> infer NSEDS --------------------------------------
    d50_um = np.asarray(_take_floats(lines[cursor]), dtype=float)
    cursor += 1
    n_class = d50_um.size

    tau_ce_dynes = np.asarray(_take_floats(lines[cursor])[:n_class], dtype=float)
    cursor += 1
    tau_cs_dynes = np.asarray(_take_floats(lines[cursor])[:n_class], dtype=float)
    cursor += 1
    settling_cm_s = np.asarray(_take_floats(lines[cursor])[:n_class], dtype=float)
    cursor += 1

    tau_ce_pa = tau_ce_dynes * DYNES_CM2_TO_PA
    tau_cs_pa = tau_cs_dynes * DYNES_CM2_TO_PA

    # --- Per-interpolant block (NSICM size) ----------------------------
    size_interpolants_um = np.asarray(_take_floats(lines[cursor])[:nsicm], dtype=float)
    cursor += 1
    taucrit_dynes = np.asarray(_take_floats(lines[cursor])[:nsicm], dtype=float)
    cursor += 1
    taucrit_per_size_pa = taucrit_dynes * DYNES_CM2_TO_PA

    # --- ENRATE table or (A, n, max_rate) tuples ----------------------
    remaining_blocks = []
    while cursor < len(lines):
        remaining_blocks.append(_take_floats(lines[cursor]))
        cursor += 1

    erate_active_table: np.ndarray | None = None
    actdep_a: np.ndarray | None = None
    actdep_n: np.ndarray | None = None
    actdep_max: np.ndarray | None = None

    # Determine NSEDFLUME mode from the trailing block layout.
    if (
        len(remaining_blocks) >= nsicm
        and all(len(row) == itbm for row in remaining_blocks[:nsicm])
    ):
        # NSEDFLUME == 1: ENRATE table NSICM × ITBM.
        nsedflume = 1
        erate_active_table = np.asarray(
            remaining_blocks[:nsicm], dtype=float
        )
    elif (
        len(remaining_blocks) >= nsicm
        and all(len(row) == 3 for row in remaining_blocks[:nsicm])
    ):
        # NSEDFLUME == 2: (A, n, max_rate) per interpolant.
        nsedflume = 2
        block = np.asarray(remaining_blocks[:nsicm], dtype=float)
        actdep_a = block[:, 0]
        actdep_n = block[:, 1]
        actdep_max = block[:, 2]
    else:
        raise ValueError(
            f"bed.sdf trailing block has unexpected layout "
            f"(NSICM={nsicm}, ITBM={itbm}, "
            f"got {len(remaining_blocks)} rows of widths "
            f"{[len(r) for r in remaining_blocks]!r})"
        )

    out: dict[str, Any] = {
        "n_layers": n_layers,
        "var_bed": var_bed,
        "icalc_bl": icalc_bl,
        "nsedflume": nsedflume,
        "zb_skin_um": zb_skin_um,
        "tau_const_pa": tau_const_pa,
        "bedload_cutoff_um": bedload_cutoff_um,
        "max_deposit_limit": max_deposit_limit,
        "d50_um": d50_um,
        "tau_ce_pa": tau_ce_pa,
        "tau_cs_pa": tau_cs_pa,
        "settling_cm_s": settling_cm_s,
        "size_interpolants_um": size_interpolants_um,
        "taucrit_per_size_pa": taucrit_per_size_pa,
        # ITBM exposed for the erate.sdf parser
        "_itbm": itbm,
        "_nsicm": nsicm,
        "_n_class": n_class,
    }
    if nsedflume == 1:
        out["erate_active_table"] = erate_active_table
    else:
        out["actdep_a"] = actdep_a
        out["actdep_n"] = actdep_n
        out["actdep_max"] = actdep_max
    return out


# ---------------------------------------------------------------------------
# erate.sdf
# ---------------------------------------------------------------------------


def parse_erate_sdf(
    path: Path | str,
    n_layers: int,
    n_class: int,
    nsedflume: int,
) -> dict:
    """Parse ``erate.sdf`` per SAND2008-5621 Figure 1 layout.

    The first data line in ``erate.sdf`` is the active-layer multiplier
    (TACTM, scalar). Then for each core (loop INCORE times):

    1. critical shear stresses per layer (1 line, KB values)
       — converted from dynes/cm² to Pa.
    2. layer thicknesses per layer (1 line, KB values, cm).
    3. bulk densities per layer (1 line, KB values, g/cm³).
    4. water density and sediment solid density (1 line, 2 values, g/cm³).
    5. particle size distribution (KB lines, NSEDS values per line —
       mass percentages 0–100).
    6. erosion-rate block, format depends on ``nsedflume``:

       * ``nsedflume == 1``: ITBM × (shear-stress line + erosion-rate
         line of KB values), interleaved.
       * ``nsedflume == 2``: KB lines of ``(A, n, max_rate)`` per layer.

    The number of cores (INCORE) is inferred from the file size: the
    parser keeps reading per-core blocks until the line stream is
    exhausted.
    """
    lines = list(_iter_data_lines(path))
    cursor = 0

    # --- TACTM (active-layer multiplier) -------------------------------
    # We read it but the bundle currently has no explicit field for it.
    tactm_tokens = _take_floats(lines[cursor])
    cursor += 1
    tactm = tactm_tokens[0] if tactm_tokens else None

    # --- Per-core blocks ------------------------------------------------
    # We collect into Python lists then stack at the end; INCORE is
    # determined dynamically.
    taucrit_per_core: list[np.ndarray] = []
    thickness_per_core: list[np.ndarray] = []
    bulkdens_per_core: list[np.ndarray] = []
    psd_per_core: list[np.ndarray] = []
    erate_per_core: list[np.ndarray] = []      # nsedflume==1
    ea_per_core: list[np.ndarray] = []         # nsedflume==2
    en_per_core: list[np.ndarray] = []         # nsedflume==2
    maxrate_per_core: list[np.ndarray] = []    # nsedflume==2
    water_density: float | None = None
    solid_density: float | None = None
    tau_levels_pa: np.ndarray | None = None

    while cursor < len(lines):
        # Per-layer τ_crit (KB values)
        tau_dynes = np.asarray(_take_floats(lines[cursor])[:n_layers], dtype=float)
        cursor += 1
        taucrit_per_core.append(tau_dynes * DYNES_CM2_TO_PA)

        # Per-layer thickness
        thickness = np.asarray(_take_floats(lines[cursor])[:n_layers], dtype=float)
        cursor += 1
        thickness_per_core.append(thickness)

        # Per-layer bulk density
        bulk = np.asarray(_take_floats(lines[cursor])[:n_layers], dtype=float)
        cursor += 1
        bulkdens_per_core.append(bulk)

        # Water density, solid density (1 line, 2 numbers)
        wd = _take_floats(lines[cursor])
        cursor += 1
        water_density = wd[0]
        solid_density = wd[1]

        # Particle size distribution: KB lines × n_class
        psd_block = np.zeros((n_layers, n_class), dtype=float)
        for k in range(n_layers):
            row = _take_floats(lines[cursor])[:n_class]
            cursor += 1
            psd_block[k, :] = row
        psd_per_core.append(psd_block)

        # Erosion-rate block ------------------------------------------
        if nsedflume == 1:
            # ITBM blocks of (shear, KB-erate)
            tau_levels: list[float] = []
            erate_block = np.zeros(
                (n_layers, 0), dtype=float
            )  # placeholder; widened below
            erate_cols: list[np.ndarray] = []
            # Read until next core's first τ_crit line — but since the
            # number of (shear, erate) pairs equals ITBM and ITBM is
            # only known via bed.sdf, we instead read pairs until we
            # detect a layout boundary. The simplest robust rule: keep
            # reading pairs as long as the "shear" line is a single
            # number AND the next line has KB tokens.
            while cursor + 1 < len(lines):
                shear_tokens = _take_floats(lines[cursor])
                if len(shear_tokens) != 1:
                    break
                next_tokens = _take_floats(lines[cursor + 1])
                if len(next_tokens) != n_layers:
                    break
                tau_levels.append(shear_tokens[0])
                erate_cols.append(np.asarray(next_tokens[:n_layers], dtype=float))
                cursor += 2
            erate_block = np.stack(erate_cols, axis=1) if erate_cols else \
                np.zeros((n_layers, 0), dtype=float)
            # Convert τ levels (Pa already per SAND2008 Figure 1
            # caption "applied shears are 0, 2, 4, 8, 10, 20 Pa").
            this_tau_levels = np.asarray(tau_levels, dtype=float)
            if tau_levels_pa is None:
                tau_levels_pa = this_tau_levels
            erate_per_core.append(erate_block)
        else:
            # nsedflume == 2: KB lines of (A, n, max_rate)
            ea = np.zeros(n_layers, dtype=float)
            en = np.zeros(n_layers, dtype=float)
            mx = np.zeros(n_layers, dtype=float)
            for k in range(n_layers):
                triple = _take_floats(lines[cursor])
                cursor += 1
                ea[k] = triple[0]
                en[k] = triple[1]
                mx[k] = triple[2]
            ea_per_core.append(ea)
            en_per_core.append(en)
            maxrate_per_core.append(mx)
            if tau_levels_pa is None:
                tau_levels_pa = np.array([0.0, 1000.0])

    n_cores = len(taucrit_per_core)
    out: dict[str, Any] = {
        "n_cores": n_cores,
        "layer_taucrit_pa": np.stack(taucrit_per_core, axis=0),
        "layer_thickness_cm": np.stack(thickness_per_core, axis=0),
        "bulk_density_g_cm3": np.stack(bulkdens_per_core, axis=0),
        "water_density_g_cm3": float(water_density) if water_density is not None else 1.0,
        "solid_density_g_cm3": float(solid_density) if solid_density is not None else 2.65,
        "particle_size_distribution_pct": np.stack(psd_per_core, axis=0),
        "tau_levels_pa": tau_levels_pa if tau_levels_pa is not None
            else np.zeros(0, dtype=float),
        "_tactm": tactm,
    }
    if nsedflume == 1:
        out["erate_per_core_cm_s"] = np.stack(erate_per_core, axis=0)
        out["ea_per_core"] = None
        out["en_per_core"] = None
        out["max_rate_per_core_cm_s"] = None
    else:
        out["erate_per_core_cm_s"] = None
        out["ea_per_core"] = np.stack(ea_per_core, axis=0)
        out["en_per_core"] = np.stack(en_per_core, axis=0)
        out["max_rate_per_core_cm_s"] = np.stack(maxrate_per_core, axis=0)
    return out


# ---------------------------------------------------------------------------
# core_field.sdf
# ---------------------------------------------------------------------------


def parse_core_field_sdf(path: Path | str) -> np.ndarray:
    """Parse ``core_field.sdf`` (DSI or SNL format auto-detected per
    s_sedic.f90:151–180).

    Returns
    -------
    np.ndarray
        Integer matrix shape ``(jc, ic)`` of core IDs.
    """
    p = Path(path)

    # --- Format detection -----------------------------------------------
    # Mirror s_sedic.f90:151–157: read the first line (up to ~120 chars)
    # and search the first 117 characters for "DSI" or "dsi". If found,
    # use the DSI standard layout; else SNL standard.
    with open(p, "r") as fh:
        first_line = fh.readline()
    is_dsi = "DSI" in first_line[:117] or "dsi" in first_line[:117]

    if is_dsi:
        return _parse_core_field_dsi(p)
    return _parse_core_field_snl(p)


def _parse_core_field_snl(path: Path) -> np.ndarray:
    """SNL standard: first non-comment line is INCORE, then row-major
    matrix of integers (jc rows × ic cols).

    Per ``s_sedic.f90:177-179`` the matrix is read top-down (J = JC..1).
    The example in SAND2008-5621 Figure 2 shows the matrix in that
    same J-descending order. We preserve that order verbatim — i.e.
    row 0 of the returned array is the first matrix row in the file.
    """
    lines = list(_iter_data_lines(path))
    incore_tokens = _take_ints(lines[0])
    incore = incore_tokens[0]  # noqa: F841 — preserved for caller validation
    # All remaining lines are matrix rows.
    matrix_rows = []
    for line in lines[1:]:
        # Tokens may be space-separated OR EFDC's "120(I1,1X)" packed
        # 1-digit + space format. Splitting on whitespace handles both
        # since each digit is followed by a blank.
        matrix_rows.append([int(t) for t in line.split()])
    if not matrix_rows:
        return np.zeros((0, 0), dtype=int)
    width = max(len(r) for r in matrix_rows)
    out = np.zeros((len(matrix_rows), width), dtype=int)
    for j, row in enumerate(matrix_rows):
        out[j, : len(row)] = row
    return out


def _parse_core_field_dsi(path: Path) -> np.ndarray:
    """DSI standard: 3 header comment lines, then INCORE line, then 2
    more header comment lines, then per-cell ``I J CORE`` triples.

    Output shape ``(jmax, imax)`` is sized to fit the maximum I and J
    seen in the file (1-indexed in the file, 0-indexed in the array).
    Cells absent from the file remain 0.
    """
    # The "header lines" the Fortran skips are read with read(20,'(A80)')
    # and discarded. They include the comment-line stars. We use
    # _iter_data_lines (which already skips comments) to find the
    # INCORE line, then read all subsequent triples as data.
    lines = list(_iter_data_lines(path))
    # First non-comment data line is INCORE.
    incore_tokens = _take_ints(lines[0])
    incore = incore_tokens[0]  # noqa: F841

    triples: list[tuple[int, int, int]] = []
    for line in lines[1:]:
        tokens = line.split()
        if len(tokens) < 3:
            continue
        try:
            i = int(tokens[0])
            j = int(tokens[1])
            core = int(tokens[2])
        except ValueError:
            continue
        triples.append((i, j, core))

    if not triples:
        return np.zeros((0, 0), dtype=int)
    imax = max(t[0] for t in triples)
    jmax = max(t[1] for t in triples)
    out = np.zeros((jmax, imax), dtype=int)
    for i, j, core in triples:
        # File is 1-indexed; numpy 0-indexed.
        out[j - 1, i - 1] = core
    return out


# ---------------------------------------------------------------------------
# Bundle assembly
# ---------------------------------------------------------------------------


def load_sedflume_bundle(
    bed_sdf: Path | str,
    erate_sdf: Path | str,
    core_field_sdf: Path | str | None = None,
) -> SedflumeBundle:
    """Convenience: parse all three SEDflume files into a single bundle."""
    bed = parse_bed_sdf(bed_sdf)
    nsedflume = bed["nsedflume"]
    n_class = bed["_n_class"]

    erate = parse_erate_sdf(
        erate_sdf,
        n_layers=bed["n_layers"],
        n_class=n_class,
        nsedflume=nsedflume,
    )

    if core_field_sdf is not None:
        core_field_ij = parse_core_field_sdf(core_field_sdf)
    else:
        core_field_ij = None

    return SedflumeBundle(
        # bed.sdf
        n_layers=bed["n_layers"],
        var_bed=bed["var_bed"],
        icalc_bl=bed["icalc_bl"],
        nsedflume=nsedflume,
        zb_skin_um=bed["zb_skin_um"],
        tau_const_pa=bed["tau_const_pa"],
        bedload_cutoff_um=bed["bedload_cutoff_um"],
        max_deposit_limit=bed["max_deposit_limit"],
        d50_um=bed["d50_um"],
        tau_ce_pa=bed["tau_ce_pa"],
        tau_cs_pa=bed["tau_cs_pa"],
        settling_cm_s=bed["settling_cm_s"],
        size_interpolants_um=bed["size_interpolants_um"],
        taucrit_per_size_pa=bed["taucrit_per_size_pa"],
        erate_active_table=bed.get("erate_active_table"),
        actdep_a=bed.get("actdep_a"),
        actdep_n=bed.get("actdep_n"),
        actdep_max=bed.get("actdep_max"),
        # erate.sdf
        n_cores=erate["n_cores"],
        layer_thickness_cm=erate["layer_thickness_cm"],
        layer_taucrit_pa=erate["layer_taucrit_pa"],
        bulk_density_g_cm3=erate["bulk_density_g_cm3"],
        water_density_g_cm3=erate["water_density_g_cm3"],
        solid_density_g_cm3=erate["solid_density_g_cm3"],
        particle_size_distribution_pct=erate["particle_size_distribution_pct"],
        tau_levels_pa=erate["tau_levels_pa"],
        erate_per_core_cm_s=erate["erate_per_core_cm_s"],
        ea_per_core=erate["ea_per_core"],
        en_per_core=erate["en_per_core"],
        max_rate_per_core_cm_s=erate["max_rate_per_core_cm_s"],
        # core_field.sdf
        core_field_ij=core_field_ij,
    )
