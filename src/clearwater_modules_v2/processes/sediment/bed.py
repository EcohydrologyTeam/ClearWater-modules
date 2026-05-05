"""Multi-layer bed state and active-layer reorganization.

The bed at each cell is a fixed-depth stack of K_B layers indexed top-down:

* Layer 1: active (sorting) layer
* Layer 2: deposition layer
* Layers 3..K_B: in-place layers (SEDflume core data)

Three layer states (LAYERACTIVE):

* 0 = absent
* 1 = active / deposited
* 2 = in-place core

Reference: SAND2008-5621; design spec §7. The borrow / promote / collapse
algorithm follows Lick (2008) and SAND2008-5621 §"S_SEDZLJ.f90"; the EFDC
Fortran source was used only as a behavioural reference, not copied.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np
import xarray as xr

from . import contracts
from .classes import SedimentClassRegistry


# Layer-state constants matching EFDC LAYERACTIVE convention.
LAYER_ABSENT = np.int8(0)
LAYER_ACTIVE = np.int8(1)
LAYER_IN_PLACE = np.int8(2)


# Mass-conservation tolerance (g/cm^2). float32 epsilon is ~1.2e-7; we
# allow a few orders of magnitude headroom for accumulated round-off
# across cells.
_MASS_CONSERVATION_TOL = 1e-5


# Type alias for the time-axis selector accepted by the BedState
# accessors: integer positional index or any label that ``.sel(time=...)``
# can match (e.g. datetime, np.datetime64, str).
TimeKey = Union[int, np.integer, np.datetime64, "object"]


def _select_time(da: xr.DataArray, t: TimeKey) -> xr.DataArray:
    """Return ``da`` selected at time ``t`` (positional if int, label otherwise)."""
    if isinstance(t, (int, np.integer)):
        return da.isel({contracts.DIM_TIME: int(t)})
    return da.sel({contracts.DIM_TIME: t})


def _assign_time(mesh: xr.Dataset, name: str, t: TimeKey, value) -> None:
    """Write ``value`` into ``mesh[name]`` at time ``t``, broadcasting/aligning.

    Uses positional ``isel`` if ``t`` is integer, otherwise label-based
    ``sel`` via ``.loc``. The underlying numpy buffer is mutated in place
    so the change is visible through any view bound to the same Dataset.
    """
    da = mesh[name]
    # Materialise the value as a plain numpy array of matching dtype.
    if isinstance(value, xr.DataArray):
        # Make sure dim ordering matches the slice we're writing into.
        slice_dims = tuple(d for d in da.dims if d != contracts.DIM_TIME)
        value = value.transpose(*slice_dims)
        arr = np.asarray(value.values, dtype=da.dtype)
    else:
        arr = np.asarray(value, dtype=da.dtype)

    if isinstance(t, (int, np.integer)):
        # Positional write into the underlying buffer.
        time_axis = da.dims.index(contracts.DIM_TIME)
        idx = [slice(None)] * da.ndim
        idx[time_axis] = int(t)
        da.values[tuple(idx)] = arr
    else:
        da.loc[{contracts.DIM_TIME: t}] = arr


@dataclass
class BedState:
    """Bound view of the bed-state DataArrays on the mesh dataset.

    Constructed once per SSM run via :func:`initialize_bed_state`; provides
    typed accessors that read from / write to the underlying xarray
    Dataset in place. The Dataset itself is the single source of truth —
    this class holds no shadow state.

    Time-axis ``t`` may be an integer positional index or a label (e.g.
    ``np.datetime64``); the accessors dispatch accordingly.
    """

    mesh: xr.Dataset
    n_layers: int
    n_classes: int
    registry: SedimentClassRegistry

    # ------------------------------------------------------------------
    # Per-layer dynamic state
    # ------------------------------------------------------------------

    def layer_mass_at(self, t: TimeKey) -> xr.DataArray:
        """TSED at time ``t`` — dims (nface, ssm_layer), units g/cm²."""
        return _select_time(self.mesh[contracts.VAR_BED_LAYER_MASS], t)

    def set_layer_mass_at(self, t: TimeKey, value) -> None:
        _assign_time(self.mesh, contracts.VAR_BED_LAYER_MASS, t, value)

    def class_fraction_at(self, t: TimeKey) -> xr.DataArray:
        """PERSED at time ``t`` — dims (nface, ssm_layer, ssm_class)."""
        return _select_time(self.mesh[contracts.VAR_BED_CLASS_FRACTION], t)

    def set_class_fraction_at(self, t: TimeKey, value) -> None:
        _assign_time(self.mesh, contracts.VAR_BED_CLASS_FRACTION, t, value)

    def layer_active_at(self, t: TimeKey) -> xr.DataArray:
        """LAYERACTIVE at time ``t`` — dims (nface, ssm_layer), int8."""
        return _select_time(self.mesh[contracts.VAR_BED_LAYER_ACTIVE], t)

    def set_layer_active_at(self, t: TimeKey, value) -> None:
        _assign_time(self.mesh, contracts.VAR_BED_LAYER_ACTIVE, t, value)

    def layer_taucrit_at(self, t: TimeKey) -> xr.DataArray:
        """TAUCOR at time ``t`` — dims (nface, ssm_layer), Pa."""
        return _select_time(self.mesh[contracts.VAR_BED_LAYER_TAUCRIT], t)

    def set_layer_taucrit_at(self, t: TimeKey, value) -> None:
        _assign_time(self.mesh, contracts.VAR_BED_LAYER_TAUCRIT, t, value)

    def layer_thickness_at(self, t: TimeKey) -> xr.DataArray:
        """HBED at time ``t`` — dims (nface, ssm_layer), m."""
        return _select_time(self.mesh[contracts.VAR_BED_LAYER_THICKNESS], t)

    def set_layer_thickness_at(self, t: TimeKey, value) -> None:
        _assign_time(self.mesh, contracts.VAR_BED_LAYER_THICKNESS, t, value)

    def layer_age_at(self, t: TimeKey) -> xr.DataArray:
        """Per-layer mean age (s) at time ``t`` — dims (nface, ssm_layer).

        Used by the consolidation model to compute the effective τ_ce.
        Older layers have larger ages and (for cohesive sediment) higher
        effective critical shear stress.
        """
        return _select_time(self.mesh[contracts.VAR_BED_LAYER_AGE], t)

    def set_layer_age_at(self, t: TimeKey, value) -> None:
        _assign_time(self.mesh, contracts.VAR_BED_LAYER_AGE, t, value)

    def bedload_mass_at(self, t: TimeKey) -> xr.DataArray:
        """CBL at time ``t`` — dims (nface, ssm_class), g/cm²."""
        return _select_time(self.mesh[contracts.VAR_BEDLOAD_MASS], t)

    def set_bedload_mass_at(self, t: TimeKey, value) -> None:
        _assign_time(self.mesh, contracts.VAR_BEDLOAD_MASS, t, value)

    # ------------------------------------------------------------------
    # Static / per-cell properties
    # ------------------------------------------------------------------

    @property
    def layer_initial_mass(self) -> xr.DataArray:
        """TSED0 — dims (nface, ssm_layer), g/cm². Time-invariant."""
        return self.mesh[contracts.VAR_BED_LAYER_INITIAL_MASS]

    @property
    def layer_bulk_density(self) -> xr.DataArray:
        """BULKDENS — dims (nface, ssm_layer), g/cm³. Time-invariant in SEDZLJ."""
        return self.mesh[contracts.VAR_BED_LAYER_BULK_DENSITY]


# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------


def _zarr_chunk_for(spec: contracts.VarSpec, sizes: dict[str, int]) -> dict[str, int]:
    """Return Zarr-friendly chunk dictionary for a VarSpec, capped by actual sizes.

    See ``contracts.DEFAULT_*_CHUNK`` for the targets.
    """
    chunks: dict[str, int] = {}
    for d in spec.dims:
        if d == contracts.DIM_TIME:
            target = contracts.DEFAULT_TIME_CHUNK
        elif d == contracts.DIM_NFACE:
            target = contracts.DEFAULT_FACE_CHUNK
        elif d == contracts.DIM_LAYER:
            target = contracts.DEFAULT_LAYER_CHUNK
        elif d == contracts.DIM_CLASS:
            target = contracts.DEFAULT_CLASS_CHUNK
        else:
            target = sizes.get(d, 1)
        chunks[d] = max(1, min(target, sizes.get(d, target)))
    return chunks


def initialize_bed_state(
    mesh: xr.Dataset,
    registry: SedimentClassRegistry,
    n_layers: int,
    initial_layer_mass: np.ndarray,        # (nface, n_layers)  g/cm^2
    initial_class_fraction: np.ndarray,    # (nface, n_layers, n_class)
    bulk_density: np.ndarray,              # (nface, n_layers) g/cm^3
    initial_layer_active: np.ndarray,      # (nface, n_layers) int8
    taucor_initial: np.ndarray,            # (nface, n_layers) Pa
    streaming_chunks: dict | None = None,
) -> BedState:
    """Allocate all bed-state DataArrays on the mesh dataset and populate ICs.

    Adds variables listed in :data:`contracts.BED_STATE_SPECS` to ``mesh``
    in place, with correct dims/dtypes and Zarr chunking hints. Returns
    a :class:`BedState` view bound to the mesh.

    The mesh is required to already carry ``time`` and ``nface``
    coordinates. ``ssm_layer`` and ``ssm_class`` coordinates are added
    here if absent (zero-based integer indices).
    """
    n_classes = len(registry)

    if contracts.DIM_TIME not in mesh.dims:
        raise ValueError(
            f"mesh must have a {contracts.DIM_TIME!r} dimension before "
            "initialize_bed_state is called"
        )
    if contracts.DIM_NFACE not in mesh.dims:
        raise ValueError(
            f"mesh must have a {contracts.DIM_NFACE!r} dimension before "
            "initialize_bed_state is called"
        )

    n_time = mesh.sizes[contracts.DIM_TIME]
    n_face = mesh.sizes[contracts.DIM_NFACE]

    # Add layer / class coordinates (idempotent).
    if contracts.DIM_LAYER not in mesh.coords:
        mesh.coords[contracts.DIM_LAYER] = np.arange(n_layers, dtype="int32")
    elif mesh.sizes[contracts.DIM_LAYER] != n_layers:
        raise ValueError(
            f"mesh already has {contracts.DIM_LAYER}="
            f"{mesh.sizes[contracts.DIM_LAYER]} but n_layers={n_layers} requested"
        )

    if contracts.DIM_CLASS not in mesh.coords:
        mesh.coords[contracts.DIM_CLASS] = np.arange(n_classes, dtype="int32")
    elif mesh.sizes[contracts.DIM_CLASS] != n_classes:
        raise ValueError(
            f"mesh already has {contracts.DIM_CLASS}="
            f"{mesh.sizes[contracts.DIM_CLASS]} but registry has {n_classes} classes"
        )

    # Validate initial-condition shapes.
    expected_layer_mass = (n_face, n_layers)
    expected_class_fraction = (n_face, n_layers, n_classes)
    if initial_layer_mass.shape != expected_layer_mass:
        raise ValueError(
            f"initial_layer_mass shape {initial_layer_mass.shape} != {expected_layer_mass}"
        )
    if initial_class_fraction.shape != expected_class_fraction:
        raise ValueError(
            f"initial_class_fraction shape {initial_class_fraction.shape} "
            f"!= {expected_class_fraction}"
        )
    if bulk_density.shape != expected_layer_mass:
        raise ValueError(
            f"bulk_density shape {bulk_density.shape} != {expected_layer_mass}"
        )
    if initial_layer_active.shape != expected_layer_mass:
        raise ValueError(
            f"initial_layer_active shape {initial_layer_active.shape} "
            f"!= {expected_layer_mass}"
        )
    if taucor_initial.shape != expected_layer_mass:
        raise ValueError(
            f"taucor_initial shape {taucor_initial.shape} != {expected_layer_mass}"
        )

    # Map dim names to sizes for shape-building below.
    dim_sizes = {
        contracts.DIM_TIME: n_time,
        contracts.DIM_NFACE: n_face,
        contracts.DIM_LAYER: n_layers,
        contracts.DIM_CLASS: n_classes,
    }

    # Per-spec allocate-zero, then assign initial conditions where applicable.
    for spec in contracts.BED_STATE_SPECS:
        shape = tuple(dim_sizes[d] for d in spec.dims)
        data = np.zeros(shape, dtype=spec.dtype)
        da = xr.DataArray(
            data,
            dims=spec.dims,
            name=spec.name,
            attrs={"units": spec.units, "description": spec.description, "role": spec.role},
        )

        # Apply optional Zarr chunking hints.
        chunk = _zarr_chunk_for(spec, dim_sizes)
        # Only chunk if dask is present in the mesh; otherwise leave as numpy.
        # (Tests run with plain numpy; chunking attrs are still recorded.)
        if streaming_chunks is not None:
            chunk.update({k: v for k, v in streaming_chunks.items() if k in da.dims})
        da.attrs["zarr_chunks"] = chunk

        mesh[spec.name] = da

    # ------------------------------------------------------------------
    # Populate initial conditions at t=0.
    # ------------------------------------------------------------------
    mesh[contracts.VAR_BED_LAYER_INITIAL_MASS].values[:] = initial_layer_mass.astype(
        mesh[contracts.VAR_BED_LAYER_INITIAL_MASS].dtype
    )
    mesh[contracts.VAR_BED_LAYER_BULK_DENSITY].values[:] = bulk_density.astype(
        mesh[contracts.VAR_BED_LAYER_BULK_DENSITY].dtype
    )

    bed = BedState(
        mesh=mesh, n_layers=n_layers, n_classes=n_classes, registry=registry
    )

    bed.set_layer_mass_at(0, initial_layer_mass)
    bed.set_class_fraction_at(0, initial_class_fraction)
    bed.set_layer_active_at(0, initial_layer_active)
    bed.set_layer_taucrit_at(0, taucor_initial)

    # Initial bed thickness consistent with TSED / BULKDENS.
    safe_dens = np.where(bulk_density > 0, bulk_density, 1.0).astype("float64")
    init_thickness = (0.01 * initial_layer_mass / safe_dens).astype("float32")
    bed.set_layer_thickness_at(0, init_thickness)
    init_total_thickness = init_thickness.sum(axis=-1).astype("float32")
    _assign_time(mesh, contracts.VAR_BED_TOTAL_THICKNESS, 0, init_total_thickness)

    return bed


# ---------------------------------------------------------------------------
# Active-layer reorganization
# ---------------------------------------------------------------------------


def _find_slln(layer_mass: np.ndarray, layer_active: np.ndarray) -> np.ndarray:
    """Per-cell index of the first non-empty layer below layer 1 (the "SLLN").

    Returns a 1-D int array of length nface. If no non-empty sub-layer
    exists the value is -1 (caller must handle).

    ``layer_mass`` and ``layer_active`` are (nface, n_layers) numpy arrays.
    Layer indices in the return are zero-based; layer 0 is the active layer
    so we search 1..n_layers-1.
    """
    n_face, n_layers = layer_mass.shape
    slln = np.full(n_face, -1, dtype=np.int32)
    # Scan layers 1..K-1 in order; first non-empty wins.
    for k in range(1, n_layers):
        candidate = (slln == -1) & (layer_active[:, k] != LAYER_ABSENT) & (layer_mass[:, k] > 0.0)
        slln = np.where(candidate, k, slln)
    return slln


def reorganize_active_layer(
    bed: BedState,
    t: TimeKey,
    tau_pa: xr.DataArray,                  # (nface,) Pa
    tau_crit_pa: xr.DataArray,             # (nface,) Pa
    d50_surface_um: xr.DataArray,          # (nface,) μm
    bulk_density_layer1: xr.DataArray,     # (nface,) g/cm^3
    tactm: float = contracts.DEFAULT_TACTM,
) -> None:
    """In-place active-layer reorganization (borrow / promote / collapse).

    Three branches per cell, applied vectorized across all faces:

    1. **Net deposition** (m_1 > T_act): push (m_1 - T_act) excess from
       layer 1 to layer 2; layer-2 PERSED becomes mass-weighted blend.
    2. **Borrow** (m_1 < T_act, τ > τ_crit(SLLN), enough mass available):
       borrow (T_act - m_1) from SLLN to layer 1; layer-1 PERSED becomes
       mass-weighted blend.
    3. **Promote / collapse** (m_1 < T_act, τ > τ_crit(SLLN), SLLN
       insufficient): merge what remains of SLLN into layer 1, mark SLLN
       absent, and identify the next non-empty layer below for the next
       step.

    T_act formula (design spec §5.6, Lick 2008):

    .. math::

        T_{\\rm act} = T_{\\rm actm}\\, D_{50,{\\rm avg}}\\,
                       \\max(1, \\tau/\\tau_{\\rm crit})\\,
                       \\rho_b / 10000   \\quad [{\\rm g/cm^2}]

    Mass conservation across all branches: ``sum(TSED, axis=layer)``
    invariant to within ``_MASS_CONSERVATION_TOL`` (asserted in
    development; can be elevated to a debug-only check for performance).
    """
    # Pull current state as numpy for vectorised in-place math.
    layer_mass = np.asarray(
        bed.layer_mass_at(t).values, dtype="float64"
    ).copy()                                                # (nface, n_layers)
    class_fraction = np.asarray(
        bed.class_fraction_at(t).values, dtype="float64"
    ).copy()                                                # (nface, n_layers, n_class)
    layer_active = np.asarray(
        bed.layer_active_at(t).values, dtype="int8"
    ).copy()                                                # (nface, n_layers)
    layer_taucrit = np.asarray(
        bed.layer_taucrit_at(t).values, dtype="float64"
    )                                                       # (nface, n_layers) Pa
    # Per-layer age (s) — for consolidation. May be all-zero if the
    # consolidation feature isn't in use; the bookkeeping is cheap so
    # we update it unconditionally.
    layer_age = np.asarray(
        bed.layer_age_at(t).values, dtype="float64"
    ).copy()                                                # (nface, n_layers)

    tau = np.asarray(tau_pa.values, dtype="float64")            # (nface,)
    tau_crit = np.asarray(tau_crit_pa.values, dtype="float64")  # (nface,)
    d50 = np.asarray(d50_surface_um.values, dtype="float64")    # (nface,)
    bd1 = np.asarray(bulk_density_layer1.values, dtype="float64")  # (nface,)

    n_face, n_layers = layer_mass.shape
    if n_layers < 2:
        # Single-layer bed has no reorganization to do.
        return

    mass_before = layer_mass.sum(axis=-1).copy()

    # T_act per cell (g/cm^2). Guard against τ_crit = 0.
    safe_taucrit = np.where(tau_crit > 0.0, tau_crit, 1.0)
    ratio = np.where(tau_crit > 0.0, tau / safe_taucrit, 1.0)
    factor = np.maximum(1.0, ratio)
    t_act = tactm * d50 * factor * bd1 / 10000.0  # (nface,)

    m1 = layer_mass[:, 0]
    p1 = class_fraction[:, 0, :]   # (nface, n_class)
    m2 = layer_mass[:, 1]
    p2 = class_fraction[:, 1, :]
    a1 = layer_age[:, 0]           # (nface,) — layer-1 mean age (s)
    a2 = layer_age[:, 1]           # (nface,) — layer-2 mean age (s)

    # ------------------------------------------------------------------
    # Branch (a): net deposition — layer 1 has more mass than T_act needs.
    # Push (m1 - T_act) into layer 2, blend its PERSED mass-weighted.
    # ------------------------------------------------------------------
    branch_a = (m1 > t_act) & (t_act > 0.0)
    if np.any(branch_a):
        excess = m1 - t_act                                    # (nface,)
        new_m2 = m2 + excess                                   # (nface,)
        # Mass-weighted blend: (p2*m2 + p1*excess) / (m2+excess)
        denom = new_m2[:, None]
        # Avoid divide-by-zero: where denom==0 the result doesn't matter
        # because new_m2 will be 0 too (no mass → fraction is moot). Set
        # denom=1 there to keep the formula well-defined.
        safe_denom = np.where(denom > 0.0, denom, 1.0)
        new_p2 = (p2 * m2[:, None] + p1 * excess[:, None]) / safe_denom

        # Age inheritance for layer 2: mass-weighted blend of its
        # existing age and the (younger) excess transferred from
        # layer 1. Layer 1's mean age is unchanged (uniform-aged
        # mass is removed from the top).
        safe_denom_age = np.where(new_m2 > 0.0, new_m2, 1.0)
        new_a2 = (a2 * m2 + a1 * excess) / safe_denom_age

        layer_mass[:, 1] = np.where(branch_a, new_m2, m2)
        layer_mass[:, 0] = np.where(branch_a, t_act, m1)
        class_fraction[:, 1, :] = np.where(branch_a[:, None], new_p2, p2)
        layer_age[:, 1] = np.where(branch_a, new_a2, a2)
        # Layer 2 is "active/deposited" once it has any mass.
        layer_active[:, 1] = np.where(
            branch_a & (new_m2 > 0.0), LAYER_ACTIVE, layer_active[:, 1]
        )

    # Refresh views after branch (a).
    m1 = layer_mass[:, 0]
    p1 = class_fraction[:, 0, :]
    a1 = layer_age[:, 0]

    # ------------------------------------------------------------------
    # Branches (b) and (c): erosion regime. Need SLLN — index of next
    # non-empty layer below layer 1.
    # ------------------------------------------------------------------
    slln = _find_slln(layer_mass, layer_active)             # (nface,) int
    has_slln = slln >= 0

    # Read τ_crit and mass at SLLN per cell. Where slln == -1 use safe
    # fallback values that exclude the cell from branches (b)/(c) anyway.
    safe_slln = np.where(has_slln, slln, 0)                  # (nface,) safe index for gather
    rows = np.arange(n_face)
    m_slln = np.where(has_slln, layer_mass[rows, safe_slln], 0.0)
    tau_slln = np.where(has_slln, layer_taucrit[rows, safe_slln], np.inf)
    p_slln = np.where(
        has_slln[:, None],
        class_fraction[rows, safe_slln, :],
        np.zeros_like(class_fraction[rows, 0, :]),
    )
    a_slln = np.where(has_slln, layer_age[rows, safe_slln], 0.0)

    erosion_regime = (m1 < t_act) & has_slln & (tau > tau_slln) & (t_act > 0.0)

    # Branch (b): sufficient mass below to top up to T_act.
    branch_b = erosion_regime & ((m1 + m_slln) > t_act)
    # Branch (c): not enough mass below — collapse SLLN entirely into layer 1.
    branch_c = erosion_regime & ~branch_b

    if np.any(branch_b):
        deficit = t_act - m1                                # (nface,)
        denom_b = np.where(t_act > 0.0, t_act, 1.0)[:, None]
        new_p1_b = (p1 * m1[:, None] + p_slln * deficit[:, None]) / denom_b
        new_m_slln = m_slln - deficit

        # Age inheritance for layer 1: mass-weighted blend of its
        # existing age and the (older) age of the donor SLLN. SLLN's
        # mean age is preserved (uniform-aged mass is removed from
        # below).
        denom_b_age = np.where(t_act > 0.0, t_act, 1.0)
        new_a1_b = (a1 * m1 + a_slln * deficit) / denom_b_age

        # Update SLLN mass.
        layer_mass[rows, safe_slln] = np.where(
            branch_b, new_m_slln, layer_mass[rows, safe_slln]
        )
        # Update layer 1.
        layer_mass[:, 0] = np.where(branch_b, t_act, layer_mass[:, 0])
        class_fraction[:, 0, :] = np.where(
            branch_b[:, None], new_p1_b, class_fraction[:, 0, :]
        )
        layer_age[:, 0] = np.where(branch_b, new_a1_b, layer_age[:, 0])

    if np.any(branch_c):
        # Mass-weighted blend of (layer 1, SLLN) into layer 1.
        new_m1 = m1 + m_slln
        denom_c = np.where(new_m1 > 0.0, new_m1, 1.0)[:, None]
        new_p1_c = (p1 * m1[:, None] + p_slln * m_slln[:, None]) / denom_c

        # Age inheritance for layer 1: mass-weighted blend of its
        # existing age and SLLN's age (now fully merged into layer 1).
        denom_c_age = np.where(new_m1 > 0.0, new_m1, 1.0)
        new_a1_c = (a1 * m1 + a_slln * m_slln) / denom_c_age

        class_fraction[:, 0, :] = np.where(
            branch_c[:, None], new_p1_c, class_fraction[:, 0, :]
        )
        layer_mass[:, 0] = np.where(branch_c, new_m1, layer_mass[:, 0])
        layer_age[:, 0] = np.where(branch_c, new_a1_c, layer_age[:, 0])
        # Zero out SLLN.
        layer_mass[rows, safe_slln] = np.where(
            branch_c, 0.0, layer_mass[rows, safe_slln]
        )
        # Mark SLLN absent (LAYERACTIVE = 0). Use a layer-mask so we don't
        # accidentally deactivate other layers.
        layer_idx = np.arange(n_layers)[None, :]                # (1, n_layers)
        slln_mask = (layer_idx == safe_slln[:, None]) & branch_c[:, None]
        layer_active = np.where(slln_mask, LAYER_ABSENT, layer_active)
        # Also clear the class fractions of the now-absent SLLN.
        class_fraction[rows, safe_slln, :] = np.where(
            branch_c[:, None],
            np.zeros_like(class_fraction[rows, safe_slln, :]),
            class_fraction[rows, safe_slln, :],
        )
        # Also reset the now-absent SLLN's age.
        layer_age[rows, safe_slln] = np.where(
            branch_c, 0.0, layer_age[rows, safe_slln]
        )

    # Mass-conservation invariant (development-time check).
    mass_after = layer_mass.sum(axis=-1)
    if not np.allclose(mass_after, mass_before, atol=_MASS_CONSERVATION_TOL):
        delta = np.abs(mass_after - mass_before).max()
        raise AssertionError(
            f"Active-layer reorganization violated mass conservation: "
            f"max |Δsum(TSED)|={delta:.3e} > tol={_MASS_CONSERVATION_TOL:.0e}"
        )

    # If layer 1 has any mass it should be marked active.
    layer_active[:, 0] = np.where(
        layer_mass[:, 0] > 0.0, LAYER_ACTIVE, layer_active[:, 0]
    )

    # Write back to mesh.
    bed.set_layer_mass_at(t, layer_mass)
    bed.set_class_fraction_at(t, class_fraction)
    bed.set_layer_active_at(t, layer_active)
    bed.set_layer_age_at(t, layer_age)


# ---------------------------------------------------------------------------
# Age dilution on deposition (used by ssm.run after deposition is added)
# ---------------------------------------------------------------------------


def dilute_layer1_age_on_deposition(
    bed: BedState,
    t: TimeKey,
    layer1_mass_before: np.ndarray,        # (nface,) g/cm^2 — before deposition added
    deposited_mass: np.ndarray,            # (nface,) g/cm^2 — Δm added to layer 1
) -> None:
    """Update layer-1 age in place via mass-weighted dilution by fresh deposit.

    .. math::

        t_{1,\\rm new} = \\frac{t_1 \\, m_1}{m_1 + \\Delta m}

    The new mass enters with age 0; the new layer-mean age is the
    existing age weighted by the existing-mass fraction of the total.
    This is the simplest and most common dilution rule in the
    depth-averaged consolidation literature (Sanford & Maa 2001).

    Where there is no existing mass and no deposit, the layer age is
    left at its current value (typically zero).
    """
    age = np.asarray(bed.layer_age_at(t).values, dtype="float64").copy()
    m1 = np.asarray(layer1_mass_before, dtype="float64")
    dm = np.asarray(deposited_mass, dtype="float64")
    new_m1 = m1 + dm
    # Where new_m1 == 0, no mass present at all → leave age unchanged
    # (still 0). Otherwise dilute by the mass-weighted formula.
    safe = np.where(new_m1 > 0.0, new_m1, 1.0)
    new_a1 = age[:, 0] * m1 / safe
    age[:, 0] = np.where(new_m1 > 0.0, new_a1, age[:, 0])
    bed.set_layer_age_at(t, age)


# ---------------------------------------------------------------------------
# Bed elevation
# ---------------------------------------------------------------------------


def update_bed_elevation(bed: BedState, t: TimeKey, dt_seconds: float = 0.0) -> None:
    """Recompute bed thickness, total thickness, and bed-change diagnostics.

    Per layer:

    .. math:: H_{bed}(K) = 0.01 \\cdot \\frac{m_K}{\\rho_{b,K}} \\quad [\\text{m}]

    Total thickness is the layer-sum; ``bed_change`` is the per-step delta
    relative to the previous time slice (zero at t=0); ``cumulative_bed_change``
    is the running sum from t=0.

    All quantities are written back to the mesh dataset in place.

    Parameters
    ----------
    bed : BedState
    t : TimeKey
        Time slice to update.
    dt_seconds : float, optional
        SSM time step (s). When ``> 0``, the per-layer age field is
        advanced by ``dt`` for every layer that currently holds mass.
        Defaults to ``0.0`` for backward compatibility (no aging).
    """
    layer_mass = np.asarray(bed.layer_mass_at(t).values, dtype="float64")
    bulk_density = np.asarray(bed.layer_bulk_density.values, dtype="float64")

    # Per-layer thickness in metres (CGS mass / density → cm; ×0.01 → m).
    safe_dens = np.where(bulk_density > 0.0, bulk_density, 1.0)
    layer_thickness = 0.01 * layer_mass / safe_dens               # (nface, n_layers)
    layer_thickness = np.where(bulk_density > 0.0, layer_thickness, 0.0)

    total_thickness = layer_thickness.sum(axis=-1)                # (nface,)

    bed.set_layer_thickness_at(t, layer_thickness.astype("float32"))
    _assign_time(
        bed.mesh, contracts.VAR_BED_TOTAL_THICKNESS, t, total_thickness.astype("float32")
    )

    # bed_change: delta vs previous time slice; 0 at t=0.
    if isinstance(t, (int, np.integer)):
        idx = int(t)
    else:
        # Resolve label → integer index along the time axis.
        idx = int(
            bed.mesh.indexes[contracts.DIM_TIME].get_loc(t)
        )

    if idx == 0:
        bed_change = np.zeros_like(total_thickness, dtype="float32")
        cumulative = bed_change.copy()
    else:
        prev_total = np.asarray(
            bed.mesh[contracts.VAR_BED_TOTAL_THICKNESS]
            .isel({contracts.DIM_TIME: idx - 1})
            .values,
            dtype="float64",
        )
        bed_change = (total_thickness - prev_total).astype("float32")

        # Cumulative is previous cumulative + this step's delta.
        prev_cum = np.asarray(
            bed.mesh[contracts.VAR_BED_CUMULATIVE_CHANGE]
            .isel({contracts.DIM_TIME: idx - 1})
            .values,
            dtype="float64",
        )
        cumulative = (prev_cum + bed_change).astype("float32")

    _assign_time(bed.mesh, contracts.VAR_BED_CHANGE, t, bed_change)
    _assign_time(bed.mesh, contracts.VAR_BED_CUMULATIVE_CHANGE, t, cumulative)

    # ------------------------------------------------------------------
    # Advance layer age by dt for every layer that currently holds mass.
    # Empty layers stay at age 0 (no consolidation clock for ghost mass).
    # ------------------------------------------------------------------
    if dt_seconds > 0.0:
        age = np.asarray(bed.layer_age_at(t).values, dtype="float64").copy()
        has_mass = layer_mass > 0.0
        age = np.where(has_mass, age + float(dt_seconds), age)
        # Empty layers: pin to zero (covers any leftover float drift).
        age = np.where(has_mass, age, 0.0)
        bed.set_layer_age_at(t, age)
