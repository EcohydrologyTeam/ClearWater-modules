"""Visually compelling SSM demonstration on a synthetic 1D shear gradient.

Unlike ``smoke_plot.py`` (which mirrors the smoke test: spatially uniform
forcing on a 49-cell mesh), this script imposes a *spatial gradient* of
bed shear stress across N synthetic cells while holding the initial bed
identical at every cell. The differential response across cells then
produces three figures suitable for the SSM ERDC Technical Note:

  gradient_taub_spatial.png   — bar plot of imposed τ across cells
  gradient_bed_mass.png       — fanned-out per-cell bed-mass time series
  gradient_armoring.png       — surface D50 vs time at selected cells

Run (no editable install required):

    python tests/sediment/figures/demonstration_plot.py

Output PNGs are written to:

    Publication-ClearWater-SSM/technical_note/latex/src/results_figures/

(the script targets the TN figures directory directly so the LaTeX
include statements pick up the new artefacts on the next build).

Design notes
------------
* Mesh: synthetic ``xarray.Dataset`` with ``(time, nface) = (N_STEPS+2,
  N_CELLS)``. No hydraulics; no edge connectivity. SSM is configured
  with ``shear_driver="external"`` and ``bedload_solver="off"`` so the
  synthetic mesh is acceptable.
* Initial bed: same multi-layer SAND2008 SEDflume bundle replicated at
  every cell (single core ID = 1). 8 size classes, ~5 in-place layers.
* Forcing: τ ramps linearly from ``TAU_LOW`` (0.3 Pa) at cell 0 to
  ``TAU_HIGH`` (5 Pa) at cell N-1 — straddles τ_ce of the SAND2008
  classes so we get qualitatively distinct erosion / armoring behaviour
  across the gradient.
* Step size: ``dt = 600 s`` × 150 steps ≈ 25 hours simulated.
"""

from __future__ import annotations

import sys
from datetime import timedelta
from pathlib import Path

import numpy as np
import xarray as xr

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm  # noqa: F401  -- kept for ScalarMappable
from matplotlib.colors import Normalize

try:
    import seaborn as sns

    sns.set_style("whitegrid")
    _HAS_SEABORN = True
except ImportError:  # pragma: no cover
    _HAS_SEABORN = False


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

_THIS_FILE = Path(__file__).resolve()
# tests/sediment/figures/demonstration_plot.py → repo root is parents[3].
_MODULES_REPO_ROOT = _THIS_FILE.parents[3]
_MODULES_SRC = _MODULES_REPO_ROOT / "src"
if str(_MODULES_SRC) not in sys.path:
    sys.path.insert(0, str(_MODULES_SRC))

_SAND2008_DATA_DIR = (
    _MODULES_REPO_ROOT / "tests" / "sediment" / "data" / "sand2008_example"
)
_BED_SDF = _SAND2008_DATA_DIR / "bed.sdf"
_ERATE_SDF = _SAND2008_DATA_DIR / "erate.sdf"
_CORE_FIELD_SDF = _SAND2008_DATA_DIR / "core_field.sdf"

# TN figures directory (sibling repo).
_TN_FIG_DIR = (
    _MODULES_REPO_ROOT.parent
    / "Publication-ClearWater-SSM"
    / "technical_note"
    / "latex"
    / "src"
    / "results_figures"
)


# ---------------------------------------------------------------------------
# Scenario parameters
# ---------------------------------------------------------------------------

N_CELLS: int = 20
N_STEPS: int = 192
DT_SECONDS: int = 3600                 # 1-hour sediment step (long horizon)
TAU_LOW_PA: float = 0.3                # below τ_ce of medium / coarse classes
TAU_HIGH_PA: float = 4.0               # above τ_ce of all classes
DPI: int = 150
FIGSIZE = (11.0, 7.0)                  # ≈ 1650×1050 at dpi=150 → resampled


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def _check_fixtures() -> None:
    missing: list[Path] = [
        p for p in (_BED_SDF, _ERATE_SDF, _CORE_FIELD_SDF) if not p.is_file()
    ]
    if missing:
        msg = "\n  ".join(str(m) for m in missing)
        print(
            "ERROR: required SAND2008 fixtures are missing:\n  " + msg,
            file=sys.stderr,
        )
        sys.exit(2)


def _build_synthetic_mesh(n_cells: int, n_steps: int) -> xr.Dataset:
    """Synthetic mesh: (time, nface) only; no hydraulics, no edges."""
    from clearwater_modules_v2.processes.sediment import contracts

    n_time = n_steps + 2  # +2 buffer slots so writes never go out of range
    return xr.Dataset(
        coords={
            contracts.DIM_TIME: np.arange(n_time, dtype="int64"),
            contracts.DIM_NFACE: np.arange(n_cells, dtype="int64"),
        }
    )


def _build_ssm(mesh: xr.Dataset, time_step: timedelta):
    from clearwater_modules_v2.processes.sediment import SSM, contracts
    from clearwater_modules_v2.processes.sediment.io.sedflume import (
        load_sedflume_bundle,
    )
    from clearwater_modules_v2.processes.sediment.ssm import (
        _build_classes_from_bundle,
    )

    bundle = load_sedflume_bundle(_BED_SDF, _ERATE_SDF, _CORE_FIELD_SDF)
    registry = _build_classes_from_bundle(bundle)
    n_face = mesh.sizes[contracts.DIM_NFACE]

    ssm = SSM(
        sediment_classes=registry,
        sedflume_bundle=bundle,
        shear_driver="external",
        shear_options={"growth_limit": 0.0},
        bedload_solver="off",
        nsedflume=bundle.nsedflume,
        biostabilization_alpha=0.0,
        time_step=time_step,
        core_id=np.ones(n_face, dtype=np.int64),
    )
    ssm.bind_mesh(mesh)
    return ssm, registry, bundle


def _seed_external_tau_ramp(
    mesh: xr.Dataset, tau_low: float, tau_high: float
) -> np.ndarray:
    """Allocate the external τ field and write a linear cell-wise ramp.

    The ramp is constant in time, so every τ slot in the mesh holds the
    same per-cell vector. Returns the per-cell τ vector for plotting.
    """
    from clearwater_modules_v2.processes.sediment import contracts

    n_time = mesh.sizes[contracts.DIM_TIME]
    n_face = mesh.sizes[contracts.DIM_NFACE]
    tau_per_cell = np.linspace(tau_low, tau_high, n_face, dtype="float32")
    arr = np.broadcast_to(
        tau_per_cell[None, :], (n_time, n_face)
    ).astype("float32").copy()
    mesh[contracts.VAR_BED_SHEAR_STRESS_INPUT] = (
        (contracts.DIM_TIME, contracts.DIM_NFACE),
        arr,
    )
    return tau_per_cell


# ---------------------------------------------------------------------------
# Per-step diagnostic capture (mirrors smoke_plot.py shape)
# ---------------------------------------------------------------------------


def _surface_d50_per_cell(
    class_fraction: np.ndarray,        # (nface, n_layer, n_class)
    layer_active: np.ndarray,          # (nface, n_layer)
    d50_um: np.ndarray,                # (n_class,)
) -> np.ndarray:
    """Mass-weighted D50 of the topmost non-absent layer at each cell."""
    n_face = class_fraction.shape[0]
    out = np.zeros(n_face, dtype="float64")
    for f in range(n_face):
        non_absent = np.flatnonzero(layer_active[f] != 0)
        k = int(non_absent[0]) if non_absent.size else 0
        frac = class_fraction[f, k]
        s = frac.sum()
        out[f] = float((frac * d50_um).sum() / s) if s > 0.0 else 0.0
    return out


def _topmost_in_place_d50_per_cell(
    class_fraction: np.ndarray,        # (nface, n_layer, n_class)
    layer_active: np.ndarray,          # (nface, n_layer) int8 {0=absent, 1=active, 2=in_place}
    d50_um: np.ndarray,                # (n_class,)
) -> np.ndarray:
    """Mass-weighted D50 of the *topmost in-place* layer per cell.

    The active layer (layer index 0/1) is constantly refreshed from the
    underlying in-place layer by the reorganize step; tracking the top
    in-place layer directly is the cleanest way to see armoring emerge
    over the simulated window. Mirrors the SAND2008 reference test
    technique (see ``test_sand2008_reference.py``::``test_in_place_layer
    _d50_increases_under_selectively_gated_shear``).
    """
    n_face, n_layer, _ = class_fraction.shape
    LAYER_IN_PLACE = 2  # bed_mod.LAYER_IN_PLACE
    out = np.zeros(n_face, dtype="float64")
    for f in range(n_face):
        ip = np.flatnonzero(layer_active[f] == LAYER_IN_PLACE)
        if ip.size == 0:
            # Fall back to topmost non-absent.
            non_absent = np.flatnonzero(layer_active[f] != 0)
            k = int(non_absent[0]) if non_absent.size else 0
        else:
            k = int(ip[0])
        frac = class_fraction[f, k]
        s = frac.sum()
        out[f] = float((frac * d50_um).sum() / s) if s > 0.0 else 0.0
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    _check_fixtures()
    _TN_FIG_DIR.mkdir(parents=True, exist_ok=True)

    from clearwater_modules_v2.processes.sediment import contracts

    mesh = _build_synthetic_mesh(N_CELLS, N_STEPS)
    ssm, registry, bundle = _build_ssm(
        mesh, time_step=timedelta(seconds=DT_SECONDS)
    )
    tau_per_cell = _seed_external_tau_ramp(mesh, TAU_LOW_PA, TAU_HIGH_PA)

    n_class = len(registry)
    d50_um_array = registry.d50_um_array
    suspended_names = [cls.suspended_var for cls in registry]

    # Capture buffers --------------------------------------------------------
    bed_total_mass = np.zeros((N_STEPS, N_CELLS), dtype="float64")
    d50_surface = np.zeros((N_STEPS, N_CELLS), dtype="float64")
    d50_in_place = np.zeros((N_STEPS, N_CELLS), dtype="float64")
    domain_bed_mass = np.zeros(N_STEPS, dtype="float64")
    domain_source_cum_mass = np.zeros(N_STEPS, dtype="float64")

    # Initial-state snapshot at t=0 (post-init, before any erosion).
    init_layer_mass = np.asarray(
        ssm._bed.layer_mass_at(0).values, dtype="float64"
    )
    init_class_fraction = np.asarray(
        ssm._bed.class_fraction_at(0).values, dtype="float64"
    )
    init_layer_active = np.asarray(
        ssm._bed.layer_active_at(0).values, dtype="int8"
    )
    initial_total_mass_per_cell = init_layer_mass.sum(axis=-1)
    initial_d50_surface = _surface_d50_per_cell(
        init_class_fraction, init_layer_active, d50_um_array
    )
    initial_d50_in_place = _topmost_in_place_d50_per_cell(
        init_class_fraction, init_layer_active, d50_um_array
    )

    # Time-step loop ---------------------------------------------------------
    for t_idx in range(N_STEPS):
        if t_idx > 0:
            ssm._bed.set_layer_mass_at(
                t_idx, ssm._bed.layer_mass_at(t_idx - 1).values
            )
            ssm._bed.set_class_fraction_at(
                t_idx, ssm._bed.class_fraction_at(t_idx - 1).values
            )
            ssm._bed.set_layer_active_at(
                t_idx, ssm._bed.layer_active_at(t_idx - 1).values
            )
            ssm._bed.set_layer_taucrit_at(
                t_idx, ssm._bed.layer_taucrit_at(t_idx - 1).values
            )

        ssm.run(time=t_idx)

        layer_mass = np.asarray(
            ssm._bed.layer_mass_at(t_idx).values, dtype="float64"
        )
        bed_total_mass[t_idx] = layer_mass.sum(axis=-1)

        class_fraction = np.asarray(
            ssm._bed.class_fraction_at(t_idx).values, dtype="float64"
        )
        layer_active = np.asarray(
            ssm._bed.layer_active_at(t_idx).values, dtype="int8"
        )
        d50_surface[t_idx] = _surface_d50_per_cell(
            class_fraction, layer_active, d50_um_array
        )
        d50_in_place[t_idx] = _topmost_in_place_d50_per_cell(
            class_fraction, layer_active, d50_um_array
        )

        domain_bed_mass[t_idx] = float(layer_mass.sum())
        step_source_g_cm2 = 0.0
        for name in suspended_names:
            src_name = f"{name}_source"
            if src_name in mesh.data_vars:
                step_source_g_cm2 += float(np.nansum(mesh[src_name].values))
        domain_source_cum_mass[t_idx] = (
            domain_source_cum_mass[t_idx - 1] if t_idx > 0 else 0.0
        ) + step_source_g_cm2

    # ---- Quick sanity checks (printed; not asserted) ----------------------
    delta_bed = domain_bed_mass[-1] - float(initial_total_mass_per_cell.sum())
    sum_check = delta_bed + domain_source_cum_mass[-1]
    print("Demonstration scenario summary")
    print("------------------------------")
    print(f"  cells              : {N_CELLS}")
    print(f"  steps              : {N_STEPS}  (dt = {DT_SECONDS} s)")
    print(f"  simulated time     : "
          f"{N_STEPS * DT_SECONDS / 3600.0:.2f} h")
    print(f"  τ range            : [{TAU_LOW_PA:.2f}, "
          f"{TAU_HIGH_PA:.2f}] Pa across cells (linear)")
    print(f"  initial bed mass   : {initial_total_mass_per_cell.sum():.4f} g/cm² (Σcells)")
    print(f"  final bed mass     : {domain_bed_mass[-1]:.4f} g/cm²")
    print(f"  Δ bed              : {delta_bed:+.4e} g/cm²")
    print(f"  cum source ΣΔt(E−D): {domain_source_cum_mass[-1]:+.4e} g/cm²")
    print(f"  Σ check (≈0)       : {sum_check:+.4e} g/cm²")

    # Per-cell qualitative metrics.
    bed_loss_per_cell = initial_total_mass_per_cell - bed_total_mass[-1]
    d50_change = d50_in_place[-1] - initial_d50_in_place
    print("\nPer-cell (selected): "
          "cell  τ(Pa)   Δbed(g/cm²)   D50_inplace_init→final (μm)")
    for c in (0, N_CELLS // 4, N_CELLS // 2, 3 * N_CELLS // 4, N_CELLS - 1):
        print(
            f"  cell {c:2d}  τ={tau_per_cell[c]:5.2f}  "
            f"Δbed={-bed_loss_per_cell[c]:+8.4e}  "
            f"D50_ip: {initial_d50_in_place[c]:7.2f} → "
            f"{d50_in_place[-1, c]:7.2f}  (Δ {d50_change[c]:+6.2f})"
        )

    # ---- Plotting ---------------------------------------------------------
    cmap = matplotlib.colormaps["viridis"]
    norm = Normalize(vmin=tau_per_cell.min(), vmax=tau_per_cell.max())
    cell_colors = [cmap(norm(t)) for t in tau_per_cell]

    time_hours = np.arange(N_STEPS) * DT_SECONDS / 3600.0

    out_paths: list[Path] = []

    # ----- 1) Spatial τ ramp ---------------------------------------------
    fig, ax = plt.subplots(figsize=FIGSIZE)
    bars = ax.bar(
        np.arange(N_CELLS), tau_per_cell,
        color=cell_colors, edgecolor="0.3", linewidth=0.5,
    )
    # τ_ce reference lines (per-class) so the reader sees which classes
    # are above/below threshold at each cell.
    tau_ce_classes = bundle.tau_ce_pa
    for s, tc in enumerate(tau_ce_classes):
        ax.axhline(
            tc, color="0.55", linestyle="--", linewidth=0.7, alpha=0.65,
        )
        ax.text(
            N_CELLS - 0.4, tc, f" τ_ce class {s} ({bundle.d50_um[s]:.0f} μm)",
            fontsize=7, va="center", ha="left", color="0.35",
        )
    ax.set_xlabel("Cell index")
    ax.set_ylabel(r"Imposed bed shear stress  $\tau_b$  (Pa)")
    ax.set_xlim(-0.6, N_CELLS - 0.4 + 4.0)  # right margin for τ_ce labels
    ax.set_title(
        f"Synthetic 1D channel: linear $\\tau_b$ ramp across {N_CELLS} cells "
        f"({TAU_LOW_PA:.1f} → {TAU_HIGH_PA:.1f} Pa). Dashed lines: per-class $\\tau_{{ce}}$."
    )
    fig.tight_layout()
    out_taub = _TN_FIG_DIR / "gradient_taub_spatial.png"
    fig.savefig(out_taub, dpi=DPI)
    plt.close(fig)
    out_paths.append(out_taub)

    # ----- 2) Fanned-out bed-mass time series ----------------------------
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for c in range(N_CELLS):
        ax.plot(
            time_hours, bed_total_mass[:, c],
            color=cell_colors[c], linewidth=1.6, alpha=0.95,
        )
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.04)
    cbar.set_label(r"Imposed bed shear stress  $\tau_b$  (Pa)")
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Total bed mass per cell  (g cm$^{-2}$, summed over layers)")
    ax.set_title(
        f"SSM demonstration: differential bed-mass response under a "
        f"$\\tau_b$ gradient ({N_CELLS} cells, "
        f"{N_STEPS * DT_SECONDS / 3600.0:.0f} h, SAND2008 multi-class bed)"
    )
    fig.tight_layout()
    out_bed = _TN_FIG_DIR / "gradient_bed_mass.png"
    fig.savefig(out_bed, dpi=DPI)
    plt.close(fig)
    out_paths.append(out_bed)

    # ----- 3) Surface D50 vs time at selected cells (armoring trace) -----
    # Selection rule: the visually interesting cells are the ones where
    # τ falls between the τ_ce of two adjacent size classes — partial-
    # gating regime, where some classes erode and some are protected.
    # We rank cells by |Δ D50| at the final step and keep the top six,
    # then add the two endpoints (cell 0 and cell N-1) for context so
    # the reader can see that very-low-τ and very-high-τ cells stay flat.
    d50_change_now = np.abs(d50_in_place[-1] - initial_d50_in_place)
    armoring_rank = np.argsort(-d50_change_now)
    top_armoring = [int(c) for c in armoring_rank[:6]]
    # Always include the endpoints for context.
    selected_cells = sorted(set(top_armoring) | {0, N_CELLS - 1})

    fig, ax = plt.subplots(figsize=FIGSIZE)
    for c in selected_cells:
        ax.plot(
            time_hours, d50_in_place[:, c],
            color=cell_colors[c], linewidth=1.8, marker="o",
            markersize=3.5, markevery=max(1, N_STEPS // 20),
            label=f"cell {c}  (τ = {tau_per_cell[c]:.2f} Pa)",
        )
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel(r"In-place layer mass-weighted $D_{50}$  ($\mu$m)")
    ax.set_title(
        "Armoring signature: $D_{50}$ of the topmost in-place layer rises as "
        "fines preferentially erode under partial-shear gating"
    )
    ax.legend(loc="best", fontsize=8, ncol=2, frameon=True)
    fig.tight_layout()
    out_d50 = _TN_FIG_DIR / "gradient_armoring.png"
    fig.savefig(out_d50, dpi=DPI)
    plt.close(fig)
    out_paths.append(out_d50)

    # ---- Report -----------------------------------------------------------
    print("\nWrote PNGs to TN figures directory:")
    for p in out_paths:
        print(f"  {p}  ({p.stat().st_size / 1024.0:.1f} KiB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
