"""Generate matplotlib visualizations of SSM behaviour from the
``plan01_10x5`` RAS-coupled smoke-test scenario.

This script replicates the setup in
``tests/sediment/test_ssm_riverine_smoke.py::test_ssm_riverine_smoke_plan01_10x5``
but runs ~10 SSM steps (vs. 5 in the smoke test) and captures per-step
diagnostics that are then plotted as four PNGs in this directory:

  01_bed_mass_timeseries.png
  02_bed_shear_spatial.png
  03_armoring_class_fractions.png
  04_mass_conservation.png

Run from any cwd (uses absolute paths):

    python tests/sediment/figures/smoke_plot.py

The script is purely additive — it does not modify SSM source or the
smoke test itself. If the sibling ClearWater-Riverine-streaming repo or
the SAND2008 SEDflume bundle is missing, it exits with a clear message.
"""

from __future__ import annotations

import sys
from datetime import timedelta
from pathlib import Path

import numpy as np
import xarray as xr

# --- matplotlib / seaborn (force Agg backend for headless environments) ----
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import seaborn as sns

    sns.set_style("whitegrid")
    _HAS_SEABORN = True
except ImportError:  # pragma: no cover
    _HAS_SEABORN = False


# ---------------------------------------------------------------------------
# Path resolution (mirrors test_ssm_riverine_smoke.py)
# ---------------------------------------------------------------------------

_THIS_FILE = Path(__file__).resolve()
# tests/sediment/figures/smoke_plot.py → repo root is parents[3].
_MODULES_REPO_ROOT = _THIS_FILE.parents[3]
_RIVERINE_REPO_ROOT = _MODULES_REPO_ROOT.parent / "ClearWater-Riverine-streaming"
_RIVERINE_FIXTURE_DIR = (
    _RIVERINE_REPO_ROOT / "tests" / "data" / "simple_test_cases"
)
_PLAN01_DIR = _RIVERINE_FIXTURE_DIR / "plan01_10x5"
_PLAN01_HDF = _PLAN01_DIR / "clearWaterTestCases.p01.hdf"
_PLAN01_IC = _PLAN01_DIR / "cwr_initial_conditions_p01.csv"
_PLAN01_BC = _PLAN01_DIR / "cwr_boundary_conditions_p01.csv"

_SAND2008_DATA_DIR = (
    _MODULES_REPO_ROOT / "tests" / "sediment" / "data" / "sand2008_example"
)
_BED_SDF = _SAND2008_DATA_DIR / "bed.sdf"
_ERATE_SDF = _SAND2008_DATA_DIR / "erate.sdf"
_CORE_FIELD_SDF = _SAND2008_DATA_DIR / "core_field.sdf"

_OUTPUT_DIR = _THIS_FILE.parent

# Make sure the modules repo's src/ is on sys.path so this script can be
# launched directly without an editable install.
_MODULES_SRC = _MODULES_REPO_ROOT / "src"
if str(_MODULES_SRC) not in sys.path:
    sys.path.insert(0, str(_MODULES_SRC))


def _import_or_die() -> None:
    """Import ClearWater-Riverine, falling back to its sibling src dir."""
    try:
        import clearwater_riverine  # noqa: F401
        return
    except ImportError:
        riverine_src = _RIVERINE_REPO_ROOT / "src"
        if riverine_src.is_dir() and str(riverine_src) not in sys.path:
            sys.path.insert(0, str(riverine_src))
        try:
            import clearwater_riverine  # noqa: F401
            return
        except ImportError as e:
            print(
                f"ERROR: cannot import clearwater_riverine even after adding "
                f"{riverine_src} to sys.path: {e}",
                file=sys.stderr,
            )
            sys.exit(2)


def _check_fixtures() -> None:
    missing: list[Path] = []
    for p in (_PLAN01_HDF, _PLAN01_IC, _PLAN01_BC,
              _BED_SDF, _ERATE_SDF, _CORE_FIELD_SDF):
        if not p.is_file():
            missing.append(p)
    if missing:
        msg = "\n  ".join(str(m) for m in missing)
        print(
            "ERROR: required fixtures are missing:\n  " + msg,
            file=sys.stderr,
        )
        sys.exit(2)


# ---------------------------------------------------------------------------
# Scenario builders (mirror test_ssm_riverine_smoke.py)
# ---------------------------------------------------------------------------


def _make_riverine(n_steps: int):
    import clearwater_riverine as cwr

    return cwr.ClearwaterRiverine(
        flow_field_file_path=str(_PLAN01_HDF),
        diffusion_coefficient_input=0.001,
        constituent_dict={
            "tracer": {
                "initial_conditions": str(_PLAN01_IC),
                "boundary_conditions": str(_PLAN01_BC),
                "units": "mg/L",
            },
        },
        datetime_range=(0, n_steps + 2),
    )


def _build_ssm(riverine, time_step: timedelta):
    from clearwater_modules_v2.processes.sediment import SSM, contracts
    from clearwater_modules_v2.processes.sediment.io.sedflume import (
        load_sedflume_bundle,
    )
    from clearwater_modules_v2.processes.sediment.ssm import (
        _build_classes_from_bundle,
    )

    bundle = load_sedflume_bundle(_BED_SDF, _ERATE_SDF, _CORE_FIELD_SDF)
    registry = _build_classes_from_bundle(bundle)
    n_face = riverine.mesh.sizes[contracts.DIM_NFACE]

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
    ssm.bind_mesh(riverine.mesh)
    return ssm, registry, bundle


def _seed_external_tau(mesh: xr.Dataset, tau_value_pa: float) -> None:
    from clearwater_modules_v2.processes.sediment import contracts

    n_time = mesh.sizes[contracts.DIM_TIME]
    n_face = mesh.sizes[contracts.DIM_NFACE]
    arr = np.full((n_time, n_face), tau_value_pa, dtype="float32")
    mesh[contracts.VAR_BED_SHEAR_STRESS_INPUT] = (
        (contracts.DIM_TIME, contracts.DIM_NFACE), arr,
    )


# ---------------------------------------------------------------------------
# Main: replicate smoke-test loop, capture diagnostics, write PNGs
# ---------------------------------------------------------------------------


def main() -> int:
    _import_or_die()
    _check_fixtures()

    from clearwater_modules_v2.processes.sediment import contracts

    n_steps = 10  # ~2x the smoke test, longer-trend visualization

    riverine = _make_riverine(n_steps)
    mesh = riverine.mesh
    ssm, registry, bundle = _build_ssm(
        riverine, time_step=timedelta(seconds=60),
    )

    _seed_external_tau(mesh, tau_value_pa=4.0)

    # Identify "real" cells (Riverine pads with ghost cells); we only plot
    # bed/τ for real cells.
    nreal = int(mesh.attrs.get("nreal", mesh.sizes[contracts.DIM_NFACE]))
    n_face_total = mesh.sizes[contracts.DIM_NFACE]
    n_class = len(registry)

    time_labels = mesh[contracts.DIM_TIME].values

    # --- Capture buffers ---------------------------------------------------
    bed_total_mass = np.zeros((n_steps, n_face_total), dtype="float64")
    bed_elevation = np.zeros((n_steps, n_face_total), dtype="float64")
    tau_b = np.zeros((n_steps, n_face_total), dtype="float64")
    d50_surface = np.zeros((n_steps, n_face_total), dtype="float64")
    surface_class_frac = np.zeros(
        (n_steps, n_face_total, n_class), dtype="float64",
    )
    domain_bed_mass = np.zeros(n_steps, dtype="float64")
    domain_suspended_mass = np.zeros(n_steps, dtype="float64")
    domain_bedload_mass = np.zeros(n_steps, dtype="float64")
    domain_source_cum_mass = np.zeros(n_steps, dtype="float64")

    d50_um_array = registry.d50_um_array
    suspended_names = [cls.suspended_var for cls in registry]
    class_labels = [cls.label for cls in registry]

    # ---- Time-step loop (mirrors smoke test, with capture) ---------------
    for t_idx in range(n_steps):
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

        riverine.update()
        ssm.run(time=time_labels[t_idx])

        # --- Capture per-cell bed mass (sum over layers) ------------------
        layer_mass = np.asarray(
            ssm._bed.layer_mass_at(t_idx).values, dtype="float64"
        )                                                  # (nface, n_layers)
        bed_total_mass[t_idx] = layer_mass.sum(axis=-1)

        # --- Bed elevation -------------------------------------------------
        if contracts.VAR_BED_ELEVATION in mesh.data_vars:
            be = mesh[contracts.VAR_BED_ELEVATION].isel(
                {contracts.DIM_TIME: t_idx}
            ).values
        else:
            # Fallback: bulk * thickness; use cumulative bed-change if absent.
            be = np.zeros(n_face_total, dtype="float64")
        bed_elevation[t_idx] = be

        # --- τ_b -----------------------------------------------------------
        if contracts.VAR_BED_SHEAR_STRESS in mesh.data_vars:
            tau_b[t_idx] = mesh[contracts.VAR_BED_SHEAR_STRESS].isel(
                {contracts.DIM_TIME: t_idx}
            ).values

        # --- Surface-layer class fractions and D50 -------------------------
        class_frac = np.asarray(
            ssm._bed.class_fraction_at(t_idx).values, dtype="float64",
        )                                                  # (nface, n_layer, n_class)
        layer_active = np.asarray(
            ssm._bed.layer_active_at(t_idx).values, dtype="int8",
        )                                                  # (nface, n_layer)
        # Find topmost non-absent layer per cell (LAYER_ABSENT = 0).
        # SSM logic: scan top-down for first nonzero entry.
        top_idx = np.zeros(n_face_total, dtype=np.int64)
        for f in range(n_face_total):
            row = layer_active[f]
            non_absent = np.flatnonzero(row != 0)
            top_idx[f] = int(non_absent[0]) if non_absent.size else 0
        # Gather (nface, n_class).
        surf = class_frac[np.arange(n_face_total), top_idx, :]
        surface_class_frac[t_idx] = surf
        # Mass-weighted mean D50: Σ(frac * d50) / Σ(frac). NaN-safe.
        sum_frac = surf.sum(axis=-1)
        with np.errstate(invalid="ignore", divide="ignore"):
            d50 = np.where(
                sum_frac > 0.0,
                (surf * d50_um_array[None, :]).sum(axis=-1) / np.where(
                    sum_frac > 0.0, sum_frac, 1.0
                ),
                0.0,
            )
        d50_surface[t_idx] = d50

        # --- Domain totals -------------------------------------------------
        # Bed mass: Σ over layers and cells (real cells only).
        # Convert g/cm² × cell area would be most accurate, but cell areas
        # are not consistently exposed; we report the cell-summed g/cm²
        # which is conservation-relevant when paired with the same units
        # for the source flux.
        domain_bed_mass[t_idx] = float(layer_mass[:nreal].sum())

        # Suspended mass (mg/L on Riverine side). Real cells only. We
        # convert mg/L * volume(m³) → grams, then to g/cm² equivalent by
        # dividing by total real-cell area (best-effort). Volumes come
        # from VAR_VOLUME if present.
        susp_total_g = 0.0
        for name in suspended_names:
            if name not in mesh.data_vars:
                continue
            da = mesh[name].isel({contracts.DIM_TIME: t_idx}).values
            # mg/L = g/m³. Multiply by cell volume in m³ to get grams.
            if contracts.VAR_VOLUME in mesh.data_vars:
                vol = mesh[contracts.VAR_VOLUME].isel(
                    {contracts.DIM_TIME: t_idx}
                ).values
            else:
                vol = np.ones_like(da)
            susp_total_g += float(
                np.nansum(da[:nreal] * vol[:nreal])
            )
        domain_suspended_mass[t_idx] = susp_total_g

        # Bedload (off in this scenario, so will be zero; capture for
        # completeness if the variable is present).
        if contracts.VAR_BEDLOAD_MASS in mesh.data_vars:
            bl = mesh[contracts.VAR_BEDLOAD_MASS].isel(
                {contracts.DIM_TIME: t_idx}
            ).values
            domain_bedload_mass[t_idx] = float(bl[:nreal].sum())

        # Cumulative source-injection mass per real cell, summed over
        # classes and steps. This is the SSM-side bookkeeping.
        step_source_g_cm2 = 0.0
        for name in suspended_names:
            src_name = f"{name}_source"
            if src_name in mesh.data_vars:
                src = mesh[src_name].values
                step_source_g_cm2 += float(np.nansum(src[:nreal]))
        domain_source_cum_mass[t_idx] = (
            domain_source_cum_mass[t_idx - 1] if t_idx > 0 else 0.0
        ) + step_source_g_cm2

    # ---- Plotting --------------------------------------------------------
    figsize = (10, 6)  # ~1000 × 600 at 100 dpi (we save at 150 dpi).
    dpi = 150

    # Restrict cell-level figures to real cells.
    real_slice = slice(0, nreal)

    # ----- 01: bed mass timeseries ---------------------------------------
    fig, ax = plt.subplots(figsize=figsize)
    n_plot = min(10, nreal)
    cell_idx_plot = np.linspace(0, nreal - 1, n_plot, dtype=int)
    for i, c in enumerate(cell_idx_plot):
        ax.plot(
            np.arange(n_steps),
            bed_total_mass[:, c],
            marker="o", markersize=4, linewidth=1.2,
            label=f"cell {c}",
        )
    ax.set_xlabel("Step index")
    ax.set_ylabel("Total bed mass per cell  (g cm$^{-2}$, summed over layers)")
    ax.set_title(
        "SSM smoke test (plan01_10x5): bed mass under sustained τ_b = 4 Pa"
    )
    ax.legend(
        ncol=2, fontsize=8, frameon=True, loc="best",
    )
    fig.tight_layout()
    out_01 = _OUTPUT_DIR / "01_bed_mass_timeseries.png"
    fig.savefig(out_01, dpi=dpi)
    plt.close(fig)

    # ----- 02: bed shear spatial -----------------------------------------
    fig, ax = plt.subplots(figsize=figsize)
    if "face_x" in mesh.coords or "face_x" in mesh.data_vars:
        fx = np.asarray(mesh["face_x"].values, dtype="float64")[real_slice]
        fy = np.asarray(mesh["face_y"].values, dtype="float64")[real_slice]
    else:
        # Fallback grid layout (10x5 = 50, plan01 has 49 real cells).
        fx = np.arange(nreal, dtype="float64")
        fy = np.zeros(nreal, dtype="float64")
    tau_final = tau_b[-1, real_slice]
    sc = ax.scatter(
        fx, fy, c=tau_final, cmap="viridis", vmin=0.0, vmax=5.0,
        s=120, edgecolor="white", linewidth=0.4,
    )
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Bed shear stress τ_b  (Pa)")
    ax.set_xlabel("face_x")
    ax.set_ylabel("face_y")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title(f"Bed shear stress (Pa) at final step (step={n_steps - 1})")
    fig.tight_layout()
    out_02 = _OUTPUT_DIR / "02_bed_shear_spatial.png"
    fig.savefig(out_02, dpi=dpi)
    plt.close(fig)

    # ----- 03: armoring class fractions ----------------------------------
    # Pick the cell with the largest cumulative bed-mass loss over the
    # window. Tie-break: lowest cell index.
    bed_loss = bed_total_mass[0, :nreal] - bed_total_mass[-1, :nreal]
    cell_focus = int(np.argmax(bed_loss))
    fracs = surface_class_frac[:, cell_focus, :]  # (steps, n_class)
    fig, ax = plt.subplots(figsize=figsize)
    ax.stackplot(
        np.arange(n_steps),
        fracs.T,
        labels=class_labels,
        alpha=0.85,
    )
    ax.set_xlabel("Step index")
    ax.set_ylabel("Surface-layer mass fraction")
    ax.set_ylim(0.0, 1.0)
    ax.set_title(
        f"Surface-layer class fractions at cell {cell_focus} "
        f"(armoring evolution; Δbed_mass = {bed_loss[cell_focus]:.3e} g/cm²)"
    )
    ax.legend(
        ncol=2, fontsize=8, loc="upper right",
        bbox_to_anchor=(1.0, -0.12),
    )
    fig.tight_layout()
    out_03 = _OUTPUT_DIR / "03_armoring_class_fractions.png"
    fig.savefig(out_03, dpi=dpi)
    plt.close(fig)

    # ----- 04: mass conservation -----------------------------------------
    # We compare three quantities, all expressed as "change since t=0"
    # in cell-summed g/cm² (the SSM-native unit), so they live on the
    # same scale and a single y-axis is readable:
    #
    #   Δ bed mass               = domain_bed_mass[t] − domain_bed_mass[0]
    #   Δ bedload mass           = domain_bedload_mass[t] − domain_bedload_mass[0]
    #   cum source ΣΔt(E − D)    = ∫ net source into Riverine over t
    #
    # Conservation invariant for "off" bedload + no across-boundary loss:
    #
    #   Δbed + Δbedload + cum_source ≈ 0
    #
    # We plot all three, plus their algebraic sum (which should hover
    # near zero at float32 precision) and a tolerance band. Suspended
    # mass on the Riverine side is reported as an annotation rather
    # than a line because its native unit (g, from mg/L × m³) is on a
    # very different scale.
    fig, ax = plt.subplots(figsize=figsize)
    steps = np.arange(n_steps)

    delta_bed = domain_bed_mass - domain_bed_mass[0]
    delta_bedload = domain_bedload_mass - domain_bedload_mass[0]
    sum_check = delta_bed + delta_bedload + domain_source_cum_mass

    ax.plot(
        steps, delta_bed,
        label="Δ bed mass (cell-summed g/cm²)",
        marker="o", color="tab:blue",
    )
    ax.plot(
        steps, delta_bedload,
        label="Δ bedload mass (g/cm²) — bedload=off",
        marker="s", linestyle=":", color="#888",
    )
    ax.plot(
        steps, domain_source_cum_mass,
        label="cum source ΣΔt(E−D) (g/cm²)",
        marker="d", color="tab:red",
    )
    ax.plot(
        steps, sum_check,
        label="Σ check: Δbed + Δbedload + cum-source",
        linestyle="--", color="black", linewidth=1.4,
    )

    # Tolerance band around zero. Use a relative tolerance pegged to
    # initial bed mass × float32 epsilon × an order-of-magnitude buffer.
    # float32 ≈ 1.2e-7 → with cell-summed mass ~4188, expect drift
    # ~ 5e-4. We shade ±max(5e-4, 1% of |sum_check|.max()) as visual aid.
    eps_band = max(5.0e-4, 0.01 * float(np.max(np.abs(sum_check)) or 0.0))
    ax.axhspan(-eps_band, eps_band, color="0.85", alpha=0.5, zorder=0,
               label=f"float32 tolerance band ±{eps_band:.2e}")

    ax.axhline(0.0, color="0.4", linewidth=0.6, zorder=0)
    ax.set_xlabel("Step index")
    ax.set_ylabel("Mass change since t=0  (cell-summed g cm$^{-2}$)")
    ax.set_title(
        "Domain total mass conservation: Δbed + cum-source should ≈ 0"
    )
    ax.legend(fontsize=8, loc="best")

    susp_note = (
        f"Suspended tracer mass (Riverine, mg/L × m³): "
        f"{domain_suspended_mass[0]:.3e} g (t=0) → "
        f"{domain_suspended_mass[-1]:.3e} g (t=last)"
    )
    ax.annotate(
        susp_note + f"\nfloat32 round-off ≈ 1e-7 × |bed mass|",
        xy=(0.5, -0.18), xycoords="axes fraction",
        ha="center", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.9),
    )
    fig.tight_layout()
    out_04 = _OUTPUT_DIR / "04_mass_conservation.png"
    fig.savefig(out_04, dpi=dpi)
    plt.close(fig)

    # ---- Report ----------------------------------------------------------
    print("Wrote PNGs:")
    for p in (out_01, out_02, out_03, out_04):
        size_kb = p.stat().st_size / 1024.0
        print(f"  {p}  ({size_kb:.1f} KiB)")
    print(f"\nScenario summary:")
    print(f"  n_steps                 = {n_steps}")
    print(f"  n_face (with ghosts)    = {n_face_total}")
    print(f"  nreal                   = {nreal}")
    print(f"  n_class                 = {n_class}")
    print(f"  τ_b (uniform)           = 4.0 Pa")
    print(
        f"  bed mass at t=0 (Σreal) = "
        f"{domain_bed_mass[0]:.4f} g/cm²"
    )
    print(
        f"  bed mass at t=last      = "
        f"{domain_bed_mass[-1]:.4f} g/cm²  "
        f"(Δ = {domain_bed_mass[-1] - domain_bed_mass[0]:+.4e})"
    )
    print(
        f"  cum source ΣΔt(E−D)     = "
        f"{domain_source_cum_mass[-1]:+.4e} g/cm²"
    )
    print(
        f"  Σ check (Δbed + cum src)= "
        f"{(domain_bed_mass[-1] - domain_bed_mass[0]) + domain_source_cum_mass[-1]:+.4e}"
    )
    print(
        f"  cell of max erosion     = "
        f"{cell_focus} (Δbed_mass = {bed_loss[cell_focus]:+.4e} g/cm²)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
