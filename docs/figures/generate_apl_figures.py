"""Generate the APL presentation figures (PNG + SVG).

Three figures, one per story beat in ``docs/APL_slide_captions.md``:

  1. apl_context          -- what the APL is (system-context diagram)
  3. apl_dependency_dag   -- the library resolves execution order
  4. apl_firing_timeline  -- per-process timesteps (multi-rate schedule)

The dependency edges in figure 2 are read from the live process classes
(``upstream_processes`` and ``output_variables``), so the diagram tracks the
code rather than a hand-drawn snapshot. The firing timeline replicates the
exact schedule rule used by ``Model.__build_process_schedule``: a process
fires at substep ``i`` when ``(i * base_step_seconds) % interval_seconds == 0``.

Run with the conda ``clearwater`` env so the process classes import:

    PYTHONPATH=src python docs/figures/generate_apl_figures.py
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
DASH = (0, (5, 3))

# ---------------------------------------------------------------------------
# Shared style
# ---------------------------------------------------------------------------
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "svg.fonttype": "none",
        "figure.dpi": 110,
    }
)

# ClearWater-ish palette
C = {
    "host": "#1f6f8b",       # transport drivers
    "apl": "#2a9d8f",        # the library
    "apl_fill": "#e7f4f1",
    "registry": "#264653",   # shared state
    "heat": "#e76f51",       # temperature / forcing
    "nutrient": "#457b9d",   # N, P
    "carbon": "#6d6875",     # C, CBOD, alkalinity, DIC
    "oxygen": "#2a9d8f",     # DOX
    "biology": "#52b788",    # algae, pathogen
    "matter": "#b08968",     # POM, N2
    "coupler": "#9d4edd",    # Riverine
    "ink": "#1d3557",
    "muted": "#6c757d",
    "grid": "#ced4da",
}

GROUP_COLOR = {
    "Temperature": C["heat"],
    "Nitrogen": C["nutrient"],
    "Phosphorus": C["nutrient"],
    "Carbon": C["carbon"],
    "CBOD": C["carbon"],
    "Alkalinity": C["carbon"],
    "DOX": C["oxygen"],
    "FloatingAlgae": C["biology"],
    "BenthicAlgae": C["biology"],
    "Pathogen": C["biology"],
    "POM": C["matter"],
    "N2": C["matter"],
    "Riverine": C["coupler"],
}

DISPLAY = {
    "Temperature": "Temperature",
    "Nitrogen": "Nitrogen",
    "Phosphorus": "Phosphorus",
    "Carbon": "Carbon",
    "CBOD": "CBOD",
    "Alkalinity": "Alkalinity",
    "DOX": "DOX",
    "FloatingAlgae": "Floating\nAlgae",
    "BenthicAlgae": "Benthic\nAlgae",
    "Pathogen": "Pathogen",
    "POM": "POM",
    "N2": "N₂",
    "Riverine": "Riverine\n(coupler)",
}


def _save(fig, name):
    for ext in ("png", "svg"):
        path = os.path.join(HERE, f"{name}.{ext}")
        fig.savefig(path, bbox_inches="tight", facecolor="white",
                    dpi=300 if ext == "png" else None)
        print(f"  wrote {os.path.relpath(path, os.path.dirname(HERE))}")
    plt.close(fig)


def _box(ax, x, y, w, h, text, fc, ec, tc="white", fs=11, weight="bold",
         r=2.0, ls="-", z=3, va="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle=f"round,pad=0.3,rounding_size={r}",
                 linewidth=1.8, edgecolor=ec, facecolor=fc, linestyle=ls, zorder=z))
    if text:
        ty = y + h / 2 if va == "center" else y + h - 1.6
        ax.text(x + w / 2, ty, text, ha="center", va=va, color=tc,
                fontsize=fs, fontweight=weight, zorder=z + 1)


def _arrow(ax, x1, y1, x2, y2, color, two=True, lw=2.0, ls="-", rad=0.0, z=5):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                 arrowstyle="<|-|>" if two else "-|>", mutation_scale=14,
                 linewidth=lw, color=color, linestyle=ls,
                 connectionstyle=f"arc3,rad={rad}", zorder=z,
                 shrinkA=1, shrinkB=1))


# ---------------------------------------------------------------------------
# Figure 1 -- APL system-context diagram
#   side-by-side APL (kinetics) and ClearWater-Riverine (transport) over the
#   single shared registry. No arrow crosses a box.
# ---------------------------------------------------------------------------
def figure_context():
    fig, ax = plt.subplots(figsize=(13, 7.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    ax.text(50, 97.5, "ClearWater Modeling System",
            ha="center", va="top", fontsize=17, fontweight="bold", color=C["ink"])

    # ---- shared registry (the seam), bottom, full width ----
    _box(ax, 5, 4, 90, 10.5,
         "ClearWater-data  ·  VariableRegistry\n"
         "the single shared, per-cell, time-indexed state  —  the seam",
         C["registry"], C["registry"], tc="white", fs=11)

    # ---- APL container (left) ----
    _box(ax, 3, 24, 58, 60, "", C["apl_fill"], C["apl"], r=3, z=2)
    ax.text(32, 80.5, "ClearWater Aquatic Processes Library (APL)",
            ha="center", va="center", fontsize=12, fontweight="bold", color=C["apl"])
    ax.text(32, 76, "computes the kinetics", ha="center", va="center",
            fontsize=9.5, style="italic", color=C["apl"])

    # Temperature
    _box(ax, 6, 40, 13, 22, "Temperature\n\n(heat budget)", "white", C["heat"],
         tc=C["heat"], fs=9)

    # WQ Constituents panel
    _box(ax, 21, 28, 27, 40, "", "white", C["nutrient"], r=1.5, z=3)
    ax.text(34.5, 64.5, "WQ Constituents", ha="center", va="center",
            fontsize=9.5, fontweight="bold", color=C["nutrient"])
    chips = ["Nitrogen", "Phosphorus", "Carbon", "DOX", "CBOD", "Alkalinity",
             "N₂", "POM", "Floating\nAlgae", "Benthic\nAlgae", "Pathogen"]
    cols, cw, ch, gx, gy = 3, 7.6, 5.4, 0.7, 1.1
    gx0, gy0 = 22.4, 54.5
    for i, name in enumerate(chips):
        r, c = divmod(i, cols)
        cx = gx0 + c * (cw + gx)
        cy = gy0 - r * (ch + gy)
        ax.add_patch(FancyBboxPatch((cx, cy), cw, ch,
                     boxstyle="round,pad=0.1,rounding_size=0.8",
                     linewidth=1.0, edgecolor=C["nutrient"], facecolor="#eef3f7",
                     zorder=4))
        ax.text(cx + cw / 2, cy + ch / 2, name, ha="center", va="center",
                fontsize=6.7, color=C["ink"], zorder=5)

    # Riverine coupler (right edge of APL)
    _box(ax, 49, 42, 10, 20, "Riverine\ncoupler", "white", C["coupler"],
         tc=C["coupler"], fs=8.5)

    # ---- host transport engine (right) ----
    _box(ax, 65, 42, 19, 22,
         "ClearWater-\nRiverine\n\nHEC-RAS-2D\ntransport\n(advection–diffusion)",
         "white", C["host"], tc=C["host"], fs=8)
    # HEC-RAS HDF input
    _box(ax, 65, 86, 19, 8.5, "HEC-RAS-2D HDF output", "white", C["muted"],
         tc=C["muted"], fs=7.6, weight="normal")
    # second verified host
    _box(ax, 87, 42, 10, 22, "ClearWater-\nHMS\n\n(verified)", "white", C["host"],
         tc=C["host"], fs=7.6)

    # ---- arrows (none cross a box; labels sit clear of the arrows) ----
    # APL <-> registry  (kinetics)
    _arrow(ax, 18, 24, 18, 14.5, C["apl"], lw=2.4)
    ax.text(21.5, 19.25, "kinetics:\nreads forcings, writes reaction updates",
            ha="left", va="center", fontsize=8.2, color=C["apl"])
    # both hosts <-> registry  (transport; share the same registry)
    _arrow(ax, 74.5, 42, 74.5, 14.5, C["host"], lw=2.4)
    _arrow(ax, 92, 42, 92, 14.5, C["host"], lw=2.4)
    ax.text(83, 28, "transport\n(shared registry)",
            ha="center", va="center", fontsize=8.2, color=C["host"])
    # coupler <-> CWR  (the bridge)
    _arrow(ax, 59, 53, 65, 53, C["coupler"], lw=2.2)
    # HEC-RAS-2D HDF -> ClearWater-Riverine
    _arrow(ax, 74.5, 86, 74.5, 64, C["muted"], two=False, lw=1.6)

    ax.text(50, 1.2,
            "The host owns transport and the clock; the APL owns the chemistry; "
            "the registry is the single shared state. One library, multiple verified "
            "hosts — ClearWater-Riverine (HEC-RAS-2D) and ClearWater-HMS.",
            ha="center", va="center", fontsize=9, style="italic", color=C["muted"])

    _save(fig, "apl_context")


# ---------------------------------------------------------------------------
# Figure 2 -- process dependency DAG (edges read from the live classes)
# ---------------------------------------------------------------------------
def _load_processes():
    """Return {name: {'upstream': (...), 'writes': [...]}} from live classes."""
    from clearwater_modules_v3 import processes as P
    names = ["Temperature", "Nitrogen", "Phosphorus", "Carbon", "DOX", "CBOD",
             "Alkalinity", "N2", "POM", "FloatingAlgae", "BenthicAlgae",
             "Pathogen", "Riverine"]
    out = {}
    for n in names:
        cls = getattr(P, n)
        out[n] = {
            "upstream": tuple(getattr(cls, "upstream_processes", ())),
            "writes": list(getattr(cls, "output_variables", []) or []),
        }
    return out


def figure_dependency_dag(procs):
    fig, ax = plt.subplots(figsize=(13, 7.8))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    ax.text(50, 97.5, "APL process dependency graph",
            ha="center", va="top", fontsize=15, fontweight="bold", color=C["ink"])
    ax.text(40, 92.8,
            "execution order resolved from each process's validated upstream_processes",
            ha="center", va="top", fontsize=9.3, style="italic", color=C["muted"])

    nw, nh = 13, 8.6

    # Layout: sources (with downstream edges) in one row, sinks below them,
    # independent processes pulled out into a side panel so no edge ever
    # passes behind an unrelated node.
    src_x = {"Nitrogen": 13, "Carbon": 28, "CBOD": 43,
             "FloatingAlgae": 58, "BenthicAlgae": 73}
    src_y = 55
    sinks = {"DOX": (34, 20), "Phosphorus": (64, 20)}
    pos = {n: (x, src_y) for n, x in src_x.items()}
    pos.update(sinks)

    def node(name, x, y, ec=None, ls="-"):
        col = ec or GROUP_COLOR.get(name, C["ink"])
        ax.add_patch(FancyBboxPatch((x - nw / 2, y - nh / 2), nw, nh,
                     boxstyle="round,pad=0.2,rounding_size=1.4",
                     linewidth=1.8, edgecolor=col, facecolor="white",
                     linestyle=ls, zorder=4))
        ax.text(x, y, DISPLAY.get(name, name), ha="center", va="center",
                fontsize=8.6, fontweight="bold", color=col, zorder=5)

    # Temperature forcing band (data-flow forcing, not an ordering edge)
    ax.add_patch(FancyBboxPatch((4, 77), 75, 11,
                 boxstyle="round,pad=0.2,rounding_size=2",
                 linewidth=0, facecolor="#fdece7", zorder=1))
    node("Temperature", 19, 82.5, ec=C["heat"])
    ax.text(73, 82.5,
            "Temperature forces every Arrhenius rate\n"
            "(read-only forcing — not an ordering edge)",
            ha="right", va="center", fontsize=8.2, color=C["heat"], style="italic")

    # ---- ordering edges (upstream -> reader), read live ----
    def edge(up, dst, rad):
        x1, y1 = pos[up]
        x2, y2 = pos[dst]
        ax.add_patch(FancyArrowPatch((x1, y1 - nh / 2), (x2, y2 + nh / 2),
                     arrowstyle="-|>", mutation_scale=13, linewidth=2.0,
                     color=C["ink"], connectionstyle=f"arc3,rad={rad}",
                     zorder=3, shrinkA=1, shrinkB=2))

    n_edges = 0
    for name, meta in procs.items():
        if name not in pos:
            continue
        for up in meta["upstream"]:
            if up not in pos:
                continue
            dx = pos[up][0] - pos[name][0]
            rad = 0.16 if dx > 0 else (-0.16 if dx < 0 else 0.0)
            edge(up, name, rad)
            n_edges += 1

    for n in src_x:
        node(n, *pos[n])
    node("DOX", *sinks["DOX"])
    node("Phosphorus", *sinks["Phosphorus"])

    # ---- independent-process side panel ----
    px, pw = 84, 14
    ax.add_patch(FancyBboxPatch((px, 29), pw, 41,
                 boxstyle="round,pad=0.3,rounding_size=1.5",
                 linewidth=1.2, edgecolor=C["grid"], facecolor="#f8f9fa", zorder=1))
    ax.text(px + pw / 2, 67, "No execution-order\nconstraint", ha="center",
            va="center", fontsize=8.2, fontweight="bold", color=C["muted"])
    ax.text(px + pw / 2, 61.8, "(run in any order)", ha="center", va="center",
            fontsize=7, style="italic", color=C["muted"])
    for i, name in enumerate(["POM", "N2", "Alkalinity", "Pathogen"]):
        cy = 55.5 - i * 7.4
        col = GROUP_COLOR[name]
        ax.add_patch(FancyBboxPatch((px + 1.6, cy - 2.6), pw - 3.2, 5.2,
                     boxstyle="round,pad=0.1,rounding_size=0.8",
                     linewidth=1.4, edgecolor=col, facecolor="white", zorder=2))
        ax.text(px + pw / 2, cy, DISPLAY[name].replace("\n", " "),
                ha="center", va="center", fontsize=7.6, fontweight="bold",
                color=col, zorder=3)

    # Riverine coupler note
    ax.add_patch(FancyBboxPatch((px, 14), pw, 11,
                 boxstyle="round,pad=0.2,rounding_size=1.2",
                 linewidth=1.6, edgecolor=C["coupler"], facecolor="white", zorder=2))
    ax.text(px + pw / 2, 21.3, "Riverine", ha="center", va="center",
            fontsize=8, fontweight="bold", color=C["coupler"])
    ax.text(px + pw / 2, 17.3, "coupler — exchanges\nall state each step",
            ha="center", va="center", fontsize=6.6, color=C["coupler"])

    # ---- legend ----
    ax.add_patch(FancyArrowPatch((6, 7), (15, 7), arrowstyle="-|>",
                 mutation_scale=13, linewidth=2.0, color=C["ink"]))
    ax.text(16.5, 7,
            "execution-order dependency  (validated  upstream_processes :  "
            "the writer must run before the reader within a substep)",
            ha="left", va="center", fontsize=8.2, color=C["ink"])
    ax.text(43, 2.2,
            f"{n_edges} ordering edges read live from the process classes; "
            "the Model validates this order at init.",
            ha="center", va="center", fontsize=8, style="italic", color=C["muted"])

    _save(fig, "apl_dependency_dag")


# ---------------------------------------------------------------------------
# Figure 3 -- multi-rate firing timeline
# ---------------------------------------------------------------------------
def figure_firing_timeline():
    base_s = 300          # base substep = 5 min (the finest clock)
    horizon_s = 4 * 3600  # 4-hour window
    n = horizon_s // base_s

    # Illustrative per-process timesteps (each an integer multiple of base_s).
    # The bundled demo hands every process the same dt; the framework
    # supports the heterogeneous rates shown here.
    intervals = [  # (process, interval_seconds), ordered fast -> slow
        ("Riverine (coupler)", 300),
        ("Temperature", 300),
        ("DOX / reaeration", 300),
        ("CBOD", 600),
        ("Nitrogen", 900),
        ("Carbon", 900),
        ("Floating Algae", 1200),
        ("Phosphorus", 1800),
        ("Alkalinity", 1800),
        ("N₂", 1800),
        ("POM", 1800),
        ("Pathogen", 3600),
        ("Benthic Algae", 3600),
    ]
    color_for = {
        "Riverine (coupler)": C["coupler"], "Temperature": C["heat"],
        "DOX / reaeration": C["oxygen"], "CBOD": C["carbon"],
        "Nitrogen": C["nutrient"], "Carbon": C["carbon"],
        "Floating Algae": C["biology"], "Phosphorus": C["nutrient"],
        "Alkalinity": C["carbon"], "N₂": C["matter"], "POM": C["matter"],
        "Pathogen": C["biology"], "Benthic Algae": C["biology"],
    }

    fig, ax = plt.subplots(figsize=(12.5, 6.6))
    rows = list(reversed(intervals))  # slowest at bottom row 0
    for r, (name, interval) in enumerate(rows):
        N = interval // base_s
        fires = [i for i in range(n + 1) if (i * base_s) % interval == 0]
        ax.hlines(r, -0.4, n + 0.4, color=C["grid"], lw=0.8, zorder=1)
        ax.scatter(fires, [r] * len(fires), s=46, color=color_for[name],
                   zorder=3, edgecolor="white", linewidth=0.6)
        ax.text(-1.2, r, name, ha="right", va="center", fontsize=9.5, color=C["ink"])
        ax.text(n + 1.4, r, f"Δt = {interval // 60} min   (N={N})",
                ha="left", va="center", fontsize=8.3, color=C["muted"])

    ax.set_xlim(-1, n + 1)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    xticks = list(range(0, n + 1, 6))  # every 30 min
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(t) for t in xticks])
    ax.set_xlabel("base substep index  (base step = 5 min, the finest clock)", fontsize=10)

    secax = ax.secondary_xaxis("top")
    secax.set_xticks(xticks)
    secax.set_xticklabels([f"{t * base_s // 3600:d}:{(t * base_s % 3600) // 60:02d}"
                           for t in xticks])
    secax.set_xlabel("elapsed time  (h:mm)", fontsize=10)

    ax.set_title("Per-process timesteps — multi-rate firing on a shared base substep",
                 fontsize=14, fontweight="bold", color=C["ink"], pad=26)
    ax.text(0.5, -0.17,
            "Each process fires every N base substeps, N = Δt ÷ base step "
            "(exact rule: fires when (i·base) mod Δt = 0). Fast kinetics step "
            "finely; slow pools update less often.  Illustrative rates — the "
            "shipped demo uses a uniform Δt.",
            transform=ax.transAxes, ha="center", va="top", fontsize=8.6,
            style="italic", color=C["muted"], wrap=True)

    fig.subplots_adjust(left=0.16, right=0.86, top=0.84, bottom=0.16)
    _save(fig, "apl_firing_timeline")


def main():
    print("Generating APL figures into docs/figures/ ...")
    figure_context()
    try:
        procs = _load_processes()
    except Exception as exc:  # pragma: no cover
        print(f"  ! could not import process classes ({exc});", file=sys.stderr)
        print("    run with PYTHONPATH=src and the clearwater env.", file=sys.stderr)
        sys.exit(1)
    figure_dependency_dag(procs)
    figure_firing_timeline()
    print("Done.")


if __name__ == "__main__":
    main()
