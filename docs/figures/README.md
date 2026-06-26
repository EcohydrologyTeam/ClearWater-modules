# APL figures

Presentation/documentation figures for the ClearWater Aquatic Processes
Library (APL). Each figure maps to a slide in
[`../APL_slide_captions.md`](../APL_slide_captions.md).

| File | Slide | Shows |
|---|---|---|
| `apl_context.{png,svg}` | 1 — *what the APL is* | Transport hosts (ClearWater-Riverine, ClearWater-HMS) ↔ the APL (Temperature + 11 WQ Constituents + Riverine coupler) ↔ the shared `VariableRegistry` state contract. |
| `apl_dependency_dag.{png,svg}` | 2 — *the library resolves execution order* | Process dependency graph. Edges are the validated `upstream_processes` ordering constraints, read live from the process classes. |
| `apl_firing_timeline.{png,svg}` | 3 — *per-process timesteps* | Multi-rate firing on a shared base substep. Replicates the exact rule in `Model.__build_process_schedule`. |

## Regenerate

```bash
PYTHONPATH=src python docs/figures/generate_apl_figures.py
```

Run in the conda `clearwater` env so the process classes import. Outputs
300-dpi PNG (slides) and SVG (vector / further editing) for each figure.

Notes:

- The DAG edges are not hand-drawn — they are pulled from each process class's
  `upstream_processes` and `output_variables`, so the figure tracks the code.
- The firing timeline uses an **illustrative** set of per-process timesteps to
  exercise the multi-rate capability; the bundled demo hands every process
  the same Δt.
