# clearwater_modules_v3

v3 is the convergence of two parallel ClearWater Modules tracks:

- **v1** (`clearwater_modules`): function-style framework with kernel
  optimization, wet-mask gating, hotstart from `xr.Dataset`, latent-heat
  unit fix, and thin-water stability guard.
- **v2** (`clearwater_modules_v2`): class-based framework with `Process`
  composition, YAML configuration via `init_from_file`, per-process
  substepping, and a chunking execution path.

v3 keeps v2's framework as the architectural baseline and adds v1's
optimization and correctness work. See:

- `design/clearwater_modules_v3_architecture_specification.md` — umbrella
- `design/clearwater_modules_v3_tsm_design_specification.md` — TSM specifics
- `design/clearwater_modules_v3_tsm_gap_analysis.md` — Phase 0 diff table

## Status

Phase 1 (scaffold). Every public symbol is a thin overlay re-export from v2.
Subsequent phases replace overlay imports with v3-native implementations:

| Phase | Adds |
|---|---|
| Phase 2 | v3-native `processes/temperature.py` (merged TSM) |
| Phase 3 | v3-native `model.py` (kernel optimization, wet-mask, hotstart, chunking) and `config/init.py` (extra YAML keys) |
| Phase 4 | Test suite ports and v2/v3 parity tests |
| Phase 5 | README updates and migration notes for downstream users |

## Migration (TSM)

| v2 import | v3 equivalent |
|---|---|
| `import clearwater_modules_v2 as cwm` | `import clearwater_modules_v3 as cwm` |
| `from clearwater_modules_v2.config import init_from_file` | `from clearwater_modules_v3.config import init_from_file` |
| `from clearwater_modules_v2.processes.temperature import Temperature` | `from clearwater_modules_v3.processes.temperature import Temperature` |
