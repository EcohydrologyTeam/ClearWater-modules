# clearwater_modules_v3

`clearwater_modules_v3` is the convergence of two parallel ClearWater Modules tracks:

- **v1** (`clearwater_modules`): function-style framework with kernel optimization, wet-mask gating, hotstart from `xr.Dataset`, latent-heat unit fix, and thin-water stability guard.
- **v2** (`clearwater_modules_v2`): class-based framework with `Process` composition, YAML configuration via `init_from_file`, per-process substepping, and a chunking execution path. Authored by LimnoTech (Paul Tomasula, Anthony Aufdenkampe, Jason Rutyna, Sarah Jordan).

v3 keeps v2's framework as the architectural baseline and folds in v1's optimization and correctness work, plus the synthesis required to make them coexist. v3 is co-developed with LimnoTech as the v2 framework authors and will eventually supersede both v1 and v2.

## Status

v3 is the consolidated, **v3-native** package: v2 has been removed entirely
(the top-level `__init__.py` no longer re-exports any v2 module). TSM and
NSM1 are merged and v3-native; the NSM1 gold-standard correctness gate
(CA-1, SCI-N1, SCI-A3, SCI-A2, CB1, DOX-F1/F2, SCI-A1) is resolved or
explicitly NSM2-deferred with regression tests.

- v3 `Temperature` (`processes/temperature.py`) is **v3-native**: it carries the v1 latent-heat unit correction, the v1 thin-water depth ramp and per-hour `dT/dt` cap, the Fortran-parity sediment thermal diffusivity, and the dynamic sediment-temperature evolution that all three Python ports (v1, v2, v3) had previously dropped relative to the canonical Fortran reference.
- v3 `Model` (`model.py`) is **v3-native**: it adds the kernel-optimization compute schedule (precomputed and timezone-independent), registry-level wet-mask gating, hotstart from `xr.Dataset`, and chunking aligned to integer step indices (immune to floating-point drift in `current_time` arithmetic).
- v3 `init_from_file` (`config/init.py`) is **v3-native**: it accepts the v2 YAML schema unchanged and adds two optional top-level keys, `hotstart` and `wet_mask`.
- The process classes (`Riverine`, `BenthicAlgae`, `FloatingAlgae`, `Nitrogen`, plus `Alkalinity`, `Carbon`, `CBOD`, `DOX`, `N2`, `Pathogen`, `POM`) and the `Process` / `ProcessFactory` base are all **v3-native in-tree** — not v2 re-exports.

The v3 multi-agent code review on 2026-05-04 closed all 10 CRITICAL findings (C1--C10). The authoritative MAJOR triage records **17 of 18 MAJOR findings resolved**; the single remainder (M4) is an explicit NSM2 deferral, guarded and documented. See `design/clearwater_modules_v3_review_findings.md` for the full triage and the NSM1 line-level review under `design/v3_nsm1_review_2026-05-15/`.

### Phase status

| Phase | Scope | Status |
|---|---|---|
| Phase 0 | Inventory and gap analysis (TSM and NSM1) | Complete |
| Phase 1 | Overlay scaffold (historical; v2 since removed) | Complete |
| Phase 2 | v3-native `processes/temperature.py` (merged TSM) | Complete |
| Phase 3 | v3-native `model.py` and `config/init.py` (extra YAML keys) | Complete |
| Phase 4 | Test suite ports and v2/v3 parity tests | Complete |
| Phase 5 | README updates and migration notes for downstream users | Complete |
| Phases 6--10 | v3-native NSM1 (pattern alignment + per-process merge) | Complete |
| Gold-standard gate | NSM1 correctness fixes (CA-1, SCI-N1, SCI-A3, SCI-A2, CB1, DOX-F1/F2) + SCI-A1 NSM2 deferral | Complete |
| Phase R-1 | Review-finding cleanup: all 10 CRITICAL findings resolved (2026-05-04) | Complete |
| Phase R-2 | Review-finding cleanup: physics correctness CRITICAL fixes (sediment-diffusivity, mixing_ratio_air, sediment-T evolution) | Complete |
| Phase R-3 | Review-finding cleanup: MAJOR findings (validation, NaN guards, ordering hazards) | Complete (17 of 18; M4 NSM2-deferred) |

## What's new in v3 (relative to v2)

The following capabilities exist only in v3. They are not in v2 and, in two cases (sediment-diffusivity Fortran parity, dynamic sediment-temperature evolution), they are not in v1 either.

- **Latent-heat unit fix.** v2 evaluated the latent-heat-of-vaporization polynomial with a Kelvin temperature where the polynomial expects Celsius. The error was approximately 26% in `Lv` at 20 °C, propagating directly into the latent heat flux. v3 evaluates the polynomial in Celsius.
- **Thin-water depth ramp + per-hour `dT/dt` cap.** v3 ramps the net heat flux down as wetted depth approaches the configured `q_net_depth_ramp_ref`, and caps `dT/dt` at the configured `dTdt_max_per_hour`. Both guards are required for stability in shallow riverine cells where transient depth approaches zero.
- **Vectorized `mixing_ratio_air` guard.** v2's guard against zero or negative `(P_air - e_air)` was scalar-valued and silently broke on multi-cell `xr.DataArray` inputs. v3's guard uses `xr.where` against `(denom <= 0.0)` and works for any shape.
- **Sediment-diffusivity Fortran parity.** The Fortran TSM reference (HEC-RAS-WQ `modGlobal.f90`) declares `alphas` in m²/day with default `0.0432`. v2 inherited a transcription error introducing the value `0.0061` with an inconsistent docstring claiming m²/s. v3 restores `0.0432 m²/day` and the consistent formula.
- **Dynamic sediment-temperature evolution.** The Fortran TSM evolves `T_sed` on each substep via `ΔT_s = α / (0.5 · h₂²) · (T_w − T_s) · dt / 86400`, paired with the sediment heat flux `q_sed` so that water and sediment exchange identical enthalpy (energy-conservative). v1 and v2 dropped this update, treating `T_sed` as a static forcing. v3 reinstates it under `evolve_sediment_temperature=True` (default).
- **Kernel optimization with precomputed schedule.** v3 `Model` precomputes a per-process firing schedule indexed by integer step number. The hot loop tests `step_index % process_step_indices == 0` rather than recomputing `start_time.timestamp()` modular arithmetic. The schedule is timezone-independent and exact under floating-point time arithmetic.
- **Registry-level wet-mask gating.** v3 `Model.__apply_wet_mask` honors `process.output_variables` (the variables the process writes) rather than `process.variables` (which includes inputs). Dry-cell forcings such as `wind_speed`, `air_temperature`, and `volume` are preserved across substeps; only outputs are NaN-masked. A `getattr` fallback to `variables` retains backward compatibility with processes that have not migrated.
- **Hotstart from `xr.Dataset`.** v3 `Model` accepts a `hotstart_dataset` and `hotstart_timestep`, seeds the registry from the dataset slice at that timestep, and offers each `Process` an optional pair of opt-in hooks `to_hotstart() -> dict` / `from_hotstart(state: dict)` for substep-internal state (e.g., v3 `Temperature.from_hotstart` clears the "skip first time step" flag so the resumed run does not skip its first substep).
- **Chunking aligned to integer step indices.** v3 `Model.__process_loop_chunked` precomputes `interior_chunk_step_indices: set[int]` from `chunk_size_seconds / time_step_seconds`. Boundary detection uses integer comparison, immune to floating-point drift in repeated `current_time +=` arithmetic and immune to timezone effects on `datetime.timestamp()`. `chunk_size` must be an integer multiple of `time_step`; otherwise `init_from_file` raises `ValueError`.
- **YAML schema extensions.** v3 accepts two optional top-level keys, `hotstart` and `wet_mask`. See the schema below.

## Migration

### Top-level package

| v2 import | v3 equivalent |
|---|---|
| `import clearwater_modules_v2 as cwm` | `import clearwater_modules_v3 as cwm` |
| `from clearwater_modules_v2 import Model` | `from clearwater_modules_v3 import Model` |

### Configuration

| v2 import | v3 equivalent |
|---|---|
| `from clearwater_modules_v2.config import init_from_file` | `from clearwater_modules_v3.config import init_from_file` |
| `from clearwater_modules_v2.config import read_config` | `from clearwater_modules_v3.config import read_config` |

### Processes

| v2 import | v3 equivalent |
|---|---|
| `from clearwater_modules_v2.processes.riverine import Riverine` | `from clearwater_modules_v3.processes.riverine import Riverine` |
| `from clearwater_modules_v2.processes.temperature import Temperature` | `from clearwater_modules_v3.processes.temperature import Temperature` |
| `from clearwater_modules_v2.processes.benthic_algae import BenthicAlgae` | `from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae` |
| `from clearwater_modules_v2.processes.floating_algae import FloatingAlgae` | `from clearwater_modules_v3.processes.floating_algae import FloatingAlgae` |
| `from clearwater_modules_v2.processes.nitrogen import Nitrogen` | `from clearwater_modules_v3.processes.nitrogen import Nitrogen` |

### Process base classes

| v2 import | v3 equivalent |
|---|---|
| `from clearwater_modules_v2.processes.base import Process` | `from clearwater_modules_v3.processes.base import Process` |
| `from clearwater_modules_v2.processes.base import ProcessFactory` | `from clearwater_modules_v3.processes.base import ProcessFactory` |

### Utility modules

| v2 import | v3 equivalent |
|---|---|
| `from clearwater_modules_v2.utils.constants import ...` | `from clearwater_modules_v3.utils.constants import ...` |
| `from clearwater_modules_v2.utils.conversions import ...` | `from clearwater_modules_v3.utils.conversions import ...` |

### YAML schema extensions (v3 only)

v3 accepts every key the v2 schema accepts, plus two optional top-level keys:

```yaml
hotstart:                       # optional; if absent, fresh-start semantics
  dataset_path: hotstart.nc     # any path xr.open_dataset can load
  timestep: '2022-05-13 12:00:00'   # the time slice to seed the registry from

wet_mask:                       # optional; if absent, no wet-mask gating
  variable: wetted_surface_area # registry variable used as the mask source
  threshold: 1.0                # cells with mask <= threshold are treated as dry

# everything else identical to the v2 schema
```

`hotstart` seeds the registry from `dataset_path` at `timestep`. v3 `Model` reads the dataset, locates the time slice, and writes each matching variable into the registry before the run starts. Each `Process` may optionally implement `to_hotstart()` and `from_hotstart()` hooks to preserve substep-internal state across the resume; processes that do not implement these hooks fall back to fresh-start semantics for their internal flags (the registry-level state is always preserved).

`wet_mask` enables registry-level wet-mask gating. v3 `Model` reads `variable` from the registry on each substep, computes the boolean mask `(variable > threshold)`, and applies it to each process's `output_variables` after that process runs. Dry-cell outputs are written as NaN; dry-cell forcings are preserved.

## Backward compatibility

Every v2 YAML configuration that does not contain a `hotstart` or `wet_mask` top-level key runs unchanged on v3. The v3 `init_from_file` returns the same `Model` API surface as v2's. The process classes are **v3-native** (not v2 re-exports) but are registered under the same names with the same `ProcessFactory` and accept the v2 YAML schema unchanged, so the migration is import-path-only for any v2 user who is not opting into hotstart or wet-mask. (The v3 NSM1 processes additionally carry the gold-standard correctness fixes; these change kinetics relative to v1/v2 by design and are documented in `src/clearwater_modules_v3/parameter_defaults_corrections.md` and the per-domain audit docs.)

The two v3-only YAML keys are strictly additive. v3 also tolerates v2's existing keys without modification.

## Testing

v3 carries:

- A v3 TSM regression suite mirroring v1's TSM tests (15 calculation tests, 4 latent-heat tests, 6 stability-ramp tests, plus hotstart-roundtrip tests).
- Parity tests pinning v3 outputs against v1 (corrected) and v2 (verified Sumwere baseline) within floating-point tolerance.
- Review-finding regression tests: every CRITICAL finding (C1--C10) has dedicated tests under `tests/v3/` that pin the resolution and prevent regression.

The end-to-end coupled TSM + Riverine demo on Sumwere Creek (4,320 timesteps) runs on v3 with no notebook code changes beyond import statements.

## References

- `design/clearwater_modules_v3_architecture_specification.md` --- umbrella architecture, package layout, v1/v2/v3 contribution table, goals and non-goals.
- `design/clearwater_modules_v3_tsm_design_specification.md` --- TSM-specific design, including the merge plan for each v1 capability into the v2 framework, YAML schema, and the migration table.
- `design/clearwater_modules_v3_tsm_gap_analysis.md` --- Phase 0 inventory of v1 vs v2 vs v3 TSM differences.
- `design/clearwater_modules_v3_review_findings.md` --- multi-agent code review findings (10 CRITICAL, 18 MAJOR, 19 MINOR, 6 observations) with per-finding triage and resolution status.
