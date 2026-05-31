# Design Spec: shared `Solid` (suspended-solids) input source

**Date:** 2026-05-30
**Component:** cross-cutting — `processes/pathogen.py`, `processes/phosphorus.py`,
`processes/floating_algae.py`, `processes/benthic_algae.py`, `processes/riverine.py`
(bridge), `examples/nsm1_demo_setup.py`.
**Status:** Proposed.
**Gates:** `clearwater_modules_v3_light_extinction.md`,
`clearwater_modules_v3_phosphorus_partitioning.md` (both name this as their
"shared prerequisite"). Also tightens the `pathogen` light term in coupled runs.

## Summary

Inorganic suspended solids (`Solid`, mg/L) feed the Beer-Lambert light
extinction coefficient (`utils.light.L`, term `lambdas * Solid`) and the
phosphorus sorption fraction (`utils.partitioning.fdp`, term `kdpo4 * Solid`).
Today `Solid` is supplied **two different ways**, and neither works for a
coupled run with spatially varying solids:

| Consumer | How it reads `Solid` today | Fallback when absent |
|---|---|---|
| `pathogen` | optional **registry** read (`_get_optional`, `pathogen.py:235`) | `0.0` + one-time warning |
| `phosphorus` | constructor **parameter** `self.Solid` (`phosphorus.py:393`) | `1.0` (DEFAULTS) |
| `floating_algae` | constructor **parameter** `self.Solid` (`floating_algae.py:526`) | `1.0` (DEFAULTS) |
| `benthic_algae` | constructor **parameter** `self.Solid` (`benthic_algae.py:422`) | `1.0` (DEFAULTS) |

Consequences:

- In a **coupled run** nothing registers `Solid`, so `pathogen` uses `0` (and
  warns every step) while the three kinetic processes silently use the toy
  scalar `1.0` — they disagree about the same physical field.
- The parameter convention is a single scalar, so the three kinetic processes
  cannot see a **per-cell, time-varying** solids field even when one exists
  (riverine mesh, or a future solids model).

This spec defines `Solid` as a single **canonical registry input** that all
four consumers read the same way, with a documented set of providers, while
preserving byte-identical behavior for runs that register no `Solid` (so it
needs **no re-baseline**).

## Canonical definition

- **Registry key:** `Solid`.
- **Units:** mg/L, inorganic suspended-solids (ISS) concentration.
- **Shape:** per-cell `xr.DataArray` (like every other registry forcing);
  scalar-broadcast is allowed.
- **Optionality:** **optional** input. It is NOT added to any
  `Process.variables` list, so it does not become a hard provider requirement
  in `Model.__init_model` Step 8 and standalone runs without a `Solid` provider
  keep working via the per-process fallback below.

## Read convention (the `wind_shelter_coefficient` pattern)

Every consumer reads `Solid` **registry-first, scalar-fallback**, exactly as
`Temperature` reads the optional `wind_shelter_coefficient` forcing
(`temperature.py:470-475`):

```text
if "Solid" in registry:
    Solid = registry.get_at_time("Solid", time)   # shared per-cell field
else:
    Solid = <per-process scalar fallback>          # legacy back-compat
```

- `pathogen` already does this (its `_get_optional` is the same idea; fallback
  `0.0`).
- `phosphorus`, `floating_algae`, `benthic_algae` change from "always use the
  `self.Solid` parameter" to "use the registry `Solid` when present, else
  `self.Solid`". The `self.Solid` parameter is retained as the fallback so
  nothing breaks when no provider is wired.

When a provider registers `Solid`, all four consumers see the **same** field —
that is the "shared source."

## Providers (the "source")

`Solid` can be supplied by any one of these, depending on the run:

1. **Standalone / demo** — register a constant or per-cell `Solid` field in the
   registry (initial conditions). This is how a user pins a realistic ISS value
   for a closed-system or synthetic-mesh run.
2. **Coupled riverine run** — bridge a Clearwater-Riverine suspended-solids mesh
   constituent to `Solid` by adding an entry to `Riverine._MESH_TO_CANONICAL`
   (`riverine.py:22`), e.g. `"<CWR-solids-name>": "Solid"`, *iff* the mesh
   carries one. If the mesh has no solids constituent, the coupling config
   registers a constant `Solid` forcing instead. **Coordinate with the
   ClearWater-riverine session** (it owns `riverine.py` and the mesh constituent
   list); this spec only fixes the canonical name and read convention so the
   bridge mapping is unambiguous.
3. **Future solids model (SSM)** — SSM now lives in its own repo
   (`ClearWater-modules-phase2-SSM`) and is not production-ready. When it
   matures, an SSM process writes `Solid` to the registry each substep and all
   four consumers pick it up with no further change. This spec deliberately
   does **not** implement any solids kinetics.

## Backward compatibility — no re-baseline

The change is read-path-only and **preserves each consumer's existing fallback
value**. The coupled demo (`build_nsm1_demo`) and the standalone test fixtures
register **no** `Solid`, so every consumer falls back to exactly the value it
uses today (`pathogen` → `0.0`; the three kinetic processes → `self.Solid =
1.0`). The 4,320-substep coupled trajectory is therefore byte-identical to the
active baseline `d530a3a` — **no re-baseline is required**. (Verified as an
explicit test below.)

## Known inconsistency — deferred, NOT fixed here

The fallbacks disagree: `pathogen` defaults missing solids to `0.0` while the
three kinetic processes default to `1.0`. In a standalone run with no `Solid`
provider, pathogen sees clear water while phosphorus/algae see 1 mg/L. Unifying
the fallback (e.g. all to `0.0`, or all to a shared `global_vars` default) would
change standalone results and **force a re-baseline**, so it is out of scope for
this step. It is flagged here as a follow-up decision: pick one canonical
no-provider default and re-baseline once. Likewise, **do not** add a `Solid`
field to the demo's default initial conditions in this step — that too would
change the trajectory and force a re-baseline; demo-as-provider is a separate,
reviewed decision.

## Companion input `kdpo4`

`utils.partitioning.fdp(use_TIP, Solid, kdpo4)` needs a partition coefficient
`kdpo4` alongside `Solid`. Wiring `kdpo4` (and setting a non-1.0 default) is the
job of `clearwater_modules_v3_phosphorus_partitioning.md`; this spec only
guarantees a consistent `Solid` field for it to multiply.

## Required changes

1. `phosphorus.py`, `floating_algae.py`, `benthic_algae.py`: read `Solid`
   registry-first with `self.Solid` fallback (small, identical edit at each
   `Solid=self.Solid` call site / forcing-read block). Keep the `self.Solid`
   parameter and its `1.0` default as the fallback.
2. `pathogen.py`: no change (already registry-first).
3. Documentation: this spec; plus a one-line note in `riverine.py`'s
   `_MESH_TO_CANONICAL` that `Solid` is a recognized canonical name awaiting a
   mesh-constituent or constant-forcing provider.

## Verification

- **Unit, per consumer:** register a nonzero per-cell `Solid` in a fixture
  registry, run the process, and assert the light-extinction / partitioning
  term reflects it (differs from the `self.Solid` fallback result).
- **Fallback:** with no `Solid` registered, assert the result equals the
  current `self.Solid`-based result (per-process default preserved).
- **Shared-field consistency:** register one `Solid` field and assert all four
  consumers compute extinction/partitioning from the same value.
- **No-re-baseline guard:** the coupled-demo parity test
  (`test_coupled_demo_parity.py`) still passes bit-identically against
  `baseline_coupled_trajectory_d530a3a.nc` (no `Solid` registered in the demo).
```
