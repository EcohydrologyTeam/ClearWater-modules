# Design Spec: Pathogen process canonical algae-name reference

**Date:** 2026-05-30
**Component:** `src/clearwater_modules_v3/processes/pathogen.py`
**Severity:** Minor. Pathogen light-attenuation term silently uses zero algae.
**Status:** Proposed.

## Summary

The Pathogen process looks up the phytoplankton state under the name `ap` for
its light-extinction / shading term. The canonical registry name for floating
algae is `algae_floating` (registered by `FloatingAlgae` and by the riverine
bridge). Because `ap` is absent, the process emits

```
Pathogen: optional registry variable 'ap' not present; treating as 0 for the
light-extinction calculation.
```

and computes its light term with zero algal shading. In a coupled run with a
dense bloom this understates the attenuation the pathogen die-off term should
see.

## Evidence

- `FloatingAlgae` and the riverine bridge register floating algae as
  `algae_floating`; no process registers `ap`.
- The Pathogen process reads `ap` (optional) and falls back to 0, printing the
  warning above on every step of a coupled run.

## Required change

In `pathogen.py`, read the floating-algae state under the canonical name
`algae_floating` (matching the other processes that consume it). If backward
compatibility with `ap` is desired during migration, prefer `algae_floating`
and fall back to `ap` only when the canonical name is absent, mirroring the
`tip` / `phosphorus_total_inorganic` precedence pattern already used in the
algae processes.

## Verification

- Unit: register `algae_floating` (nonzero) in a fixture registry, run the
  Pathogen process, and assert the light term reflects the algal shading (no
  fallback warning, result differs from the zero-algae case).
- Confirm no remaining `'ap'`-not-present warnings in a full coupled run.
