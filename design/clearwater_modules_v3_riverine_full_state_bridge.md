# Design Spec: Riverine process full-state bridge (Tier 2)

**Date:** 2026-05-30
**Component:** `src/clearwater_modules_v3/processes/riverine.py`
**Severity:** Feature gap. Blocks a full NSM-I run through `init_from_file`; the
manual `build_v3_modules` path is the current workaround.
**Status:** Proposed. Depends on the MeshView item-access fix
(`clearwater_modules_v3_riverine_process_meshview_compat.md`, Tier 1).

## Summary

After the Tier 1 fix, `Riverine.init_process` bridges exactly five constituents
from the transport mesh into the registry (`Ap, NH4, NO3, TIP, DOX`). Every
other NSM1 process state (`organic_nitrogen, n2, organic_phosphorus, poc, doc,
dic, cbod, pom, benthic_algae, alkalinity, pathogen`) is **not** registered from
transport, so a config that enables the full NSM1 process set fails at the first
process that reads an unbridged state, and constituents the riverine model is
transporting are not coupled back to the kinetics. Tier 2 replaces the hardcoded
five-state bridge with a **config-driven, all-state bridge**.

This is the gap that forced the Report 2 nutrient re-run onto the manual
`build_v3_modules` coupling. With Tier 2 in place, `init_from_file` supports the
same full-NSM-I state coverage as `build_v3_modules`.

## Problem

- The bridge set is hardcoded (`Ap, NH4, NO3, TIP, DOX`) and gated on
  `model.has_process("FloatingAlgae")` (Tier 1 widens the gate to mesh presence,
  but still only for those five names).
- The riverine config already lists the constituents it transports. The bridge
  should be derived from that list rather than hardcoded, so any constituent the
  user configures for transport is exchanged with the kinetics under its
  canonical name.

## Required change

1. **Canonical map for the full NSM-I state set** (CW-Riverine fork name ->
   v3 canonical registry name):

   ```python
   _MESH_TO_CANONICAL = {
       "Ap": "algae_floating",   "Ab": "benthic_algae",
       "NH4": "ammonium",        "NO3": "nitrate",
       "OrgN": "organic_nitrogen", "N2": "n2",
       "TIP": "tip",             "OrgP": "organic_phosphorus",
       "POC": "poc",             "DOC": "doc",            "DIC": "dic",
       "CBOD": "cbod",           "POM": "pom",
       "Alk": "alkalinity",      "PX": "pathogen",
   }
   ```

2. **Bridge by mesh presence.** Register the canonical alias for every entry of
   the map that the transport mesh actually carries (item access, shared buffer
   via `copy(deep=False)` so the coupling stays two-way), exactly as the Tier 1
   loop does for the five states.

3. **Two-way vs transport-only.** Optionally honor a per-constituent flag in the
   riverine config (default two-way) so a constituent can be transported without
   being reacted, or reacted without overwriting transport. The default
   (two-way for every bridged constituent) reproduces `build_v3_modules`.

4. **Forcings unchanged.** `depth`, `volume`, `wetted_surface_area`, and the
   met forcings continue to be provided as today (Tier 1 keeps the `depth`
   placeholder guard).

## Verification

- Unit: a registry/`MeshView` fixture carrying all 16 fork-named constituents;
  assert every canonical alias is registered and that mutating a process state
  is visible on the corresponding mesh array (shared-buffer check).
- Integration: `init_from_file` with a config listing all 16 constituents and
  the full NSM1 process set advances `model.run()` to completion and reproduces,
  within tolerance, the `build_v3_modules` trajectory on the same inputs.
- Regression: a config listing only five constituents still bridges only those.

## Coordination note

This is the change that lets the Report 2 (and future) coupled runs use the
canonical config-driven path instead of a hand-rolled stepping loop. It does not
change results, only the assembly path. Reproducing the coupled run also needs a
single environment carrying ClearWater-Riverine, clearwater-modules (v3), and
clearwater-data at compatible versions (the Report 2 run used the riverine `dev`
pixi env with modules-v3 on `PYTHONPATH`).
