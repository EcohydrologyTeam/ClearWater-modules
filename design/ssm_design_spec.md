# SSM design spec — moved

The Sediment Simulation Module (SSM) has been extracted to its own
repository:

**`ClearWater-modules-phase2-SSM`**
(local path: `../ClearWater-modules-phase2-SSM/`)

The full algorithmic specification now lives at
`design/ssm_design_spec.md` inside that repository.

See also `design/ssm_improvement_plan.md` (defect inventory and
migration plan), `design/ssm_consolidation.md` (Sanford-Maa cohesive-bed
consolidation), and `design/ssm_bedload_functions.md` (bedload closure
review) — all in the new repo.

## Why it moved

SSM is a research prototype with documented CRITICAL defects. It was
decoupled from `ClearWater-modules-streaming` so v3 TSM and NSM kinetic
work can proceed on its own cadence, and so SSM can sit alongside ESM
and the planned BSM as a sibling repository.
