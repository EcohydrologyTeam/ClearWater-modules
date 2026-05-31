# Design Spec: Light extinction coefficient computation (NSM1-I port completion)

**Date:** 2026-05-30
**Component:** `src/clearwater_modules_v3/processes/floating_algae.py` (and the
shared optical inputs it would read)
**Severity:** Incomplete Fortran -> Python port. Affects algal light limitation
and therefore the full algae / dissolved-oxygen response.
**Status:** IMPLEMENTED 2026-05-30. `FloatingAlgae` now computes lambda each
step via `utils.light.L` (default; `use_computed_light_extinction=True`), with
the scalar `light_attenuation_coefficient` retained as an override
(`use_computed_light_extinction=False`). Optical defaults (`lambda0/1/2`,
`lambdas`, `lambdam`, `fcom`) added as `_LIGHT_EXT_DEFAULTS`, verified against the
NSM1-I Fortran (see "Implementation note" below). Tests:
`tests/v3/nsm1/test_light_extinction_v3.py`. Trajectory-perturbing → coupled
baseline re-captured. **BenthicAlgae was wired in a follow-up (2026-05-30):** its
own `limit_light` (`exp(-lambda*depth)`, light reaching the bed) now uses the same
computed lambda (`utils.light.L`) from the water-column constituents — the
floating-algae `Ap`, suspended solids, and POC that attenuate light before it
reaches the benthos. The benthic biomass does not enter lambda (NSM1-I computes
one global lambda from the water column). Same `use_computed_light_extinction`
override; trajectory-perturbing → a second baseline re-capture. The Report 2
Willowbend re-run / Methods update remains a follow-up in the Report 2 tracker.

## Implementation note (2026-05-30, verified against the Fortran)

The annotation's call to "resolve against the util, not the prior report text"
was confirmed directly against the NSM1-I Fortran the user provided
(`/Users/todd/GitHub/ecohydrology/ClearWater/modules/fortran/NSM1`):

- `02_global/nsmi_global_params.f90:421-427` assembles lambda as
  `lambda0 + lambdas*Solid + (use_POC) lambdam*POC/focm + (use_Algae) lambda1*Ap
  + lambda2*Ap**0.66667` — **identical** to `utils.light.L`.
- `01_extraneous/nsmi_main.f90` parameter table gives the NSM1-I defaults
  (last value before the `.false.` flag): `lambda0=0.02`, `lambdas=0.052`,
  `lambdam=0.174`, `lambda1=0.0088` (linear, `m⁻¹·(µg-Chla/L)⁻¹`),
  `lambda2=0.054` (nonlinear, `^-2/3`). These match `_LIGHT_EXT_DEFAULTS` /
  Pathogen's `_LIGHT_DEFAULTS` exactly.

So the prior report's lambda1/lambda2/lambdam role inconsistency is resolved:
lambda1/lambda2 are the linear/nonlinear **algal** terms, lambdam is the
**organic-matter** (POC/fcom) term, lambdas is the **suspended-solids** term.

---

**Original proposal (superseded by the implementation above):** This is a
port-completion item, not a new feature.

> **Review note (2026-05-30, code-verified):** The Summary's claim that "no
> optical-component lambda computation exists anywhere in v3" is **incorrect** —
> the search was scoped to `processes/` and missed `utils/`. `utils/light.py:L(...)`
> already computes the NSM1-I optical lambda from `lambda0 / lambda1 / lambda2 /
> lambdas / lambdam, Solid, POC, Ap`, and `processes/pathogen.py` already calls it
> for its light-decay term. The coefficient roles the spec proposes to "resolve
> from the Fortran" are already settled in `utils/light.py` — resolve against the
> util, not the prior report text.
>
> **Corrected scope:** wire `FloatingAlgae` to `utils.light.L` (replacing the
> constant `light_attenuation_coefficient`) and add the coefficient defaults to
> `ALGAE_DEFAULTS` / the shared optical params. This is *wiring an existing util*,
> not implementing an optical model — substantially smaller and lower-risk than
> written. **Shared prerequisite:** a `Solid` (suspended-solids) input source,
> which is currently undefined for a coupled run (a constant default param in
> `benthic_algae`, an optional registry read in `pathogen`). See the matching note
> in `clearwater_modules_v3_phosphorus_partitioning.md`.

## Summary

In NSM1-I (Fortran), the light extinction coefficient lambda is **computed** from
the optical constituents of the water column (background attenuation plus
contributions from inorganic suspended solids, particulate organic carbon, and
algal self-shading). In v3, `FloatingAlgae` instead takes lambda as a scalar
parameter, `light_attenuation_coefficient` (default 1.0 / m), with an explicit
code comment acknowledging the gap:

```python
# floating_algae.py:347-348
# ``light_attenuation_coefficient`` (lambda) is not in ALGAE_DEFAULTS
# (Fortran/v1 compute lambda from the POM/Chla sum in modGlobalParam).
```

No optical-component lambda computation exists anywhere in v3 (confirmed by
search across `processes/`). The light-limitation factor `limit_light`
(`floating_algae.py:815`) consumes this constant lambda directly, so with the
default the model attenuates light at a fixed 1.0 / m regardless of the simulated
POC, solids, and chlorophyll fields. This is the gap that forced the Report 2
demonstration to run at a constant lambda and forced the report to document the
parameterized form rather than the optical model.

## Required change

Implement the NSM1-I optical-component computation of lambda, evaluated each step
from the current state, following \ct{zhang2016aquatic} / the Fortran
`modGlobalParam` formulation:

```
lambda = lambda0
       + lambda_s * S                      (inorganic suspended solids)
       + lambda_1 * POC                     (linear POC attenuation)
       + lambda_2 * (f_com * POC)           (nonlinear POC term)
       + lambda_m * Ap^(2/3)                (algal self-shading)
```

Confirm the exact coefficient roles and the algal exponent against the Fortran
source before implementing (the prior report text had the lambda_1 / lambda_2 /
lambda_m roles inconsistent; resolve from the Fortran, not the prior doc). Add
the coefficients (`lambda0`, `lambda_s`, `lambda_1`, `lambda_2`, `lambda_m`,
`f_com`) to `ALGAE_DEFAULTS` (or the shared global optical parameters) with the
NSM1-I default values. Retain the ability to override lambda with a directly
specified scalar for tests and didactic runs, but make the computed form the
default behavior so a coupled run reproduces NSM1-I optics.

Terms should be conditionally included based on which constituents are present
(POC term only when carbon is active, algal term only when algae are active),
matching the Fortran conditional structure.

## Verification

- Unit: with known POC / solids / Ap fields, assert the computed lambda matches
  the Fortran/NSM1-I value to tolerance; assert each term activates/deactivates
  with its constituent.
- Regression vs Fortran NSM1-I: a single-cell light-limitation comparison
  against the legacy code for a range of chlorophyll and POC values.
- Integration: confirm the coupled Willowbend run with computed lambda produces
  physically sensible algal light limitation (the oxbow, with high Ap, should
  show stronger self-shading than the constant-lambda run).

## Report follow-up (tracked in the Report 2 repo)

Once implemented:
- Re-run the Willowbend demonstration so the results reflect the computed lambda
  (this may change the algal-bloom decline rate and the DO response).
- Restore the optical-component light-extinction equation to the Methods chapter
  (Section 2.2.6, "Light Extinction"), replacing the parameterized-coefficient
  paragraph currently in place, with coefficient roles matching the corrected
  implementation.

## Coordination note

This is the highest-priority NSM1-I port-completion item for Report 2 because it
feeds algal growth and therefore the headline algae and dissolved-oxygen results.
Tracked in the Report 2 repository's code-completion tracking document.
