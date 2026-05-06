# Phase 9.F.4 — POM `h2` physical-role research

## Summary

- **Physical meaning:** `h2` is the **active sediment layer thickness (m)** — i.e., the
  thickness of the bed-sediment compartment into which water-column particulates
  settle, are buried, and from which dissolved species can flux back to the water.
  Both the legacy Fortran `modGlobalParam.f90` declaration and every v1 NSM1
  Python docstring state this explicitly. The Phase 0 audit's "unclear physical
  role" `FIXME` is unjustified by the source code; the meaning is documented in
  v1 and Fortran but was lost in the v3 parameter file.
- **Cross-model analog:** `h2 = 0.1 m` is exactly the Di Toro (2001) /
  QUAL2K convention for the **lower anaerobic sediment layer thickness `H₂`**
  (10 cm). QUAL2K p. 68 (Section 5.6, SOD/Nutrient Flux Model) describes a
  two-layer Di Toro sediment with "a thin (≈ 1 mm) surface aerobic layer
  underlain by a thicker (10 cm) lower anaerobic layer," and Eq. 214
  (p. 70) defines `H₂` = "the thickness of the anaerobic layer [m]" used as
  the volumetric divisor in the bed-POM mass balance. NSM1's POM compartment
  represents bed POM (Fortran state variable `POM2`, the "2" suffix denoting
  the lower / anaerobic layer in Di Toro notation), so `h2` is the canonical
  `H₂` thickness, not a misnamed parameter or unit-conversion factor.
- **Recommended v3 disposition:** **Document, do not rename or revalue.** The
  parameter name (`h2`), value (`0.1` m), and physical meaning are all correct
  and traceable to Di Toro (2001) and QUAL2K. The fix is editorial only:
  replace the misleading `FIXME(phase1-audit)` comment in
  `src/clearwater_modules_v3/parameters/pom.py` with the correct description
  ("active sediment layer thickness, m; Di Toro (2001) / QUAL2K H₂ anaerobic
  layer convention") and add a one-paragraph docstring note in
  `src/clearwater_modules_v3/processes/pom.py` explaining that POM in NSM1
  represents the bed-sediment compartment (Fortran `POM2`), and that `h2`
  is the divisor that converts areal fluxes (m·mg/L/d) from the water column
  into volumetric concentration changes (mg/L/d) in that bed layer. A
  separate, **higher-priority** finding is unrelated to `h2` itself: see
  the "Related finding" section below on the v1 `vb` units mismatch.

## Trace through v3 / v1 / Fortran code paths

### v3 `pom.py`

File: [`src/clearwater_modules_v3/processes/pom.py`](../src/clearwater_modules_v3/processes/pom.py)

`h2` is consumed in three rate terms in `POM.rate(...)` (lines 261–270):

| Term | Formula | Dimensional check |
|---|---|---|
| `rate_burial` (line 265) | `vb * pom / h2` | `[m/d] · [mg/L] / [m] = [mg/L/d]` ✓ |
| `rate_poc_settling` (line 270) | `vsoc * poc / h2 / fcom` | `[m/d] · [mg/L] / [m] / [-] = [mg/L/d]` ✓ |
| `rate_algal_settling` / `rate_benthic_mortality` | read from upstream Process caches that already include `/ h2` | both already in `[mg/L/d]` ✓ |

Module docstring (lines 26–27) already states correctly: "POM is normalized
by `h2` (active sediment layer thickness, m), not by water depth." The
parameter file [`src/clearwater_modules_v3/parameters/pom.py`](../src/clearwater_modules_v3/parameters/pom.py)
line 11 is the only place still carrying the stale `FIXME` comment:

```python
'h2': 0.1,   # FIXME(phase1-audit): m; sediment burial/sedimentation depth denominator, unclear physical role
```

This `FIXME` contradicts the docstring two files away and should be cleared.

### v1 `processes.py` POM block (lines 2185–2313)

File: [`src/clearwater_modules/nsm1/processes.py`](../src/clearwater_modules/nsm1/processes.py)

Every `h2` consumer in v1 documents the parameter consistently as
**"active sediment layer thickness (m)"**:

- `POM_algal_settling` (L2200–2218): `vsap * Ap * rda / h2`, docstring L2213
- `POM_POC_settling` (L2236–2254): `vsoc * POC / h2 / fcom`, docstring L2248
- `POM_benthic_algae_mortality` (L2257–2277): `Ab * kdb_tc * Fb * (1-Fw) / h2`, docstring L2272
- `POM_burial` (L2281–2293): `vb * POM / h2`, docstring L2291

The static-variable registry [`src/clearwater_modules/nsm1/static_variables.py:921–927`](../src/clearwater_modules/nsm1/static_variables.py)
is the canonical declaration:

```python
Variable(
    name='h2',
    long_name='active sediment layer thickness',
    units='m',
    description='active sediment layer thickness',
    use='static'
)
```

Default value: `h2 = 0.1` m (`src/clearwater_modules/nsm1/constants.py`).

### Fortran `modPOM.f90` and `modGlobalParam.f90`

Files (legacy Fortran NSM1):
- `/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/modPOM.f90`
- `/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/modGlobalParam.f90`

`modGlobalParam.f90` line 38 declares `h2` with an inline comment:

```fortran
real(R8), allocatable, dimension(:) :: h2   ! active Sediment layer thickness (m)
```

with default `h2 = 0.1` (line 134). This is the single authoritative
upstream definition.

`modPOM.f90` operates on the bed-POM state variable `POM2` (lines 6, 22, 113),
not water-column POM. The four consumer terms are:

```fortran
ApSettling_POM2  = vsap(r) * Ap * rda(r) / h2(r)              ! L98
AbDeath_POM2     = AbDeath * Fb(r) * (1.0 - Fw(r)) / h2(r)    ! L103
POCSettling_POM2 = vsoc(r) * POC / focm(r) / h2(r)            ! L108
POM2_Burial      = vb(r) / 365.0 * POM2 / h2(r)               ! L114
```

Dimensional analysis confirms `h2` is the bed-layer thickness used to
convert areal sources to bed volumetric concentrations.

**Note on the `2` subscript:** `POM2` denotes the bed compartment, matching
the Di Toro `H₂` / "layer 2" convention. The water-column counterpart in
some Di Toro derivatives is `POM1` (or surface aerobic layer); NSM1 only
implements `POM2` as a state variable, with water-column algal/POC sources
treated as fluxes into it. v1's rename of `POM2` → `POM` (and v3's adoption
of that rename) drops this subscripting clue; the conceptual identity as
the bed/anaerobic compartment must therefore be carried in docstrings.

## Cross-reference with QUAL2K / Di Toro

QUAL2K v2.11 documentation (Chapra, Pelletier, Tao; December 2008), pp. 68–75,
Section 5.6 "SOD/Nutrient Flux Model":

- p. 68: "Sediment nutrient fluxes and sediment oxygen demand (SOD) are
  based on a model developed by Di Toro (Di Toro et al. 1991, Di Toro and
  Fitzpatrick 1993, Di Toro 2001)... The sediments are divided into 2
  layers: a thin (≈ 1 mm) surface aerobic layer underlain by a thicker
  (10 cm) lower anaerobic layer."
- p. 70, Eq. 214 (POC mass balance, anaerobic layer):
  `H₂ · dPOC₂,G1/dt = J_POC,G1 − k_POC,G1 · θ^(T−20) · H₂ · POC₂,G1 − w₂ · POC₂,G1`
  with `H₂ = the thickness of the anaerobic layer [m]`.
- pp. 71–75: `H₂` is reused as the divisor / volumetric scaling factor in the
  ammonium (Eq. 219), nitrate (Eq. 227), and inorganic phosphorus (Eq. 243)
  anaerobic-layer mass balances. `H₁` (aerobic layer, ≈ 1 mm) plays the
  analogous role for the surface layer (Eqs. 218, 226, 242).

The NSM1 `h2 = 0.1` m is therefore identical in name, value, and physical
meaning to the QUAL2K / Di Toro `H₂` (lower anaerobic sediment layer
thickness). NSM1 omits the `H₁` aerobic layer because it does not implement
the two-layer Di Toro flux model — it carries only the bed POM compartment
as a state variable, with first-order burial/dissolution kinetics rather
than the full Di Toro diagenesis. (Implementing the full two-layer flux
model is the future Phase 4.x "NSM2 sediment diagenesis" scope per the
project memory note.)

## Recommendation and rationale

**Recommendation:** Editorial documentation fix only — no rename, no value
change.

1. In [`src/clearwater_modules_v3/parameters/pom.py`](../src/clearwater_modules_v3/parameters/pom.py)
   line 11, replace the stale `FIXME(phase1-audit)` comment with:
   `# m; active sediment layer thickness (Di Toro 2001 / QUAL2K H₂ anaerobic-layer convention, ≈ 10 cm)`.
2. In [`src/clearwater_modules_v3/processes/pom.py`](../src/clearwater_modules_v3/processes/pom.py)
   module docstring, add a one-paragraph note that NSM1's POM state
   variable represents the **bed-sediment** POM compartment (Fortran
   `POM2`), and that `h2` converts areal water-column fluxes (m · mg/L/d)
   into bed volumetric concentration changes (mg/L/d). Cross-reference
   QUAL2K Section 5.6 / Di Toro (2001).
3. Optional: register `h2` in a v3 static-variables / parameter-metadata
   registry analogous to v1's `static_variables.py:921` so the long_name,
   units, and description are queryable, not buried in a comment.

**Rationale:** The parameter is correctly named, correctly defaulted, and
dimensionally consistent with the canonical Di Toro / QUAL2K framework.
The Phase 0 audit comment "unclear physical role" reflects a documentation
gap (the description was not carried over from v1's `static_variables.py`
into the v3 `parameters/pom.py` dict), not a substantive scientific issue.
Renaming or revaluing would create a diff against v1 and the Fortran
reference for no scientific benefit.

## Related finding (out of scope for Phase 9.F.4 but flagged)

While tracing `h2`, I noticed a separate latent issue worth recording for a
future audit phase: the **`vb` burial-velocity units convention differs
between v1/v3 and Fortran**.

- Fortran (`modGlobalParam.f90:39`): `vb` units are `m/a` (m/yr), default
  `0.0025 m/yr`. The Fortran `POM2_Burial = vb / 365.0 * POM2 / h2` term
  performs the m/yr → m/d conversion inline.
- v1 (`static_variables.py:235`): `vb` units relabeled to `m/d`, default
  `0.01 m/d` (`constants.py:325`), and the burial term is
  `vb * POM / h2` with the inline comment "removed 365 from FORTRAN"
  (`processes.py:2293`).
- v3 (`processes/pom.py:265`): inherits the v1 form `self.vb * pom / self.h2`
  with default `vb = 0.01 m/d` (inline `_POM_GLOBAL_DEFAULTS`).

If a user supplied the Fortran-default value `0.0025` interpreting it as
m/d (per v1/v3 docstring) without the `/365` conversion, the bed-POM burial
flux would be ~365× faster than the legacy Fortran model. Conversely, the
v1/v3 default of `vb = 0.01 m/d` corresponds to `3.65 m/yr` — over 1000×
the Fortran default of `0.0025 m/yr`. The v1 author's "removed 365 from
FORTRAN" change effectively redefined the parameter; it is unclear whether
the new default `0.01 m/d` was deliberately chosen or back-derived from
typical literature `vb` values. **This belongs in a separate audit ticket
covering parameter-default provenance**, not in the `h2` task. No action
recommended here.

## Sources

- [Legacy Fortran `modPOM.f90`](file:///Users/todd/Downloads/NSM_comparison/NSM1/Source%20Files/modPOM.f90)
- [Legacy Fortran `modGlobalParam.f90`](file:///Users/todd/Downloads/NSM_comparison/NSM1/Source%20Files/modGlobalParam.f90) (lines 38, 134, 197 for `h2`)
- [v1 `processes.py` POM block](../src/clearwater_modules/nsm1/processes.py) (lines 2185–2313)
- [v1 `static_variables.py` `h2` declaration](../src/clearwater_modules/nsm1/static_variables.py) (lines 921–927)
- [v3 `processes/pom.py`](../src/clearwater_modules_v3/processes/pom.py)
- [v3 `parameters/pom.py`](../src/clearwater_modules_v3/parameters/pom.py)
- QUAL2K Documentation (Chapra, Pelletier, Tao, December 2008), Section 5.6, pp. 68–75 — local copy at `~/.claude/projects/-Users-todd-GitHub-ecohydrology-ClearWater-modules-streaming/683eb018-76f5-4fd6-a559-c4dcbad5bf26/tool-results/webfetch-1778027776950-krf3re.pdf`
- Di Toro, D.M., 2001. *Sediment Flux Modeling*. Wiley-Interscience. (referenced by QUAL2K as the source of the two-layer `H₁`/`H₂` framework)
- Di Toro, D.M., Paquin, P.R., Subburamu, K., Gruber, D.A., 1991. Sediment oxygen demand model: methane and ammonia oxidation. *J. Environ. Eng.* 116(5).
- Di Toro, D.M., Fitzpatrick, J.J., 1993. Chesapeake Bay sediment flux model. Tech. Rep. EL-93-2, U.S. Army Eng. Waterways Exp. Stn., Vicksburg, MS.
