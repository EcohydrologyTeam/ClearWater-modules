# SSM bedload transport functions — design memo

**Status:** Implemented (Stage-1 menu, seven functions; Stage-2 solver wiring)
**Date:** 2026-05-02
**Scope:** Pluggable bedload transport-rate closure for the Sediment
Simulation Module (SSM) in `clearwater_modules_v2.processes.sediment`.
**Files:** `bedload.py`, `contracts.py`, `ssm.py`, `tests/sediment/test_bedload.py`.
**Audience:** Modelers selecting a transport closure for production runs;
SSM developers adding or extending closures.

---

## 1. Motivation

The original SSM port (Stage 1, before this memo) shipped exactly one
bedload transport function: van Rijn (1984a) Part II. That is faithful
to the EFDC SEDZLJ algorithm (`s_bedload.f90`), but it is *not* what the
mainstream sediment-transport competitors offer:

| Code               | Bedload-formula menu                                                           |
|--------------------|--------------------------------------------------------------------------------|
| **SRH-2D** (USBR)  | van Rijn, Wu, Wilcock-Crowe, Parker, Engelund-Hansen, Yang, Meyer-Peter-Müller |
| **TELEMAC + Gaia** | Meyer-Peter-Müller, Einstein-Brown, Engelund-Hansen, van Rijn, Wilcock-Crowe   |
| **MIKE 21**        | Engelund-Hansen, van Rijn, Meyer-Peter-Müller, Yang, Wu, Garcia-Parker         |
| **Delft3D-MOR**    | van Rijn, Engelund-Hansen, Meyer-Peter-Müller, Wilcock-Crowe, Soulsby-Van Rijn |

The single-formula SSM cannot compete on this axis. This memo documents
the abstraction we added to put SSM on equal footing: a Protocol-based
pluggable closure with seven peer-reviewed transport functions ready to
select via YAML, plus a registry that future formulas can extend without
touching the solver code.

The change is fully backwards-compatible: the default
(`transport_function: van_rijn`) reproduces the prior behaviour to
within float64 round-off (`test_van_rijn_class_matches_helpers`).

---

## 2. Abstraction

### 2.1 Protocol

```python
@runtime_checkable
class BedloadTransportFunction(Protocol):
    name: str  # e.g. "wilcock_crowe", "parker", "yang", ...

    def transport_rate(
        self,
        tau_pa: xr.DataArray,        # (nface,) bed shear stress, Pa
        d50_um: float,               # class median grain size
        tau_ce_pa: float,            # class critical shear for erosion
        velocity_m_s: xr.DataArray,  # (nface,) depth-averaged velocity
        depth_m: xr.DataArray,       # (nface,) hydraulic depth
        slope: xr.DataArray | float, # energy slope
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray: ...           # q_b in g cm⁻¹ s⁻¹
```

### 2.2 Output convention

Every closure returns the per-cell bedload mass discharge per unit
channel width in **g cm⁻¹ s⁻¹**. This is van Rijn's natural unit and
matches the existing CGS internals of SSM (CBL is g/cm², bedload
velocities are cm/s, the saltation height is cm).

To convert:

* Volumetric rate per unit width (cm² s⁻¹): `q_b_vol = q_b / ρ_s`.
* SI mass rate per unit width (kg m⁻¹ s⁻¹): `q_b_si = q_b / 10`.
* Dimensionless Einstein parameter:
  `q_b* = q_b / (ρ_s · √((s − 1) g D₅₀³))`.

### 2.3 The `registry_context` argument

Two of the seven formulas need information beyond a single grain size:

| Formula        | Required context keys                                              |
|----------------|---------------------------------------------------------------------|
| Wilcock-Crowe  | `surface_sand_fraction` (Fs), `surface_class_fraction`, `surface_geometric_mean_um` |
| Wu (2000)      | `pe_ph_ratio` (exposed/hidden probability ratio), `surface_class_fraction` |

The orchestrator builds this dict from the bed surface composition
(`bed.layer_class_fraction[surface_layer]`) and hands it to the closure
on each `transport_rate` call. When the orchestrator omits the dict
(e.g. unit tests, single-class runs), the closures fall back to benign
defaults — Wilcock-Crowe defaults to `Fs = 0.15` (typical gravel-bed
condition), Wu defaults to `pe/ph = 1` (no hiding/exposure correction).

### 2.4 Selection via YAML

```yaml
sediment:
  bedload:
    solver: standalone               # or riverine, off
    transport_function: van_rijn     # or wilcock_crowe, parker, yang,
                                     # wu, engelund_hansen, toffaleti
```

`from_config` validates the name eagerly against
`contracts.BEDLOAD_FUNCTIONS` and raises `ValueError` on a typo.

---

## 3. The seven supported closures

Each is implemented as a class in `bedload.py` and registered in
`BEDLOAD_TRANSPORT_FUNCTIONS`:

| YAML name          | Class                                       | Citation                  | Best for                   |
|--------------------|---------------------------------------------|---------------------------|----------------------------|
| `van_rijn`         | `VanRijn1984TransportFunction`              | van Rijn 1984a            | Sand, sand-gravel mixtures |
| `wilcock_crowe`    | `WilcockCrowe2003TransportFunction`         | Wilcock & Crowe 2003      | Sand-gravel mixed beds     |
| `parker`           | `Parker1990TransportFunction`               | Parker 1990               | Gravel rivers              |
| `yang`             | `YangTransportFunction`                     | Yang 1973, 1979           | Sand (total load)          |
| `wu`               | `Wu2000TransportFunction`                   | Wu, Wang & Jia 2000       | Non-uniform sediment       |
| `engelund_hansen`  | `EngelundHansen1967TransportFunction`       | Engelund & Hansen 1967    | Sand-bed total load        |
| `toffaleti`        | `Toffaleti1968TransportFunction`            | Toffaleti 1968            | Sand-bed depth-integrated  |

### 3.1 van Rijn (1984)

Wraps the existing `van_rijn_bedload_velocity_cm_s`,
`van_rijn_bedload_height_cm`, and `van_rijn_equilibrium_concentration`
helpers via `q_b = u_BL · δ_BL · C_eq`. Mathematically equivalent to
van Rijn's eq. 23 to within the documented 5 % precision; verified
by `test_van_rijn_class_matches_helpers`.

**Domain of applicability:** sand and fine gravel (D₅₀ 0.2–2 mm).
**Strengths:** consistent unit set with the existing SSM CGS internals.
**Weaknesses:** sand-only; over-predicts gravel transport.

### 3.2 Wilcock & Crowe (2003)

Surface-based bedload for sand-gravel mixtures. Uses the surface sand
fraction `Fs` to set the reference Shields stress for the geometric-mean
grain size, then applies a hiding-and-exposure correction
`b_i = 0.67/(1 + exp(1.5 − d_i/d_sg))` to obtain a per-class reference.
Two-regime closure on the stress ratio `φ = τ/τ_ri`:
`W*_i = 0.002·φ^7.5` if `φ < 1.35`, else `14·(1 − 0.894/√φ)^4.5`.

**Domain:** sand-gravel mixtures, particularly those with `Fs` in 0.05–0.45.
**Strengths:** State of the art for *mixed-bed* transport.
**Weaknesses:** Needs surface composition (provided via `registry_context`);
ill-defined on uniform-sand or pure-gravel beds.

### 3.3 Parker (1990)

Surface-based gravel formula with similarity collapse onto a single
universal function `W*(φ_50)`. Three regimes: low (`φ < 0.95`,
`W* = 0.00218·φ^14.2`), mid (exponential blow-up
`exp(14.2·x − 9.28·x²)` for `x = φ − 1`), high (asymptotic
`(1 − 0.853/φ)^4.5`).

**Domain:** Gravel-bed rivers (D₅₀ ≥ 4 mm); validated on Oak Creek and
the Boise River.
**Strengths:** The standard reference for pure-gravel transport.
**Weaknesses:** Sensitive to φ in the mid regime (exponential); sand-bed
predictions are unreliable.

### 3.4 Yang (1973, 1979)

Total-load formula based on **unit stream power** `V·S`. Sand formula
(`D < 2 mm`):

```
log Ct = 5.435 − 0.286·log(ws·d/ν) − 0.457·log(u*/ws)
       + (1.799 − 0.409·log(ws·d/ν) − 0.314·log(u*/ws))·log((V·S − Vcr·S)/ws)
```

with the incipient-motion velocity `Vcr` from Yang 1973 eq. 8 and the
settling velocity `ws` from Cheng 1997. Above 2 mm we switch to Yang
1984 gravel coefficients.

**Domain:** Sand-bed rivers with D₅₀ in 0.062–10 mm; validated on 1289
USGS data points.
**Strengths:** Robust *total*-load predictor; widely used in the U.S.
**Weaknesses:** Fits a *total*-load relation, not a bedload-only one.
SSM treats Yang as a bedload proxy by feeding the entire `q_t` through
the bedload window; the orchestrator can scale by a separate suspended
fraction if needed.

### 3.5 Wu, Wang & Jia (2000)

Non-uniform-sediment formula with hiding-and-exposure correction. Per-
class critical Shields stress `τ_ci = 0.03·(ρ_s − ρ_w)·g·d_i ·
(p_e/p_h)^(-0.6)`. Bedload component:
`Φ_b,i = 0.0053·(τ_b/τ_ci − 1)^2.2`.
Returns `q_b,i = Φ_b · p_i · ρ_s · √((s−1) g d_i³)`.

**Domain:** Non-uniform sand-gravel beds (D₅₀ 0.06–128 mm); validated on
70+ flume datasets.
**Strengths:** Single formula spans sand-gravel transition cleanly.
**Weaknesses:** Hiding/exposure ratio `p_e/p_h` requires the full bed
size distribution; defaults to 1 (uniform-bed limit) when SSM provides
no context.

### 3.6 Engelund & Hansen (1967)

Total-load formula for sand-bed rivers:
`q_t = 0.05·V⁵ / (√g · C³ · Δ² · d)`,
with Chézy `C = V/√(R·S)` (`R ≈ h` in wide channels) and submerged
specific gravity `Δ = (ρ_s − ρ_w)/ρ_w`.

**Domain:** Sand-bed rivers, fine to medium sand (D₅₀ 0.15–0.93 mm).
Validated on Danish flume data.
**Strengths:** Closed-form, no calibration constants beyond the
prefactor `0.05`.
**Weaknesses:** Treats *total* load (bed + suspended); fails on gravel;
pathological at zero slope (C → ∞). SSM guards against `S = 0` by
returning zero.

### 3.7 Toffaleti (1968)

Depth-integrated total-load procedure used by USACE channel-stability
studies. Toffaleti's full procedure decomposes the suspended-load
profile into four vertical zones (lower, middle, upper, surface) using
power-law exponents `z_i = w_s/(κ u_*)`. **SSM ships a single-zone
reduction** of the BR-1 procedure; the multi-zone integral is queued
for a phase-3 enhancement when SSM gains a coupled vertical concentration
profile.

The single-zone form: `q_t = M·V^n_v · d^0.33 / (0.00058)^(n_v − 1)` in
tons day⁻¹ ft⁻¹, with `M = 0.6` and `n_v = 1.5` (medium sand at 60 °F).
Converted to g cm⁻¹ s⁻¹ via the factor `(907.185 / 86400) / 0.3048 · 10`.

**Domain:** Sand-bed channels (D₅₀ 0.062–4 mm).
**Strengths:** Federal-agency standard (USACE).
**Weaknesses:** Single-zone reduction loses the depth-integration
fidelity Toffaleti is known for; the constants `M = 0.6, n_v = 1.5`
are pinned to medium sand at 60 °F. **Use with caution** until the
multi-zone form lands.

---

## 4. Selection guidance

| Situation                                | Recommended closure                         |
|------------------------------------------|---------------------------------------------|
| Sand bed, simple geometry                | `van_rijn` (default) or `engelund_hansen`   |
| Sand-gravel mixed bed                    | `wilcock_crowe`                             |
| Pure gravel bed                          | `parker`                                    |
| Need a U.S. regulatory standard          | `yang` (sand) or `toffaleti` (USACE BR-1)   |
| Non-uniform bed with hiding/exposure     | `wu`                                        |
| Comparing against legacy SEDZLJ output   | `van_rijn` (preserves prior behaviour)      |

When in doubt, run two closures and report the spread as an
uncertainty band — the seven formulas commonly disagree by an order of
magnitude, which is a fair representation of bedload predictability.

---

## 5. Implementation notes

### 5.1 Solver vs. closure separation

The two solvers (`BedloadStandaloneExplicit`,
`BedloadRiverineConstituent`) now consume the configured
`BedloadTransportFunction` on every step. See §6 ("Solver wiring")
below for the full mechanism, the `q_b → u_eff` derivation, and the
backwards-compatibility argument. The closure is also exposed on the
SSM instance (`ssm._bedload_transport_function`) for external code that
wants to compute `q_b` directly — for instance, a mass-budget reporter
or a future implicit bed-evolution step.

### 5.2 Why the closures take SI inputs but emit CGS

The Protocol takes Pa, m, m/s, m²/s for the input fields (matching the
RAS hydraulics layer and the ESM coupling contract) but returns g cm⁻¹
s⁻¹ for the output (matching the existing SSM CGS bed state). This
asymmetry is deliberate: it minimizes the unit-conversion footprint at
the boundary (one conversion per closure call) and keeps the bed-state
math in CGS.

### 5.3 Backwards compatibility

The default value of `bedload_transport_function` in the SSM
constructor is `"van_rijn"` (via
`contracts.DEFAULT_BEDLOAD_TRANSPORT_FUNCTION`). Existing
configurations and test fixtures continue to work without modification.

### 5.4 Adding a new closure

To extend the menu (e.g. Meyer-Peter-Müller):

1. Add the class to `bedload.py` with `name: str = "mpm"` and a
   `transport_rate(...)` method returning g cm⁻¹ s⁻¹.
2. Add an entry to `BEDLOAD_TRANSPORT_FUNCTIONS` in `bedload.py`.
3. Add `"mpm"` to `contracts.BEDLOAD_FUNCTIONS`.
4. Add a `TestMeyerPeterMuller` class to `tests/sediment/test_bedload.py`
   with at least one hand-derived reference value.

No changes are needed to `ssm.py` or `from_config` — the registry
lookup handles everything else.

---

## 6. Solver wiring (Stage-2)

### 6.1 Why the solvers must consume the closure

Stage-1 shipped a `BedloadTransportFunction` registry but the two
solvers continued to advect CBL using van Rijn `u_BL`. That left users
of, e.g., `transport_function: wilcock_crowe` with a config that
*selected* Wilcock-Crowe but *ran* van Rijn — the closure object was
bound on the SSM instance and ignored by the solver. Stage-2 closes
this gap.

### 6.2 How the closure is bound

`SSM._instantiate_drivers` instantiates the configured closure once
(`bedload_mod.get_transport_function(name)`) and passes the same
instance into the solver constructor:

```python
self._bedload_transport_function = bedload_mod.get_transport_function(
    self.bedload_transport_function_name
)
self._bedload_solver = bedload_mod.BedloadStandaloneExplicit(
    self.registry_classes,
    bedload_cutoff_um=...,
    transport_function=self._bedload_transport_function,
)
```

The solver stores the closure on `self.transport_function` and consumes
it on every `step` call. The default value (omitting the kwarg) is
`VanRijn1984TransportFunction()`, so existing callers — including every
unit-test fixture — keep working bit-for-bit.

### 6.3 The `q_b → u_eff` adapter

Each closure returns a per-cell, per-class **bedload transport rate**
`q_b` in g cm⁻¹ s⁻¹ that already incorporates whatever combination of
velocity, height, concentration, or stream power the underlying formula
uses. The solvers, however, advect a **mass field** (CBL, g/cm²) on the
mesh using an upwind face-flux step that needs an **effective bedload
velocity** `u_eff` (cm/s). The adapter is

$$u_{\rm eff} = \frac{q_b}{\delta_{BL} \cdot C_{eq}}$$

where `δ_BL` and `C_eq` are van Rijn's saltation-layer height and
equilibrium concentration. Implementation:
`bedload._qb_to_effective_velocity_cm_s`. This is the **standard
practice in production codes** that mix transport functions with a
fixed bedload-layer model (Delft3D-MOR, MIKE 21). Two important
properties:

1. **Van Rijn identity.** When `q_b` itself comes from van Rijn,
   `q_b = u_BL · δ_BL · C_eq` algebraically, so `u_eff ≡ u_BL` to
   float64 round-off — this is why
   `test_explicit_van_rijn_matches_default` passes at `rtol=1e-12` and
   the prior `test_one_step_conserves_mass_in_closed_domain` is
   unchanged.
2. **Mass conservation.** The upwind face-flux step is unchanged
   (same edge / face geometry, same gather/scatter pattern); only the
   per-cell magnitude of `u_eff` changes. Closed-domain conservation
   on a 3-cell mesh is verified for Wilcock-Crowe in
   `test_wilcock_crowe_closed_domain_mass_conservation`.

The Riverine-constituent solver writes the per-edge mean of `u_eff`
(in m/s) into `mesh[advection_coef_var_name(class.label)]` — the same
field Stage-1 wrote, just sourced from the configured closure. Its
`long_name` attribute now records the closure name
("effective bedload velocity (edge mean, wilcock_crowe) — gravel_5mm")
so downstream inspection can identify the closure that produced it.

### 6.4 Surface-composition wiring

`SSM.run` builds a registry-wide `registry_context` dict before each
solver `step` and forwards it through the Protocol's new
`registry_context` argument. The dict carries everything the seven
shipped closures might need from the bed surface:

| Key                          | Shape           | Consumed by                |
|------------------------------|-----------------|----------------------------|
| `surface_class_fraction`     | `(nface, ssm_class)` | Wilcock-Crowe, Wu      |
| `surface_sand_fraction`      | `(nface,)`      | Wilcock-Crowe (`F_s`)      |
| `surface_geometric_mean_um`  | `(nface,)`      | Wilcock-Crowe (`d_sg`)     |
| `pe_ph_ratio`                | scalar (1.0)    | Wu (uniform-bed default)   |
| `registry`                   | `SedimentClassRegistry` | (reserved for future use) |

The orchestrator computes these from the **post-erosion-and-deposition**
surface fractions (i.e. consistent with the rest of the step's bed
update). The sand fraction `F_s` follows the Wentworth classification
(D50 ∈ [62.5, 2000] μm). The geometric mean `d_sg` uses the standard
mass-weighted log-mean (`exp(Σ F_i · log d_i)`).

The closures see a **single-class** view of this dict via
`bedload._class_context`, which slices the per-class fields down to
the requested class index (averaging across cells to a scalar, since
the existing W-C / Wu signatures take a scalar `F_i`). Closures that
don't need surface composition (van Rijn, Engelund-Hansen, Toffaleti,
Yang) ignore the dict entirely.

When `registry_context` is `None` (e.g. unit tests that call
`solver.step` directly), each closure falls back to its
documented default — Wilcock-Crowe to `F_s = 0.15` (typical gravel
bed), Wu to `p_e/p_h = 1` (uniform bed). This preserves the Stage-1
behaviour for tests that pre-date the wiring.

### 6.5 Parity tests

`test_wilcock_crowe_yields_matching_advection_coefficient` asserts that
the Riverine-constituent solver's per-edge advection coefficient field
equals the standalone solver's per-edge mean of `u_eff` (within float32
storage tolerance), under Wilcock-Crowe. This validates the
equivalence between the two solver paths under any chosen function: as
long as the two solvers agree at the cell-centred `u_eff`, their
respective advection steps will agree (the standalone path runs an
explicit upwind step on the mesh, the Riverine path delegates an
implicit step to Riverine's linear solver, both reading the same
velocity field).

The previously-skipped `test_standalone_vs_riverine_parity_on_1d_channel`
remains skipped because it requires Batch-C
(`linalg.LHS.update_values`) before Riverine can actually run the
implicit step — when that lands, the parity test should run with
`transport_function=WilcockCrowe2003()` (or any other) and the standalone
and Riverine CBL fields should agree to within numerical tolerance.

---

## 7. Limitations and future work

* **Toffaleti single-zone reduction.** The current implementation pins
  `M = 0.6, n_v = 1.5` to medium sand at 60 °F. A full four-zone
  integration is queued for phase 3 once SSM has a coupled vertical
  concentration profile.
* **Yang as bedload proxy.** Yang is a *total*-load formula; SSM treats
  the entire `q_t` as bedload. A cleaner split using van Rijn's PSUS
  (suspended-fraction) factor is queued.
* **Per-cell vs. scalar surface fraction.** The Wilcock-Crowe and Wu
  closures currently consume a *scalar* `F_i` (averaged across cells
  in `_class_context`). A future enhancement will broadcast the full
  `(nface,)` field through the closures so spatial heterogeneity in
  surface composition propagates into `q_b` per cell. (For uniform
  fixtures the scalar reduction is exact, so the Stage-2 tests are
  unaffected.)
* **`p_e/p_h` is pinned to 1.0.** Wu's hiding/exposure correction is
  disabled by default. A future enhancement will derive
  `p_e/p_h` from the full surface size-distribution per Wu, Wang & Jia
  (2000) eqs. 3-4.
* **No Meyer-Peter-Müller, Soulsby-Van Rijn, Garcia-Parker, or
  Einstein-Brown.** These are obvious next additions but were out of
  scope for the Stage-1 menu (which prioritized formulas that
  competitors uniquely offer).

---

## 8. Tests

`tests/sediment/test_bedload.py` ships 36 tests (plus 1 skipped parity
test queued for Batch C):

* `TestTransportFunctionRegistry` — validates the seven-function
  registry, factory lookup, and a smoke-loop over all closures.
* `TestVanRijnClassWrapper` — verifies the wrapper agrees with the
  standalone helpers to float64 round-off.
* `TestWilcockCrowe2003`, `TestParker1990`, `TestYang`, `TestWu2000`,
  `TestEngelundHansen1967`, `TestToffaleti1968` — per-formula tests
  with at least one hand-derived reference value (annotated with the
  derivation in the docstring) and threshold-behaviour checks.
* `TestStandaloneSolverWithCustomTransportFunction` (Stage-2) — solver
  consumes the `transport_function=` argument; Wilcock-Crowe and Wu
  smoke tests validate that the per-cell `q_b` matches the standalone
  unit-test value for the same conditions; closed-domain mass
  conservation holds for non-van-Rijn closures too.
* `TestStandaloneRiverineParityNonVanRijn` (Stage-2) — both solvers
  agree on the per-edge effective velocity under Wilcock-Crowe.

The smoke loop in `test_menu_smoke_all_functions_return_nonnegative`
guards against regressions when new closures are added.

---

## 9. References

* Engelund, F., and Hansen, E. (1967). *A monograph on sediment
  transport in alluvial streams.* Teknisk Forlag, Copenhagen.
* Parker, G. (1990). "Surface-based bedload transport relation for
  gravel rivers." *J. Hydraul. Res.* 28(4), 417–436.
  DOI: 10.1080/00221689009499058.
* Toffaleti, F. B. (1968). "A procedure for computation of the total
  river sand discharge and detailed distribution, bed to surface."
  Tech. Report No. 5, Committee on Channel Stabilization,
  U.S. Army Corps of Engineers.
* van Rijn, L. C. (1984). "Sediment transport, Part I: Bed load
  transport." *J. Hydraul. Eng.* 110(10), 1431–1456.
* Wilcock, P. R., and Crowe, J. C. (2003). "Surface-based transport
  model for mixed-size sediment." *J. Hydraul. Eng.* 129(2), 120–128.
  DOI: 10.1061/(ASCE)0733-9429(2003)129:2(120).
* Wu, W., Wang, S. S. Y., and Jia, Y. (2000). "Nonuniform sediment
  transport in alluvial rivers." *J. Hydraul. Res.* 38(6), 427–434.
  DOI: 10.1080/00221680009498296.
* Yang, C. T. (1973). "Incipient motion and sediment transport."
  *J. Hydraul. Div. ASCE* 99(HY10), 1679–1704.
* Yang, C. T. (1979). "Unit stream power equations for total load."
  *J. Hydrology* 40(1–2), 123–138.
  DOI: 10.1016/0022-1694(79)90092-1.
