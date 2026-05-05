# SSM — Sediment Simulation Module: Design Specification

**Status:** Draft for review
**Author:** Generated 2026-05-02 from EFDC SEDZLJ source review (commit `dsi-llc/EFDCPlus_Stable@main`, 2026-04-02), SAND2008-5621 (Thanh, Grace & James 2008), and architectural mapping of ClearWater-Riverine-streaming, ClearWater-modules-streaming (v2), and ClearWater-modules-phase2-ESM-streaming.
**Target:** New module `clearwater_modules_v2.processes.sediment` (SSM) implementing EFDC SEDZLJ's multi-class, multi-layer, SEDflume-driven cohesive + non-cohesive sediment transport, decoupled from HEC-RAS hydrodynamics, with a documented coupling contract to ESM.

---

## 1. Executive summary

SSM ports the EFDC SEDZLJ algorithm (Jones & Lick 2001; Ziegler 2002; James, Jones, Roberts & Hayter — implemented in `dsi-llc/EFDCPlus_Stable/EFDC/SedTran-SEDZLJ/`) into the ClearWater Python stack as a v2-pattern process. It introduces:

1. **N suspended sediment classes** as transported constituents in ClearWater-Riverine (each is a scalar PDE solved by the existing implicit-upwind machinery), with settling and bed-exchange terms supplied by SSM as source/sink arrays.
2. **A multi-layer bed state** owned by SSM and stored on the mesh dataset as additional `(nface, n_layers, n_class)` DataArrays — never advected, only updated by SSM each sediment time step.
3. **A shear-stress driver** that can either compute τ\_b internally from RAS face velocities (Parker 2004 log-law + optional Christoffersen–Jonsson wave–current combination) or accept τ\_b as an external input field.
4. **Bedload transport** on the mesh as a thin saltation-layer mass balance (van Rijn 1984), advected separately from suspended load.
5. **A SEDflume input loader** that consumes the original SEDZLJ `bed.sdf`, `erate.sdf`, and `core_field.sdf` files unchanged, plus a CSV/YAML alternative for new datasets.
6. **An ESM coupling contract** that (a) supplies ESM with `bed_change`, `bed_elevation`, `sediment_concentration`, and class composition (`d50`, % sand/silt/clay) for burial mortality, scour mortality, light extinction, and habitat suitability; and (b) accepts vegetation feedbacks (composite Manning's n, biostabilization factor, root cohesion) to modulate critical shear stress.

The module follows the `clearwater_modules_v2` Process pattern (`base.py:14`, `temperature.py:15`) — pure xarray vectorization, factory-registered, YAML-configurable, no Numba/JIT in the first cut. Performance optimization is a phase-2 concern.

---

## 2. Scope, goals, and non-goals

### 2.1 Goals

| # | Goal |
|---|---|
| G1 | Faithfully reproduce SEDZLJ erosion/deposition/bed-armoring physics in Python on an unstructured RAS mesh |
| G2 | Preserve `bed.sdf`/`erate.sdf`/`core_field.sdf` input formats so existing SEDflume datasets are reusable verbatim |
| G3 | Decouple from HEC-RAS hydrodynamics: τ\_b is computed from RAS face-velocity output OR supplied externally; no FORTRAN linkage to HEC sediment solver |
| G4 | Provide a single, documented xarray contract between SSM and ESM that supports vegetation–sediment two-way feedback |
| G5 | Conform to the `clearwater_modules_v2` Process pattern so SSM composes cleanly with TSM v2, NSM v2, ESM in `Model.run()` |
| G6 | Government-public-domain compatible: re-implement from algorithms (SAND2008-5621 + Jones & Lick 2001 + van Rijn 1984 + Cheng 1997 + Soulsby 1997 + Christoffersen & Jonsson 1985) rather than direct line-by-line GPL port. Treat the EFDC source as a reference implementation, not a derivation source. |

### 2.2 Non-goals (initial release)

| # | Non-goal | Rationale |
|---|---|---|
| NG1 | Cohesive-bed consolidation with **time-varying porosity / bulk density** (Gibson-class finite-strain self-weight consolidation) | A simpler age-dependent τ_ce model (Sanford & Maa 2001) is delivered in §5.10; full Gibson-style porosity evolution is deferred to a follow-on release |
| NG2 | 3-D vertical sediment profiles | ClearWater-Riverine is depth-averaged (KC=1 single layer in EFDC parlance); 3-D is out of scope |
| NG3 | Wave forcing from STWAVE/SWAN | Stub out the API; fetch and STWAVE wave-shear coupling can be added later behind the same `compute_shear()` interface |
| NG4 | Toxics partitioning | ChemFate is its own module; SSM will expose hooks but not own toxics state |
| NG5 | Propwash erosion | EFDC's propwash extension (NSEDS2 "fast" classes) is out of scope; the API leaves room for it |
| NG6 | Morphology feedback to RAS | One-way coupling only: SSM updates bed elevation locally, but does not edit RAS geometry. Same loose-coupling stance as ESM (see ESM `README.md:33–39`). |

### 2.3 Constraints (from user)

- **License path TBD.** Default plan: re-implement from published equations + SAND2008-5621 manual (clean-room style) and treat the GPL-2.0 EFDCPlus source as documentation, not as code to copy. Variable names and file formats are fair game; algorithmic structure must come from the published literature.
- **Bed-layer state must be designed carefully** (§7).
- **Shear-stress driver must support both RAS face-velocity computation and external import** (§8).
- **SEDflume input format must be preserved** (§9).

---

## 3. Where SSM lives in the ClearWater stack

```
┌─────────────────────────────────────────────────────────────────────┐
│ HEC-RAS 2D HDF5 (offline)                                          │
│   face_velocity, water_surface_elev, volume, mannings_n             │
└────────────────┬────────────────────────────────────────────────────┘
                 │ cwr.read_ras() (clearwater_riverine/io/hdf.py)
                 ▼
┌─────────────────────────────────────────────────────────────────────┐
│ ClearWater-Riverine (transport.py)                                  │
│   - xarray mesh Dataset (nface, nedge, time)                        │
│   - Implicit upwind A·c = b per constituent per step                │
│   - Streaming Zarr output via _flush/_release_to_stream             │
└────────────────┬────────────────────────────────────────────────────┘
                 │ constituents (suspended sediment classes)
                 │ + update_concentration source/sink hook
                 ▼
┌─────────────────────────────────────────────────────────────────────┐
│ ClearWater-modules-v2 Model loop                                    │
│   - Riverine process (wraps cwr)                                    │
│   - Temperature process (TSM)                                       │
│   - Nitrogen process (NSM)                                          │
│   - **SSM** (new)  ────────────┐                                    │
│   - ESM process (phase 2)      │                                    │
└────────────────────────────────┼────────────────────────────────────┘
                                 │
                                 ▼
                      ┌──────────────────────┐
                      │  SSM internals       │
                      │  ──────────────────  │
                      │  shear.py            │
                      │  bed.py              │
                      │  erosion.py          │
                      │  deposition.py       │
                      │  bedload.py          │
                      │  classes.py          │
                      │  io/sedflume.py      │
                      └──────────────────────┘
```

The package lives under `ClearWater-modules-streaming/src/clearwater_modules_v2/processes/sediment/`. The bed-state DataArrays live on the same mesh `xr.Dataset` that all v2 processes share, under a `ssm_*` namespace (§6).

---

## 4. Module layout

Mirroring `clearwater_modules_v2/processes/temperature.py` and `nutrients/`:

```
src/clearwater_modules_v2/processes/sediment/
├── __init__.py              # exports SSM
├── ssm.py                   # SSM(Process) — top-level driver class (analogous to Temperature)
├── classes.py               # SedimentClass dataclass; class registry
├── shear.py                 # bed shear from face velocities (Parker log-law); wave/current combo stub
├── settling.py              # Cheng (1997) settling velocities; class-level cache
├── deposition.py            # Gessler / Krone / van Rijn deposition probabilities
├── erosion.py               # SEDflume table interpolation (NSEDFLUME=1) + power-law (NSEDFLUME=2)
├── bed.py                   # multi-layer bed state container + active-layer reorganization
├── bedload.py               # van Rijn bedload velocity, height, mass-balance step
├── armoring.py              # D50_avg, TAUCRIT(D50_avg) interpolation, mass-fraction sorting
├── coupling.py              # ESM contract: vegetation → critical-shear, biostabilization → erosion-rate scaling
└── io/
    ├── sedflume.py          # bed.sdf / erate.sdf / core_field.sdf parsers
    ├── csv_loader.py        # alternative CSV/YAML format for new datasets
    └── hotstart.py          # SEDBED_HOT.SDF analog as NetCDF
```

Tests under `tests/sediment/` per the v2 convention (one file per submodule plus integration).

---

## 5. Governing equations and algorithms

All references are to SAND2008-5621 unless noted. EFDC source citations are given as `s_<file>.f90:<line>` for verification; **no code is copied**, only the published algorithm is reproduced.

### 5.1 Sediment classes

A class `s ∈ {1..N_s}` is defined by:

| Symbol | Variable | Units | Source |
|---|---|---|---|
| D₅₀,s | `d50` | μm | user |
| τ\_ce,s | `tau_ce` | dynes/cm² | user (Soulsby 1997 if blank) |
| τ\_cs,s | `tau_cs` | dynes/cm² | user (van Rijn 1984 eqs 8–9 if blank) |
| w\_s,s | `settling_velocity` | cm/s | user, or computed from D₅₀ via Cheng (1997) |
| ρ\_s,s | `solid_density` | g/cm³ | typically 2.65 |

**Cheng (1997) settling velocity** (SAND2008-5621 eq. 3, also `s_sedic.f90:444`):

$$d_*  = D_{50}\left[\frac{(s_s-1)g}{\nu^2}\right]^{1/3}, \qquad
w_s = \frac{\nu}{D_{50}}\left(\sqrt{25 + 1.2\,d_*^{2}} - 5\right)^{1.5}$$

with ν = 0.01 cm²/s, g = 980 cm/s², s\_s = ρ\_s/ρ\_w.

### 5.2 Bed shear stress τ\_b (`shear.py`)

Three modes selected at runtime via YAML:

**Mode A — `external`**: τ\_b is read directly from the mesh dataset variable `bed_shear_stress` (Pa, time-varying per cell). For users who already have τ from a coupled wave model or measurements.

**Mode B — `current_only`** (default for HEC-RAS 2D users, equivalent to EFDC `ISWAVE=0`):

$$\tau_b = \rho_w\, f_c\, |U|^2,
\qquad f_c = \left(\frac{0.42}{\ln\!\bigl(11\,h/(2\,k_n)\bigr)}\right)^{2}$$

(Parker 2004 log-law; `s_shear.f90:261`). Velocity magnitude |U| comes from RAS face velocities averaged to cell centroids using the existing `edge_to_face` map in `clearwater_riverine.utilities`. Roughness `k_n = max(D50_avg, ZBSKIN)` in meters; `D50_avg` is the bed-surface mean grain diameter (§5.6).

**Mode C — `wave_current`** (Christoffersen & Jonsson 1985; `s_shear.f90:285–306`): combined wave–current friction with iterative solution for σ\_w, M\_w, J\_w, δ\_w, k'. Wave parameters (orbital velocity U\_δ, frequency ω, direction θ\_w) are read from the mesh dataset variables `wave_orbital_velocity`, `wave_frequency`, `wave_direction`. Stub in v1; full implementation in phase 2.

A 10 %-per-step growth limiter is applied to τ to mirror SEDZLJ's stability device (`s_shear.f90:315`).

### 5.3 Critical shear stresses

- **τ\_ce (erosion threshold)** — class-level constant, supplied or computed via Soulsby (1997).
- **τ\_cs (suspension threshold)** — class-level, supplied or via van Rijn (1984) eqs 8–9.
- **τ\_crit(D₅₀\_avg)** for the active/deposited layer — interpolated from the SEDflume table `(SCND, TAUCRITE)` using the layer's mass-weighted mean D₅₀ (`s_sedzlj.f90:265`).
- **τ\_crit(layer K)** for in-place layers — read directly from `bed.sdf` per layer per core (`s_sedic.f90:254`).

### 5.4 Erosion rate (`erosion.py`)

Two formulations, switched by `nsedflume` config flag:

**`nsedflume=1` (table interpolation)** — bilinear in (τ, layer fractional mass remaining) of the SEDflume `ERATE(K, L, M)` table, log-space in the depth dimension and linear in τ (`s_sedzlj.f90:498` and `:535`):

$$E_{\text{rate}}(K,\tau) = \mathrm{interp}_{\log}\bigl[\,\tau \in (\tau_{\!\downarrow},\tau_{\!\uparrow}),\; m_K/m_{K,0} \in (0,1)\,\bigr]$$

returning cm/s, then multiplied by the layer's dry bulk density to give g/cm²/s.

**`nsedflume=2` (power law)** — per core, per layer:

$$E_{\text{rate}} = A\,(\tau\,[\text{Pa}])^{n}, \qquad E_{\text{rate}} \le E_{\max}$$

with (A, n, E\_max) per layer (`s_sedzlj.f90:508`, `:539`).

Per-class erosion mass over Δt is `ELAY(s) = PERSED(s,K) × ERATEMOD × Δt`, gated by `τ ≥ τ_ce(s)` (`s_sedzlj.f90:572`). The class-by-class gate is what produces emergent armoring.

### 5.5 Deposition rate (`deposition.py`)

Per class:

$$D_s = P_{\!d,s}\, C_{b,s}\, w_{s,s}\, \Delta t$$

with bottom-cell concentration `C_b` estimated from a single-layer exponential profile if needed (`s_sedzlj.f90:113`). Probability `P_d`:

- **Sand-class (D₅₀ ≥ `bedload_cutoff`, default 64 μm)** — Gessler (1965) erfc form (`s_sedzlj.f90:139`):
  $$P_y = \frac{1}{0.57}\!\left(\frac{\tau_{cs}}{\tau} - 1\right), \quad
    P_d = \tfrac{1}{2}\,\mathrm{erfc}(-P_y/\sqrt{2})$$
  implemented via the Abramowitz & Stegun rational approximation already used in `s_sedzlj.f90:141–148`.
- **Mud-class (D₅₀ < 64 μm)** — Krone (`s_sedzlj.f90:151`): $P_d = \max(1 - \tau/\tau_{cs},\, 0)$.

Each cell's deposition is capped at the available mass in the bottom water layer (`MAXDEPLIMIT × C × h`).

### 5.6 Active-layer reorganization (`bed.py`, `armoring.py`)

Required active-layer thickness (Lick 2008; `s_sedzlj.f90:273`):

$$T_{\text{act}} = T_{\text{actm}}\, D_{50,\text{avg}}\,
  \max\!\left(1,\frac{\tau}{\tau_{\text{crit}}}\right)\,\frac{\rho_b}{10000}$$

(units: g/cm² when D₅₀ in μm, ρ\_b in g/cm³). Default `T_actm = 2`.

Reorganization rules (executed before erosion each step):

1. If `m_1 > T_act` (net deposition that step), excess is pushed to layer 2 (deposition layer), preserving mass-weighted PERSED.
2. If `m_1 < T_act` and τ > τ\_crit(SLLN), borrow from the next-non-empty layer SLLN to top up to T\_act.
3. If borrowed mass insufficient, promote the next layer up and zero the empty one.

D₅₀\_avg is the mass-weighted mean of class D₅₀s. Armoring emerges automatically because erosion fractionates per class but borrowing from below preserves that layer's composition.

### 5.7 Bedload transport (`bedload.py`)

Only for classes with D₅₀ ≥ `bedload_cutoff`. State variable: `cbl[s]` = bedload mass per unit bed area (g/cm²) on each cell.

Transport parameter (van Rijn 1984; `s_bedload.f90:147`):

$$T_R = \max\!\left(\frac{\tau - \tau_{ce}}{\tau_{ce}},\, 0\right)$$

Bedload velocity (eq. 20a):

$$u_{BL} = 1.5\, T_R^{0.6}\, \sqrt{(s_s-1)\,g\,D_{50}}\quad[\mathrm{cm/s}]$$

Bedload (saltation) layer height (eq. 20b):

$$\delta_{BL} = 0.3\, D_{50}\, d_*^{0.7}\,\sqrt{T_R}\quad[\mathrm{cm}]$$

Suspended fraction (PSUS) — log-interpolation between `√τ_cs/w_s` and shear-velocity ratio of 4 (`s_bedload.f90:126`).

Equilibrium concentration for deposition probability (van Rijn 1981 eq. 21; `s_sedzlj.f90:191`):

$$C_{eq,s} = 0.117\,\rho_s\,\frac{T_R}{d_*}$$

The bedload mass balance is a separate explicit upwind step on the mesh (face-by-face flux), updated each sediment time step. **Implementation note**: this is mathematically a small first-order PDE on (nface,) per class. The cleanest reuse of ClearWater-Riverine machinery is to register `cbl[s]` as additional constituents with a custom advection coefficient = `u_BL` instead of `edge_velocity`. To be evaluated in §11.

**Pluggable transport-rate closure (Stage-1 menu, Stage-2 solver wiring).** The per-cell, per-class bedload transport rate `q_b` is now selected by name from a registry of seven peer-reviewed closures (`bedload.BEDLOAD_TRANSPORT_FUNCTIONS`), and the standalone / Riverine-constituent solvers consume the configured closure on every step (no longer hard-wired to van Rijn `u_BL`). See `ssm_bedload_functions.md` §6 ("Solver wiring") for the `q_b → u_eff` derivation, the surface-composition wiring (`registry_context`), and the parity tests:

| YAML name          | Citation                          | Best for                   |
|--------------------|-----------------------------------|----------------------------|
| `van_rijn`         | van Rijn 1984a (default)          | Sand                       |
| `wilcock_crowe`    | Wilcock & Crowe 2003              | Sand-gravel mixed beds     |
| `parker`           | Parker 1990                       | Gravel rivers              |
| `yang`             | Yang 1973, 1979                   | Sand (total load)          |
| `wu`               | Wu, Wang & Jia 2000               | Non-uniform sediment       |
| `engelund_hansen`  | Engelund & Hansen 1967            | Sand-bed total load        |
| `toffaleti`        | Toffaleti 1968 (single-zone)      | USACE BR-1 sand            |

Each closure conforms to a `BedloadTransportFunction` Protocol and returns `q_b` in g cm⁻¹ s⁻¹. Selection is via YAML:

```yaml
sediment:
  bedload:
    solver: standalone               # or riverine, off
    transport_function: van_rijn     # or wilcock_crowe, parker, yang,
                                     # wu, engelund_hansen, toffaleti
```

The default (`van_rijn`) preserves backwards compatibility with the original SEDZLJ-port behaviour. See `ssm_bedload_functions.md` for the full design memo, per-formula domains of applicability, and selection guidance.

### 5.8 Bed-elevation update

After each sediment step, per cell:

$$h_{bed}(L) = \sum_{K=1}^{KB} \frac{m_K}{\rho_{b,K}} \times 0.01 \quad[\mathrm{m}]$$

`bed_change` is the per-step delta; `cumulative_bed_change` is the running sum. Both are written back to the mesh dataset for ESM consumption (§6.3).

### 5.10 Cohesive-bed consolidation (`consolidation.py`)

SEDZLJ holds the per-layer critical shear stress τ_ce constant once initialised (`s_sedzlj.f90:707` "SEDZLJ DOES NOT HAVE CONSOLIDATION"). For freshly-deposited cohesive sediment this under-predicts the strength gain that occurs over hours-to-weeks of self-weight consolidation, biasing erosion fluxes high relative to MIKE 21, Delft3D, and TELEMAC, which all carry some form of consolidation model.

SSM closes this gap with an **opt-in** Sanford & Maa (2001) single-mode age-dependent τ_ce:

$$\tau_{ce}^{\rm eff}(t_{\rm age}) = \tau_{ce,\infty} - (\tau_{ce,\infty} - \tau_{ce,0})\,\exp(-t_{\rm age}/T_c)$$

where:

* $\tau_{ce,0}$ — freshly-deposited critical shear stress (lower bound; default 0.10 Pa)
* $\tau_{ce,\infty}$ — fully-consolidated critical shear stress (upper bound; default 0.50 Pa, 5× the lower bound per Sanford & Maa 2001 Fig. 4)
* $T_c$ — consolidation time scale (default 7 days, the typical e-folding time observed in flume experiments on estuarine mud)
* $t_{\rm age}$ — layer's mass-weighted mean age, advanced by `dt` each step in `update_bed_elevation` for every layer that holds mass

Consolidation is applied **only to cohesive classes** (D₅₀ < `bedload_cutoff`, default 64 μm); non-cohesive (sand) classes retain the static τ_ce.

**Per-layer age tracking** lives on the mesh under `ssm_bed_layer_age` (s; (time, nface, ssm_layer)). Three propagation rules govern age evolution:

1. **Time advancement:** every step, `update_bed_elevation` adds `dt` to every layer with non-zero mass. Empty layers are pinned to age 0.
2. **Age dilution on deposition:** when fresh mass Δm enters layer 1, the new layer-mean age is $t_{1,\rm new} = t_1 \cdot m_1 / (m_1 + \Delta m)$ (the deposit enters with age 0).
3. **Age inheritance on borrow / promote / collapse:** mass-weighted blending of the donor and recipient layers' ages, mirroring the existing PERSED blend logic in `reorganize_active_layer`.

**Configuration (opt-in):**

```yaml
sediment:
  consolidation:
    enabled: true                  # default false (SEDZLJ-equivalent behaviour)
    model: sanford_maa
    tau_ce_zero_pa: 0.10
    tau_ce_inf_pa: 0.50
    consolidation_time_s: 604800   # 7 days
```

**Limitations.** This first release scopes consolidation to the τ_ce(age) formulation only. A complete model would additionally evolve time-varying porosity / bulk density (Gibson, England & Hussey 1967), gel-point dynamics (Toorman 1999), and finite-strain self-weight consolidation. Those are deferred to a follow-on release. See `design/ssm_consolidation.md` for the full design memo, calibration discussion, and limitations relative to a Gibson-class model.

References:
* Sanford, L. P., and Maa, J. P.-Y. (2001). "A unified erosion formulation for fine sediments." *Marine Geology* 179(1–2), 9–23. DOI: 10.1016/S0025-3227(01)00201-8.
* Mehta, A. J., and Partheniades, E. (1975). "An investigation of the depositional properties of flocculated fine sediments." *J. Hydraul. Res.* 13(4), 361–381. DOI: 10.1080/00221687509499694.

---

## 6. Data contracts

### 6.1 Inputs SSM reads from the mesh dataset

| Variable | Dims | Units | Source | Required |
|---|---|---|---|---|
| `volume` | (time, nface) | m³ | RAS HDF5 | yes |
| `water_surface_elev` | (time, nface) | m | RAS HDF5 | yes |
| `face_hydraulic_depth` | (time, nface) | m | RAS HDF5 or derived | yes |
| `edge_velocity` | (time, nedge) | m/s | RAS HDF5 | yes (Mode B) |
| `mannings_n` | (nface,) | s/m^(1/3) | RAS or ESM | yes (Mode B) |
| `bed_shear_stress` | (time, nface) | Pa | external | only Mode A |
| `wave_orbital_velocity`, `wave_frequency`, `wave_direction` | (time, nface) | m/s, rad/s, rad | wave model | only Mode C |
| `vegetation_biostabilization` | (time, nface) | dimensionless [0,1] | ESM (optional) | optional |
| `vegetation_root_cohesion` | (time, nface) | Pa | ESM (optional) | optional |

`composite_manning_n` (from ESM) is read in place of static `mannings_n` when ESM is in the run, allowing vegetation roughness to feed back into shear computation.

### 6.2 Suspended-sediment constituents (transported by Riverine)

Each class is registered as a Riverine constituent with name `ssm_suspended_<label>` (e.g., `ssm_suspended_silt_fine`, `ssm_suspended_sand_medium`):

```python
constituent_dict = {
    f"ssm_suspended_{cls.label}": {
        "initial_conditions": ic_path,   # CSV: Cell_Index, Concentration
        "boundary_conditions": bc_path,  # CSV: RAS2D_TS_Name, Datetime, Concentration
        "units": "mg/L",
        "decay_rate": 0.0,               # not used; SSM owns sources/sinks
    }
    for cls in ssm.classes
}
```

SSM injects bed-exchange flux (erosion − deposition) each step via `ClearwaterRiverine.update(update_concentration={...})`.

### 6.3 Bed-state DataArrays SSM owns on the mesh dataset

| Variable | Dims | Units | Description |
|---|---|---|---|
| `ssm_bed_layer_mass` | (time, nface, n_layers) | g/cm² | per-layer dry mass (TSED) |
| `ssm_bed_layer_initial_mass` | (nface, n_layers) | g/cm² | TSED0, for SEDflume depth interp |
| `ssm_bed_class_fraction` | (time, nface, n_layers, n_class) | dimensionless | PERSED |
| `ssm_bed_layer_active` | (time, nface, n_layers) | int8 {0,1,2} | LAYERACTIVE |
| `ssm_bed_layer_taucrit` | (time, nface, n_layers) | Pa | TAUCOR |
| `ssm_bed_layer_bulk_density` | (nface, n_layers) | g/cm³ | BULKDENS (constant per SEDZLJ) |
| `ssm_bed_thickness` | (time, nface, n_layers) | m | HBED |
| `ssm_bed_layer_age` | (time, nface, n_layers) | s | per-layer mean age since deposition (consolidation, §5.10) |
| `ssm_bed_total_thickness` | (time, nface) | m | sum over layers |
| `ssm_bed_d50_surface` | (time, nface) | μm | mass-weighted mean D₅₀ of surface layer |
| `ssm_bed_elevation` | (time, nface) | m | absolute bed elevation |
| `ssm_bed_change` | (time, nface) | m/step | per-step delta |
| `ssm_bed_cumulative_change` | (time, nface) | m | running sum |
| `ssm_bedload_mass` | (time, nface, n_class) | g/cm² | CBL |
| `ssm_bed_shear_stress` | (time, nface) | Pa | computed/imported τ\_b |
| `ssm_bed_critical_shear_stress` | (time, nface) | Pa | τ\_crit(D₅₀\_avg) at surface |

Bed-state arrays use the streaming pattern from `transport.py:850–912` — flushed to Zarr every `streaming_interval` steps and replaced with NaN-padded buffers — so a 92-day, 587k-cell run does not blow up RAM. The NaN-sentinel + `released_time_range` attribute machinery already exists; SSM must register its arrays so `_release_to_stream` knows about them.

### 6.4 Outputs to ESM (consumed by `esm.io.clearwater_interface`)

ESM's existing variable names (from `esm/model.py:468–471`) are kept verbatim; SSM populates them as aliases of its own:

| ESM-side name | SSM source | Notes |
|---|---|---|
| `bed_change` | `ssm_bed_change` | unchanged |
| `cumulative_bed_change` | `ssm_bed_cumulative_change` | unchanged |
| `current_bed_elevation` | `ssm_bed_elevation` | unchanged |
| `sediment_concentration` | sum over `ssm_suspended_*` | total TSS, kg/m³ |
| `bed_d50_surface` (new) | `ssm_bed_d50_surface` | enables d₅₀-dependent habitat in ESM |
| `bed_shear_stress` (new) | `ssm_bed_shear_stress` | replaces velocity-only scour mortality |
| `bed_critical_shear_stress` (new) | `ssm_bed_critical_shear_stress` | for excess-shear scour formulation |

Add a new `ClearWaterInterface.get_sediment_state()` method to `esm/io/clearwater_interface.py` that returns a typed dataclass; existing ESM consumers continue to work via the legacy field names.

### 6.5 Inputs from ESM (vegetation feedback)

| Variable | Dims | Units | ESM source | SSM use |
|---|---|---|---|---|
| `composite_manning_n` | (time, nface) | s/m^(1/3) | `esm/processes/roughness/composite.py` | replaces static `mannings_n` in shear computation |
| `vegetation_biostabilization` | (time, nface) | [0,1] | new ESM output | multiplies τ\_ce: `tau_ce_eff = tau_ce × (1 + α·B)` with α calibrated |
| `vegetation_root_cohesion` | (time, nface) | Pa | new ESM output | additive to τ\_crit in deeper layers |
| `vegetation_frontal_area` | (time, nface) | m²/m² | `esm/io/clearwater_interface.py:411` (already declared, not yet wired) | optional skin-vs-form drag partitioning in shear |

These are all opt-in; if absent, SSM uses defaults and ESM coupling reverts to one-way.

---

## 7. Bed-layer state design

This is the most subtle piece and the user flagged it for careful design.

### 7.1 Layout

The bed at each cell is represented by **a fixed-depth stack of K\_B layers**, indexed top-down (`K=1` is the surface):

```
K=1  ┌─ Active layer (dynamic; sorting layer)         ─┐
K=2  ├─ Deposition layer (dynamic; new sediment)        │
K=3  ├─ In-place layer 1 (SEDflume core data)           │ K_B layers total
...  ├─ ...                                             │
K=KB ├─ In-place layer N (deepest core data)           ─┘
```

Three layer states (`LAYERACTIVE`):
- **0 = absent** — no mass; layer fully eroded
- **1 = active/deposited** — top two layers; properties evolve with sorting and deposition
- **2 = in-place** — original SEDflume core data; fixed bulk density, fixed per-layer τ\_crit, depth-interpolated erosion rate

`K_B` is a configuration constant per run (typically 5–10), set in the SEDflume `bed.sdf` file and validated against any hotstart.

### 7.2 Storage

Two design alternatives were considered:

**A. Wide xarray DataArray on the mesh dataset** (chosen):
```python
mesh["ssm_bed_class_fraction"]  # dims: (time, nface, n_layers, n_class)
mesh["ssm_bed_layer_mass"]      # dims: (time, nface, n_layers)
```
Pros: integrates with existing streaming/Zarr machinery; no parallel state object to keep in sync; ESM can read it directly.
Cons: nontrivial memory footprint; n\_class × n\_layers × nface × ntime can grow large.

**B. Separate `BedState` object owned by SSM with its own xarray Dataset**.
Pros: cleaner separation; can be checkpointed independently.
Cons: requires building a parallel streaming/release machinery.

→ **Decision: A**, with two mitigations:
1. The bed state has its own configurable streaming interval (default = 10× the suspended-sediment interval) since it changes more slowly.
2. `ssm_bed_class_fraction` and `ssm_bed_layer_mass` are written to Zarr as `int8` (PERSED scaled to 0–255) and `float32` respectively, halving on-disk size vs. float64.

### 7.3 Memory budget (worked example)

For a 587k-cell, 92-day run with `dt_sediment = 1 hr` and `n_layers = 8`, `n_class = 6`:

- Live in RAM (one timestep): `587k × 8 × 6 × 4 B = 113 MB` for class fractions; `587k × 8 × 4 B = 19 MB` for layer mass. Total bed live state ≈ **150 MB**.
- Suspended live state: `587k × 6 × 4 B = 14 MB` (one timestep) plus the t/t+1 buffer the implicit solver needs ≈ **30 MB**.
- Streaming writes flush the rest to Zarr and replace with NaN-padded buffers.

This fits comfortably within the existing streaming budget (the design memo `cw_riverine_streaming_in_memory_release.md` targets 5 GB total for the same case).

### 7.4 Active-layer reorganization order

Per cell per sediment step, in order:

1. **Compute τ\_b** (`shear.py`).
2. **Compute D₅₀\_avg of surface layer** → look up τ\_crit (active) or read τ\_crit (in-place).
3. **Compute T\_act** (§5.6 formula).
4. **Reorganize**: if τ\_b > τ\_crit, ensure layer 1 is full to T\_act by borrowing from layer 2 then SLLN; if depositional, push excess from layer 1 to layer 2.
5. **Compute deposition** D\_s for each class (Gessler/Krone gate).
6. **Compute erosion** E\_s for each class (per-class τ\_ce gate, per-layer mass cap).
7. **Update layer mass and PERSED** with mass conservation (`s_sedzlj.f90:621–660` algorithm; the `RATIOMASS` rescaling protects float32 precision).
8. **Update bed elevation** and write `bed_change`.
9. **Inject net flux** (E − D) per class into the suspended-sediment source/sink array for Riverine's next step.

This ordering is critical: D₅₀\_avg and τ\_crit must be evaluated **before** erosion strips fines, otherwise the active-layer reorganization can permit non-conservative mass loss when τ approaches τ\_crit.

---

## 8. Shear-stress driver design

Both options live behind a single interface:

```python
class ShearStressDriver(Protocol):
    def compute(self, mesh: xr.Dataset, t: int, ssm_state: SSMState) -> xr.DataArray:
        """Return τ_b (Pa) on (nface,) for time index t."""
```

Concrete implementations:

| Class | Mode | Inputs |
|---|---|---|
| `ExternalShearDriver` | A | `mesh["bed_shear_stress"]` |
| `CurrentOnlyShearDriver` | B (default) | RAS face velocity → cell centroid; Manning's n; depth; D₅₀\_avg from SSM |
| `WaveCurrentShearDriver` | C (phase 2) | Mode B inputs + wave parameters; runs Christoffersen–Jonsson 1985 iteration |

Mode B implementation notes:

- **Velocity reconstruction**: ClearWater-Riverine stores velocity on edges (`edge_velocity`, m/s). For shear we need cell-centroid magnitude. Reuse the existing edge-to-face averaging in `clearwater_riverine.utilities` (lines ~220–235) where the same operation is already done for diffusion-coefficient computation.
- **Roughness**: `k_n = max(D50_surface_meters, ZBSKIN_meters)`, where `ZBSKIN` defaults to 1.5 mm (typical SEDZLJ value, per `bed.sdf` example). This separates skin friction (used for sediment) from total form drag (Manning's n, used for hydrodynamics).
- **Manning's n vs k\_n**: Mode B uses k\_n directly; Manning's n is *not* used in the shear formula (Parker log-law gives f\_c from k\_n alone). However, when ESM supplies `composite_manning_n`, we offer an alternative formulation `f_c = g·n²/h^(1/3)` (selected by `shear.formulation: "manning"` in YAML) so users who trust their calibrated Manning's n can use it.
- **Stability limiter**: Apply the 10 %-per-step growth limiter from `s_shear.f90:315` to suppress shock-induced oscillations.

The driver is selected in YAML:

```yaml
processes:
  - sediment:
      shear:
        mode: current_only       # external | current_only | wave_current
        formulation: log_law     # log_law | manning  (Mode B only)
        zb_skin: 0.0015          # m, skin roughness fallback
        growth_limit: 0.10
```

---

## 9. SEDflume input file compatibility

### 9.1 Files preserved verbatim

`io/sedflume.py` parses the three native SEDZLJ files unchanged:

- **`bed.sdf`** — global parameters (KB, NSICM, ZBSKIN, TAUCONST, BEDLOAD\_CUTOFF, etc.), per-class properties (D₅₀, TCRE, TCRSUS, settling), per-interpolant properties (SCND, TAUCRITE), and either (a) ENRATE table for `nsedflume=1` or (b) (A, n, max\_rate) tuples for `nsedflume=2`. Format and card structure exactly per SAND2008-5621 Figure 3.
- **`erate.sdf`** — per core: per-layer thicknesses, critical shear stresses, bulk densities, water/sediment density, particle size distribution (mass percent per layer per class), and per-shear-level erosion rates. Format per SAND2008-5621 Figure 1.
- **`core_field.sdf`** — integer matrix mapping each cell (i, j) to a core ID. Both DSI standard (one row per cell with i, j, core) and SNL standard (raster matrix) formats are supported per `s_sedic.f90:151–180`.

A reference dataset (the SEDflume example from SAND2008-5621 §"BED.SDF") is included in `tests/data/sedflume_example/` and round-trip parsed in `test_io_sedflume.py`.

### 9.2 Mesh-mapping for unstructured meshes

The original `core_field.sdf` is structured-grid (i, j) addressed. For the unstructured RAS mesh we need a `core_id` per `nface`. Two paths:

1. **Spatial mapping file**: `core_field_unstructured.csv` with columns `Cell_Index, Core_ID`. SSM's `ConfigLoader` reads this when the mesh is unstructured.
2. **Polygon overlay**: a GIS polygon shapefile of core extents; cells are assigned by point-in-polygon at load time. Phase 2.

### 9.3 Optional CSV/YAML alternative

For new projects without legacy SEDZLJ inputs, a YAML-equivalent config:

```yaml
sediment_classes:
  - label: silt_fine
    d50_um: 32
    tau_ce_pa: 0.15
    tau_cs_pa: 0.20
    settling_cm_s: -1   # -1 = compute via Cheng (1997)
  - label: sand_medium
    d50_um: 250
    ...

bed_layers:
  n_layers: 8
  layer_thickness_cm: [0, 0, 5, 5, 5, 10, 20, 50]   # active and deposition start at 0
  bulk_density_g_cm3: [1.6, 1.6, 1.6, 1.7, 1.7, 1.8, 1.9, 1.9]
  initial_class_mass_fraction_pct:
    - layer: 3
      values: [10, 20, 30, 25, 10, 5]   # one per class

cores:
  - id: 1
    cells: [0, 1, 2, ...]   # nface indices
    erate_table:            # SEDflume data
      shear_levels_pa: [0, 0.2, 0.4, 0.8, 1.0, 2.0]
      erate_cm_s:             # per layer, per shear level
        - [0, 6.6e-5, 4.66e-4, 3.29e-3, 6.17e-3, 4.36e-2]
        - ...
```

`io/csv_loader.py` translates this to the same internal data structures as `sedflume.py`.

### 9.4 Hotstart / restart

`SEDBED_HOT.SDF` is replaced by NetCDF (`ssm_state_<timestamp>.nc`) consistent with the existing checkpoint pattern in `tests/test_checkpoint.py`. Backwards-compatibility shim: `io/hotstart.py` can read the legacy ASCII format and emit the NetCDF equivalent.

---

## 10. ESM coupling

### 10.1 Replacing ESM's bed-mortality inputs

Currently (`esm/model.py:1140–1180`), burial mortality reads `bed_change` and `bed_elevation` directly from RAS HDF5 via `HECRASSediment` (`hecras_reader.py:46–58`). With SSM in the run, those same variables are populated by SSM each ESM time step instead. The `ClearWaterInterface` already auto-detects fields; add a sediment branch.

### 10.2 Upgrading scour mortality

ESM's `scour_death_type` (constants.py:361–377) currently exposes velocity-based modes (1, 3, 5) and erosion-based modes (2, 3, 5). Adding mode 6 = **excess-shear**:

```python
if scour_death_type == 6:
    excess_shear = max(0, ssm_bed_shear_stress - ssm_bed_critical_shear_stress)
    scour_rate = k_scour × excess_shear / tau_ref
```

This becomes the recommended mode when SSM is present.

### 10.3 Light extinction from SSM-resolved TSS

ESM's `compute_light_extinction_from_sediment` (`light.py:131–150`) takes a single TSS scalar. With SSM:

```python
TSS = sum(mesh[f"ssm_suspended_{cls.label}"] for cls in ssm.classes)
k_d = k_0 + sum(a_class[cls] × mesh[f"ssm_suspended_{cls.label}"] for cls in ssm.classes)
```

per-class extinction coefficients `a_class` (typically larger for fines) give more physical light attenuation than a single constant.

### 10.4 Vegetation-biostabilization feedback (new ESM output)

ESM does not currently compute biostabilization. Add to `esm/processes/`:

```python
# esm/processes/biostabilization.py
def compute_biostabilization(
    biomass: xr.DataArray,        # g/m²
    root_density: xr.DataArray,   # g/m²
    species: xr.DataArray,
    params: dict,
) -> xr.DataArray:
    """Return biostabilization factor B ∈ [0,1] per cell."""
    # Scaled e.g. by Le Hir et al. (2007) or empirical species lookup
    ...
```

SSM consumes via `tau_ce_eff = tau_ce × (1 + α·B)`. The α coefficient is calibration-grade; default 0.5 (50 % τ\_ce increase at full vegetation cover).

### 10.5 Coupling cadence

Per ESM `README.md:40–49`, the orchestrator pattern is: hourly transport/kinetics → daily ESM update. SSM follows the transport cadence (typically `dt_sediment = dt_transport`), and exposes daily-aggregated bed-state fields to ESM via `streaming_helpers.py`. ESM-to-SSM vegetation-feedback fields update at ESM's cadence and SSM picks them up on the next sediment step (within-day step-function constancy).

---

## 11. Performance considerations

The first cut prioritizes correctness and clarity (pure xarray, like TSM v2). Optimization guidance for phase 2:

1. **Hot loop is `bed.py:reorganize_active_layer`** — runs per cell per step. The Fortran does this in a per-cell loop with branching (`s_sedzlj.f90:300–342`). Vectorizing across cells with `np.where` masks for each branch is feasible; the layer-by-layer recursion is small (≤ KB iterations).
2. **Erosion-rate interpolation** is a 2-D bilinear lookup per (cell, class). This is straightforwardly vectorizable; precompute the τ-bracket indices in a single sort per step.
3. **Bedload advection** — registering `cbl[s]` as a Riverine constituent reuses the sparse implicit solver, but the bedload velocity is *class-specific and shear-driven*, not the same as the suspended-sediment advection coefficient. This requires a per-constituent advection coefficient slot in `linalg.py` — currently there's only one. Two options:
   - **Phase 1**: implement bedload as an explicit upwind step in `bedload.py` (a few lines of NumPy per face), not as a Riverine constituent. Simpler, no Riverine modification needed.
   - **Phase 2**: extend Riverine to accept per-constituent advection coefficients. Bigger lift but better consistency.
   → **Recommend phase 1** for the initial release.
4. **Numba/JAX**: defer until profiling shows a real hotspot. The xarray-pure approach in TSM v2 is fast enough on the 587k-cell Albany case.
5. **Streaming**: bed state has its own `streaming_interval`, default 10× longer than suspended (bed evolves slowly).

---

## 12. Testing strategy

| Tier | What | Reference |
|---|---|---|
| Unit | Per-equation: Cheng settling, Gessler probability, Krone probability, van Rijn bedload velocity, Soulsby τ\_ce, log-law f\_c | `tests/sediment/test_*.py` |
| Unit | SEDflume erate.sdf round-trip parse | `tests/sediment/test_io_sedflume.py` |
| Unit | Active-layer reorganization on hand-built 3-layer cells | `tests/sediment/test_bed.py` |
| Integration | Mass conservation: total mass (suspended + bed + bedload) constant within float64 epsilon over N steps in a closed domain | `tests/sediment/test_conservation.py` |
| Integration | Armoring emergence: 6-class bed under steady τ should preferentially erode fines and coarsen surface D₅₀ over time | `tests/sediment/test_armoring.py` |
| Integration | Idealized channel: uniform-slope flume, single sand class, equilibrium bedload concentration matches van Rijn equilibrium analytically | `tests/sediment/test_van_rijn_equilibrium.py` |
| Integration | Reference SEDflume dataset from SAND2008-5621 Figure 3 → compare per-layer erosion rates to published values | `tests/sediment/test_sand2008_reference.py` |
| Regression | Replicate one of the published SEDZLJ application cases (e.g., the Lower Duwamish or Kalamazoo River cases) at low resolution and compare bed-evolution time series to EFDC+ output | `tests/sediment/test_efdc_comparison.py` |
| End-to-end | SSM + Riverine + TSM v2 + ESM in a coupled run on the existing Albany test mesh (587k cells, 92 days) | `tests/integration/test_ssm_esm_albany.py` |

The EFDC-comparison regression test is the **primary correctness gate**: build a small EFDC+ test case, run both EFDC+ and SSM with the same `bed.sdf`/`erate.sdf`/`core_field.sdf`, and compare bed-state and TSS time series at sentinel cells. This is the closest thing to a published "answer key" for SEDZLJ.

---

## 13. Phased implementation plan

**Phase 1 — single-class, single-layer, current-only (4–6 weeks)**
- `classes.py`, `settling.py` (Cheng), `shear.py` Mode B current-only
- Single-class suspended sediment as a Riverine constituent
- Single-layer bed (active layer only) — no reorganization, no in-place layers
- Erosion via power law (`nsedflume=2`)
- Krone deposition only (mud)
- Mass conservation test passing
- Goal: end-to-end pipeline working, validates the data contracts

**Phase 2 — multi-class, multi-layer, SEDflume input (6–8 weeks)**
- `bed.py` with full active + deposition + in-place layer machinery
- `armoring.py` with D₅₀\_avg and τ\_crit interpolation
- `erosion.py` with table interpolation (`nsedflume=1`)
- Gessler deposition (sand)
- `io/sedflume.py` round-trip
- SAND2008-5621 reference dataset test passing
- ESM output contract wired (bed\_change, bed\_elevation, sediment\_concentration via `ClearWaterInterface`)

**Phase 3 — bedload + ESM feedback (4–6 weeks)**
- `bedload.py` van Rijn velocity, height, mass balance
- ESM: `processes/biostabilization.py` + scour-mortality mode 6
- Vegetation → critical-shear feedback
- van Rijn equilibrium test passing

**Phase 4 — production polish (3–4 weeks)**
- Hotstart NetCDF; legacy `SEDBED_HOT.SDF` shim
- Streaming integration for bed state
- Optimization where profiling indicates
- Albany 587k-cell coupled run

**Phase 5 — wave-current shear, propwash hooks (deferred)**
- Mode C `WaveCurrentShearDriver`
- Stub propwash erosion API

Total effort estimate: **17–24 weeks of focused work** to reach phase 4. This excludes calibration data acquisition and case-study validation.

---

## 14. License & provenance strategy

- The EFDC+ source (GPL-2.0) is treated as a **reference implementation and documentation source**, not as code to copy. Algorithms are re-derived from:
  - SAND2008-5621 (public release, government work)
  - Jones & Lick (2001), van Rijn (1984, 1981), Cheng (1997), Soulsby (1997), Christoffersen & Jonsson (1985), Gessler (1965), Krone (1962), Lick (2008) — all peer-reviewed and not under EFDC's GPL
  - Variable names from the SEDZLJ Fortran (TSED, PERSED, TAUCRIT, etc.) — names are not copyrightable
  - File formats (bed.sdf, erate.sdf) — also not copyrightable
- ClearWater is government-funded and intended for public release; SSM should match (CC0 or US Government Public Domain). This is compatible with the algorithm-only re-implementation path.
- The EFDC+ source can be cited in docstrings as a verification reference — *"Cross-checked against EFDCPlus_Stable s_sedzlj.f90 commit X for behavioral equivalence"* — without creating a derivative-work tie.

---

## 15. Open questions for review

1. **Bedload as Riverine constituent vs. standalone explicit step?** Recommendation: standalone explicit (phase 1 simpler), revisit if bedload turns out to need implicit stability.
2. **Bed-state as wide DataArray vs. separate Dataset?** Recommendation: wide DataArray (chosen), monitor memory.
3. **Manning's n vs. log-law for shear?** Both supported; default to log-law (consistent with SEDZLJ). Should we expose both as user-facing options or just log-law?
4. **ESM scour-mortality mode 6 (excess shear)** — make it the default when SSM is present, or leave the existing modes as default with mode 6 opt-in?
5. **Calibration coefficient α for biostabilization** — start at 0.5 (placeholder), but we will need a published source to defend the default. Le Hir et al. (2007) is one candidate.
6. **Phase 1 scope** — is single-class/single-layer too thin to be useful as a milestone, or is it the right MVP? Could merge phase 1 + phase 2 for a meatier first release at the cost of a longer wait.
7. **Wave coupling** — anyone with a near-term need for wave–current interaction? If not, deferring Mode C to phase 5 is safe.
8. **Naming convention** — `ssm_suspended_<label>` for constituents is verbose. Alternatives: `sed_<label>`, `ssm_<label>`. Defer to the ClearWater team's established naming convention.

---

## 16. Appendix — key file:line references

EFDC SEDZLJ source (verification only):
- Top driver — `EFDC/SedTran-SEDZLJ/s_main.f90:10` (`SEDZLJ_MAIN`)
- Per-cell core — `s_sedzlj.f90:9` (`SEDZLJ(L)`)
- Initialization — `s_sedic.f90:9` (`SEDIC`)
- Shear — `s_shear.f90:9` (`SEDZLJ_SHEAR`); Parker log-law `:261`; Christoffersen-Jonsson `:285–306`
- Bedload — `s_bedload.f90:9` (`BEDLOADJ`); van Rijn equations `:147–150`
- Slope — `s_slope.f90:9`
- Toxics linkage — `EFDC/ChemFate/caltoxb.f90`

ClearWater-Riverine architecture:
- Solver — `transport.py:548` (`update`)
- Constituent registration — `transport.py:76, 188–201`; `constituents.py:271–368`
- Streaming — `transport.py:850–912`; design memo `design/cw_riverine_streaming_in_memory_release.md`
- HDF ingest — `io/hdf.py:48–86`

ClearWater-modules-v2 pattern:
- Process ABC — `processes/base.py:14`
- TSM template — `processes/temperature.py:15`
- Factory — `processes/base.py:61`

ESM phase 2 sediment touchpoints:
- HEC-RAS sediment reader — `esm/io/hecras_reader.py:46–58`
- Light extinction from TSS — `esm/processes/light.py:131–150`
- Burial mortality — `esm/model.py:1140–1180`
- Scour-mortality modes — `esm/constants.py:361–377`
- Future-feature stub — `esm/io/clearwater_interface.py:406–421`

References (algorithms):
- Thanh, Grace & James (2008) SAND2008-5621 — primary algorithmic reference
- Jones & Lick (2001) — SEDZLJ erosion-rate formulation
- van Rijn (1984a, 1984b) — bedload, suspension threshold
- Cheng (1997) — settling velocity
- Soulsby (1997) — critical shear stress
- Christoffersen & Jonsson (1985) — wave–current friction
- Gessler (1965), Krone (1962) — deposition probabilities
- Lick (2008) — active-layer thickness, slope effects
- Le Hir et al. (2007) — biostabilization (proposed)
