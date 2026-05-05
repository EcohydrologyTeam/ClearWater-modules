# SSM — Cohesive-bed consolidation: design memo

**Status:** Draft for review (initial release scope)
**Author:** Generated 2026-05-02 from a literature scan and an algorithmic mapping into the existing SSM bed-state model.
**Audience:** SSM developers; reviewers familiar with EFDC SEDZLJ, MIKE 21, Delft3D, or TELEMAC cohesive-sediment treatments.
**Companion document:** `ssm_design_spec.md` §5.10.

---

## 1. Motivation

EFDC SEDZLJ — the algorithmic basis for SSM — explicitly does not model cohesive-bed consolidation. The Fortran source flags this:

```
! s_sedzlj.f90:707
! NOTE: SEDZLJ DOES NOT HAVE CONSOLIDATION
```

This is a known limitation. For freshly-deposited cohesive sediment, the *in situ* critical shear stress for erosion (τ_ce) increases over hours to weeks as the bed dewaters and grain-grain bonds tighten. Holding τ_ce constant at the unconsolidated value:

* over-predicts erosion of fresh deposits during subsequent flow events,
* makes erosion fluxes a strong (and often unphysical) function of the time-discretisation of deposition events,
* puts SSM at a behavioural disadvantage relative to MIKE 21, Delft3D, and TELEMAC, all of which carry some form of consolidation model.

The first SSM release closes this gap with the simplest defensible formulation: an opt-in, age-dependent τ_ce based on Sanford & Maa (2001).

## 2. Formulation

For each cohesive class (D₅₀ < `bedload_cutoff`, default 64 μm), the layer-effective critical shear stress at age $t_{\rm age}$ is:

$$
\tau_{ce}^{\rm eff}(t_{\rm age}) =
  \tau_{ce,\infty} - (\tau_{ce,\infty} - \tau_{ce,0})\,\exp(-t_{\rm age}/T_c)
$$

where:

| Symbol | Default | Description |
|---|---|---|
| $\tau_{ce,0}$ | 0.10 Pa | Freshly-deposited critical shear stress (lower bound; $t_{\rm age} \to 0$ asymptote) |
| $\tau_{ce,\infty}$ | 0.50 Pa | Fully-consolidated critical shear stress (upper bound; $t_{\rm age} \to \infty$ asymptote) |
| $T_c$ | 7 days (604 800 s) | Consolidation e-folding time |
| $t_{\rm age}$ | computed | Layer's mass-weighted mean age, advanced by `dt` each step |

The formulation has three virtues:

1. **Empirically grounded.** It reproduces the qualitative behaviour observed in flume experiments by Mehta & Partheniades (1975) and Sanford & Maa (2001) on estuarine mud: the strength gain is rapid in the first few days and asymptotes within a few weeks.
2. **Mathematically simple.** A three-parameter closed form, no iteration, no auxiliary state beyond the per-layer age field.
3. **Calibration tractable.** The three parameters are physically interpretable: $\tau_{ce,0}$ reflects sediment composition (organic content, clay mineralogy), $\tau_{ce,\infty}$ reflects the fully-consolidated bulk density, and $T_c$ reflects the dewatering time scale (a function of permeability and bed thickness).

For non-cohesive (sand) classes the effective τ_ce remains the static value; consolidation is a clay/silt phenomenon and the SEDflume tabular erosion-rate model already captures sand-bed sorting via the per-(D₅₀, τ) interpolation.

### 2.1 Behavioural checks

| $t_{\rm age}$ | $\tau_{ce}^{\rm eff}$ | Comment |
|---|---|---|
| 0 | $\tau_{ce,0}$ | Fresh deposit — minimum strength |
| $T_c$ | $\tau_{ce,0} + (1 - 1/e)(\tau_{ce,\infty} - \tau_{ce,0}) \approx \tau_{ce,0} + 0.632\,\Delta\tau$ | One e-folding of recovery |
| $5 T_c$ | $\tau_{ce,0} + 0.993\,\Delta\tau$ | Effectively fully consolidated |
| $\to \infty$ | $\tau_{ce,\infty}$ | Asymptote |

Unit tests in `tests/sediment/test_consolidation.py` lock these checkpoints in: `test_age_zero_returns_tau_ce_zero`, `test_age_at_tc_matches_one_minus_one_over_e`, `test_age_large_approaches_tau_ce_inf`.

## 3. Per-layer age tracking

The age field lives on the mesh as `ssm_bed_layer_age` (s; dims `(time, nface, ssm_layer)`, dtype `float32`). It is allocated by `bed.initialize_bed_state` to zero everywhere at $t=0$, mirroring the assumption that the IC corresponds to a fully-consolidated bed (worst case for erosion fluxes — an analyst who knows the bed is fresh can override the IC).

### 3.1 Time advancement

Every step, `bed.update_bed_elevation(dt_seconds=dt)` adds `dt` to every layer that holds non-zero mass, and pins empty layers to age 0. Backward compatibility: passing `dt_seconds=0` (the default) leaves the age field untouched, so legacy callers that don't yet propagate dt see the same behaviour as before.

### 3.2 Age dilution on deposition

When fresh mass $\Delta m$ is deposited into layer 1 (existing mass $m_1$, existing age $t_1$):

$$
t_{1,\rm new} = t_1 \cdot \frac{m_1}{m_1 + \Delta m}
$$

The deposit enters with age 0; the new layer-mean age is the existing age weighted by the existing-mass fraction of the total. This is the simplest and most common dilution rule in the depth-averaged consolidation literature and matches the mass-weighted blending already used elsewhere in `bed.py` for class fractions.

The rule is exposed as the standalone helper `bed.dilute_layer1_age_on_deposition`, called by `ssm.run` immediately after the per-step deposition mass is added to layer 1. When no deposit occurred (cell-wise), the helper is a no-op.

### 3.3 Age inheritance during reorganization

The active-layer reorganization in `bed.reorganize_active_layer` has three branches; each propagates ages by mass-weighted blending, mirroring the existing PERSED logic:

* **Branch (a) — net deposition (m_1 > T_act):** Excess mass $(m_1 - T_{\rm act})$ is pushed from layer 1 to layer 2. Layer 1's mean age is unchanged (uniform-aged mass is removed off the top); layer 2's new age is the mass-weighted blend $(t_2 m_2 + t_1 (m_1 - T_{\rm act})) / (m_2 + (m_1 - T_{\rm act}))$.
* **Branch (b) — borrow from SLLN:** Deficit $(T_{\rm act} - m_1)$ is borrowed from the next non-empty sub-layer (SLLN). Layer 1's new age is $(t_1 m_1 + t_{\rm SLLN} (T_{\rm act} - m_1)) / T_{\rm act}$; SLLN's mean age is unchanged.
* **Branch (c) — collapse:** SLLN is fully merged into layer 1. New layer 1 age is $(t_1 m_1 + t_{\rm SLLN} m_{\rm SLLN}) / (m_1 + m_{\rm SLLN})$; SLLN's age is reset to 0 (empty).

Tests in `test_consolidation.py::TestAgeInheritanceInReorganization` and `test_bed.py::TestLayerAgeTracking::test_age_inheritance_branch_b_borrow` lock the arithmetic in.

## 4. Calibration parameters and defaults

The package ships three constants in `contracts.py`:

```python
DEFAULT_CONSOLIDATION_TAU_CE_ZERO_PA = 0.10   # Pa
DEFAULT_CONSOLIDATION_TAU_CE_INF_PA  = 0.50   # Pa  (5× the lower bound)
DEFAULT_CONSOLIDATION_TIME_S         = 604_800  # 7 days
```

Provenance:

* The 0.10 Pa lower bound sits in the middle of the 0.05–0.15 Pa range typical for freshly-deposited estuarine mud reported by Sanford & Maa (2001). Sediments rich in organic matter or with high CEC clays may sit toward the high end; coarse silts toward the low end.
* The 5× ratio for $\tau_{ce,\infty}/\tau_{ce,0}$ matches the Sanford & Maa (2001) Fig. 4 Chesapeake Bay mud calibration; ratios of 3–5× are typical, with sandy mixed-grain beds at the low end and pure clay beds at the high end.
* The 7-day e-folding time is the canonical value from Mehta & Partheniades (1975) flume work; field-scale dewatering can be slower (weeks) for thick deposits or faster (hours) for thin layers in agitated environments. Site-specific calibration is recommended for any application with bed-evolution claims.

## 5. Configuration

The model is **opt-in**. By default, SSM behaves as today (no consolidation, τ_ce constant per layer). To enable:

```yaml
sediment:
  consolidation:
    enabled: true
    model: sanford_maa
    tau_ce_zero_pa: 0.10               # optional; defaults from contracts.py
    tau_ce_inf_pa: 0.50                # optional
    consolidation_time_s: 604800       # optional; 7 days
```

`from_config` parses this block into a `SanfordMaaConsolidation` instance and threads it through to the erosion model and the bed reorganization. The model is also accepted as a constructor kwarg (`SSM(..., consolidation_model=...)`) for programmatic instantiation.

A `model: sanford_maa` selector is the only currently supported value. Future single-mode variants (e.g. logistic recovery, two-mode with a fast and slow exponent) can register additional names without breaking the schema.

## 6. Wiring into the run loop

In `ssm.run`, after the per-class effective τ_ce (post-vegetation) has been computed:

1. Broadcast `tau_ce_eff` (per-class) across the layer dimension to obtain a baseline `tau_ce_layer_class[face, layer, class]` of shape `(nface, n_layers, n_class)`.
2. If a consolidation model is configured, call `erosion.apply_consolidation` to overwrite the cohesive-class entries with the age-adjusted values, broadcasting the per-`(face, layer)` aged value across the cohesive class indices.
3. Use the per-`(face, layer, class)` τ_ce in the per-layer erosion gate:

   ```python
   gate_k = (tau_arr[:, None] >= tau_ce_layer_class[:, k, :])
   ```

   replacing the previous per-`(face, class)` gate.

When no consolidation model is configured, the broadcast in step 1 is the entire effect: `tau_ce_layer_class[face, k, c] == tau_ce_eff[face, c]` for all `k`, recovering the previous gating logic exactly.

## 7. Limitations and what's deferred

This first release is intentionally narrow. A complete cohesive-bed consolidation model would additionally evolve:

| Phenomenon | Rationale for deferral |
|---|---|
| Time-varying porosity / bulk density (Gibson, England & Hussey 1967) | Requires a 1-D continuity / momentum solve in the bed column at each cell each step; a non-trivial expansion of the bed-state vector and the time-stepping loop. Defer to v2. |
| Gel point and floc-network mechanics (Toorman 1999) | Adds an upper-bound constraint on solids volume fraction; tractable but not needed for a first release where the bulk density is held constant. |
| Finite-strain (large-deformation) self-weight consolidation | The Gibson PDE in its non-linearised form requires moving-mesh handling in the bed column; a substantial implementation lift. |
| Multi-mode consolidation (fast + slow exponents) | The single-mode form is adequate for most applications; multi-mode can be added by composing two `SanfordMaaConsolidation` instances and summing the effects. |
| Bioturbation / biostabilization of consolidated beds | ESM coupling exists for biostabilization (`coupling.read_vegetation_feedback`); the consolidation model treats bioturbation as zero-order (no rejuvenation of layer ages). |
| Layer-structured τ_ce profile in fresh deposits | The current model treats the full layer as having a single mass-weighted age. Real fresh deposits may have a measurable depth-gradient (the bottom of the layer is older than the top). Mitigated in practice by the multi-layer K_B stack. |

When the time comes to add Gibson-style porosity evolution, the natural extension is to add a parallel `ssm_bed_layer_porosity` field on the mesh, advance it via a per-cell column solve in `update_bed_elevation`, and let `consolidation.SanfordMaaConsolidation` consume the porosity as a fourth input alongside the age. The current API is forward-compatible with that extension: the Protocol takes one DataArray now and can be extended to take more later without breaking existing callers.

## 8. Testing strategy

Unit tests in `tests/sediment/test_consolidation.py` cover:

* Asymptotes ($t \to 0$, $t = T_c$, $t \to \infty$) and monotonicity.
* Validation rejects non-positive $T_c$ or $\tau_{ce,\infty} < \tau_{ce,0}$.
* `apply_consolidation_per_class` leaves non-cohesive classes untouched.
* Layer-age advances by `dt` each call to `update_bed_elevation`; empty layers pinned to 0; `dt=0` is a no-op (backward compatibility).
* Age dilution on deposition gives the mass-weighted result; deposit into an empty layer yields age 0; zero-deposit is a no-op.
* Borrow / promote / collapse correctly propagate ages.
* Integration: a single-cell single-class run with sustained low shear plus periodic deposition shows τ_ce rising over time, monotonically and approaching the asymptote.

Tests in `tests/sediment/test_bed.py::TestLayerAgeTracking` exercise the bed-side mechanics independently of the consolidation closure.

## 9. References

* Gibson, R. E., England, G. L., and Hussey, M. J. L. (1967). "The theory of one-dimensional consolidation of saturated clays." *Géotechnique* 17(3), 261–273.
* Mehta, A. J., and Partheniades, E. (1975). "An investigation of the depositional properties of flocculated fine sediments." *J. Hydraul. Res.* 13(4), 361–381. DOI: 10.1080/00221687509499694.
* Sanford, L. P., and Maa, J. P.-Y. (2001). "A unified erosion formulation for fine sediments." *Marine Geology* 179(1–2), 9–23. DOI: 10.1016/S0025-3227(01)00201-8.
* Toorman, E. A. (1999). "Sedimentation and self-weight consolidation: constitutive equations and numerical modelling." *Géotechnique* 49(6), 709–726.

EFDC SEDZLJ source for behavioural cross-reference (verification only, not derivation): `EFDC/SedTran-SEDZLJ/s_sedzlj.f90:707`.
