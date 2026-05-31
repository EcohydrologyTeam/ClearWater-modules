# Design Spec: Inorganic phosphorus sorption partitioning (NSM1-I port completion)

**Date:** 2026-05-30
**Component:** `src/clearwater_modules_v3/processes/phosphorus.py`,
`src/clearwater_modules_v3/parameters/phosphorus.py`
**Severity:** Incomplete Fortran -> Python port. Affects TIP settling and algal
phosphorus limitation.
**Status:** COMPLETE 2026-05-30 (no code change beyond the shared `Solid` wiring).
The `fdp` machinery was already fully wired in v3 — `phosphorus.py` uses
`(1-fdp)·TIP` for TIP settling, and `floating_algae`/`benthic_algae` pass
`fdp(use_TIP, Solid, kdpo4)` into `limit_phosphorus` as the dissolved fraction
(`FP = fdp·TIP/(KsP+fdp·TIP)`). The shared-`Solid`-input work made all three read
the canonical `Solid` field. **The only open item was the `kdpo4` default, and
the research conclusion is: keep `kdpo4 = 0.0`.**

## kdpo4 research (2026-05-30 — NSM1 Fortran + HEC-RAS-WQ + v1 + CE-QUAL-W2)

The water-column inorganic-P partition coefficient default is **0.0 (L/kg)** in
every reference implementation:

| Source | water-column `kdpo4` default |
|---|---|
| NSM1 Fortran (`fortran/NSM1/01_extraneous/nsmi_main.f90:1101`) | **0.0** (typ. 0–200, max 0–80000) |
| HEC-RAS-WQ NSMII (`modMain.f90:408`) | **0.0** |
| v1 legacy Python (`clearwater_modules/nsm1/constants.py:238`) | **0.0** |
| CE-QUAL-W2 `PARTP` | application-set, **off (0) by default**; redox-coupled form |

The only non-zero default (`kdpo42 = 20.0 L/kg`, `modSedFlux.f90:479`) is the
**sediment-layer** partition coefficient — sediment-flux / NSM2 territory, not
the water-column `kdpo4`. That is the distinction the v3 `parameters/phosphorus.py`
"NSM2 territory" comment was reaching for.

**Conclusion:** there is no defensible non-zero NSM1-I water-column default to
adopt; `kdpo4 = 0.0` (→ `fdp = 1.0`, fully dissolved) **is** the NSM1-I default
and matches all four references. Sorption is enabled per-application when site
data exist. The spec's earlier "set a non-1.0 `kdpo4` default" wording was
imprecise. No trajectory change, no re-baseline. (Original proposal below
retained for the algebra/reference.)

---

**Original proposal:** Proposed. Confirmed against the Fortran NSM1-I source.

> **Review note (2026-05-30, code-verified):** `utils/partitioning.py:fdp(...)`
> **already exists** and its docstring shows it was cross-checked against the
> Fortran (`modGlobalParam.f90:228`). It uses the **dissolved-fraction** convention
> `fdp = 1 / (1 + kdpo4*Solid*1e-6)` (<= 1), consistent with `phosphorus.py`'s
> `(1 - fdp)*TIP` settling term. The Fortran snippet quoted below uses the
> reciprocal **denominator** convention `fdp = 1 + sum(kdpo4*Solid/1e6)` (>= 1);
> with that convention the quoted settling `(1 - fdp)*TIP` would be **negative**,
> so the snippet as written is internally inconsistent. Reconcile toward v3's
> existing util (dissolved-fraction), not the snippet.
>
> **Corrected scope:** wire the existing `utils.partitioning.fdp` into the algal /
> benthic P-limitation factors and TIP settling, supply the `Solid` / `kdpo4`
> inputs, and set a non-1.0 `kdpo4` default — not implement `fdp` from scratch.
> **Shared prerequisite:** a `Solid` input source (see the matching note in
> `clearwater_modules_v3_light_extinction.md`).

## Summary

v3 sets the dissolved-phosphorus fraction `fdp = 1.0` (TIP fully dissolved), with
a `FIXME(phase1-audit)` noting the simplification. The Fortran NSM1-I
**computes** `fdp` from suspended solids and a sorption partition coefficient.
This is a port gap, not a modeling choice.

## Fortran NSM1-I reference

`fortran/NSM1/04_phosphorus/nsmi_phosphorus.f90`:

```fortran
fdp = 1.0
do i = 1, nGS
  fdp = fdp + kdpo4(i,r) * Solid(i) / 1.0E6
end do
DIP = TIP / fdp
TIP_Settling = vs(r) / depth * (1.0 - fdp) * TIP   ! particulate (sorbed) fraction settles
```

and the algal P-limitation factors use the partitioned dissolved fraction
(`fortran/NSM1/03_dissolved_oxygen/nsmi_algae.f90` and `nsmi_benthic_algae.f90`):

```fortran
FP  = fdp * TIP / (KsP(r)  + fdp * TIP)
FPb = fdp * TIP / (KsPb(r) + fdp * TIP)
```

So `fdp` (a) reduces the bioavailable dissolved inorganic P seen by algae, and
(b) sets the sorbed fraction `(1 - fdp)`-equivalent that settles to the bed.

## Required change

1. Add the sorption partition coefficient `kdpo4` (per solids class) and the
   suspended-solids field `Solid` as inputs (suspended solids are already needed
   by the light-extinction work; coordinate the shared `Solid` input).
2. Compute `fdp` from `kdpo4` and `Solid` each step, matching the Fortran
   accumulation, and use it in:
   - the algal and benthic-algal P-limitation factors (`FP`, `FPb`), replacing
     the bare `TIP` with `fdp * TIP`; and
   - TIP settling of the particulate (sorbed) fraction.
3. Default `kdpo4` to the NSM1-I value(s). With `kdpo4 = 0` (or no solids), `fdp`
   reduces to 1.0, recovering the current fully-dissolved behavior, so the change
   is backward-compatible when sorption is off.

Resolve the exact `fdp` / `DIP` / settling algebra against the Fortran (the
`fdp = 1 + sum(kdpo4*Solid)` convention with `DIP = TIP/fdp` is unusual; reproduce
it faithfully rather than re-deriving).

## Verification

- Unit: with known `Solid` and `kdpo4`, assert `fdp`, `DIP`, the P-limitation
  factors, and TIP settling match the Fortran NSM1-I to tolerance; assert
  `kdpo4 = 0` recovers `fdp = 1.0`.
- Integration: the coupled Willowbend run with partitioning active shows reduced
  bioavailable P and nonzero TIP settling relative to the fully-dissolved run.

## Report follow-up (Report 2 repo)

After the fix: re-run the demonstration (TIP and algal-P-limitation results
change) and document inorganic-P partitioning in the Methods phosphorus section,
removing the fully-dissolved assumption note. Tracked as item A3.
