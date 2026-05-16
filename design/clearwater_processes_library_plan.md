# ClearWater Processes Library: Architecture, Naming, and Cutover Plan

> **Status:** Plan of record. Naming and architecture decided; execution gated (§7).
> **Date:** 2026-05-16
> **Author:** Todd Steissberg (ERDC)
> **Scope:** Identity, architectural role, embedding model, and rename/cutover plan for ClearWater's core water-column process library.

## 1. Summary

ClearWater's core reaction engine consists of the merged limnological water-column processes that simulate the heat budget (temperature), constituent kinetics (e.g., nutrient cycling and dissolved oxygen dynamics), algae, related processes (e.g., light limitation), and other features (e.g., pathogens). The project currently stores these processes in `src/clearwater_modules_v3/`. The team is developing this code as a standalone library that inter-links with the other components of the ClearWater modeling system. The project will describe and use the library as follows:

- **Repository:** `ClearWater-processes`
- **Package / import:** `from clearwater import processes`
- **Identity (brand):** the ClearWater Aquatic Processes Library (APL)

The repository and import use plain, convention-consistent words (the address). The project will use "APL" as the spoken and written identity in documentation, citation, presentations, and when the library is embedded in third-party host models. The rename is deferred to the v3 1.0 cutover (§6--§7).

## 2. The ClearWater taxonomy

ClearWater is a modeling system. Its prior acronym, Corps Library for Environmental Analysis of Watersheds, will be modified to better describe what it has evolved into: a Computational Laboratory for Environmental Analysis and Restoration of Watersheds.

| Component | Role | Kind |
|---|---|---|
| ClearWater-Riverine | Transport operator (advection-diffusion) | Core engine |
| ClearWater-processes (APL) | Water-column physical/chemical/biological reaction core | Core engine |
| BSM | Aqueous geochemistry / speciation (PHREEQC core) | Coupled module |
| ESM | Vegetation / ecohydrology | Coupled module |
| SSM | Sediment / solids | Coupled module |
| Sediment diagenesis | Di Toro multi-G benthic flux model | Coupled module (synchronously flux-coupled; §7) |

Definitions used throughout:

- **Core** -- an in-place operator on the shared transport state vector. The processes library is the reaction core, the sibling of ClearWater-Riverine's transport core. "Core" is an internal architectural term, not the public name.
- **Module** -- an encapsulated component with its own domain, state, and cadence, coupled through explicit state exchange via the orchestrator. The taxonomy reserves this term for BSM/ESM/SSM and the sediment-diagenesis module. The processes library is not a module within ClearWater.
- **Library** -- the frame-invariant description of the processes engine: callable, embeddable, driver-agnostic code. It is a core within ClearWater and, when embedded in a host model, a library/module of that host (§5).

## 3. The processes library

ClearWater-processes is the Aquatic Processes Library (APL), the core limnological algorithms for the water column: heat budget, light, gas exchange, nutrient cycling (N, P, C), dissolved oxygen, alkalinity, algae, benthic algae, organic matter, pathogens, and CBOD (the physical, chemical, and biological processes of inland waters).

- **Scope is general.** The library targets rivers, reservoirs, and watershed/runoff applications; it is not river-specific.
- **Positioning.** The v3 code is effectively new, derived from the primary water-quality literature, not a translation of legacy Fortran/Python kernels. The authoritative validation is by literature derivation and real-world case studies (e.g., Santiam-Salem, Willamette). Legacy v1-parity tests are internal regression scaffolding that retire with v1; they are not the external validation basis. (Open item: formalize the validation report -- §7.)

## 4. Naming and identity

The naming uses two deliberately separated layers.

**Address -- plain, convention-consistent:**

- Repository `ClearWater-processes` matches `ClearWater-Riverine` and `ClearWater-data`: plain words that say what the repository is, and the name corrects the taxonomically false `ClearWater-modules` (this is the core library, not a module collection).
- Import `from clearwater import processes` parallels `riverine` and `data`.

**Identity / brand:** the ClearWater Aquatic Processes Library (APL), used in the repository description, README, citations, presentations, and whenever the library is embedded standalone in host models.

Rationale:

- **Address vs. identity.** The address is intentionally generic and descriptive; the brand (APL) carries identity and discoverability. ClearWater has no captive distribution channel, so the project must carry the APL brand deliberately and consistently in all documentation and citation. The brand is the load-bearing identity layer, not optional polish.
- **"Aquatic" is intentional and frame-relative.** Under the ClearWater umbrella, "Aquatic" is mildly redundant (ClearWater already implies water). Standalone, embedded in a host model such as AdH or HEC-RAS-1D, there is no "ClearWater" prefix, and "Aquatic Processes Library" is self-describing and not redundant. The project chose the expansion for the embedded context, where the standalone name carries the weight.

## 5. Embedding and dual-topology coupling

The v3 code already separates reaction (each `Process.run(time, registry)`) from transport/hydraulics ingestion (the `Model` loads data sources into the registry) from orchestration (the `Model` owns the schedule and time loop). This is operator splitting realized in code, so one core serves two topologies:

- **ClearWater-driven:** ClearWater-Riverine owns transport and the loop and calls the reaction operator.
- **Host-driven:** a host model (AdH, HEC-RAS-1D) owns the mesh, transport, and loop and calls the same reaction operator as an embedded library, replacing, for example, the legacy C++ NSM linked into AdH.

Design commitments:

- The library exposes a thin reaction-operator API (`init / set_state / react(dt) / get_state`), batched over all active cells per step, and ClearWater-Riverine calls the library through that same API so the embedding contract is the only contract.
- The design follows the BMI / PhreeqcRM pattern (a reaction engine embedded in external transport hosts).

Caveats for any embedding work: batched calls are mandatory (a per-cell Python call inside a host solver dominates runtime); a constituent/units adapter maps the host's transported set to the library's state; the host transports temperature while the library supplies the energy-budget source term; sediment-diagenesis flux coupling is stiff and must live inside the kernel call, not across the host boundary; cross-topology results must agree within operator-split and Δt tolerance; the team must not target bit-parity with the legacy C++ NSM (validate against the literature-grounded library with a documented behavior diff).

Prerequisite: a kernel-purity audit. The team must verify and then enforce as a tested invariant that no process reaches outside the registry and that the registry exposes array-level get/set efficient enough to drive from a host adapter. "Zero hidden transport/IO/scheduling dependencies" must be an explicit, tested invariant.

## 6. Cutover sequence

The rename is deferred and executed as one discrete step within the v3 1.0 integration:

1. Close the open v3 NSM1 review CRITICALs.
2. Retire v1 (`clearwater_modules`): delete the package; migrate parity tests to frozen reference values.
3. Merge `streaming` → upstream `EcohydrologyTeam/ClearWater-modules` `main` (the v3 1.0 integration).
4. In one discrete step at that point:
   - Rename the repository `ClearWater-modules` to `ClearWater-processes`. GitHub redirects preserve existing clones, links, and citation URLs.
   - Rename the package `clearwater_modules*` to `clearwater.processes`. No redirect exists for an import path, so this is the breaking change. Grep the downstream ecosystem (ESM, Riverine, notebooks) and announce the import change ahead of the cut.
   - Set the repository description and README to the APL identity.

## 7. Status and open items

**Decided:** the taxonomy (§2); the processes-library positioning (§3); the address/identity naming and the `ClearWater-processes` / `from clearwater import processes` / APL set (§4); the dual-topology direction and embedding API approach (§5); the gated cutover sequence (§6).

**Open, resolve before or during execution:**

- **Packaging convention.** `from clearwater import processes` implies a `clearwater` namespace package. Existing ecosystem packages (`clearwater_data`, `clearwater_riverine`) are flat. The team must decide whether the ecosystem adopts the `clearwater.` namespace or reconciles this another way, and align Riverine/data accordingly.
- **Sediment-diagenesis coupling interface:** separate synchronously flux-coupled module vs. bound sub-model (the per-step DO/N/P exchange is stiff).
- **ESM co-location:** whether and when ESM lives in-repo as its own module package.
- **Embedding delivery format:** embed CPython vs. compiled native kernel vs. out-of-process coupling.
- **Kernel-purity audit** (§5): the team must run and enforce it before committing to the host-embedding path.
- **Validation report:** formalize the literature-derivation and case-study validation that establishes the library's authority independent of legacy parity.

**Execution gate:** nothing in §6 begins before step 1 (CRITICALs closed). The rename lands as part of the upstream-merge step, not piecemeal on the `streaming` fork beforehand.
