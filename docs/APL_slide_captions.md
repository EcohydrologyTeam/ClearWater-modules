# APL Slide Captions

Each slide gives a short **on-screen caption** (one line, for the slide itself), the **on-screen bullets** (the slide body, kept deliberately sparse so nothing crowds), and **speaker notes** (the fuller narration, for the notes pane below the slide). Slides are ordered: what the APL is → why it is a library, not a module → per-process timesteps.

---

## Slide 1 — What the APL is

**On-screen caption:**
The Aquatic Processes Library (APL) — heat and water-quality kinetics for ClearWater host models, now running in ClearWater-HMS.

**On-screen bullets:**
- **Temperature** — water-column heat budget + sediment temperature
- **WQ Constituents (11):**
  - Nitrogen · Phosphorus · Carbon · Dissolved Oxygen
  - CBOD · Alkalinity · N₂ · POM
  - Floating Algae · Benthic Algae · Pathogen
- **Riverine** — coupler that exchanges state with ClearWater-Riverine
- **Implemented, tested, and verified in ClearWater-HMS**

**Speaker notes:**
The APL is the library of water-column heat and water-quality processes that supply the kinetics for ClearWater's host transport models. It is built around the Temperature process — the full water-column heat budget, sediment-temperature evolution, the latent-heat unit fix, and a thin-water stability guard — and the eleven WQ Constituent processes: Nitrogen (ammonium, nitrate, organic N; nitrification / denitrification), Phosphorus (soluble reactive P, organic P, sediment–water exchange), Carbon (dissolved inorganic and organic pools), Dissolved Oxygen (re-aeration, SOD, BOD, photosynthesis / respiration), CBOD, Alkalinity (with pH coupling), N₂ (the denitrification product), POM, Floating Algae (single phytoplankton group in the release version, multi-group in dev), Benthic Algae (periphyton / attached algae), and Pathogen decay. The Riverine coupler bridges state to and from ClearWater-Riverine on HEC-RAS-2D. The library has been implemented with ClearWater-HMS and tested and verified there; the same process code couples to other transport drivers without modification. Throughout, the library computes the reaction and heat-exchange terms while the host model owns transport and advances the clock.

---

## Slide 2 — A library, not a module

**On-screen caption:**
A library, not a module — composable processes a host model *builds with*, not a monolith it *runs*.

**On-screen bullets:**
- **"Module" is legacy.** The original Fortran release packaged these capabilities as monolithic Simulation Modules — one block, run end to end.
- **The APL is a Python library.** Host models import it and call in; it never runs standalone.
- **Processes are composable objects.** Pick the set you want; the library resolves their order and computes the kinetics.
- **Kinetics here, transport and clock in the host.** That split is what makes it reusable across drivers.

**Speaker notes:**
The word "module" is inherited from the original Fortran code, where each capability was a self-contained Simulation Module you ran end to end. The modernized code is structured differently: it is a Python library of composable process objects. A host model constructs the library's `Model`, hands it the set of processes it wants, and the library resolves the dependency order and computes the kinetics each step. Nothing in the library runs on its own. It is infrastructure that ClearWater-HMS, ClearWater-Riverine, and other drivers build on. The one-line takeaway: a module is something you run; a library is something you build with, and the APL is the latter.

---

## Slide 3 — Per-process timesteps

**On-screen caption:**
Each process runs on its own clock — multi-rate execution on a shared base substep.

**On-screen bullets:**
- **Each process carries its own timestep** and integrates its own kinetics over it.
- **The base substep is the finest clock.** A process fires every *N* base steps (*N* = its step ÷ base step).
- **Rule:** each process step is an integer multiple of the base step, validated at startup.
- **Payoff:** fast kinetics step finely; slow pools update less often — compute where the dynamics need it.

**Speaker notes:**
In the current framework, each process is constructed with its own timestep and integrates its kinetics over that interval, rather than every process sharing a single global timestep as in the earlier framework. The model runs on a base substep — the finest clock in the simulation — and the orchestrator precomputes, for each step index, exactly which processes fire. A process fires every *N* base substeps, where *N* is its own timestep divided by the base step. The one rule is that each process's timestep must be an integer multiple of the base step; the model validates this at startup and raises if it is not, so the base step is the fastest any process can run and every other rate is a clean multiple of it. The payoff is multi-rate execution: fast-responding kinetics such as re-aeration and temperature can be stepped finely while slow pools such as benthic algae and sediment exchange update less frequently, putting compute where the dynamics demand it. (The bundled demo hands every process the same timestep for simplicity; the framework supports heterogeneous per-process rates.)
