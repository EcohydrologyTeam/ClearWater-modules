"""Pytest fixtures and helpers for v3 NSM1 Tier 1 closed-system tests.

The Tier 1 contract (design spec Section 9, Section 14 resolved Q7):

    Closed system + no boundaries + no settling + balanced source/sink
    pairs --> total mass of N, P, C, O2-equivalents, Alk constant to
    floating-point roundoff AND ``diagnostics.clip_events == {}``.

This module exposes the building blocks every Tier 1 conservation test
needs:

* ``ClosedSystemConfig`` — a dataclass collecting the parameter overrides
  that turn off boundaries, settling, sediment exchange, and atmospheric
  exchange. Phase 2-6 tests pass this through to whatever Process they
  exercise.
* ``InMemoryRegistry`` — a stripped-down ``VariableRegistry`` stand-in
  that supports ``register`` / ``get`` / ``get_at_time`` / ``set_at_time``
  in-place (no time axis). The Tier 1 tests don't exercise the real
  registry's chunked time-axis behavior; they need a fast in-memory
  store for state.
* ``closed_system_config`` (fixture) — the canonical zero-flux config.
* ``initial_state_5cell`` (fixture) — five-cell initial conditions for
  every state variable that a Tier 1 test might inspect. Values are
  positive, modest in magnitude, and chosen so that source/sink balances
  are non-trivial.
* ``in_memory_registry`` (fixture) — registers ``initial_state_5cell``
  into a fresh ``InMemoryRegistry``.
* ``total_n`` / ``total_p`` / ``total_c`` / ``total_o2_equivalents`` /
  ``total_alkalinity`` — helper functions that sum mass-equivalent
  contributions across all reservoirs. Stoichiometric coefficients are
  taken from ``clearwater_modules_v3.parameters`` (algae N:Chl-a etc.).

Phase 1.4 scaffolding note: the helpers cover every state variable named
in design-spec Section 4 / Appendix A. Reservoirs whose Process classes
don't yet exist (CBOD, OrgN, TIP, OrgP, POC, DOC, DIC, POM, DOX, Alk,
N2, Pathogen) are summed if present in the registry and silently
omitted if absent. The helpers therefore work today with only NH4 +
NO3 + algae, and remain correct as Phases 2-6 land more state.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Any, Iterable

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.parameters.algae import DEFAULTS as ALGAE_DEFAULTS
from clearwater_modules_v3.parameters.balgae import DEFAULTS as BALGAE_DEFAULTS


# ---------------------------------------------------------------------------
# Stoichiometric ratios (sourced from parameter library)
# ---------------------------------------------------------------------------
# Floating algae per-Chla mass ratios (mg-X / ug-Chla) from algae.DEFAULTS:
#   AWn  -- N per Chla
#   AWp  -- P per Chla
#   AWc  -- C per Chla
#   AWd  -- dry weight per Chla (used to convert biomass back to org-X)
# Benthic algae are reported in g-D / m^2 with BWn / BWp / BWc per g-D, so
# the conversion to mg-X is BWn * Ab (g-D/m^2) * 1000 (mg/g).

AP_N_PER_CHLA: float = float(ALGAE_DEFAULTS["AWn"])     # mg-N per ug-Chla
AP_P_PER_CHLA: float = float(ALGAE_DEFAULTS["AWp"])     # mg-P per ug-Chla
AP_C_PER_CHLA: float = float(ALGAE_DEFAULTS["AWc"])     # mg-C per ug-Chla
AB_N_PER_GD: float = float(BALGAE_DEFAULTS["BWn"])      # mg-N per g-D
AB_P_PER_GD: float = float(BALGAE_DEFAULTS["BWp"])      # mg-P per g-D
AB_C_PER_GD: float = float(BALGAE_DEFAULTS["BWc"])      # mg-C per g-D

# C:N (Redfield-ish) used for converting POM / CBOD to N-equivalents in
# Tier 1 closed-system tests. v3 NSM1 sets per-process ratios; Phase 3+
# tests should override these via fixture parameters when their Process
# class declares its own ratio.
POM_C_TO_N: float = 5.68    # mg-C per mg-N (Redfield 106:16)
POM_C_TO_P: float = 41.1    # mg-C per mg-P (Redfield 106:1 by atoms => 41.1 by mass)

# CBOD-to-O2 coupling: 1 mg CBOD == 1 mg-O2 by definition (BOD is an
# oxygen-demand quantity).
CBOD_O2_PER_MASS: float = 1.0

# DOC-to-O2 oxidation stoichiometry: roughly 32/12 = 2.67 mg-O2 per mg-C.
DOC_O2_PER_C: float = 32.0 / 12.0


# ---------------------------------------------------------------------------
# Closed-system config
# ---------------------------------------------------------------------------


@dataclass
class ClosedSystemConfig:
    """Parameter overrides that put the v3 NSM1 system in closed-system mode.

    Phase 2-6 tests construct Process classes with these overrides so the
    Tier 1 conservation invariant holds. Every keyword corresponds to a
    physical pathway that, if non-zero, would either remove or add mass
    to the closed system:

    * Settling velocities (``vsap``, ``vsbp``, ``vs``, ``vsop``, ``vsom``,
      ``vsoc``) — bulk transport of particulates out of the water column.
    * Sediment fluxes (``NH4fromBed`` / ``rnh4_20``, ``DIPfromBed``,
      ``NO3_BedDenit`` / ``vno3_20``, ``DIC_sed_release``, ``SOD_20``) —
      mass exchange with the sediment compartment.
    * Atmospheric exchange (``kah_20_user``, ``kaw_20_user``) —
      reaeration / DIC outgassing.
    * Boundary inflows / outflows — handled at the Model orchestration
      level (zero by construction in the closed-system fixture).
    """

    # Settling velocities (m/d). Each turns off settling for one reservoir.
    vsap: float = 0.0
    vsbp: float = 0.0
    vs: float = 0.0
    vsop: float = 0.0
    vsom: float = 0.0
    vsoc: float = 0.0

    # Sediment fluxes. Names match the v1 NSM1 constants and the v3
    # nitrogen.DEFAULTS / phosphorus.DEFAULTS / carbon.DEFAULTS keys
    # (post-Phase 1 audit may rename; the ``# FIXME`` markers in those
    # files are tracked separately).
    rnh4_20: float = 0.0           # NH4 from bed; nitrogen.DEFAULTS key
    DIPfromBed: float = 0.0
    vno3_20: float = 0.0           # NO3 bed denitrification
    DIC_sed_release: float = 0.0
    SOD_20: float = 0.0
    SOD_theta: float = 1.0         # benign even if SOD_20 == 0

    # Atmospheric exchange. ``-1`` selectors in the v1 menu mean
    # "user-supplied"; setting to 0 disables atmospheric flux entirely.
    kah_20_user: float = 0.0
    kaw_20_user: float = 0.0
    hydraulic_reaeration_option: int = -1
    wind_reaeration_option: int = -1

    # Hotstart-friendly placeholder for future fluxes; tests can extend.
    extra: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        """Flat dict suitable for passing as ``**kwargs`` to a Process."""
        out = {
            "vsap": self.vsap,
            "vsbp": self.vsbp,
            "vs": self.vs,
            "vsop": self.vsop,
            "vsom": self.vsom,
            "vsoc": self.vsoc,
            "rnh4_20": self.rnh4_20,
            "DIPfromBed": self.DIPfromBed,
            "vno3_20": self.vno3_20,
            "DIC_sed_release": self.DIC_sed_release,
            "SOD_20": self.SOD_20,
            "SOD_theta": self.SOD_theta,
            "kah_20_user": self.kah_20_user,
            "kaw_20_user": self.kaw_20_user,
            "hydraulic_reaeration_option": self.hydraulic_reaeration_option,
            "wind_reaeration_option": self.wind_reaeration_option,
        }
        out.update(self.extra)
        return out


# ---------------------------------------------------------------------------
# In-memory registry
# ---------------------------------------------------------------------------


class InMemoryRegistry:
    """Minimal ``VariableRegistry`` stand-in for Tier 1 closed-system tests.

    The real ``clearwater_data.variables.VariableRegistry`` carries a time
    axis and chunk-store hooks. Tier 1 closed-system conservation tests
    need neither: they read pre-step state, write post-step state, and
    sum across cells. ``InMemoryRegistry`` matches the subset of the
    interface that v3 ``Process.run`` and the v3 ``Model`` actually call
    (``register``, ``get``, ``get_at_time``, ``set_at_time``,
    ``get_variable``, ``__contains__``).
    """

    def __init__(self) -> None:
        self._data: dict[str, xr.DataArray] = {}

    def register(self, name: str, value: xr.DataArray) -> None:
        self._data[name] = value

    def get(self, name: str) -> xr.DataArray:
        return self._data[name]

    def get_at_time(self, name: str, time: datetime) -> xr.DataArray:
        return self._data[name]

    def set_at_time(self, name: str, time: datetime, value: xr.DataArray) -> None:
        self._data[name] = value

    def get_variable(self, name: str):
        raise KeyError(name)

    def __contains__(self, name: str) -> bool:
        return name in self._data

    def keys(self) -> Iterable[str]:
        return self._data.keys()

    def snapshot(self) -> dict[str, xr.DataArray]:
        """Shallow copy for end-of-run conservation comparison."""
        return {k: v.copy() for k, v in self._data.items()}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="function")
def n_cells() -> int:
    """Number of cells in the synthetic mesh."""
    return 5


@pytest.fixture(scope="function")
def closed_system_config() -> ClosedSystemConfig:
    """Default closed-system parameter overrides (everything zero)."""
    return ClosedSystemConfig()


@pytest.fixture(scope="function")
def initial_state_5cell(n_cells: int) -> dict[str, xr.DataArray]:
    """Physically reasonable non-zero initial conditions on a 5-cell mesh.

    Units (matching v1 NSM1 conventions):
    * NH4, NO3, OrgN, TIP, OrgP, POC, DOC, DIC, POM, CBOD, DOX -- mg/L
    * Ap (floating algae) -- ug-Chla / L
    * Ab (benthic algae)  -- g-D / m^2
    * N2 (TDG)            -- mg-N / L
    * Alk                 -- mg-CaCO3 / L
    * Pathogen (PX)       -- count / 100 mL (units ignored for Tier 1)
    * water_temperature   -- degC
    * depth, volume, surface_area -- m, m^3, m^2

    Values are spread across cells to avoid trivial uniform-state tests.
    """

    def _da(values: list[float]) -> xr.DataArray:
        return xr.DataArray(np.array(values, dtype=float), dims="cell")

    return {
        # Nitrogen reservoirs
        "ammonium": _da([0.05, 0.10, 0.15, 0.20, 0.30]),       # NH4, mg-N/L
        "nitrate": _da([1.0, 2.0, 3.0, 4.0, 5.0]),             # NO3, mg-N/L
        "organic_nitrogen": _da([0.20, 0.25, 0.30, 0.35, 0.40]),
        "n2": _da([10.0, 10.5, 11.0, 11.5, 12.0]),
        # Phosphorus reservoirs
        "tip": _da([0.10, 0.12, 0.14, 0.16, 0.18]),
        "organic_phosphorus": _da([0.05, 0.06, 0.07, 0.08, 0.09]),
        # Carbon reservoirs
        "poc": _da([1.0, 1.2, 1.4, 1.6, 1.8]),
        "doc": _da([2.0, 2.2, 2.4, 2.6, 2.8]),
        "dic": _da([5.0, 5.5, 6.0, 6.5, 7.0]),
        # Particulate organic matter
        "pom": _da([3.0, 3.5, 4.0, 4.5, 5.0]),
        # CBOD (single group; multi-group adds 'cbod_2', 'cbod_3', ... in Phase 3)
        "cbod": _da([2.0, 2.5, 3.0, 3.5, 4.0]),
        # DOX
        "oxygen_dissolved": _da([8.0, 8.5, 9.0, 9.5, 10.0]),
        # Alkalinity
        "alkalinity": _da([100.0, 105.0, 110.0, 115.0, 120.0]),
        # Pathogen
        "pathogen": _da([1e3, 5e3, 1e4, 5e4, 1e5]),
        # Algae
        "ap": _da([5.0, 6.0, 7.0, 8.0, 10.0]),                 # ug-Chla/L
        "ab": _da([1.0, 1.5, 2.0, 2.5, 3.0]),                  # g-D/m^2
        # Forcings (unchanged through a closed-system test)
        "water_temperature": _da([20.0, 20.0, 20.0, 20.0, 20.0]),
        "depth": _da([1.0, 1.5, 2.0, 2.5, 3.0]),
        "volume": _da([1.0, 1.0, 1.0, 1.0, 1.0]),
        "surface_area": _da([1.0, 1.0, 1.0, 1.0, 1.0]),
    }


@pytest.fixture(scope="function")
def in_memory_registry(initial_state_5cell: dict[str, xr.DataArray]) -> InMemoryRegistry:
    """Fresh ``InMemoryRegistry`` pre-loaded with ``initial_state_5cell``."""
    registry = InMemoryRegistry()
    for name, value in initial_state_5cell.items():
        registry.register(name, value)
    return registry


@pytest.fixture(scope="function")
def closed_system_time_window() -> tuple[datetime, datetime, timedelta]:
    """``(start_time, end_time, time_step)`` for a 100-step run.

    Hundred 5-minute substeps => 8h 20m of simulated time, a duration
    over which closed-system kinetic balances should expose any
    integrator drift while still being short enough for tests to stay
    snappy.
    """
    start = datetime(2026, 1, 1, 0, 0, 0)
    time_step = timedelta(minutes=5)
    end = start + 100 * time_step
    return start, end, time_step


# ---------------------------------------------------------------------------
# Conservation helper functions
# ---------------------------------------------------------------------------
# Each helper returns the closed-system total of a conserved quantity as a
# 0-d DataArray (sum across cells). Tier 1 tests compare ``helper(initial)``
# against ``helper(final)`` with ``np.testing.assert_allclose(rtol=1e-12)``.
#
# For Phase 1.4 the helpers are intentionally permissive: variables not
# present in the registry are silently skipped. As Phases 2-6 land
# additional Process classes, the same helper functions sum more
# reservoirs without re-coding.


def _get(registry, name: str) -> xr.DataArray | None:
    """Return ``registry[name]`` or ``None`` if not present.

    Both ``InMemoryRegistry`` and the real ``VariableRegistry`` support
    ``__contains__`` and ``.get(name)``.
    """
    if name not in registry:
        return None
    return registry.get(name)


def _sum_over_cells(da: xr.DataArray) -> xr.DataArray:
    """Sum across the cell dimension, returning a 0-d DataArray."""
    if "cell" in da.dims:
        return da.sum(dim="cell")
    return da.sum()


def total_n(registry: Any) -> xr.DataArray:
    """Total nitrogen mass (mg-N) summed across all reservoirs.

    Reservoirs included (when present):
    * NH4, NO3, OrgN — direct mg-N/L pools
    * 2 * N2 — N2 gas reported as mg-N/L (each molecule is 2 N atoms,
      so molar count carries 2 N atoms; v1 reports N2 already as
      mg-N/L so the factor is 1; this helper assumes the v1
      convention)
    * Ap * AWn — algal Chl-a converted to N at the per-Chla ratio
    * Ab * BWn — benthic algae dry-weight converted to N
    * POM * (1/POM_C_to_N) * (1/(C-to-N)) -- POM as N-equivalent;
      treated as zero contribution if POM is in dry-weight units
      that don't trivially convert to N (Phase 3 will revisit when
      POM Process lands and the units convention is final).
    * CBOD does NOT contribute N directly; CBOD is an oxygen-demand
      quantity, not a nitrogen reservoir.
    * DOC, POC contribute via their N-equivalent through the C:N ratio
      *only if* the user wants that included. By default this helper
      treats POC/DOC as carbon reservoirs only (no N coupling).
    """
    pieces: list[xr.DataArray] = []
    for name in ("ammonium", "nitrate", "organic_nitrogen", "n2"):
        da = _get(registry, name)
        if da is not None:
            pieces.append(_sum_over_cells(da))

    # Algae N-equivalent
    ap = _get(registry, "ap")
    if ap is not None:
        pieces.append(_sum_over_cells(ap * AP_N_PER_CHLA))
    ab = _get(registry, "ab")
    if ab is not None:
        # g-D/m^2 -> mg-N per cell (per unit area; the Tier 1 test sums
        # cell-sums, so the area factor is constant by construction).
        pieces.append(_sum_over_cells(ab * AB_N_PER_GD))

    if not pieces:
        return xr.DataArray(0.0)
    return sum(pieces[1:], pieces[0])


def total_p(registry: Any) -> xr.DataArray:
    """Total phosphorus mass (mg-P) summed across all reservoirs.

    Reservoirs included (when present):
    * TIP, OrgP — direct mg-P/L pools
    * Ap * AWp, Ab * BWp — algal P-equivalent
    """
    pieces: list[xr.DataArray] = []
    for name in ("tip", "organic_phosphorus"):
        da = _get(registry, name)
        if da is not None:
            pieces.append(_sum_over_cells(da))

    ap = _get(registry, "ap")
    if ap is not None:
        pieces.append(_sum_over_cells(ap * AP_P_PER_CHLA))
    ab = _get(registry, "ab")
    if ab is not None:
        pieces.append(_sum_over_cells(ab * AB_P_PER_GD))

    if not pieces:
        return xr.DataArray(0.0)
    return sum(pieces[1:], pieces[0])


def total_c(registry: Any) -> xr.DataArray:
    """Total carbon mass (mg-C) summed across all reservoirs.

    Reservoirs included (when present):
    * POC, DOC, DIC -- direct mg-C/L pools
    * Ap * AWc, Ab * BWc -- algal C-equivalent
    * POM contributes if a C:N or C:dry-weight ratio is defined; for
      Phase 1.4 POM contributes zero (Phase 3's POM Process will pin
      the convention).
    * CBOD contributes via 1 mg-CBOD == (1/DOC_O2_PER_C) mg-C; this is
      conservative under the closed-system test where CBOD oxidation
      converts CBOD to DIC at exactly that ratio.
    """
    pieces: list[xr.DataArray] = []
    for name in ("poc", "doc", "dic"):
        da = _get(registry, name)
        if da is not None:
            pieces.append(_sum_over_cells(da))

    ap = _get(registry, "ap")
    if ap is not None:
        pieces.append(_sum_over_cells(ap * AP_C_PER_CHLA))
    ab = _get(registry, "ab")
    if ab is not None:
        pieces.append(_sum_over_cells(ab * AB_C_PER_GD))

    cbod = _get(registry, "cbod")
    if cbod is not None:
        pieces.append(_sum_over_cells(cbod / DOC_O2_PER_C))

    if not pieces:
        return xr.DataArray(0.0)
    return sum(pieces[1:], pieces[0])


def total_o2_equivalents(registry: Any) -> xr.DataArray:
    """Total oxygen-equivalent mass (mg-O2) summed across all reservoirs.

    Reservoirs included (when present):
    * DOX -- direct mg-O2/L pool
    * CBOD -- 1 mg-CBOD == 1 mg-O2 by definition
    * DOC -- 1 mg-C * 32/12 mg-O2/mg-C == DOC's stoichiometric oxygen
      demand under DOC -> DIC oxidation
    * POC -- treated identically to DOC for closed-system bookkeeping
    """
    pieces: list[xr.DataArray] = []
    dox = _get(registry, "oxygen_dissolved")
    if dox is not None:
        pieces.append(_sum_over_cells(dox))

    cbod = _get(registry, "cbod")
    if cbod is not None:
        pieces.append(_sum_over_cells(cbod * CBOD_O2_PER_MASS))

    for c_name in ("doc", "poc"):
        da = _get(registry, c_name)
        if da is not None:
            pieces.append(_sum_over_cells(da * DOC_O2_PER_C))

    if not pieces:
        return xr.DataArray(0.0)
    return sum(pieces[1:], pieces[0])


def total_alkalinity(registry: Any) -> xr.DataArray:
    """Total alkalinity (mg-CaCO3) summed across cells.

    v3 1.0.0 treats Alk as a simple tracer with source/sink terms (no
    carbonate solver). Closed-system Tier 1 alkalinity invariance is
    therefore "Alk pool sum is constant"; no derived contributions
    enter at this stage.
    """
    alk = _get(registry, "alkalinity")
    if alk is None:
        return xr.DataArray(0.0)
    return _sum_over_cells(alk)


__all__ = [
    "ClosedSystemConfig",
    "InMemoryRegistry",
    "AP_N_PER_CHLA",
    "AP_P_PER_CHLA",
    "AP_C_PER_CHLA",
    "AB_N_PER_GD",
    "AB_P_PER_GD",
    "AB_C_PER_GD",
    "POM_C_TO_N",
    "POM_C_TO_P",
    "CBOD_O2_PER_MASS",
    "DOC_O2_PER_C",
    "total_n",
    "total_p",
    "total_c",
    "total_o2_equivalents",
    "total_alkalinity",
]
