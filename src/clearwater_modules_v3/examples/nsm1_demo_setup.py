"""Headless NSM1 demo construction helper.

Phase 8.A (v3 NSM1 implementation plan): a programmatic builder that
constructs all 11 v3 NSM1 Process classes against a synthetic in-memory
registry, wires them together via a minimal ``Model`` stand-in, and
returns a runnable handle. The helper closes the gap between the v3
``init_from_file`` YAML entry point (which requires a Riverine mesh or
fully-specified ``data_sources`` mapping) and the NSM1 unit tests
(which build their own ad-hoc fixtures). It is the recommended
"NSM1-only, no transport, no I/O" entry point.

Why this exists
---------------

The Phase 7.C demo notebook (`examples/V3/04_Example_NSM1.ipynb`)
manually constructed a 5-cell synthetic mesh, instantiated each
Process class with parameter overrides, set five legacy attributes on
the v2 ``Nitrogen`` overlay that the YAML config path normally
injects, mirrored the Phosphorus state into the legacy
``phosphorus_total_inorganic`` registry name, and built a private
``_DemoModel`` stand-in to satisfy the ``init_process(model, registry)``
hook. That construction was correct but fragile — every consumer had
to repeat the same five workarounds.

Phase 8.A (this module) consolidates the construction into a single
helper. The legacy-attribute and state-name workarounds are no longer
needed because:

* v2 ``Nitrogen.__init__`` now seeds the five legacy attributes from
  the v3 ``ALGAE_DEFAULTS`` / ``BALGAE_DEFAULTS`` (``Fix 1``).
* v2 ``FloatingAlgae.run`` and ``BenthicAlgae.run`` now read the
  canonical v3 inorganic-P state name ``tip`` (``Fix 2``).

Public API
----------

* :class:`InMemoryRegistry` — minimal ``VariableRegistry`` stand-in
  supporting the subset of the API the Process classes call
  (``register``, ``get``, ``get_at_time``, ``set_at_time``,
  ``__contains__``, ``keys``).
* :func:`default_initial_conditions` — returns a fresh dict mapping
  every NSM1 state variable + forcing to a 5-cell ``xr.DataArray`` with
  physically reasonable mesotrophic-stream values.
* :func:`default_process_parameters` — returns a fresh dict mapping
  each Process class name to its parameter-override dict. Every
  Process gets ``{}`` (use v3 defaults) except ``DOX``, which enables
  hydraulic reaeration so the kinetics run in physical equilibrium.
* :func:`build_nsm1_demo` — constructs the registry, the 11 Process
  instances, the firing schedule, and a minimal Model stand-in;
  returns a :class:`Nsm1Demo` handle with a ``step()`` method.

Usage
-----

::

    from datetime import datetime, timedelta
    from clearwater_modules_v3.examples import build_nsm1_demo

    demo = build_nsm1_demo(time_step=timedelta(minutes=5))
    t = datetime(2026, 1, 1)
    for _ in range(12):  # 1 hour at 5-min substep
        demo.step(t)
        t += demo.time_step

    # Inspect any state variable on the registry.
    final_dox = demo.registry.get("oxygen_dissolved")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Any, Iterable

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.alkalinity import Alkalinity
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae
from clearwater_modules_v3.processes.carbon import Carbon
from clearwater_modules_v3.processes.cbod import CBOD
from clearwater_modules_v3.processes.dox import DOX
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.n2 import N2
from clearwater_modules_v3.processes.nitrogen import Nitrogen
from clearwater_modules_v3.processes.pathogen import Pathogen
from clearwater_modules_v3.processes.phosphorus import Phosphorus
from clearwater_modules_v3.processes.pom import POM
from clearwater_modules_v3.utils.numerics import Diagnostics

# Default mesh dimension for the synthetic demo. Five cells matches the
# Tier 1 conservation-test fixture in tests/v3/nsm1/conftest.py so the
# demo and the test suite produce comparable trajectories.
N_CELLS_DEFAULT = 5


# ---------------------------------------------------------------------------
# In-memory registry stand-in
# ---------------------------------------------------------------------------


class InMemoryRegistry:
    """Minimal ``VariableRegistry`` stand-in for headless demos.

    Implements just the subset of the registry API that the Process
    classes call: ``register``, ``get``, ``get_at_time``, ``set_at_time``,
    ``__contains__``, ``keys``. Mirrors the fixture in
    ``tests/v3/nsm1/conftest.py`` so demos and unit tests use the same
    shape.
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

    def get_variable(self, name: str) -> Any:
        # The Tier 1 conftest raises KeyError; the v3 Model only calls
        # ``get_variable`` for output-variable bookkeeping, which is
        # not exercised by this demo.
        raise KeyError(name)

    def __contains__(self, name: str) -> bool:
        return name in self._data

    def keys(self) -> Iterable[str]:
        return self._data.keys()

    def snapshot(self) -> dict[str, xr.DataArray]:
        """Return a deep copy of every registry variable.

        Useful for time-series tracking in demo loops without having to
        re-derive each variable's name list.
        """
        return {k: v.copy() for k, v in self._data.items()}


# ---------------------------------------------------------------------------
# Default initial conditions and parameters
# ---------------------------------------------------------------------------


def _da(values: list[float], dim: str = "cell") -> xr.DataArray:
    return xr.DataArray(np.asarray(values, dtype=float), dims=dim)


def default_initial_conditions(n_cells: int = N_CELLS_DEFAULT) -> dict[str, xr.DataArray]:
    """Return a fresh dict of NSM1 initial conditions on an ``n_cells`` mesh.

    Values are physically reasonable for a mesotrophic stream and
    overlap with the Tier 1 conservation-test fixture in
    ``tests/v3/nsm1/conftest.py`` so demo runs and unit tests are
    cross-comparable. Forcings are encoded as cell-aligned constants so
    the demo can reuse the same registry plumbing as the state
    variables.
    """

    def const(value: float) -> xr.DataArray:
        return _da([value] * n_cells)

    return {
        # --- Nitrogen ---
        "ammonium":           const(0.05),    # mg-N/L
        "nitrate":            const(5.0),     # mg-N/L
        "organic_nitrogen":   const(1.726),   # mg-N/L
        "n2":                 const(1.0),     # mg-N/L
        # --- Phosphorus (canonical v3 name "tip") ---
        "tip":                const(0.07),    # mg-P/L
        "organic_phosphorus": const(0.24),    # mg-P/L
        # --- Carbon ---
        "poc":                const(4.0),     # mg-C/L
        "doc":                const(1.0),     # mg-C/L
        "dic":                const(1.0),     # mg-C/L
        # --- Particulate organic matter ---
        "pom":                const(10.0),    # mg/L
        # --- CBOD (single group) ---
        "cbod":               const(5.0),     # mg-O2/L
        # --- DOX ---
        "oxygen_dissolved":   const(8.0),     # mg-O2/L
        # --- Pathogen ---
        "pathogen":           const(1.0),     # count/100mL (relative)
        # --- Alkalinity ---
        "alkalinity":         const(1.0),     # mg-CaCO3/L (small for visibility)
        # --- Algae ---
        "algae_floating":     const(40.0),    # ug-Chla/L
        "benthic_algae":      const(24.0),    # g-D/m^2
        # --- Forcings ---
        "water_temperature":   const(20.0),    # deg C
        "depth":               const(1.5),     # m
        "solar_radiation":     const(250.0),   # W/m^2
        "atmospheric_pressure": const(1013.25), # mb
    }


def default_process_parameters() -> dict[str, dict[str, Any]]:
    """Return per-Process parameter overrides for the headless demo.

    All Process classes default to v3 ``DEFAULTS`` (empty override
    dict). The exception is ``DOX``, which enables hydraulic
    reaeration with ``kah_20_user=20.0 1/d`` so air-water exchange can
    balance algal photosynthesis at the demo's solar forcing. Without
    that override, DOX would accumulate non-physically because the v3
    default ``kah_20_user=0.0`` disables reaeration.
    """
    return {
        "FloatingAlgae": {},
        "BenthicAlgae":  {},
        "Nitrogen":      {},
        "Phosphorus":    {},
        "Carbon":        {},
        "POM":           {},
        "CBOD":          {},
        "DOX": {
            "pressure_mb": 1013.25,
            "kah_20_user": 20.0,
            "hydraulic_reaeration_option": 1,
        },
        "N2": {
            "pressure_mb": 1013.25,
        },
        "Pathogen":   {},
        "Alkalinity": {},
    }


# ---------------------------------------------------------------------------
# Minimal Model stand-in
# ---------------------------------------------------------------------------


class _DemoModel:
    """Minimal ``Model`` stand-in for ``init_process(model, registry)``.

    The real v3 ``Model`` carries a process schedule, registry, output
    store, etc. The Process ``init_process`` hooks query just two
    methods on the model: ``has_process`` and ``get_process``. Plus a
    ``diagnostics`` attribute for clip-with-log routing. This stand-in
    matches the same surface so demos avoid the YAML data-source
    plumbing.
    """

    def __init__(self, processes: list) -> None:
        self._processes = list(processes)
        self.diagnostics = Diagnostics()

    def has_process(self, name_or_type) -> bool:
        if isinstance(name_or_type, str):
            return any(type(p).__name__ == name_or_type for p in self._processes)
        return any(isinstance(p, name_or_type) for p in self._processes)

    def get_process(self, name_or_type):
        if isinstance(name_or_type, str):
            return next(p for p in self._processes if type(p).__name__ == name_or_type)
        return next(p for p in self._processes if isinstance(p, name_or_type))


# ---------------------------------------------------------------------------
# Public builder
# ---------------------------------------------------------------------------


@dataclass
class Nsm1Demo:
    """Handle returned by :func:`build_nsm1_demo`.

    Attributes:
        registry: The in-memory registry holding state and forcings.
        processes: The 11 Process instances in firing order. Producers
            (FloatingAlgae, BenthicAlgae, Nitrogen, POM, CBOD) precede
            consumers (Carbon, Phosphorus, N2, DOX, Pathogen,
            Alkalinity) so step-scoped rate caches are populated
            before they are read.
        model: The minimal ``_DemoModel`` stand-in passed to each
            Process's ``init_process`` hook.
        time_step: Per-Process substep length (every Process fires
            every step in this demo).
    """

    registry: InMemoryRegistry
    processes: list
    model: _DemoModel
    time_step: timedelta

    def step(self, time: datetime) -> None:
        """Advance every Process by one substep at ``time``."""
        for process in self.processes:
            process.run(time, self.registry)

    def run(
        self,
        start: datetime,
        n_steps: int,
        snapshot_every: int | None = None,
    ) -> list[dict[str, xr.DataArray]]:
        """Advance the demo by ``n_steps`` substeps.

        Args:
            start: Initial wall-clock time. Mostly cosmetic for the
                in-memory registry (``get_at_time`` ignores its
                argument), but the value is forwarded to each
                Process's ``run`` for any process that observes it.
            n_steps: Number of substeps to execute.
            snapshot_every: If given, capture a deep-copy of the
                registry every ``snapshot_every`` substeps and return
                the list. Defaults to ``None`` (no snapshots; returns
                an empty list).

        Returns:
            A list of registry snapshots (each a dict mapping name to
            ``xr.DataArray``). Empty list when ``snapshot_every`` is
            ``None``.
        """
        snapshots: list[dict[str, xr.DataArray]] = []
        t = start
        for i in range(n_steps):
            self.step(t)
            t += self.time_step
            if snapshot_every is not None and (i + 1) % snapshot_every == 0:
                snapshots.append(self.registry.snapshot())
        return snapshots


def build_nsm1_demo(
    time_step: timedelta = timedelta(minutes=5),
    n_cells: int = N_CELLS_DEFAULT,
    initial_conditions: dict[str, xr.DataArray] | None = None,
    process_parameters: dict[str, dict[str, Any]] | None = None,
) -> Nsm1Demo:
    """Build a runnable headless v3 NSM1 demo.

    Constructs a 5-cell synthetic mesh by default, instantiates all 11
    NSM1 Process classes with v3 ``DEFAULTS`` (plus the DOX reaeration
    override), wires them through the ``init_process(model, registry)``
    hook on a minimal ``_DemoModel`` stand-in, and returns a
    :class:`Nsm1Demo` handle. Callers can step the demo via
    ``demo.step(t)`` or ``demo.run(start, n_steps)``.

    Args:
        time_step: Per-Process substep length. Every Process fires
            every step. Default: 5 minutes.
        n_cells: Number of mesh cells in the default initial
            conditions. Ignored when ``initial_conditions`` is
            supplied. Default: 5.
        initial_conditions: Optional override dict mapping registry
            variable names to ``xr.DataArray`` values. When supplied,
            ``n_cells`` is ignored and the caller is responsible for
            providing every state and forcing variable the 11
            Processes read.
        process_parameters: Optional override dict mapping Process
            class names (``"FloatingAlgae"``, ``"Nitrogen"``, etc.) to
            parameter-override dicts. Missing class names fall back to
            :func:`default_process_parameters`. Pass ``{"DOX": {}}``
            to disable the hydraulic-reaeration override.

    Returns:
        :class:`Nsm1Demo` with the populated registry, the 11 Process
        instances in firing order, the ``_DemoModel`` stand-in, and
        the time step.
    """
    # --- Registry ---
    registry = InMemoryRegistry()
    ic = (
        dict(initial_conditions)
        if initial_conditions is not None
        else default_initial_conditions(n_cells=n_cells)
    )
    for name, da in ic.items():
        registry.register(name, da)

    # --- Process parameters ---
    base_params = default_process_parameters()
    user_params = dict(process_parameters) if process_parameters else {}
    for k, v in user_params.items():
        # Caller-supplied Process keys override the defaults wholesale
        # (no nested merge); this matches the v3 ``Process.DEFAULTS``
        # semantics where the caller passes a complete parameter dict.
        base_params[k] = v

    # --- Processes ---
    floating_algae = FloatingAlgae(parameters=base_params["FloatingAlgae"], time_step=time_step)
    benthic_algae  = BenthicAlgae(parameters=base_params["BenthicAlgae"],   time_step=time_step)
    nitrogen       = Nitrogen(parameters=base_params["Nitrogen"],           time_step=time_step)
    phosphorus     = Phosphorus(parameters=base_params["Phosphorus"],       time_step=time_step)
    carbon         = Carbon(parameters=base_params["Carbon"],               time_step=time_step)
    pom            = POM(parameters=base_params["POM"],                     time_step=time_step)
    cbod           = CBOD(parameters=base_params["CBOD"],                   time_step=time_step)
    dox            = DOX(parameters=base_params["DOX"],                     time_step=time_step)
    n2             = N2(parameters=base_params["N2"],                       time_step=time_step)
    pathogen       = Pathogen(parameters=base_params["Pathogen"],           time_step=time_step)
    alkalinity     = Alkalinity(parameters=base_params["Alkalinity"],       time_step=time_step)

    # --- Firing order ---
    # Producers fire before consumers so step-scoped rate caches are
    # populated before they are read.
    process_order = [
        floating_algae,    # produces algal_*_rate caches
        benthic_algae,     # produces balgae_*_rate caches
        nitrogen,          # produces nitrification_flux_rate, denitrification_flux_rate
        pom,               # consumes algae mortality
        cbod,              # produces cbod_oxidation_rate
        carbon,            # consumes pom + algae rates
        phosphorus,        # consumes algae rates
        n2,                # consumes nitrogen denitrification rate
        dox,               # consumes everything
        pathogen,          # standalone
        alkalinity,        # consumes nitrogen + algae rates
    ]

    # --- Model stand-in + init_process wiring ---
    model = _DemoModel(process_order)
    for process in process_order:
        process.init_process(model, registry)

    return Nsm1Demo(
        registry=registry,
        processes=process_order,
        model=model,
        time_step=time_step,
    )


__all__ = [
    "InMemoryRegistry",
    "Nsm1Demo",
    "build_nsm1_demo",
    "default_initial_conditions",
    "default_process_parameters",
    "N_CELLS_DEFAULT",
]
