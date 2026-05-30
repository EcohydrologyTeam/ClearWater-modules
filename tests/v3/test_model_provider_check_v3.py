"""v3 ``Model`` provider-coverage check (Change C).

Design: ``design/clearwater_modules_v3_riverine_process_meshview_compat.md``
(section "Change C").

``Model.__init_model`` runs a one-shot coverage check after every provider
has had its chance to populate the registry — data sources (step 1),
hotstart seed (step 2), ``init_process`` incl. the riverine bridge
(step 3), and ``from_hotstart`` (step 4). Any ``process.variables`` entry
still absent from the registry has no provider and would otherwise surface
as a latent runtime ``KeyError`` mid-substep. The check turns that into a
clear init-time ``KeyError`` naming the variable and the declaring
process(es).

These tests use the same lightweight ``_StubRegistry`` / ``_StubProcess``
construction as ``test_model_robustness_v3.py`` so no real ``Model`` data
source or transport solve is needed.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.model import Model
from clearwater_modules_v3.processes.base import Process


# ---------------------------------------------------------------------------
# Stubs (mirror test_model_robustness_v3.py)
# ---------------------------------------------------------------------------


class _StubRegistry:
    """Minimal VariableRegistry stand-in.

    Exposes ``_data`` (the stub convention the provider check falls back
    to when ``_registry`` is absent) and ``__contains__`` so the coverage
    check can probe membership.
    """

    def __init__(self) -> None:
        self._data: dict[str, object] = {}

    def register(self, name: str, value, overwrite: bool = False) -> None:
        self._data[name] = value

    def get(self, name: str):
        return self._data[name]

    def get_at_time(self, name: str, time):
        if name not in self._data:
            raise KeyError(name)
        return self._data[name]

    def set_at_time(self, name: str, time, value) -> None:
        self._data[name] = value

    def get_variable(self, name: str):
        raise KeyError(name)

    def __contains__(self, name: str) -> bool:
        return name in self._data


class _StubProcess(Process):
    """Process that declares ``variables`` and seeds nothing in init."""

    def __init__(
        self,
        name: str = "stub",
        variables: list[str] | None = None,
        time_step: timedelta = timedelta(minutes=5),
    ) -> None:
        super().__init__(time_step)
        self._name = name
        self.variables = list(variables or [])

    def run(self, time, registry) -> None:  # pragma: no cover - not exercised
        pass

    def process_name(self) -> str:
        return self._name


class _SeedingProcess(_StubProcess):
    """Process that registers its declared inputs in ``init_process``.

    Stands in for a provider (e.g. the riverine bridge) that populates the
    registry during step 3 so the coverage check passes.
    """

    def init_process(self, model, registry) -> None:
        for variable in self.variables:
            registry.register(variable, xr.DataArray(np.zeros(3)))


def _build_model(
    processes,
    *,
    variable_registry=None,
    variable_data_sources=None,
) -> Model:
    return Model(
        processes=tuple(processes),
        variable_registry=variable_registry or _StubRegistry(),
        variable_data_sources=variable_data_sources or {},
        start_time=datetime(2026, 1, 1, 0, 0, 0),
        end_time=datetime(2026, 1, 1, 0, 30, 0),
        time_step=timedelta(minutes=5),
        output_variables=[],
        simulation_directory=None,
        chunk_size=None,
        wet_mask_variable=None,
        wet_mask_threshold=0.0,
        hotstart_dataset=None,
        hotstart_timestep=None,
    )


# ---------------------------------------------------------------------------
# Missing provider -> clear init-time error
# ---------------------------------------------------------------------------


def test_unprovided_declared_input_raises_at_init():
    """A process declaring an input that no provider supplies must raise a
    clear ``KeyError`` at ``init_model`` time naming the variable."""
    p = _StubProcess("phosphorus", variables=["tip", "depth"])
    model = _build_model(processes=[p])
    with pytest.raises(KeyError) as exc_info:
        model.init_model()
    message = str(exc_info.value)
    # Both unprovided declared inputs are named.
    assert "tip" in message
    assert "depth" in message
    # The declaring process is named.
    assert "phosphorus" in message
    # The remediation hint mentions where a provider can come from.
    assert "data source" in message
    assert "bridge" in message


def test_error_names_only_the_missing_variable():
    """When one declared input is provided and another is not, only the
    missing one is reported."""
    seeder = _SeedingProcess("seeder", variables=["water_temperature"])
    consumer = _StubProcess("consumer", variables=["water_temperature", "depth"])
    model = _build_model(processes=[seeder, consumer])
    with pytest.raises(KeyError) as exc_info:
        model.init_model()
    message = str(exc_info.value)
    assert "depth" in message
    # water_temperature was seeded by the seeder -> must NOT be flagged.
    assert "water_temperature" not in message.split("Registered variables:")[0]


# ---------------------------------------------------------------------------
# Fully-provided model -> no false positive
# ---------------------------------------------------------------------------


def test_fully_provided_model_initializes_cleanly():
    """A Model whose declared inputs are all provided initializes without
    error (no false positive from the coverage check)."""
    consumer = _StubProcess("consumer", variables=["tip", "depth"])
    model = _build_model(processes=[consumer])
    # Pre-register the declared inputs (a pre-registration provider).
    model._Model__registry.register("tip", xr.DataArray(np.zeros(3)))
    model._Model__registry.register("depth", xr.DataArray(np.zeros(3)))
    model.init_model()
    assert model._Model__init_complete is True


def test_process_with_no_declared_inputs_initializes_cleanly():
    """A process that declares no inputs never trips the coverage check."""
    p = _StubProcess("noinputs", variables=[])
    model = _build_model(processes=[p])
    model.init_model()
    assert model._Model__init_complete is True


# ---------------------------------------------------------------------------
# Bridge-supplied depth -> covered by the check
# ---------------------------------------------------------------------------


def test_bridge_supplied_depth_passes_check():
    """A process that supplies a declared input during ``init_process``
    (as the riverine bridge supplies ``depth``) satisfies the coverage
    check for a downstream consumer that declares the same input."""
    # The "bridge" seeds depth + tip at init_process time (step 3).
    bridge = _SeedingProcess("riverine", variables=["depth", "tip"])
    # A consumer declares them as inputs; both are provided by the bridge.
    consumer = _StubProcess("phosphorus", variables=["tip", "depth"])
    model = _build_model(processes=[bridge, consumer])
    model.init_model()
    assert model._Model__init_complete is True
