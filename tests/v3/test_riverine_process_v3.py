"""Unit tests for the v3 ``Riverine`` process MeshView bridge.

Covers the chunk-safe re-bridge added in the riverine MeshView-compat
change (``design/clearwater_modules_v3_riverine_process_meshview_compat.md``):

- MeshView item access (``mesh["Ap"]``) replaces attribute access.
- The constituent gate is widened: no ``has_process("FloatingAlgae")``
  guard; present constituents bridge unconditionally.
- The inorganic-P name is reconciled to ``tip`` (was
  ``phosphorus_total_inorganic``).
- ``depth`` is bridged from ``mesh["coupling_depth"]`` (the resolved
  cell mean water-column depth, a length), not ``wetted_surface_area``
  (area). ``init_process`` first calls
  ``riverine_instance.enable_coupling_depth()`` to turn that resolved
  depth on for the coupled run.
- The bridge is re-applied each substep so chunk reloads (which
  re-register FRESH DataArrays) are picked up rather than stranded on
  the previous chunk's buffers.

These tests build a real ``MeshView`` over a hand-built
``VariableRegistry`` and a stub ``riverine_instance`` -- no HEC-RAS mesh
or transport solve is needed.
"""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_data.variables import VariableRegistry, DataArrayVariable
from clearwater_riverine.fork_compat import MeshView

from clearwater_modules_v3.processes.riverine import Riverine


# Distinguishable magnitudes per constituent so an accidental cross-wire
# (e.g. depth aliasing wetted_surface_area) is caught, not masked.
_CHUNK1_VALUES = {
    "Ap": [1.0, 1.1],
    "NH4": [2.0, 2.1],
    "NO3": [3.0, 3.1],
    "TIP": [4.0, 4.1],
    "DOX": [5.0, 5.1],
    "coupling_depth": [0.5, 0.6],
    "volume": [100.0, 110.0],
    "water_temperature": [20.0, 21.0],
    "wetted_surface_area": [200.0, 220.0],  # area: deliberately != depth
}


def _da(values: list[float]) -> xr.DataArray:
    """A 1-D nface DataArray (no time dim) for a stub constituent."""
    return xr.DataArray(np.asarray(values, dtype=float), dims=["nface"])


def _build_registry(names: list[str]) -> VariableRegistry:
    registry = VariableRegistry()
    for name in names:
        registry.register(name, DataArrayVariable(_da(_CHUNK1_VALUES[name])))
    return registry


class _StubRiverineInstance:
    """Minimal stand-in for ``ClearwaterRiverine``.

    Only ``.mesh`` (a real MeshView), ``.is_chunked`` and
    ``enable_coupling_depth()`` are exercised by these unit tests;
    ``update()`` / ``finalize()`` are not called.

    The real ``enable_coupling_depth()`` enables + seed-computes the
    resolved coupling depth and registers it under ``'coupling_depth'``.
    The stub's registry already carries ``coupling_depth``, so the stub
    only records that the coupling hook was invoked.
    """

    def __init__(self, registry: VariableRegistry) -> None:
        self.mesh = MeshView(registry)
        self.is_chunked = False
        self.enable_coupling_depth_calls = 0

    def enable_coupling_depth(self) -> None:
        self.enable_coupling_depth_calls += 1


class _DummyModel:
    """Stub model. ``init_process`` keeps the ``model`` param for interface
    compatibility but no longer reads it; this just satisfies the call."""

    def has_process(self, process_type) -> bool:  # pragma: no cover - defensive
        return True


def _make_process(registry: VariableRegistry) -> Riverine:
    return Riverine(_StubRiverineInstance(registry))


def _read(registry: VariableRegistry, name: str) -> np.ndarray:
    return registry.get_variable(name).get().values


# ---------------------------------------------------------------------------
# Unit: full constituent set bridges to canonical names with correct depth.
# ---------------------------------------------------------------------------


def test_init_process_bridges_all_constituents_and_depth():
    registry = _build_registry(
        [
            "Ap",
            "NH4",
            "NO3",
            "TIP",
            "DOX",
            "coupling_depth",
            "volume",
            "water_temperature",
            "wetted_surface_area",
        ]
    )
    process = _make_process(registry)

    process.init_process(_DummyModel(), registry)

    # init_process turns on resolved-depth computation for the coupled run.
    assert process.riverine_instance.enable_coupling_depth_calls == 1

    # Canonical aliases registered.
    for canonical in (
        "algae_floating",
        "ammonium",
        "nitrate",
        "tip",
        "oxygen_dissolved",
        "depth",
    ):
        assert canonical in registry, f"{canonical!r} not bridged"

    # Inorganic-P reconciled to ``tip``; old name absent.
    assert "tip" in registry
    assert "phosphorus_total_inorganic" not in registry

    # Values track the source mesh constituents.
    np.testing.assert_array_equal(_read(registry, "algae_floating"), [1.0, 1.1])
    np.testing.assert_array_equal(_read(registry, "tip"), [4.0, 4.1])

    # depth aliases coupling_depth (length), NOT wetted_surface_area (area).
    np.testing.assert_array_equal(_read(registry, "depth"), [0.5, 0.6])
    np.testing.assert_array_equal(
        _read(registry, "depth"), _read(registry, "coupling_depth")
    )
    assert not np.array_equal(
        _read(registry, "depth"), _read(registry, "wetted_surface_area")
    )


# ---------------------------------------------------------------------------
# Subset: only constituents present in the mesh bridge.
# ---------------------------------------------------------------------------


def test_init_process_bridges_only_present_constituents():
    # Omit Ap.
    registry = _build_registry(
        [
            "NH4",
            "NO3",
            "TIP",
            "DOX",
            "coupling_depth",
            "volume",
            "water_temperature",
            "wetted_surface_area",
        ]
    )
    process = _make_process(registry)

    process.init_process(_DummyModel(), registry)

    assert "algae_floating" not in registry  # Ap absent -> no bridge
    assert "ammonium" in registry
    assert "nitrate" in registry
    assert "tip" in registry
    assert "oxygen_dissolved" in registry
    assert "depth" in registry


# ---------------------------------------------------------------------------
# Missing depth: fail loud.
# ---------------------------------------------------------------------------


def test_init_process_raises_when_coupling_depth_missing():
    registry = _build_registry(
        [
            "Ap",
            "NH4",
            "NO3",
            "TIP",
            "DOX",
            "volume",
            "water_temperature",
            "wetted_surface_area",
        ]
    )
    process = _make_process(registry)

    with pytest.raises(KeyError, match="coupling_depth"):
        process.init_process(_DummyModel(), registry)


# ---------------------------------------------------------------------------
# Chunk-safety: the key regression test for the stale-after-chunk-1 bug.
# ---------------------------------------------------------------------------


def test_rebridge_picks_up_new_chunk_objects():
    registry = _build_registry(
        [
            "Ap",
            "NH4",
            "NO3",
            "TIP",
            "DOX",
            "coupling_depth",
            "volume",
            "water_temperature",
            "wetted_surface_area",
        ]
    )
    process = _make_process(registry)

    # Seed t0 (chunk-1).
    process.init_process(_DummyModel(), registry)
    np.testing.assert_array_equal(_read(registry, "algae_floating"), [1.0, 1.1])
    np.testing.assert_array_equal(_read(registry, "depth"), [0.5, 0.6])

    # Simulate a chunk reload: riverine's per-chunk refresh re-registers
    # FRESH DataArray objects for the constituents and coupling_depth.
    registry.register(
        "Ap", DataArrayVariable(_da([91.0, 92.0])), overwrite=True
    )
    registry.register(
        "coupling_depth", DataArrayVariable(_da([9.5, 9.6])), overwrite=True
    )

    # Without the re-bridge the canonical aliases would still point at the
    # stranded chunk-1 buffers.
    np.testing.assert_array_equal(_read(registry, "algae_floating"), [1.0, 1.1])
    np.testing.assert_array_equal(_read(registry, "depth"), [0.5, 0.6])

    # run() re-bridges after update(); exercise the helper directly.
    process._bridge_mesh_to_registry(registry)

    np.testing.assert_array_equal(_read(registry, "algae_floating"), [91.0, 92.0])
    np.testing.assert_array_equal(_read(registry, "depth"), [9.5, 9.6])


# ---------------------------------------------------------------------------
# End-to-end integration: out of scope (depends on Change A + a riverine
# mesh fixture). Stub kept to document the dependency.
# ---------------------------------------------------------------------------


@pytest.mark.skip(reason="needs Change A + riverine mesh fixture")
def test_init_from_file_coupled_integration():  # pragma: no cover
    """init_from_file with Ap/NH4/NO3/TIP/DOX + water_temperature, Phosphorus
    + FloatingAlgae enabled, riverine first; assert model.run() advances and
    Phosphorus reads finite ``tip`` and a length-scale ``depth``. Also covers
    the chunked-integration variant across a real chunk boundary."""
