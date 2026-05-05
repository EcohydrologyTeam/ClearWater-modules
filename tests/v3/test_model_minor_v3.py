"""Phase R-4/R-5 minor cleanups for the v3 ``Model`` orchestrator and
``processes`` package.

Covers M16, m11, m12, m15, and m18 (MINOR) from
``design/clearwater_modules_v3_review_findings.md``:

- **M16.** The aspirational 5-step "integrator-pattern contract" in the
  ``processes/base.py`` module docstring is demoted to a guideline.
  The numbered contract wording is removed because real v3 processes
  (notably ``Temperature``) follow a per-substep ``delta_state``
  pattern instead.
- **m11.** The unused private helper ``Model.__step_index`` is deleted
  (dead code).
- **m12.** Empty-string and ``Path("")`` arguments to
  ``simulation_directory`` are now treated explicitly as "not
  provided" and fall back to ``Path(".")``; every other value
  (including ``Path(".")`` itself) is preserved.
- **m15.** The ``processes`` parameter annotation on ``Model.__init__``
  is widened from ``tuple[Process, ...]`` to ``Iterable[Process]`` to
  reflect the constructor's actual accept-any-iterable behavior.
- **m18 (MINOR).** The unused ``RUN_ORDER`` constant is removed from
  the ``clearwater_modules_v3.processes`` package surface.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import pytest

import clearwater_modules_v3.processes as v3_processes
from clearwater_modules_v3.model import Model
from clearwater_modules_v3.processes import base as v3_process_base
from clearwater_modules_v3.processes.base import Process


# ---------------------------------------------------------------------------
# Local stubs (mirror the patterns used in test_model_orchestration_v3.py)
# ---------------------------------------------------------------------------


class _StubRegistry:
    def __init__(self) -> None:
        self._data: dict[str, object] = {}

    def register(self, name: str, value) -> None:
        self._data[name] = value

    def get(self, name: str):
        return self._data[name]

    def get_at_time(self, name: str, time):
        return self._data[name]

    def set_at_time(self, name: str, time, value) -> None:
        self._data[name] = value

    def get_variable(self, name: str):
        raise KeyError(name)

    def __contains__(self, name: str) -> bool:
        return name in self._data


class _StubProcess(Process):
    """Minimal process for orchestration smoke-tests."""

    variables: list[str] = []

    def __init__(self, name: str = "stub", time_step: timedelta = timedelta(minutes=5)) -> None:
        super().__init__(time_step)
        self._name = name

    def run(self, time, registry) -> None:  # pragma: no cover - not exercised
        pass

    def process_name(self) -> str:
        return self._name


def _build_model(
    processes,
    *,
    simulation_directory=None,
) -> Model:
    return Model(
        processes=processes,
        variable_registry=_StubRegistry(),
        variable_data_sources={},
        start_time=datetime(2026, 1, 1, 0, 0, 0),
        end_time=datetime(2026, 1, 1, 0, 30, 0),
        time_step=timedelta(minutes=5),
        output_variables=[],
        simulation_directory=simulation_directory,
        chunk_size=None,
    )


# ---------------------------------------------------------------------------
# M16: integrator-pattern contract demoted to guideline
# ---------------------------------------------------------------------------


def test_m16_process_base_docstring_no_longer_lists_5_step_contract():
    """The numbered 5-step "integrator-pattern contract" wording in the
    ``processes/base.py`` module docstring is removed; the docstring
    now frames Forward Euler and per-substep ``delta_state`` as
    co-equal patterns instead."""
    doc = v3_process_base.__doc__ or ""
    # The original numbered contract used "Forward Euler" inside an
    # enumerated list. Verify the *prescriptive* phrasing is gone.
    assert "integrator-pattern contract" not in doc
    # The replaced docstring frames both options; verify the new wording
    # is present so a future regression that re-introduces the contract
    # is caught.
    assert "Two integration patterns are common" in doc
    assert "delta_state" in doc
    # The M5 hotstart-ordering contract must remain intact.
    assert "M5 ordering contract" in doc
    assert "from_hotstart" in doc


# ---------------------------------------------------------------------------
# m11: __step_index helper deleted
# ---------------------------------------------------------------------------


def test_m11_step_index_helper_is_deleted():
    """The unused private ``__step_index`` helper has been removed from
    the Model class. Verify by checking the name-mangled attribute is
    not present at the class level."""
    assert not hasattr(Model, "_Model__step_index")


# ---------------------------------------------------------------------------
# m12: empty-string and Path("") simulation_directory handling
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "sim_dir",
    [None, "", Path("")],
    ids=["none", "empty_str", "empty_path"],
)
def test_m12_empty_simulation_directory_falls_back_to_cwd(sim_dir):
    """``None``, ``""``, and ``Path("")`` are all treated as "not
    provided" and resolve to ``Path(".")``."""
    model = _build_model(processes=[], simulation_directory=sim_dir)
    assert isinstance(model._Model__simulation_directory, Path)
    assert model._Model__simulation_directory == Path(".")


def test_m12_explicit_dot_path_is_preserved():
    """``Path(".")`` is a valid explicit user choice (not falsy under
    the new predicate) and must be preserved."""
    model = _build_model(processes=[], simulation_directory=Path("."))
    assert isinstance(model._Model__simulation_directory, Path)
    assert model._Model__simulation_directory == Path(".")


def test_m12_explicit_nonempty_path_is_preserved():
    """Sanity: a real user path is not rewritten."""
    model = _build_model(
        processes=[], simulation_directory=Path("/tmp/v3_m12_check")
    )
    assert model._Model__simulation_directory == Path("/tmp/v3_m12_check")


# ---------------------------------------------------------------------------
# m15: processes annotation widened to Iterable[Process]
# ---------------------------------------------------------------------------


def test_m15_processes_accepts_plain_list_and_stores_tuple():
    """Pass a plain ``list`` (not a ``tuple``) for ``processes``;
    construction succeeds and the internal attribute is normalized to
    a ``tuple``."""
    stub_a = _StubProcess("a")
    stub_b = _StubProcess("b")
    model = _build_model(processes=[stub_a, stub_b])  # plain list, not tuple
    assert isinstance(model._Model__processes, tuple)
    assert model._Model__processes == (stub_a, stub_b)


def test_m15_processes_accepts_generator():
    """Generators are also iterables; construction must succeed."""
    stubs = [_StubProcess("g0"), _StubProcess("g1")]
    model = _build_model(processes=(p for p in stubs))
    assert isinstance(model._Model__processes, tuple)
    assert model._Model__processes == tuple(stubs)


# ---------------------------------------------------------------------------
# m18 (MINOR): RUN_ORDER no longer exported from processes package
# ---------------------------------------------------------------------------


def test_m18_minor_run_order_not_exposed():
    """``RUN_ORDER`` is removed from the ``clearwater_modules_v3.processes``
    package surface (constant deleted; ``__all__`` entry removed)."""
    assert not hasattr(v3_processes, "RUN_ORDER")
    assert "RUN_ORDER" not in getattr(v3_processes, "__all__", [])
