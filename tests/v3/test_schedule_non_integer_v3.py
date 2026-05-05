"""M18 non-integer-second timestep test (review-findings 2026-05-04).

M18 was the deferred concern that
``__build_process_schedule``'s ``current_seconds % interval == 0.0``
test is exact for integer-second intervals but is exposed to
floating-point drift when ``time_step`` is a sub-second
``timedelta``. The C6 fix added an up-front "process.time_step_seconds
must be an integer multiple of model time_step_seconds" validation,
which closes M18 in the integer-second regime. This file pins the
actual behavior at sub-second cadences and documents the residual
gap.

Stub patterns mirror those in
``tests/v3/test_model_orchestration_v3.py``; copies are inlined.

Scenarios:

1. Integer-second non-divisor: ``time_step=5min``, process=7s.
   ``init_model()`` raises ``ValueError`` mentioning "integer multiple".
2. Sub-second 3:1 ratio: ``time_step=0.1s``, process=0.3s. Mathematically
   an integer multiple, BUT in IEEE-754 ``0.3 % 0.1`` is
   ``0.09999999999999998 != 0``. Documents the residual M18 gap: the
   current implementation incorrectly rejects this configuration.
3. Sub-second 2.5:1 ratio: ``time_step=0.1s``, process=0.25s. Genuinely
   not an integer multiple; rejected as expected.
4. Edge: ``time_step=0s``. The current implementation raises
   ``ZeroDivisionError`` (not ``ValueError``) because
   ``__count_substeps`` divides by ``time_step.total_seconds()`` before
   the cadence validation runs. Pins this behavior so a future change
   that adds an explicit zero-timestep guard surfaces visibly.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import pytest

from clearwater_modules_v3.model import Model
from clearwater_modules_v3.processes.base import Process


# ---------------------------------------------------------------------------
# Stubs (inlined from test_model_orchestration_v3.py)
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
    """Minimal Process: configurable per-process ``time_step``."""

    variables: list[str] = []

    def __init__(
        self, name: str, time_step: timedelta = timedelta(minutes=5)
    ) -> None:
        super().__init__(time_step)
        self._name = name

    def run(self, time, registry) -> None:
        pass

    def process_name(self) -> str:
        return self._name


def _build_model(
    processes,
    *,
    time_step: timedelta,
    end_time: datetime,
    start_time: datetime = datetime(2026, 1, 1, 0, 0, 0),
) -> Model:
    return Model(
        processes=tuple(processes),
        variable_registry=_StubRegistry(),
        variable_data_sources={},
        start_time=start_time,
        end_time=end_time,
        time_step=time_step,
        output_variables=[],
    )


# ---------------------------------------------------------------------------
# Scenario 1: integer-second non-divisor (5min model, 7s process)
# ---------------------------------------------------------------------------


def test_integer_second_non_divisor_rejected_by_init_model():
    """``time_step=5min`` (300s) and a process with ``time_step=7s``
    (420s/60 = 7s; 7 is not a divisor of 300). ``init_model()``
    invokes ``__build_process_schedule``, which raises ``ValueError``
    mentioning "integer multiple"."""
    bad = _StubProcess("p_7sec", time_step=timedelta(seconds=7))
    model = _build_model(
        processes=[bad],
        time_step=timedelta(minutes=5),
        end_time=datetime(2026, 1, 1, 1, 0, 0),
    )
    with pytest.raises(ValueError, match="integer multiple"):
        model.init_model()


def test_integer_second_non_divisor_rejected_by_build_schedule_directly():
    """Same configuration as above, but exercising the private
    ``__build_process_schedule`` directly. Mirrors the pattern used by
    ``test_c6_schedule_validation_rejects_non_divisor_cadence`` in
    ``tests/v3/test_model_orchestration_v3.py``."""
    bad = _StubProcess("p_7sec", time_step=timedelta(seconds=7))
    model = _build_model(
        processes=[bad],
        time_step=timedelta(minutes=5),
        end_time=datetime(2026, 1, 1, 1, 0, 0),
    )
    with pytest.raises(ValueError, match="integer multiple"):
        model._Model__build_process_schedule()


# ---------------------------------------------------------------------------
# Scenario 2: sub-second mathematically-integer ratio (0.1s model, 0.3s proc)
# ---------------------------------------------------------------------------


def test_sub_second_three_to_one_ratio_documents_m18_residual_gap():
    """``time_step=0.1s`` and process ``time_step=0.3s``: 0.3 / 0.1 == 3
    in real arithmetic, so this should be a valid 3:1 cadence. In
    IEEE-754 doubles, however, ``0.3 % 0.1 == 0.09999999999999998 !=
    0``, and the cadence-multiple validation in
    ``__build_process_schedule`` rejects it.

    This test pins the residual M18 gap: at sub-second time steps the
    bare ``%`` modulo check misclassifies mathematically-integer
    multiples as non-divisors. The fix is non-trivial (round to a
    tolerance, or refactor to integer microseconds); deferred to a
    future revision. The test is written to fail-loudly if the
    behavior changes — at which point the assertion below should be
    updated to confirm the schedule builds correctly with a 3:1 firing
    cadence.

    If a future fix lands, the expected post-fix behavior is:
    ``schedule = model._Model__build_process_schedule()`` succeeds, and
    the process fires at substeps 0, 3, 6, ... within the run window.
    """
    proc_03 = _StubProcess("p_0.3sec", time_step=timedelta(seconds=0.3))
    model = _build_model(
        processes=[proc_03],
        time_step=timedelta(seconds=0.1),
        end_time=datetime(2026, 1, 1, 0, 0, 1),  # 1s window -> 10 substeps
    )
    # Pin current (M18-residual) behavior: validation rejects.
    with pytest.raises(ValueError, match="integer multiple"):
        model._Model__build_process_schedule()


# ---------------------------------------------------------------------------
# Scenario 3: sub-second genuine non-divisor (0.1s model, 0.25s proc)
# ---------------------------------------------------------------------------


def test_sub_second_two_point_five_to_one_ratio_rejected():
    """``time_step=0.1s`` and process ``time_step=0.25s``: 0.25 / 0.1 ==
    2.5 in real arithmetic, genuinely not an integer multiple. The
    cadence-multiple validation correctly rejects this configuration
    (regardless of float-precision concerns; the rejection is correct
    here even though the residual gap above produces a false positive
    at 3:1)."""
    proc_025 = _StubProcess("p_0.25sec", time_step=timedelta(seconds=0.25))
    model = _build_model(
        processes=[proc_025],
        time_step=timedelta(seconds=0.1),
        end_time=datetime(2026, 1, 1, 0, 0, 1),
    )
    with pytest.raises(ValueError, match="integer multiple"):
        model._Model__build_process_schedule()


# ---------------------------------------------------------------------------
# Scenario 4: zero time_step edge case
# ---------------------------------------------------------------------------


def test_zero_time_step_raises_zero_division_error():
    """``time_step=timedelta(seconds=0)`` is pathological. The
    constructor (``Model.__init__``) does not currently validate that
    ``time_step > 0``; ``validate()`` checks only that
    ``start_time < end_time``. The first place the bad value hits an
    arithmetic op is ``__count_substeps``, which computes
    ``delta_seconds // 0.0`` and raises ``ZeroDivisionError``.

    This test pins that behavior so a future change that adds an
    explicit zero-timestep guard (``ValueError`` would be more
    appropriate than ``ZeroDivisionError`` for the user-facing API)
    surfaces as a deliberate test update rather than a silent change.
    """
    p = _StubProcess("p", time_step=timedelta(minutes=5))
    model = _build_model(
        processes=[p],
        time_step=timedelta(seconds=0),
        end_time=datetime(2026, 1, 1, 0, 0, 1),
    )
    # Either ZeroDivisionError or ValueError is acceptable. Pin the
    # current path (ZeroDivisionError from __count_substeps) and
    # accept ValueError as a future-friendly alternative so the test
    # does not fail merely because someone added a clean guard.
    with pytest.raises((ZeroDivisionError, ValueError)):
        model._Model__build_process_schedule()


def test_zero_time_step_constructor_does_not_validate():
    """Document the gap: ``Model.__init__`` accepts ``time_step=0``
    silently. The error only surfaces later when arithmetic is
    performed. A future improvement is to validate ``time_step >
    timedelta(0)`` in ``__init__`` or ``validate()``."""
    p = _StubProcess("p", time_step=timedelta(minutes=5))
    # No exception at construction time.
    model = _build_model(
        processes=[p],
        time_step=timedelta(seconds=0),
        end_time=datetime(2026, 1, 1, 0, 0, 1),
    )
    # ``validate()`` checks only start_time < end_time; it does not
    # catch the zero time_step.
    model.validate()  # does not raise (gap)
