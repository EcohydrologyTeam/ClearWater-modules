"""Phase 10 pattern conformance — single test file that iterates every
v3 NSM1 Process class and asserts the canonical pattern A–J shape from
``design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md``
§3.

Verifies that the refactor template established in Phase 2 (Carbon)
and applied through Phase 9 (Alkalinity) is present uniformly. Failing
a clause here means a future code change drifted away from the canon.

Retained indefinitely (not deleted in any cleanup phase). This is the
load-bearing scan that future v3 NSM2 / v3 1.1+ work checks against.
"""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import pytest

from clearwater_modules_v3.processes import (
    Alkalinity,
    BenthicAlgae,
    CBOD,
    Carbon,
    DOX,
    FloatingAlgae,
    N2,
    Nitrogen,
    POM,
    Pathogen,
    Phosphorus,
)
from clearwater_modules_v3.utils.numerics import Diagnostics


# The 11 v3 NSM1 Process classes that the pattern-alignment spec
# covers. Riverine (transport) and Temperature (TSM canonical
# exemplar) are NOT in this list; they have separate contracts.
ALL_NSM1_PROCESSES = [
    Alkalinity,
    BenthicAlgae,
    CBOD,
    Carbon,
    DOX,
    FloatingAlgae,
    N2,
    Nitrogen,
    POM,
    Pathogen,
    Phosphorus,
]


# Pathogen uses _rate_with_components (rate-form integrator per spec
# §10 Q5); all other Processes use _change_with_components.
RATE_FORM_PROCESSES = {Pathogen}


def _canonical_helper_name(cls) -> str:
    return (
        "_rate_with_components"
        if cls in RATE_FORM_PROCESSES
        else "_change_with_components"
    )


def _legacy_shadow_name(cls) -> str:
    return (
        "_rate_legacy_inline"
        if cls in RATE_FORM_PROCESSES
        else "_change_legacy_inline"
    )


# ---------------------------------------------------------------------------
# Pattern B — fused helper
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cls", ALL_NSM1_PROCESSES, ids=lambda c: c.__name__)
def test_canonical_helper_exists(cls) -> None:
    """Every Process exposes a fused
    ``_change_with_components`` (or ``_rate_with_components`` for
    rate-form integrators) helper. Spec §3 pattern B."""
    name = _canonical_helper_name(cls)
    assert hasattr(cls, name), (
        f"{cls.__name__} missing canonical helper {name!r} (spec §3 pattern B)"
    )
    helper = getattr(cls, name)
    assert callable(helper), f"{cls.__name__}.{name} is not callable"


# ---------------------------------------------------------------------------
# REGISTRY_DIAGNOSTICS class attribute (pattern G)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cls", ALL_NSM1_PROCESSES, ids=lambda c: c.__name__)
def test_registry_diagnostics_present_and_nonempty(cls) -> None:
    """Every Process declares ``REGISTRY_DIAGNOSTICS: tuple[str, ...]``
    as a non-empty class attribute. Spec §3 pattern G."""
    assert hasattr(cls, "REGISTRY_DIAGNOSTICS"), (
        f"{cls.__name__} missing REGISTRY_DIAGNOSTICS class attribute"
    )
    rd = cls.REGISTRY_DIAGNOSTICS
    assert isinstance(rd, tuple), (
        f"{cls.__name__}.REGISTRY_DIAGNOSTICS is {type(rd).__name__}, "
        "expected tuple"
    )
    assert len(rd) > 0, f"{cls.__name__}.REGISTRY_DIAGNOSTICS is empty"
    for name in rd:
        assert isinstance(name, str), (
            f"{cls.__name__}.REGISTRY_DIAGNOSTICS contains non-string {name!r}"
        )


@pytest.mark.parametrize("cls", ALL_NSM1_PROCESSES, ids=lambda c: c.__name__)
def test_registry_diagnostics_names_unique(cls) -> None:
    """Within a Process, REGISTRY_DIAGNOSTICS names are unique."""
    rd = cls.REGISTRY_DIAGNOSTICS
    assert len(rd) == len(set(rd)), (
        f"{cls.__name__}.REGISTRY_DIAGNOSTICS has duplicate names: {rd}"
    )


# ---------------------------------------------------------------------------
# Pattern J — init_process captures diagnostics
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cls", ALL_NSM1_PROCESSES, ids=lambda c: c.__name__)
def test_init_process_method_exists(cls) -> None:
    """Every Process implements ``init_process(model, registry)`` for
    sibling discovery + diagnostics capture (spec §3 pattern J)."""
    assert hasattr(cls, "init_process"), (
        f"{cls.__name__} missing init_process method"
    )


@pytest.mark.parametrize("cls", ALL_NSM1_PROCESSES, ids=lambda c: c.__name__)
def test_init_process_captures_diagnostics(cls) -> None:
    """``init_process`` reads ``model.diagnostics`` and assigns to
    ``self.diagnostics``. We verify by source inspection rather than
    runtime mocking so a refactor that changes the attribute name
    (e.g., ``self._diag``) is caught."""
    init_process_src = inspect.getsource(cls.init_process)
    # Two acceptable patterns in v3:
    #   1. ``self.diagnostics = model_diagnostics`` after a getattr
    #   2. ``self.diagnostics = getattr(model, "diagnostics", ...)``
    assert "self.diagnostics" in init_process_src, (
        f"{cls.__name__}.init_process does not assign to "
        "self.diagnostics (spec §3 pattern J)"
    )


# ---------------------------------------------------------------------------
# Pattern D — unconditional clip-with-log
# ---------------------------------------------------------------------------


_GUARD_PATTERNS = [
    re.compile(r"isinstance\([^)]*xr\.DataArray\)[^:]*self\.diagnostics"),
    re.compile(r"self\.diagnostics is not None"),
]


@pytest.mark.parametrize("cls", ALL_NSM1_PROCESSES, ids=lambda c: c.__name__)
def test_no_clip_with_log_guard_branches(cls) -> None:
    """No ``isinstance(... DataArray) and self.diagnostics is not None``
    guard branches remain in any Process source. Pattern D unconditional
    clip-with-log was harmonised in Phase 1 / Phase 6."""
    src = inspect.getsource(inspect.getmodule(cls))
    for pat in _GUARD_PATTERNS:
        match = pat.search(src)
        assert match is None, (
            f"{cls.__name__} (file {inspect.getsourcefile(cls)}) contains "
            f"a clip-with-log guard branch: {match.group(0)!r}"
        )


@pytest.mark.parametrize("cls", ALL_NSM1_PROCESSES, ids=lambda c: c.__name__)
def test_no_step_zero_placeholder_in_clip_calls(cls) -> None:
    """No ``clip_negative_state(..., step=0)`` placeholders remain.
    Step attribution comes from ``diagnostics.current_step`` (Phase 0.6
    Q1). The placeholder was removed across the per-Process phases."""
    src = inspect.getsource(inspect.getmodule(cls))
    bad_pattern = re.compile(r"clip_negative_state\([^)]*step\s*=\s*0[^)]*\)")
    matches = bad_pattern.findall(src)
    assert not matches, (
        f"{cls.__name__} (file {inspect.getsourcefile(cls)}) contains "
        f"clip_negative_state(..., step=0) placeholder(s): {matches}"
    )


# ---------------------------------------------------------------------------
# Sanity check: each Process has its REGISTRY_DIAGNOSTICS exercised
# by a registry-diagnostics test file (single source of truth).
# ---------------------------------------------------------------------------


_TESTS_DIR = Path(__file__).resolve().parent


_PROCESS_TO_TEST_FILE = {
    Alkalinity: "test_alkalinity_registry_diagnostics.py",
    BenthicAlgae: "test_benthic_algae_registry_diagnostics.py",
    CBOD: "test_cbod_registry_diagnostics.py",
    Carbon: "test_carbon_registry_diagnostics.py",
    DOX: "test_dox_registry_diagnostics.py",
    FloatingAlgae: "test_floating_algae_registry_diagnostics.py",
    N2: "test_n2_registry_diagnostics.py",
    Nitrogen: "test_nitrogen_registry_diagnostics.py",
    POM: "test_pom_registry_diagnostics.py",
    Pathogen: "test_pathogen_registry_diagnostics.py",
    Phosphorus: "test_phosphorus_registry_diagnostics.py",
}


@pytest.mark.parametrize("cls", ALL_NSM1_PROCESSES, ids=lambda c: c.__name__)
def test_process_has_registry_diagnostics_test_file(cls) -> None:
    """Every Process owns a ``test_<process>_registry_diagnostics.py``
    file exercising pattern G. Pinned so a future PR removing the test
    file by mistake is caught."""
    fname = _PROCESS_TO_TEST_FILE[cls]
    test_file = _TESTS_DIR / fname
    assert test_file.is_file(), (
        f"{cls.__name__} is missing its registry-diagnostics test file "
        f"({test_file})"
    )
