"""Shared pytest fixtures.

Path shim: prepend the streaming-repo ``src/`` to ``sys.path`` so the
local ``clearwater_modules_v2`` (with Phase 2.A/2.B fixes) and
``clearwater_modules_v3`` are preferred over the conda env's editable
install of ``clearwater_modules``, which points at a vendor copy
elsewhere on disk that does NOT carry the streaming-repo work. This
must run before any ``clearwater_modules_v2`` / ``clearwater_modules_v3``
import in tests; placing it at the top of the root ``tests/conftest.py``
ensures the path is in place before pytest collects any test file that
imports from those packages.
"""
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import pytest
import xarray as xr
from clearwater_modules.shared.types import (
    Variable,
    Process,
)


@pytest.fixture(scope='session')
def initial_array():
    """Return a 10x10 xarray.DataArray."""
    return xr.DataArray(
        data=1.0,
        dims=['y', 'x'],
        coords={'x': range(10), 'y': range(10)},
        attrs={
            'long_name': 'Initial Array',
            'units': 'm',
            'description': 'An initial array.',
        }
    )


def mock_equation_0(a: float, b: float) -> float:
    return a + b


def mock_equation_1(a: float, b: float, dynamic_0: float) -> float:
    return a * b * dynamic_0


def mock_equation_2(a: float, b: float, dynamic_1: float) -> float:
    return a / b / dynamic_1


@pytest.fixture(scope='session')
def process_functions() -> list[Process]:
    """Return a list of process functions."""
    return [
        mock_equation_0,
        mock_equation_1,
        mock_equation_2,
    ]


@pytest.fixture(scope='session')
def static_variables() -> list[Variable]:
    """Return a static Variable."""
    out_vars: list[Variable] = []
    for i in ['a', 'b']:
        out_vars.append(Variable(
            name=f'{i}',
            long_name='Static Variable 0',
            units='m',
            description='A static variable.',
            use='static',
        ))
    return out_vars


@pytest.fixture(scope='session')
def dynamic_variables(process_functions) -> list[Variable]:
    """Return a list of dynamic Variables."""
    vars: list[Variable] = []
    for i, func in enumerate(process_functions):
        vars.append(Variable(
            name=f'dynamic_{i}',
            long_name=f'Dynamic Variable {i}',
            units='m',
            description='A dynamic variable.',
            use='dynamic',
            process=func,
        ))
    return vars
