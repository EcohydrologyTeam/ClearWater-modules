"""Shared pytest fixtures."""
import pytest
import xarray as xr
from clearwater_modules_python.base import (
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

