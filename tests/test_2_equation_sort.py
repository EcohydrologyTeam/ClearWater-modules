"""Tests the processes_sort module."""
import pytest
import typing
from clearwater_modules_python.base import (
    Variable,
    Process,
    VariableTypes,
)
from clearwater_modules_python.sorter import (
    split_variables,
    get_process_args,
    sort_dynamic_variables,
)


def mock_equation_0(a: float, b: float) -> float:
    return a + b


def mock_equation_1(a: float, b: float, dynamic_0: float) -> float:
    return a * b * mock_equation_1


def mock_equation_2(a: float, b: float, dynamic_1: float) -> float:
    return a / b / mock_equation_1


@pytest.fixture
def process_functions() -> list[Process]:
    """Return a list of process functions."""
    return [
        mock_equation_0,
        mock_equation_1,
        mock_equation_2,
    ]


@pytest.fixture
def static_vars() -> list[Variable]:
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


@pytest.fixture
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


@pytest.fixture
def all_variables(static_vars, dynamic_variables) -> list[Variable]:
    """Return a list of all Variables."""
    return static_vars + dynamic_variables


def test_split_variables(static_vars, dynamic_variables, all_variables) -> None:
    """Test the split_variables function."""
    split_dict = split_variables(all_variables)

    # check keys
    for i in typing.get_args(VariableTypes):
        assert i in split_dict.keys()

    # check length
    assert len(split_dict['static']) == len(static_vars)
    assert len(split_dict['dynamic']) == len(dynamic_variables)

    # check order
    for i in range(len(static_vars)):
        assert split_dict['static'][i].name == static_vars[i].name
    for i in range(len(dynamic_variables)):
        assert split_dict['dynamic'][i].name == dynamic_variables[i].name


def test_get_process_args() -> None:
    """Test the get_process_args function."""
    args: list[str] = get_process_args(mock_equation_1)
    assert isinstance(args, list)
    assert args == ['a', 'b', 'dynamic_0']


def test_sort_dynamic_variables(dynamic_variables, all_variables) -> None:
    """Test the sort_dynamic_variables function."""
    split_dict = split_variables(all_variables)
    sorted_vars = sort_dynamic_variables(split_dict)

    # check types
    assert isinstance(sorted_vars, list)
    assert isinstance(sorted_vars[0], Variable)

    # check length
    assert len(sorted_vars) == len(dynamic_variables)

    # check order
    assert sorted_vars[0].name == 'dynamic_0'
    assert sorted_vars[1].name == 'dynamic_1'
    assert sorted_vars[2].name == 'dynamic_2'
