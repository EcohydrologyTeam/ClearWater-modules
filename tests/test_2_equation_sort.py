"""Tests the processes_sort module."""
import pytest
import typing
from clearwater_modules_python.shared.types import (
    Process,
    Variable,
    VariableTypes,
)
from clearwater_modules_python.sorter import (
    split_variables,
    get_process_args,
    sort_dynamic_variables,
)


@pytest.fixture
def all_variables(static_variables, dynamic_variables) -> list[Variable]:
    """Return a list of all Variables."""
    return static_variables + dynamic_variables


def test_split_variables(
    static_variables: list[Variable],
    dynamic_variables: list[Variable],
    all_variables: list[Variable],
) -> None:
    """Test the split_variables function."""
    split_dict = split_variables(all_variables)

    # check keys
    for i in typing.get_args(VariableTypes):
        assert i in split_dict.keys()

    # check length
    assert len(split_dict['static']) == len(static_variables)
    assert len(split_dict['dynamic']) == len(dynamic_variables)

    # check order
    for i in range(len(static_variables)):
        assert split_dict['static'][i].name == static_variables[i].name
    for i in range(len(dynamic_variables)):
        assert split_dict['dynamic'][i].name == dynamic_variables[i].name


def test_get_process_args(process_functions: list[Process]) -> None:
    """Test the get_process_args function."""
    args: list[str] = get_process_args(process_functions[1])
    assert isinstance(args, list)
    assert args == ['a', 'b', 'dynamic_0']


def test_sort_dynamic_variables(
    dynamic_variables: list[Variable], 
    all_variables: list[Variable],
) -> None:
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
