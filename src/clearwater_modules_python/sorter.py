"""Sorts and validates a stack of processes based on their annotations.

The idea here is to make the manual ordering of processes unnecessary, allowing for improved flexibility and maintainability.
Importantly this assumes processes are given names that match the required arguments of other processes.
"""
import numba
from typing import (
    TypedDict,
)
from clearwater_modules_python.base import (
    Process,
    Variable,
)


class SplitVariablesDict(TypedDict):
    """A dict containing all variables split by type.

    Attributes:
        static: A list of static variables (out).
        dynamic: A list of dynamic variables (in).
        state: A list of state variables (in/out).
    """
    static: list[Variable]
    dynamic: list[Variable]
    state: list[Variable]


def split_variables(
    variables: list[Variable],
) -> SplitVariablesDict:
    """Return a tuple of static and dynamic variables."""
    split_dict: SplitVariablesDict = {
        'static': [],
        'dynamic': [],
        'state': [],
    }
    for variable in variables:
        if variable.use in split_dict.keys():
            split_dict[variable.use].append(variable)
    return split_dict


def get_process_args(equation: Process) -> list[str]:
    """Return a tuple with the equation and its required argument names."""
    args: list[str] = list(equation.__annotations__.keys())
    args.remove('return')
    return args


@numba.njit
def __rapid_sort(
    static_vars: list[str],
    variable_args_dict: dict[str, tuple[Variable, list[str]]],
) -> list[Variable]:
    """Sorts dynamic variables based on their required arguments."""
    ordered_vars: list[Variable] = []
    previous_len: int = 0
    while len(variable_args_dict) > 0:
        for var_name, item in variable_args_dict.items():
            var, args = item
            if all(arg in ordered_vars + static_vars for arg in args):
                ordered_vars.append(var)
                del variable_args_dict[var_name]
        if len(variable_args_dict) == previous_len:
            raise ValueError(
                f'Circular dependency detected in dynamic/state variables! '
                f'Variables remaining: {list(variable_args_dict.keys())}'
            )
    return ordered_vars


def sort_dynamic_variables(variables_dict: SplitVariablesDict) -> None:
    """Sorts dynamic variables based on their required arguments."""
    static_vars: list[str] = []
    for static_var in variables_dict['static']:
        static_vars.append(static_var.name)

    variable_args: dict[str, tuple[Variable, list[str]]] = {}
    for var in variables_dict['dynamic'] + variables_dict['state']:
        if var.process is None:
            raise ValueError(
                f'Dynamic/state variable {var.name} must be calculated by a process.'
            )
        variable_args[var.name] = var, get_process_args(var.process)

    return __rapid_sort(static_vars, variable_args)
