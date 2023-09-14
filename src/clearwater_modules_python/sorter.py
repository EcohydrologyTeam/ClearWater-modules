"""Sorts and validates a stack of processes based on their annotations.

The idea here is to make the manual ordering of processes unnecessary, allowing for improved flexibility and maintainability.
Importantly this assumes processes are given names that match the required arguments of other processes.
"""
from clearwater_modules_python.shared.types import (
    Process,
    Variable,
    SplitVariablesDict,
)


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
    if 'return' in args: args.remove('return')
    return args


def __rapid_sort(
    static_vars: list[str],
    variable_args_dict: dict[str, tuple[Variable, list[str]]],
) -> list[Variable]:
    """Sorts dynamic variables based on their required arguments."""
    ordered_vars: list[Variable] = []
    previous_len: int = len(variable_args_dict)
    while len(variable_args_dict) > 0:
        drop_keys: list[str] = []
        for var_name, item in variable_args_dict.items():
            var, args = item
            ordered_names: list[str] = [var.name for var in ordered_vars]
            if all(arg in ordered_names + static_vars for arg in args):
                ordered_vars.append(var)
                drop_keys.append(var_name)
        for key in drop_keys:
            variable_args_dict.pop(key)
        if len(variable_args_dict) == previous_len:
            raise ValueError(
                f'Circular dependency detected in dynamic/state variables! '
                f'Variables remaining: {list(variable_args_dict.keys())}'
            )
        else:
            previous_len = len(variable_args_dict)
    return ordered_vars


def sort_variables_for_computation(variables_dict: SplitVariablesDict) -> list[Variable]:
    """Sorts all non-static variables based on their required arguments."""
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
