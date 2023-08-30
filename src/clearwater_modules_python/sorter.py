"""Sorts and validates a stack of processes based on their annotations.

The idea here is to make the manual ordering of processes unnecessary, allowing for improved flexibility and maintainability.
Importantly this assumes processes are given names that match the required arguments of other processes.
"""

from enum import EnumMeta
from clearwater_modules_python.shared_types import (
    Process,
    ParametersDict,
)


def get_equation_args(equation: Process) -> tuple[Process, list[str]]:
    """Return a tuple with the equation and its required argument names."""
    args: list[str] = list(equation.__annotations__.keys())
    args.remove('return')
    return equation, args

# TODO: Fix this, work in progress


def get_required_parameters(
    processes_args: list[list[str]],
    parameters: list[ParametersDict | EnumMeta],
) -> dict[str, float | int | bool]:
    """Return a dict of parameters required by the processes.

    Args:
        processes_args: A list of lists of equation arguments.
        parameters: A list of entities (Enums or TypedDicts) containing parameters.
    """
    ...


def order_processes(
    processes: list[Process],
    parameters: list[str],
) -> list[Process]:
    """Return a list of processes sorted by their required arguments."""
    sorted_processes: list[Process] = []
    while len(processes) > 0:
        for equation in processes:
            equation, args = get_equation_agrs(equation)
            if all(arg in sorted_processes for arg in args):
                sorted_processes.append(equation)
                processes.remove(equation)
    return sorted_processes
