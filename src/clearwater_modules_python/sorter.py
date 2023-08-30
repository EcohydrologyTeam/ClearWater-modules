"""Sorts and validates a stack of equations based on their annotations.

The idea here is to make the manual ordering of equations unnecessary, allowing for improved flexibility and maintainability.
Importantly this assumes equations are given names that match the required arguments of other equations.
"""

from enum import EnumMeta
from clearwater_modules_python.shared_types import (
    Process,
    ConstantsDict,
)


def get_equation_args(equation: Process) -> tuple[Process, list[str]]:
    """Return a tuple with the equation and its required argument names."""
    args: list[str] = list(equation.__annotations__.keys())
    args.remove('return')
    return equation, args

# TODO: Fix this, work in progress


def get_required_constants(
    equations_args: list[list[str]],
    constants: list[ConstantsDict | EnumMeta],
) -> dict[str, float | int | bool]:
    """Return a dict of constants required by the equations.

    Args:
        equations_args: A list of lists of equation arguments.
        constants: A list of entities (Enums or TypedDicts) containing constants.
    """
    ...


def order_equations(
    equations: list[Process],
    constants: list[str],
) -> list[Process]:
    """Return a list of equations sorted by their required arguments."""
    sorted_equations: list[Process] = []
    while len(equations) > 0:
        for equation in equations:
            equation, args = get_equation_agrs(equation)
            if all(arg in sorted_equations for arg in args):
                sorted_equations.append(equation)
                equations.remove(equation)
    return sorted_equations
