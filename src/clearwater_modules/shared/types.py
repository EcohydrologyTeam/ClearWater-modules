from dataclasses import dataclass
from typing import (
    TypedDict,
    Callable,
    Optional,
    Literal,
)

Process = Callable[..., float]
InitialVariablesDict = dict[str, float | int | bool]
VariableTypes = Literal['static', 'dynamic', 'state']

@dataclass(slots=True, frozen=True)
class Variable:
    """Variable type."""
    name: str
    long_name: str
    units: str
    description: str
    use: VariableTypes
    process: Optional[Process] = None


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



