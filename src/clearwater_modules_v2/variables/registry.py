from variables.base import Variable
from datetime import datetime


class VariableRegistry:
    """
    A simple registry for storing variables by name.
    """

    def __init__(self):
        self._registry = {}

    def register(self, name: str, value: Variable, overwrite: bool = False):
        """
        Register a variable with a given name.
        """
        if name in self._registry and not overwrite:
            raise ValueError(f"Variable {name} already registered.")
        self._registry[name] = value

    def unregister(self, name: str):
        """
        Remove a variable from the registry.
        """
        if name in self._registry:
            del self._registry[name]

    def get(self, name: str) -> object:
        """
        Retrieve a variable by name.
        """
        variable: Variable | None = None
        try:
            variable = self._registry[name]
        except KeyError:
            raise ValueError(f"Variable {name} not found.")

        return variable.get()

    def get_at_time(self, name: str, time: datetime) -> object:
        """
        Retrieve a variable by name and time.
        """
        variable: Variable | None = None
        try:
            variable = self._registry[name]
        except KeyError:
            raise ValueError(f"Variable {name} not found.")

        return variable.get_at_time(time)
