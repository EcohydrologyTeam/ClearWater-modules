class VariableRegistry:
    """
    A simple registry for storing variables by name.
    """
    def __init__(self):
        self._registry = {}

    def register(self, name, value):
        """
        Register a variable with a given name.
        """
        self._registry[name] = value

    def get(self, name):
        """
        Retrieve a variable by name. Returns None if not found.
        """
        return self._registry.get(name)

    def update(self, name, value):
        """
        Update the value of an existing variable.
        Raises KeyError if the variable does not exist.
        """
        if name in self._registry:
            self._registry[name] = value
        else:
            raise KeyError(f"{name} not found in registry.")

    def unregister(self, name):
        """
        Remove a variable from the registry.
        Raises KeyError if the variable does not exist.
        """
        if name in self._registry:
            del self._registry[name]
        else:
            raise KeyError(f"{name} not found in registry.")
