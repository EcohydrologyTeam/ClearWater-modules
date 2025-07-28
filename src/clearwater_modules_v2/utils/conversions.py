import functools

from custom_types import ArrayLike
# TODO: type hints


# @functools.lru_cache(maxsize=2)
def celsius_to_kelvin(celsius: ArrayLike) -> ArrayLike:
    """
    Convert Celsius to Kelvin
    """
    return celsius + 273.15
