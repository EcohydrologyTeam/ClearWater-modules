from clearwater_data import ArrayLike


# @functools.lru_cache(maxsize=2)
def celsius_to_kelvin(celsius: ArrayLike) -> ArrayLike:
    """
    Convert Celsius to Kelvin
    """
    return celsius + 273.15
