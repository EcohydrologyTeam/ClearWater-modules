from clearwater_data import ArrayLike


# @functools.lru_cache(maxsize=2)
def celsius_to_kelvin(celsius: ArrayLike) -> ArrayLike:
    """
    Convert Celsius to Kelvin
    """
    #return celsius + 273.15
    return celsius + 273.16 #for testing consistency with v1


def arrhenius_correction(
    water_temperature: ArrayLike,
    reaction_kinetics: ArrayLike,
    theta: ArrayLike,
) -> ArrayLike:
    """
    Computes an adjusted kinetics reaction rate coefficient for the specified water
    temperature using the van't Hoff form of the Arrhenius equation

    Parameters
    ----------
    water_temperature : ArrayLike
        Water temperature in degrees Celsius
    reaction_kinetics : ArrayLike
        Kinetics reaction (decay) coefficient at 20 degrees Celsius
    theta : ArrayLike
        Temperature correction factor

    Returns
    ----------
    ArrayLike
        Adjusted kinetics rate for the specified water temperature
    """
    return reaction_kinetics * theta ** (water_temperature - 20.0)
