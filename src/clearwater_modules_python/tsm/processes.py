import numba


@numba.jit
def air_mixing_ratio(
    eair_mb: float,
    pressure_mb: float,
) -> float:
    """Calculate air mixing ratio.

    Args:
        eair_mb: Vapour pressure of air (mb)
        pressure_mb: Atmospheric pressure (mb)

    Returns:
        Air mixing ratio (unitless)

    """
    return 0.622 * eair_mb / (pressure_mb - eair_mb)


@numba.jit
def air_density(
    pressure_mb: float,
    air_temp_k: float,
    air_mixing_ratio: float,
) -> float:
    """Calculate air density.

    Args:
        pressure_mb: Atmospheric pressure (mb)
        air_temp_k: Air temperature (K)
        air_mixing_ratio: Air mixing ratio (unitless)

    Returns:
        Air density (kg/m^3)

    """
    return (
        0.348 *
        (pressure_mb / air_temp_k) *
        (1.0 + air_mixing_ratio) / (1.0 + 1.61 * air_mixing_ratio)
    )
