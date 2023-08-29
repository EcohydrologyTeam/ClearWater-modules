import numba
from typing import TypedDict


class ProcessOutput(TypedDict):
    key: str
    name: str
    value: float
    units: str


@numba.njit
def air_mixing_ratio(
    eair_mb: float,
    pressure_mb: float,
) -> ProcessOutput:
    """Calculate air mixing ratio.

    Args:
        eair_mb: Vapour pressure of air (mb)
        pressure_mb: Atmospheric pressure (mb)

    Returns:
        Air mixing ratio (unitless)

    """
    value: float = 0.622 * eair_mb / (pressure_mb - eair_mb)
    return ProcessOutput(
        key='air_mixing_ratio',
        name='Air Mixing Ratio',
        value=value,
        units='unitless',
    )


@numba.njit
def air_density(
    pressure_mb: float,
    air_temp_k: float,
    air_mixing_ratio: float,
) -> ProcessOutput:
    """Calculate air density.

    Args:
        pressure_mb: Atmospheric pressure (mb)
        air_temp_k: Air temperature (K)
        air_mixing_ratio: Air mixing ratio (unitless)

    Returns:
        A ProcessOutput instance of air density (kg/m^3).

    """
    air_density: float = (
        0.348 *
        (pressure_mb / air_temp_k) *
        (1.0 + air_mixing_ratio) / (1.0 + 1.61 * air_mixing_ratio)
    )
    return ProcessOutput(
        key='air_density',
        name='Air Density',
        value=air_density,
        units="kg/m^3",
    )
