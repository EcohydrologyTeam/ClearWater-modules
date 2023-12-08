"""JIT compiled processes for the heat model."""
import warnings
import numba
import numpy as np
import xarray as xr
from clearwater_modules.shared.processes import (
    celsius_to_kelvin,
)


@numba.njit
def air_temp_k(
    air_temp_c: xr.DataArray,
) -> xr.DataArray:
    """Calculate air temperature (K).

    Args:
        air_temp_c: Air temperature (C)
    """
    return celsius_to_kelvin(air_temp_c)


@numba.njit
def water_temp_k(
    water_temp_c: xr.DataArray,
) -> xr.DataArray:
    """Calculate water temperature (K).

    Args:
        water_temp_c: Water temperature (C)
    """
    return celsius_to_kelvin(water_temp_c)


@numba.njit
def mixing_ratio_air(
    eair_mb: xr.DataArray,
    pressure_mb: xr.DataArray,
) -> xr.DataArray:
    """Calculate air mixing ratio (unitless).

    Args:
        eair_mb: Vapour pressure of air (mb)
        pressure_mb: Atmospheric pressure (mb)
    """
    return 0.622 * eair_mb / (pressure_mb - eair_mb)


@numba.njit
def density_air(
    pressure_mb: xr.DataArray,
    air_temp_k: xr.DataArray,
    mixing_ratio_air: xr.DataArray,
) -> xr.DataArray:
    """Calculate air density (kg/m^3).

    Args:
        pressure_mb: Atmospheric pressure (mb)
        air_temp_k: Air temperature (K)
        mixing_ratio_air: Air mixing ratio (unitless)
    """
    return (
        0.348 *
        (pressure_mb / air_temp_k) *
        (1.0 + mixing_ratio_air) / (1.0 + 1.61 * mixing_ratio_air)

    )


@numba.njit
def emissivity_air(
    air_temp_k: xr.DataArray,
) -> xr.DataArray:
    """Calculate air emissivity (unitless).

    Args:
        air_temp_k: Air temperature (K)
    """
    return 0.00000937 * air_temp_k**2.0


@numba.njit
def wind_function(
    ri_function: xr.DataArray,
    wind_a: xr.DataArray,
    wind_b: xr.DataArray,
    wind_c: xr.DataArray,
    wind_speed: xr.DataArray,
) -> xr.DataArray:
    """Calculate wind function (unitless) for latent and sensible heat.

    Args:
        ri_function: Richardson function (unitless)
        wind_a: Wind function coefficient (unitless)
        wind_b: Wind function coefficient (unitless)
        wind_c: Wind function coefficient (unitless)
        wind_speed: Wind speed (m/s)
    """
    return (
        ri_function * (
            (wind_a / 1000000.0) +
            (wind_b / 1000000.0) *
            (wind_speed**wind_c)
        )
    )


@numba.njit
def q_latent(
    pressure_mb: xr.DataArray,
    density_water: xr.DataArray,
    lv: xr.DataArray,
    wind_function: xr.DataArray,
    esat_mb: xr.DataArray,
    eair_mb: xr.DataArray,
) -> xr.DataArray:
    """Latent heat flux (W/m^2).

    Args:
        pressure_mb: Atmospheric pressure (mb)
        density_water: Water density (kg/m^3)
        lv: Latent heat of vaporization (J/kg)
        wind_function: Wind function (unitless)
        esat_mb: Saturation vapour pressure (mb)
        eair_mb: Vapour pressure of air (mb)
    """
    return (
        (0.622 / pressure_mb) *
        lv *
        density_water *
        wind_function *
        (esat_mb - eair_mb)
    )


@numba.njit
def q_sensible(
    wind_kh_kw: xr.DataArray,
    cp_air: xr.DataArray,
    density_water: xr.DataArray,
    wind_function: xr.DataArray,
    air_temp_k: xr.DataArray,
    water_temp_k: xr.DataArray,
) -> xr.DataArray:
    # TODO: check if the return units are correct
    """Sensible heat flux (W/m2).

    Args:
        wind_kh_kw: Diffusivity ratio (unitless)
        ri_function: Richardson number (unitless)
        cp_air: Specific heat of air (J/kg/K)
        density_water: Water density (kg/m^3)
        wind_function: Wind function (unitless)
        air_temp_k: Air temperature (K)
        water_temp_k: Water temperature (K)
    """
    return (
        wind_kh_kw *
        cp_air *
        density_water *
        wind_function *
        (air_temp_k - water_temp_k)
    )


@numba.njit
def q_sediment(
    use_sed_temp: xr.DataArray,
    pb: xr.DataArray,
    cps: xr.DataArray,
    alphas: xr.DataArray,
    h2: xr.DataArray,
    sed_temp_c: xr.DataArray,
    water_temp_c: xr.DataArray,
) -> xr.DataArray:
    """Sediment heat flux (W/m^2).

    Args:
        use_sed_temp: Whether to calculate sed temp or not (boolean)
        pb: Sediment bulk density (kg/m^3)
        cps: Sediment specific heat (J/kg/K)
        alphas: Sediment thermal diffusivity (m^2/s)
        h2: Sediment active layer thickness (m)
        sed_temp_c: Sediment temperature (C)
        water_temp_c: Water temperature (C)
    """
    # 86400 converts the sediment thermal diffusivity from units of m^2/d to m^2/s
    return np.where(
        use_sed_temp,
        (
            pb * cps * alphas / 0.5 / h2 *
            (sed_temp_c - water_temp_c) / 86400.0
        ),
        0.0,
    )


@numba.njit
def dTdt_sediment_c(
    use_sed_temp: xr.DataArray,
    alphas: xr.DataArray,
    h2: xr.DataArray,
    water_temp_c: xr.DataArray,
    sed_temp_c: xr.DataArray,
) -> xr.DataArray:
    """Sediments temperature change (C).

    Args:
        use_sed_temp: Whether to calculate sed temp or not (boolean)
        alphas: Sediment thermal diffusivity (m^2/s)
        h2: Sediment active layer thickness (m)
        water_temp_c: Water temperature (C)
        sed_temp_c: Sediment temperature (C)
    """
    return np.where(
        use_sed_temp,
        alphas / (0.5 * h2 * h2) * (water_temp_c - sed_temp_c),
        0.0,
    )


@numba.njit
def mf_d_esat_dT(
    water_temp_k: xr.DataArray,
    a1: xr.DataArray,
    a2: xr.DataArray,
    a3: xr.DataArray,
    a4: xr.DataArray,
    a5: xr.DataArray,
    a6: xr.DataArray,
) -> xr.DataArray:
    """
    Compute the derivative of function computing saturation vapor pressure 
    as a function of water temperature (Kelvin)

    Fitting parameters for vapor pressure:
    Brutsaert (1982) Evaporation into the Atmosphere, p42
    """
    return (
        a1 +
        (2.0 * a2 * water_temp_k) +
        (3.0 * a3 * water_temp_k**2.0) +
        (4.0 * a4 * water_temp_k**3.0) +
        (5.0 * a5 * water_temp_k**4.0) +
        (6.0 * a6 * water_temp_k**5.0)
    )


# -----------------------------------------------------------------------------------
# Define functions to be used in the latent heat formulation
# -----------------------------------------------------------------------------------

@numba.njit
def mf_q_longwave_down(
    air_temp_k: xr.DataArray,
    emissivity_air: xr.DataArray,
    cloudiness: xr.DataArray,
    stefan_boltzmann: xr.DataArray,
) -> xr.DataArray:
    """
    Compute downwelling longwave radiation (W/m2)

    Parameters:
        air_temp_k (float):              Air temperature (Kelvin)
        emissivity_air (float):     Emissivity of air (unitless)
        cloudiness (float):         Cloudiness (fraction)

    Returns:
        Downwelling longwave radiation (W/m2, float)
    """

    return (1.0 + 0.17 * cloudiness**2) * emissivity_air * stefan_boltzmann * air_temp_k**4.0


@numba.njit
def mf_q_longwave_up(
    water_temp_k: xr.DataArray,
    emissivity_water: xr.DataArray,
    stefan_boltzmann: xr.DataArray,
) -> xr.DataArray:
    """
    Compute upwelling longwave radiation (W/m2) as a function of water temperature (Kelvin)
    """

    # logger.debug(f'mf_q_longwave_up({water_temp_k:.2f})')

    return emissivity_water * stefan_boltzmann * water_temp_k**4.0


@numba.njit
def mf_esat_mb(
    water_temp_k: xr.DataArray,
    a0: xr.DataArray,
    a1: xr.DataArray,
    a2: xr.DataArray,
    a3: xr.DataArray,
    a4: xr.DataArray,
    a5: xr.DataArray,
    a6: xr.DataArray,
) -> xr.DataArray:
    """
    Compute the saturation vapor pressure as a function of water temperature (Kelvin)

    Fitting parameters for vapor pressure are defined in:
    Brutsaert (1982) Evaporation into the Atmosphere, p42.
    """

    return (
        a0 +
        water_temp_k *
        (
            a1 +
            water_temp_k *
            (
                a2 +
                water_temp_k * (
                    a3 +
                    water_temp_k *
                    (a4 + water_temp_k * (a5 + water_temp_k * a6))
                )
            )
        )
    )


# Temperature conversion functions

@numba.njit
def ri_number(
    gravity: xr.DataArray,
    density_air: xr.DataArray,
    density_air_sat: xr.DataArray,
    wind_speed: xr.DataArray,
) -> xr.DataArray:
    """Calculates the Richardson Number.

    Args:
        gravity: Gravity (m/s2)
        density_air: Density of air (kg/m3)
        density_air_sat: Saturation density of air (kg/m3)
        wind_speed: Wind speed (m/s)
    """
    return (
        gravity *
        (density_air - density_air_sat) *
        2.0 / (density_air * (wind_speed**2.0))
    )


def ri_function(ri_number: xr.DataArray) -> np.ndarray:
    """Calculates the Richardson Function from the Richardson Number.
    Richardson Number:
        Unstable: 0.01 >= ri_function
        Stable: 0.01 <= ri_function < 2
        Neutral: -0.01 <  ri_function < 0.01
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    # NOTE: conditions are evaluated in the order of the original code
    # TODO: refactor into a single set of conditions, once testing is in place
    ri_number_bounded: np.ndarray = np.select(
        condlist=[
            ri_number > 2.0,
            ri_number < -1.0,
        ],
        choicelist=[
            2.0,
            -1.0,
        ],
        default=ri_number,
    )
    out: np.ndarray = np.select(
        condlist=[
            (ri_number_bounded < 0.0) & (ri_number_bounded >= -0.01),  # neutral
            (ri_number_bounded < 0.0) & (ri_number_bounded < -0.01),  # unstable
            (ri_number_bounded >= 0.0) & (ri_number_bounded <= 0.01),  # neutral
            (ri_number_bounded >= 0.0) & (ri_number_bounded > 0.01),  # stable
        ],
        choicelist=[
            1.0,                                         # neutral
            (1.0 - 22.0 * ri_number_bounded) ** 0.80,    # unstable
            1.0,                                         # neutral
            (1.0 + 34.0 * ri_number_bounded) ** (-0.80),  # stable
        ],
    )
    warnings.filterwarnings("default", category=RuntimeWarning)
    return out

    # old method with where, same logic
    # da: xr.DataArray = xr.where(ri_number > 2.0, (1.0 + 34.0*2.0)**(-0.80),
    #                   xr.where(ri_number < -1.0,  (1.0 - 22.0*-1.0)**(-0.80),
    #                   xr.where((ri_number < 0.0) & (ri_number >= -0.01), 1.0,
    #                   xr.where(ri_number < -0.01, (1.0 - 22.0 * ri_number)**0.80,
    #                   xr.where((ri_number >= 0.0) & (ri_number <= 0.01), 1.0, (1.0 + 34.0 * ri_number)**(-0.80)
    # )))))


@numba.njit
def mf_latent_heat_vaporization(water_temp_k: xr.DataArray) -> xr.DataArray:
    """
    Compute the latent heat of vaporization (W/m2) as a function of water temperature (Kelvin)
    """

    return 2499999 - 2385.74 * water_temp_k


@numba.njit
def mf_density_water(water_temp_c: xr.DataArray) -> xr.DataArray:
    """
    Compute density of water (kg/m3) as a function of water temperature (Celsius)
    """

    # logger.debug(f'mf_density_water({water_temp_c:.2f})')

    return (
        999.973 *
        (1.0 - (
            (
                (water_temp_c - 3.9863) *
                (water_temp_c - 3.9863) *
                (water_temp_c + 288.9414)
            ) /
            (
                508929.2 *
                (water_temp_c + 68.12963)
            )
        ))
    )


@numba.njit
def mf_density_air_sat(water_temp_k: xr.DataArray, esat_mb: float, pressure_mb: float) -> xr.DataArray:
    """
    Compute the density of saturated air at water surface temperature.

    Parameters:
        water_temp_k (float):        Water temperature (Kelvin)
        esat_mb (float):        Saturation vapor pressure in millibars
        pressure_mb (float):    Air pressure in millibars

    Returns:
        Density of saturated air at water surface temperature (kg/m3, float)
    """

    # logger.debug(f'mf_density_air_sat({water_temp_k:.2f}, {esat_mb:.2f}, {pressure_mb:.2f})')

    mixing_ratio_sat = 0.622 * esat_mb / (pressure_mb - esat_mb)
    return 0.348 * (pressure_mb / water_temp_k) * (1.0 + mixing_ratio_sat) / (1.0 + 1.61 * mixing_ratio_sat)


def mf_cp_water(water_temp_c: xr.DataArray) -> xr.DataArray:
    """
    Compute the specific heat of water (J/kg/K) as a function of water temperature (Celsius).
    This is used in computing the source/sink term.
    """
    # polynomial logic, this produced nearly identical outputs for the range of temperatures tested
    # See discussion about the two methods: https://github.com/EcohydrologyTeam/ClearWater-modules/pull/44
    # return (
    #    (4.65e-6 * water_temp_c - 0.001) *
    #    (water_temp_c + 0.085858) *
    #    (water_temp_c - 3.1331) *
    #    water_temp_c +
    #    4219.793
    # )
    return np.select(
        condlist=[
            water_temp_c <= 0.0,
            water_temp_c <= 5.0,
            water_temp_c <= 10.0,
            water_temp_c <= 15.0,
            water_temp_c <= 20.0,
            water_temp_c <= 25.0,
        ],
        choicelist=[
            4218.0,
            4202.0,
            4192.0,
            4186.0,
            4182.0,
            4180.0,
        ],
        default=4178.0,
    )


@numba.njit
def q_net(
    q_sensible: xr.DataArray,
    q_latent: xr.DataArray,
    q_longwave_up: xr.DataArray,
    q_longwave_down: xr.DataArray,
    q_solar: xr.DataArray,
    q_sediment: xr.DataArray,
) -> xr.DataArray:
    """Net heat flux (W/m^2).

    Args:
        q_sensible: Sensible heat flux (W/m^2)
        q_latent: Latent heat flux (W/m^2)
        q_longwave_up: Upward longwave radiation (W/m^2)
        q_longwave_down: Downward longwave radiation (W/m^2)
        q_solar: Solar radiation (W/m^2)
        q_sediment: Sediment heat flux (W/m^2)
    """
    return (
        q_sensible +
        q_solar +
        q_sediment +
        q_longwave_down -
        q_longwave_up -
        q_latent
    )


@numba.njit
def dTdt_water_c(
    q_net: xr.DataArray,
    surface_area: xr.DataArray,
    volume: xr.DataArray,
    density_water: xr.DataArray,
    cp_water: xr.DataArray,
) -> xr.DataArray:
    """Water temperature change (C).

    Args:
        q_net: Net heat flux (W/m^2)
        surface_area: Surface area (m^2)
        volume: Volume (m^3)
        density_water: Water density (kg/m^3)
        cp_water: Water specific heat (J/kg/K)
    """
    return (
        q_net *
        surface_area /
        (volume * density_water * cp_water)
    )


@numba.njit
def t_water_c(
    water_temp_c: xr.DataArray,
    dTdt_water_c: xr.DataArray,
) -> xr.DataArray:
    """Water temperature (C).

    Args:
        t_water_c: Water temperature (C)
        dt_water_c: Water temperature change (C)
    """
    return water_temp_c + dTdt_water_c
