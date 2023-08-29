import numba


@numba.njit
def mixing_ratio_air(
    eair_mb: float,
    pressure_mb: float,
) -> float:
    """Calculate air mixing ratio (unitless).

    Args:
        eair_mb: Vapour pressure of air (mb)
        pressure_mb: Atmospheric pressure (mb)
    """
    return 0.622 * eair_mb / (pressure_mb - eair_mb)


@numba.njit
def density_air(
    pressure_mb: float,
    air_temp_k: float,
    air_mixing_ratio: float,
) -> float:
    """Calculate air density (kg/m^3).

    Args:
        pressure_mb: Atmospheric pressure (mb)
        air_temp_k: Air temperature (K)
        air_mixing_ratio: Air mixing ratio (unitless)
    """
    return (
        0.348 *
        (pressure_mb / air_temp_k) *
        (1.0 + air_mixing_ratio) / (1.0 + 1.61 * air_mixing_ratio)
    )


@numba.njit
def emissivity_air(
    air_temp_k: float,
) -> float:
    """Calculate air emissivity (unitless).

    Args:
        air_temp_k: Air temperature (K)
    """
    return 0.00000937 * air_temp_k**2.0


@numba.njit
def wind_function(
    wind_a: float,
    wind_b: float,
    wind_c: float,
    wind_speed: float,
) -> float:
    """Calculate wind function (unitless) for latent and sensible heat.

    Args:
        wind_a: Wind function coefficient (unitless)
        wind_b: Wind function coefficient (unitless)
        wind_c: Wind function coefficient (unitless)
        wind_speed: Wind speed (m/s)
    """
    return wind_a / 1000000.0 + wind_b / 1000000.0 * wind_speed**wind_c


@numba.njit
def q_latent(
    Ri_function: float,
    pressure_mb: float,
    density_water: float,
    Lv: float,
    wind_function: float,
    esat_mb: float,
    eair_mb: float,
) -> float:
    """Latent heat flux (W/m^2).

    Args:
        pressure_mb: Atmospheric pressure (mb)
        density_water: Water density (kg/m^3)
        Lv: Latent heat of vaporization (J/kg)
        wind_function: Wind function (unitless)
        esat_mb: Saturation vapour pressure (mb)
        eair_mb: Vapour pressure of air (mb)
    """
    return (
        Ri_function *
        (0.622 / pressure_mb) *
        Lv * density_water * wind_function *
        (esat_mb - eair_mb)
    )


@numba.njit
def q_sensible(
    wind_kh_kw: float,
    Ri_function: float,
    Cp_air: float,
    density_water: float,
    wind_function: float,
    air_temp_k: float,
    water_temp_k: float,
) -> float:
    # TODO: check if the return units are correct
    """Sensible heat flux (W/m2).

    Args:
        wind_kh_kw: Diffusivity ratio (unitless)
        Ri_function: Richardson number (unitless)
        Cp_air: Specific heat of air (J/kg/K)
        density_water: Water density (kg/m^3)
        wind_function: Wind function (unitless)
        air_temp_k: Air temperature (K)
        water_temp_k: Water temperature (K)
    """
    return (
        wind_kh_kw *
        Ri_function *
        Cp_air * density_water * wind_function *
        (air_temp_k - water_temp_k)
    )


@numba.njit
def q_sediment(
    pb: float,
    Cps: float,
    alphas: float,
    h2: float,
    TsedC: float,
    TwaterC: float,
) -> float:
    """Sediment heat flux (W/m^2).

    Args:
        pb: Sediment bulk density (kg/m^3)
        Cps: Sediment specific heat (J/kg/K)
        alphas: Sediment thermal diffusivity (m^2/s)
        h2: Sediment active layer thickness (m)
        TsedC: Sediment temperature (C)
        TwaterC: Water temperature (C)
    """
    # 86400 converts the sediment thermal diffusivity from units of m^2/d to m^2/s

    return (
        pb * Cps * alphas / 0.5 / h2 *
        (TsedC - TwaterC) / 86400.0
    )


@numba.njit
def dTdt_sediment_c(
    alphas: float,
    h2: float,
    TwaterC: float,
    TsedC: float,
) -> float:
    """Sediments temperature change (C).

    Args:
        alphas: Sediment thermal diffusivity (m^2/s)
        h2: Sediment active layer thickness (m)
        TwaterC: Water temperature (C)
        TsedC: Sediment temperature (C)
    """
    return (
        alphas / (0.5 * h2 * h2) *
        (TwaterC - TsedC)
    )


@numba.njit
def q_net(
    q_sensible: float,
    q_latent: float,
    q_longwave_up: float,
    q_longwave_down: float,
    q_solar: float,
    q_sediment: float,
) -> float:
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
        q_sensible -
        q_latent -
        q_longwave_up -
        q_longwave_down +
        q_solar +
        q_sediment
    )


@numba.njit
def dTdt_water_c(
    q_net: float,
    surface_area: float,
    volume: float,
    density_water: float,
    Cp_water: float,
) -> float:
    """Water temperature change (C).

    Args:
        q_net: Net heat flux (W/m^2)
        surface_area: Surface area (m^2)
        volume: Volume (m^3)
        density_water: Water density (kg/m^3)
        Cp_water: Water specific heat (J/kg/K)
    """
    return (
        q_net *
        surface_area /
        (volume * density_water * Cp_water)
    )
