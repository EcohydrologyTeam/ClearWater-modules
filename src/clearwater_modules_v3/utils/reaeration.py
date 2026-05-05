"""v3 NSM1 oxygen reaeration utilities.

Stateless reaeration primitives ported from v1 NSM1's shared library
(``clearwater_modules/shared/processes.py``). Provides hydraulic and
wind-driven reaeration rate selectors at 20 deg C and the combined
temperature-corrected effective reaeration coefficient. All functions
operate on ``xarray.DataArray`` inputs and return ``xarray.DataArray``.
"""

import numpy as np
import xarray as xr

from clearwater_modules_v3.utils.conversions import arrhenius_correction


def kah_20(
    kah_20_user: xr.DataArray,
    hydraulic_reaeration_option: xr.DataArray,
    velocity: xr.DataArray,
    depth: xr.DataArray,
    flow: xr.DataArray,
    topwidth: xr.DataArray,
    slope: xr.DataArray,
    shear_velocity: xr.DataArray,
) -> xr.DataArray:
    """Hydraulic oxygen reaeration coefficient at 20 deg C.

    Args:
        kah_20_user | 1/d | user-defined hydraulic reaeration rate at 20 deg C
            (used when ``hydraulic_reaeration_option == 1``).
        hydraulic_reaeration_option | int | selector for the empirical formula,
            1-9 (see Notes).
        velocity | m/s | average water velocity in cell.
        depth | m | average water depth in cell.
        flow | m^3/s | average flow rate in cell.
        topwidth | m | average top width of cell.
        slope | dimensionless | average bottom slope.
        shear_velocity | m/s | average bottom shear velocity.

    Returns:
        1/d | hydraulic reaeration coefficient at 20 deg C.

    Notes:
        Selector options:
            1. User-defined (``kah_20_user``).
            2. Covar (1976): O'Connor-Dobbins regime.
            3. Owens-Gibbs (1964).
            4. Churchill (1962).
            5. Tsivoglou-Wallace (1972): depth-piecewise blend of options 2-4.
            6. Padden-Gloyna (1971): flow-piecewise.
            7. USGS pool/riffle (Melching-Flores 1999): flow-piecewise.
            8. Thackston-Krenkel (1969): flow-piecewise.
            9. Langbien-Durum (1967): Froude-corrected shear-velocity form.
    """
    Uw_x_S = velocity * slope
    sqrt_g_h = (9.81 * depth) ** 0.5
    return xr.DataArray(
        np.select(
            condlist=[
                hydraulic_reaeration_option == 1,
                hydraulic_reaeration_option == 2,
                hydraulic_reaeration_option == 3,
                hydraulic_reaeration_option == 4,
                (hydraulic_reaeration_option == 5) & (depth < 0.61),
                (hydraulic_reaeration_option == 5) & (depth > 0.61),
                (hydraulic_reaeration_option == 5) & (depth == 0.61),
                (hydraulic_reaeration_option == 6) & (flow < 0.556),
                (hydraulic_reaeration_option == 6) & (flow >= 0.556),
                (hydraulic_reaeration_option == 7) & (flow < 0.556),
                (hydraulic_reaeration_option == 7) & (flow >= 0.556),
                (hydraulic_reaeration_option == 8) & (flow < 0.425),
                (hydraulic_reaeration_option == 8) & (flow >= 0.425),
                hydraulic_reaeration_option == 9,
            ],
            choicelist=[
                kah_20_user,
                (3.93 * velocity ** 0.5) / (depth ** 1.5),
                (5.32 * velocity ** 0.67) / (depth ** 1.85),
                (5.026 * velocity) / (depth ** 1.67),
                (5.32 * velocity ** 0.67) / (depth ** 1.85),
                (3.93 * velocity ** 0.5) / (depth ** 1.5),
                (5.026 * velocity) / (depth ** 1.67),
                517 * Uw_x_S ** 0.524 * flow ** -0.242,
                596 * Uw_x_S ** 0.528 * flow ** -0.136,
                88 * Uw_x_S ** 0.313 * depth ** -0.353,
                142 * Uw_x_S ** 0.333 * depth ** -0.66 * topwidth ** -0.243,
                31183 * velocity * slope,
                15308 * velocity * slope,
                2.16 * (1 + 9 * (velocity / sqrt_g_h) ** 0.25) * shear_velocity / depth,
            ],
            default=kah_20_user,
        )
    )


def kaw_20(
    kaw_20_user: xr.DataArray,
    wind_speed: xr.DataArray,
    wind_reaeration_option: xr.DataArray,
) -> xr.DataArray:
    """Wind-driven oxygen reaeration velocity at 20 deg C.

    Args:
        kaw_20_user | m/d | user-defined wind reaeration velocity at 20 deg C
            (used when ``wind_reaeration_option == 1``).
        wind_speed | m/s | wind speed at 10 m above the water surface
            (the input is referenced to 2 m and rescaled to 10 m via the
            standard 1/7 power law internally).
        wind_reaeration_option | int | selector for the empirical formula,
            1-13 (see Notes).

    Returns:
        m/d | wind-driven reaeration velocity at 20 deg C.

    Notes:
        Selector options:
            1.  User-defined (``kaw_20_user``).
            2.  Banks (1975).
            3.  Banks-Herrera (1977): piecewise at ``Uw10 = 3.5 m/s``.
            4.  Wanninkhof et al. (1991).
            5.  Chen-Kanwisher (1980).
            6.  Cole-Buchak (1995).
            7.  Banks-Herrera-Banks blend: piecewise at ``Uw10 = 5.5 m/s``.
            8.  Liss (1973).
            9.  Downing-Truesdale (1955): piecewise at ``Uw10 = 4.1 m/s``.
            10. Kanwisher (1963).
            11. Yu et al. (1977).
            12. Weiler (1974).
            13. Atkinson (1995): piecewise at ``Uw10 = 1.6 m/s``.
    """
    Uw10 = wind_speed * (10.0 / 2.0) ** 0.143
    return xr.DataArray(
        np.select(
            condlist=[
                wind_reaeration_option == 1,
                wind_reaeration_option == 2,
                (wind_reaeration_option == 3) & (Uw10 <= 3.5),
                (wind_reaeration_option == 3) & (Uw10 > 3.5),
                wind_reaeration_option == 4,
                wind_reaeration_option == 5,
                wind_reaeration_option == 6,
                (wind_reaeration_option == 7) & (Uw10 <= 5.5),
                (wind_reaeration_option == 7) & (Uw10 > 5.5),
                wind_reaeration_option == 8,
                (wind_reaeration_option == 9) & (Uw10 <= 4.1),
                (wind_reaeration_option == 9) & (Uw10 > 4.1),
                wind_reaeration_option == 10,
                wind_reaeration_option == 11,
                wind_reaeration_option == 12,
                (wind_reaeration_option == 13) & (Uw10 < 1.6),
                (wind_reaeration_option == 13) & (Uw10 >= 1.6),
            ],
            choicelist=[
                kaw_20_user,
                0.864 * Uw10,
                0.2 * Uw10,
                0.057 * Uw10 ** 2,
                0.728 * Uw10 ** 0.5 - 0.317 * Uw10 + 0.0372 * Uw10 ** 2,
                0.0986 * Uw10 ** 1.64,
                0.5 + 0.05 * Uw10 ** 2,
                0.362 * Uw10 ** 0.5,
                0.0277 * Uw10 ** 2,
                0.64 + 0.128 * Uw10 ** 2,
                0.156 * Uw10 ** 0.63,
                0.0269 * Uw10 ** 1.9,
                0.0276 * Uw10 ** 2,
                0.0432 * Uw10 ** 2,
                0.319 * Uw10,
                0.398,
                0.155 * Uw10 ** 2,
            ],
            default=kaw_20_user,
        )
    )


def ka_tc(
    kah_20: xr.DataArray,
    kaw_20: xr.DataArray,
    kah_theta: xr.DataArray,
    kaw_theta: xr.DataArray,
    T_water_C: xr.DataArray,
    depth: xr.DataArray,
) -> xr.DataArray:
    """Effective oxygen reaeration coefficient, temperature-corrected.

    Combines the temperature-corrected hydraulic rate (1/d) and the
    temperature-corrected wind-driven velocity (m/d, divided by depth to
    yield 1/d) into a single overall reaeration coefficient. Both
    components are corrected via the van't Hoff form of the Arrhenius
    equation against their respective ``theta`` factors.

    Args:
        kah_20 | 1/d | hydraulic reaeration rate at 20 deg C.
        kaw_20 | m/d | wind reaeration velocity at 20 deg C.
        kah_theta | dimensionless | Arrhenius theta for ``kah``.
        kaw_theta | dimensionless | Arrhenius theta for ``kaw``.
        T_water_C | deg C | water temperature.
        depth | m | average water depth in cell.

    Returns:
        1/d | overall temperature-corrected reaeration coefficient.
    """
    kah_tc = arrhenius_correction(T_water_C, kah_20, kah_theta)
    kaw_tc = arrhenius_correction(T_water_C, kaw_20, kaw_theta)
    return kaw_tc / depth + kah_tc
