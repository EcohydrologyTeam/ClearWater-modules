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


def _first_dataarray(*args) -> xr.DataArray | None:
    """Return the first ``xr.DataArray`` argument, or ``None`` if all are
    scalars/ndarrays. Used to recover dim/coord metadata after a
    ``np.select`` call, which strips xarray attributes from its inputs."""
    for arg in args:
        if isinstance(arg, xr.DataArray):
            return arg
    return None


def kah_20(
    kah_20_user: xr.DataArray,
    hydraulic_reaeration_option: xr.DataArray,
    velocity: xr.DataArray,
    depth: xr.DataArray,
    flow: xr.DataArray,
    topwidth: xr.DataArray,
    slope: xr.DataArray,
    shear_velocity: xr.DataArray,
    min_depth: float = 0.0,
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
        min_depth | m | OPT-IN floor on ``depth`` before the inverse-depth
            empirical formulas evaluate. ``0.0`` = OFF (DEFAULT;
            byte-identical to the unfloored expression). >0 clamps
            ``depth = max(depth, min_depth)`` so the ``depth**-1.85`` /
            ``depth**-0.66`` / ``depth**-1.5`` terms cannot blow up at a
            sub-physical (newly-wet, ~0 m) cell. Coupled HEC-RAS-2D
            transport delivers a physically faithful per-cell mean depth
            that is legitimately ~1e-6 m (or exactly 0) at a wetting front;
            without a floor those cells produce ``ka_tc`` ~ 1e11/d and a
            Forward-Euler blow-up. Mirrors the opt-in CE-QUAL-W2-style
            ``min_reaeration_ka`` (a floor on ``ka_tc``) and the TSM
            ``q_net_depth_skip_threshold`` guards.

    Returns:
        1/d | hydraulic reaeration coefficient at 20 deg C.

    Notes:
        Selector options (author attributions match Fortran
        ``modGlobalParam.f90:268-339`` inline comments):
            1. User-defined (``kah_20_user``).
            2. O'Connor and Dobbins (1958): for depths 0.3-9.1 m and velocities
               0.15-0.49 m/s.
            3. Owens et al. (1964): for depths 0.12-0.73 m and velocities
               0.03-0.55 m/s.
            4. Churchill et al. (1962): for depths 0.61-3.36 m and velocities
               0.55-1.5 m/s.
            5. Cover (1976): depth-piecewise blend of options 2-4
               (Owens for depth<0.61, O'Connor for depth>0.61, Churchill at
               depth=0.61).
            6. Melching and Flores (1999) - pool-and-riffle streams:
               flow-piecewise.
            7. Melching and Flores (1999) - channel-controlled streams:
               flow-piecewise.
            8. Tsivoglou and Neal (1976): flow-piecewise on velocity * slope.
            9. Thackston and Dawson (2001): Froude-corrected shear-velocity
               form.
    """
    # Opt-in depth floor (default OFF). Clamp BEFORE any inverse-depth
    # term evaluates so a sub-physical newly-wet depth cannot drive the
    # hydraulic coefficient to ~1e11/d. See ``min_depth`` in Args.
    if min_depth > 0.0:
        depth = np.maximum(depth, min_depth)
    Uw_x_S = velocity * slope
    sqrt_g_h = (9.81 * depth) ** 0.5
    # ``np.select`` returns a bare ndarray whose dim labels are lost when
    # wrapped via ``xr.DataArray(arr)`` (xarray invents anonymous ``dim_0``).
    # Reattach the per-cell dims/coords from ``depth`` so downstream
    # broadcasting against ``oxygen_dissolved`` etc. operates per cell
    # rather than producing a spurious ``cell × dim_0`` result.
    result = np.select(
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
    template = _first_dataarray(
        depth,
        velocity,
        flow,
        topwidth,
        slope,
        shear_velocity,
        kah_20_user,
        hydraulic_reaeration_option,
    )
    if template is None:
        return xr.DataArray(result)
    return xr.DataArray(result, coords=template.coords, dims=template.dims)


def kaw_20(
    kaw_20_user: xr.DataArray,
    wind_speed: xr.DataArray,
    wind_reaeration_option: xr.DataArray,
    wind_input_height: float = 2.0,
    wind_shelter: "xr.DataArray | float" = 1.0,
) -> xr.DataArray:
    """Wind-driven oxygen reaeration velocity at 20 deg C.

    Args:
        kaw_20_user | m/d | user-defined wind reaeration velocity at 20 deg C
            (used when ``wind_reaeration_option == 1``).
        wind_speed | m/s | observed wind speed at ``wind_input_height``
            above the water surface. The empirical formulas below are
            calibrated against the 10 m reference height (Uw10); this
            function converts the input to Uw10 via the standard
            ``(10 / wind_input_height) ** 0.143`` power law.
        wind_reaeration_option | int | selector for the empirical formula,
            1-13 (see Notes).
        wind_input_height | m | height (m) above the water surface at
            which ``wind_speed`` was observed. Default 2.0 m matches the
            pre-Phase-H behaviour and the v1/Fortran NSM1 inheritance.
            For NOAA ASOS / METAR / GridMET / NLDAS records (standard
            anemometer at 10 m), pass ``wind_input_height=10.0`` so the
            internal Uw10 calculation does not double-count the height
            correction. Phase H-8 (2026-05-21): added to align with the
            v3 TSM ``Temperature.wind_input_height`` configurable.
        wind_shelter | dimensionless | per-cell wind-shelter coefficient
            (``wind_shelter_coefficient``), the same forcing the TSM wind
            function consumes. Applied to the raw wind BEFORE the height
            rescale, mirroring ``Temperature._compute_effective_wind``
            (``raw * shelter * height_factor``) and CE-QUAL-W2
            ``w2_4_unix.f90:480``. Default ``1.0`` (no shelter) preserves
            prior numerical output for callers that don't pass it. Accepts
            a scalar or an ``xr.DataArray`` aligned with ``wind_speed``.

    Returns:
        m/d | wind-driven reaeration velocity at 20 deg C.

    Notes:
        Selector options (author attributions match Fortran
        ``modGlobalParam.f90:341-414`` inline comments):
            1.  User-defined (``kaw_20_user``).
            2.  Broecker et al. (1978).
            3.  Gelda et al. (1996): piecewise at ``Uw10 = 3.5 m/s``.
            4.  Banks and Herrera (1977).
            5.  Wanninkhof (1991).
            6.  Cole and Buchak (1993).
            7.  Banks (1975): piecewise at ``Uw10 = 5.5 m/s``.
            8.  Smith (1978).
            9.  Liss (1973): piecewise at ``Uw10 = 4.1 m/s``.
            10. Downing and Truesdale (1955).
            11. Kanwisher (1963).
            12. Yu et al. (1977).
            13. Weiler (1974): piecewise at ``Uw10 = 1.6 m/s``.
    """
    # Phase H-8 (2026-05-21): use the caller-supplied
    # ``wind_input_height`` (default 2.0 m for legacy parity) rather
    # than the hard-coded ``2.0`` rescale base. With Phase G-3 wiring
    # registry-driven wind into DOX/N2/Carbon and Phase F bumping the
    # canonical runner's ``--wind-input-height`` default to 10.0 to
    # match the KSLE ASOS anemometer, the previous hard-coded
    # ``(10/2)**0.143`` rescale would have applied a 1.35x factor on
    # top of an already-10-m measurement. ``wind_input_height == 10``
    # now produces the identity factor here, matching the Temperature
    # module's log-law contract for the same height value.
    # Wind shelter (wind_shelter_coefficient) is applied to the raw wind
    # BEFORE the height rescale, matching the composition order in
    # Temperature._compute_effective_wind (raw * shelter * height_factor)
    # and CE-QUAL-W2 w2_4_unix.f90:480. Default wind_shelter=1.0 is a no-op.
    sheltered_wind = wind_speed * wind_shelter
    if abs(wind_input_height - 10.0) < 1e-12:
        Uw10 = sheltered_wind
    else:
        Uw10 = sheltered_wind * (10.0 / wind_input_height) ** 0.143
    # See ``kah_20`` for the same ``np.select`` dim-stripping fix; reattach
    # ``wind_speed`` dims/coords to the result.
    result = np.select(
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
    template = _first_dataarray(wind_speed, kaw_20_user, wind_reaeration_option)
    if template is None:
        return xr.DataArray(result)
    return xr.DataArray(result, coords=template.coords, dims=template.dims)


def ka_tc(
    kah_20: xr.DataArray,
    kaw_20: xr.DataArray,
    kah_theta: xr.DataArray,
    kaw_theta: xr.DataArray,
    T_water_C: xr.DataArray,
    depth: xr.DataArray,
    min_depth: float = 0.0,
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
        min_depth | m | OPT-IN floor on ``depth`` before the wind term
            ``kaw_tc / depth`` divides. ``0.0`` = OFF (DEFAULT;
            byte-identical to the unfloored expression). >0 clamps
            ``depth = max(depth, min_depth)``. Pass the SAME ``min_depth``
            to ``kah_20`` so both the hydraulic and wind components see the
            floored depth; flooring here alone does not protect the
            (already computed) ``kah_20`` value. See ``kah_20`` for the
            newly-wet-cell rationale.

    Returns:
        1/d | overall temperature-corrected reaeration coefficient.
    """
    # Opt-in depth floor (default OFF), applied to the wind term's 1/depth.
    if min_depth > 0.0:
        depth = np.maximum(depth, min_depth)
    kah_tc = arrhenius_correction(T_water_C, kah_20, kah_theta)
    kaw_tc = arrhenius_correction(T_water_C, kaw_20, kaw_theta)
    return kaw_tc / depth + kah_tc
