"""v3 Temperature process (Temperature Simulation Module).

This is the v3-native merged TSM. It combines:

- v2's class-based ``Process`` framework, ``ProcessFactory`` registration,
  per-process ``time_step``, and YAML-driven ``init_from_file`` configuration.
- v1's latent-heat unit fix (Kelvin polynomial → Celsius polynomial).
- v1's thin-water stability guard (depth ramp on net flux + per-hour rate
  cap on dT/dt).
- v2's ``mixing_ratio_air`` edge guard, ``use_sediment_temperature`` flag,
  and ``__skip_first_time_step`` v1-coupling-compat logic.

Default values for the new stability parameters (``q_net_depth_ramp_ref=0.3``
and ``dTdt_max_per_hour=5.0``) match v1's hardened defaults. Both can be
disabled to reproduce the unguarded v2 behavior:

- ``q_net_depth_ramp_ref=0.0`` disables the depth ramp.
- ``dTdt_max_per_hour=float("inf")`` disables the rate cap.

References:
    https://erdc-library.erdc.dren.mil/server/api/core/bitstreams/81b728f8-87a7-4ef8-e053-411ac80adeb3/content
"""

from datetime import datetime, timedelta
import warnings

import numpy as np
import xarray as xr

from clearwater_data.custom_types import ArrayLike
from clearwater_data.variables import VariableRegistry

from clearwater_modules_v3.processes.base import Process, ProcessFactory
from clearwater_modules_v3.utils import constants, conversions


class Temperature(Process):
    """Energy-balance water temperature kinetics.

    Implements the methodology in the ERDC Water Temperature Simulation
    Module (NSM/TSM). The energy budget integrates atmospheric longwave,
    upwelling longwave, latent heat, sensible heat, sediment heat exchange,
    and incident solar radiation; the resulting net flux drives a forward
    Euler update of water temperature with the v1 thin-water guards
    applied.
    """

    # Brutsaert (1982) saturation-vapor-pressure polynomial coefficients.
    # See "Evaporation into the Atmosphere", p42.
    __A0 = 6984.505294
    __A1 = -188.903931
    __A2 = 2.133357675
    __A3 = -1.288580973e-2
    __A4 = 4.393587233e-5
    __A5 = -8.023923082e-8
    __A6 = 6.136820929e-11

    variables = [
        "water_temperature",
        "wetted_surface_area",
        "volume",
        "cloudiness",
        "air_temperature",
        "solar_radiation",
        "wind_speed",
        "atmospheric_pressure",
        "atmospheric_vapor_pressure",
        "sediment_temperature",
        "sediment_thickness",
    ]

    # C5 fix (review-findings 2026-05-04): variables this process *writes*
    # back to the registry. The orchestration-layer wet-mask in
    # ``Model.__apply_wet_mask`` reads this list and masks NaN on dry
    # cells **only for outputs**, leaving forcing inputs (wind_speed,
    # air_temperature, solar_radiation, cloudiness, atmospheric_pressure,
    # atmospheric_vapor_pressure, wetted_surface_area, volume,
    # sediment_thickness) untouched. ``sediment_temperature`` is in the
    # output list because v3 evolves it dynamically (see C10 / Fortran
    # parity); when ``evolve_sediment_temperature=False`` the kernel
    # writes the unchanged value, so masking dry cells of
    # ``sediment_temperature`` is still correct (dry cells didn't have
    # heat exchange anyway).
    output_variables = [
        "water_temperature",
        "sediment_temperature",
    ]

    def __init__(
        self,
        wind_a: float = 0.3,
        wind_b: float = 1.5,
        wind_c: float = 2.0,
        wind_input_height: float = 2.0,
        surface_z0: float = 0.001,
        wind_shelter: float = 1.0,
        sediment_density: ArrayLike = 1600.0,
        sediment_specific_heat: float = 1673.0,
        air_diffusivity_ratio: float = 1.0,
        sediment_diffusivity: float = 0.0432,
        time_step: timedelta = timedelta(minutes=5),
        use_sediment_temperature: bool = True,
        evolve_sediment_temperature: bool = True,
        q_net_depth_ramp_ref: float = 0.3,
        dTdt_max_per_hour: float = 5.0,
    ) -> None:
        """Initialize the temperature process.

        Args:
            wind_a, wind_b, wind_c: Wind-function parameters in the
                Edinger, Brady & Geyer (1974) form
                ``f(W) = (a + b * W^c) / 1e6`` (multiplied internally
                by the Richardson stability function in
                :py:meth:`wind_function`). Defaults are
                ``a = 0.3, b = 1.5, c = 2.0``.

                The exponent ``c = 2.0`` (quadratic in wind) matches the
                explicit default in the CE-QUAL-W2 user manual
                (``AFW / BFW / CFW = 9.2 / 0.46 / 2.0``) and the
                QUAL2K Brady-Graves-Geyer default
                (``19.0 / 0.95 / 2.0``, different unit system, same
                exponent). All seven example case studies shipped with
                CE-QUAL-W2 — including Spokane River, Columbia Slough
                Estuary, Detroit, DeGray, Long Lake, Bonneville Dam,
                and Berlin Milton — use ``CFW = 2.0`` regardless of
                waterbody type. River-specific tuning in W2 is done
                via the per-segment wind shelter coefficient
                (``WSC(I)``), not via different wind-function
                exponents.

                The v3 magnitude coefficients ``a = 0.3, b = 1.5`` are
                inherited from v1's ``clearwater_modules.tsm.constants``
                and use the v3 ``/1e6`` normalisation. They are **not**
                directly interchangeable with W2's SI values
                (``9.2 / 0.46``) — the units bake the normalisation
                into the coefficients. A focused calibration study to
                revisit ``a`` and ``b`` against observation is tracked
                as future work.

                ``wind_c`` is validated: values not in ``{1.0, 2.0}``
                emit a ``UserWarning``, and values outside
                ``(0.0, 3.0]`` raise ``ValueError``. ``c = 3.0`` is
                allowed at the upper bound for back-compat with
                explicit opt-ins from runs that have already been
                calibrated against the prior v3 default.

                Pass any subset to override per-instance; YAML configs
                may also override via ``wind_a / wind_b / wind_c`` keys
                at ``init_from_file`` time.

                References:

                * Edinger, J.E., D.K. Brady, and J.C. Geyer (1974),
                  *Heat exchange and transport in the environment*,
                  Report 14, Cooling Water Discharge Research Project
                  (RP-49), Electric Power Research Institute, Palo
                  Alto, CA, 125 pp.
                * Brady, D.K., W.L. Graves, and J.C. Geyer (1969),
                  *Surface heat exchange at power plant cooling
                  lakes*, Cooling Water Discharge Research Project
                  Report 5, Edison Electric Institute, New York.
                * CE-QUAL-W2 User Manual, AFW/BFW/CFW default entries.
                * Chapra, S.C. (2008), *QUAL2K User Manual*,
                  §4.1.4 wind function ``f(Uw)``.
                * ``design/clearwater_modules_v3_tsm_wind_function_specification.md``
                  in this repo (full spec with W2 source citations and
                  case-study sensitivity sweep).
            wind_input_height: Height (m) above the water surface at
                which the application's ``wind_speed`` registry
                variable was measured. Default ``2.0`` (no correction),
                matching the Edinger 1974 / CE-QUAL-W2 convention. Set
                to ``10.0`` when registering raw ASOS / METAR /
                GridMET / NLDAS wind without external pre-correction.
                The log-law conversion uses ``surface_z0`` (default
                ``0.001`` m, typical for open water) per
                ``U_2 / U_z = ln(2 / z0) / ln(z / z0)``. For ASOS
                over a flat airfield (i.e., not over the water
                surface itself) the user should additionally apply an
                airfield-to-water roughness correction at the
                application layer; the in-module log-law assumes the
                wind was measured above the same surface as the
                receiving water.
            surface_z0: Roughness length (m) of the water surface,
                used by the in-module log-law height correction when
                ``wind_input_height != 2.0``. Default ``0.001`` m.
                Larger values (e.g., ``0.003`` m) are sometimes used
                for smoother-water regimes; sensitivity is small.
            wind_shelter: Scalar wind shelter coefficient applied to
                ``wind_speed`` before the Edinger formula and before
                any height correction. Default ``1.0`` (no shelter).
                Values < 1 represent canopy or topographic wind
                reduction; typical W2 values are ``1.0`` (open lake),
                ``0.85-0.90`` (open river / reservoir), ``0.5-0.85``
                (narrow channel with riparian canopy), ``0.3-0.5``
                (heavily shaded backwater). Mirrors CE-QUAL-W2's
                ``WSC(I)`` per-segment shelter coefficient (see
                ``w2_4_unix.f90:480``). For per-cell shelter,
                register ``wind_shelter_coefficient`` as an optional
                forcing variable in the registry; the registry value
                overrides this constructor scalar when present.
            sediment_density: Sediment bulk density (kg/m^3). Fortran
                default ``pb = 1600``; matches v1.
            sediment_specific_heat: Sediment specific heat (J/kg/C).
                Fortran default ``Cps = 1673``; matches v1.
            air_diffusivity_ratio: Sensible-heat diffusivity ratio.
            sediment_diffusivity: Sediment thermal diffusivity in
                **m^2/day**. Default ``0.0432`` matches the Fortran
                ``alphas`` default in
                ``HEC-RAS-WQ/RAS-1D-WQ/Kinetics Libraries/Temperature*/Source files/modTemperature.f90``
                and v1 ``clearwater_modules.tsm.constants.alphas``. The
                ``flux_sediment`` formula divides by 86400 to convert to
                W/m^2 internally; supplying an m^2/s value will produce a
                86400x-too-small flux.
            time_step: Substep length.
            use_sediment_temperature: If False, all sediment heat exchange
                is disabled (no flux, no sediment temperature evolution).
                Matches the Fortran ``use_SedTemp`` flag.
            evolve_sediment_temperature: If True (default), sediment
                temperature evolves each substep per the Fortran formula
                ``dTsed/dt = alphas / (0.5 * h2^2) * (T_water - T_sed)``.
                The water-sediment heat exchange is energy-conservative
                in this mode. If False, sediment temperature stays at
                its registered (or hotstart-seeded) value forever — this
                reproduces the v1/v2 Python ports' behavior, which is
                **not** energy-conservative and biases sediment heat
                damping under sustained warm or cold forcing. Has no
                effect when ``use_sediment_temperature`` is False.
            q_net_depth_ramp_ref: Reference depth (m) for the thin-water
                flux ramp. The net flux is multiplied by
                ``min(1, depth / q_net_depth_ramp_ref)``. Set to ``0.0``
                to disable (legacy behavior).
            dTdt_max_per_hour: Maximum |dT/dt| (K/hr). Per-substep delta T
                is clipped to ``+/- dTdt_max_per_hour * dt_hours``. Set to
                ``float("inf")`` to disable.
        """
        # Wind-function parameters (Edinger, Brady & Geyer 1974). v3
        # defaults are ``0.3 / 1.5 / 2.0`` per the constructor docstring;
        # ``c = 2.0`` matches the CE-QUAL-W2 manual default and the
        # QUAL2K Brady-Graves-Geyer default. Callers may override
        # per-instance or via YAML at config time. See
        # design/clearwater_modules_v3_tsm_wind_function_specification.md.
        self.wind_a = wind_a
        self.wind_b = wind_b
        self.wind_c = wind_c
        # Wind input transforms: log-law height correction +
        # wind-shelter scalar. Mirrors CE-QUAL-W2's
        # ``WIND2(I) = WIND(JW) * WSC(I) * log_law`` pre-correction
        # at ``w2_4_unix.f90:480``. The composition
        # ``raw_wind * shelter * height_factor`` is applied in
        # :py:meth:`wind_function` (and in the equilibrium-temperature
        # Newton-Raphson iteration which also goes through that
        # method). Default values (``input_height=2.0``,
        # ``z0=0.001``, ``shelter=1.0``) make the transforms a no-op.
        self.wind_input_height = wind_input_height
        self.surface_z0 = surface_z0
        self.wind_shelter = wind_shelter
        # Per-cell shelter cache populated in :py:meth:`run` if the
        # registry has a ``wind_shelter_coefficient`` forcing variable;
        # ``None`` means "fall back to the scalar ``self.wind_shelter``".
        self._cached_shelter: ArrayLike | None = None
        self.sediment_density = sediment_density
        self.sediment_specific_heat = sediment_specific_heat
        self.air_diffusivity_ratio = air_diffusivity_ratio
        self.sediment_diffusivity = sediment_diffusivity
        self.use_sediment_temperature = use_sediment_temperature
        self.evolve_sediment_temperature = evolve_sediment_temperature
        self.q_net_depth_ramp_ref = q_net_depth_ramp_ref
        self.dTdt_max_per_hour = dTdt_max_per_hour

        # M1 fix (review-findings 2026-05-04): validate stability params.
        # ``q_net_depth_ramp_ref`` must be finite and >= 0 (0 disables the
        # depth ramp). NaN, +inf, and negative values are rejected so a
        # silent-disable cannot occur from a typo or bad config.
        if not (np.isfinite(q_net_depth_ramp_ref) and q_net_depth_ramp_ref >= 0.0):
            raise ValueError(
                f"q_net_depth_ramp_ref must be >= 0.0 and finite (set 0.0 to disable); "
                f"got {q_net_depth_ramp_ref!r}"
            )
        # ``dTdt_max_per_hour`` must be strictly > 0. ``+inf`` is the
        # documented disable value and passes ``> 0.0``. Zero would freeze
        # the temperature field; negatives produce a constant-pegged
        # field; NaN propagates everywhere. Reject all three.
        # Note: ``np.isnan(+inf)`` is False, so the ``> 0.0`` predicate
        # admits ``+inf`` and rejects NaN (NaN comparisons are False).
        if not (dTdt_max_per_hour > 0.0):
            raise ValueError(
                f"dTdt_max_per_hour must be > 0.0 (set float('inf') to disable); "
                f"got {dTdt_max_per_hour!r}"
            )

        # Validate ``wind_c`` against the literature reference family.
        # CE-QUAL-W2's manual default is ``CFW = 2.0`` and the source
        # supports ``CFW in {1.0, 2.0}`` (``heat-exchange.f90:78``
        # comment "CFW not determined for other values of CFW"). v3
        # defaults to ``2.0`` and warns outside ``{1.0, 2.0}``; values
        # outside ``(0.0, 3.0]`` are physically indefensible and are
        # rejected. ``c = 3.0`` is allowed at the upper bound for
        # back-compat with explicit opt-ins from runs calibrated
        # against the prior v3 default.
        if not (0.0 < wind_c <= 3.0):
            raise ValueError(
                f"wind_c must be in (0.0, 3.0]; got {wind_c!r}"
            )
        if wind_c not in (1.0, 2.0):
            warnings.warn(
                f"wind_c = {wind_c} is outside the values supported by "
                f"the Edinger family of wind-function parameterisations. "
                f"CE-QUAL-W2 explicitly defaults to CFW = 2.0 and "
                f"supports CFW = 1.0; the v3 default is c = 2.0. "
                f"QUAL2K's Brady-Graves-Geyer default is also c = 2.0. "
                f"Other values are flagged as 'CFW not determined' "
                f"(W2 heat-exchange.f90:78). Coefficient `b` is "
                f"unit-coupled to `c`; using a non-standard exponent "
                f"without re-calibrating `b` will produce unphysical "
                f"heat fluxes.",
                UserWarning,
                stacklevel=2,
            )

        # Validate ``wind_input_height`` and ``surface_z0``. The log-law
        # correction ``log(2 / z0) / log(input_height / z0)`` requires
        # both ``z0 > 0`` and ``input_height > z0`` (otherwise the
        # logarithm denominators are non-positive or zero). Default
        # ``input_height = 2.0`` makes the correction a no-op.
        if not (np.isfinite(wind_input_height) and wind_input_height > 0.0):
            raise ValueError(
                f"wind_input_height must be > 0.0 and finite (m); "
                f"got {wind_input_height!r}"
            )
        if not (np.isfinite(surface_z0) and surface_z0 > 0.0):
            raise ValueError(
                f"surface_z0 must be > 0.0 and finite (m); "
                f"got {surface_z0!r}"
            )
        if surface_z0 >= wind_input_height:
            raise ValueError(
                f"surface_z0 must be strictly less than wind_input_height; "
                f"got z0 = {surface_z0!r}, input_height = "
                f"{wind_input_height!r}. The log-law correction "
                f"log(2/z0) / log(input_height/z0) requires z0 < "
                f"input_height for both logarithms to be finite and "
                f"positive."
            )

        # Validate ``wind_shelter``. Must be > 0 (zero shuts off all
        # wind-driven flux including the wind-independent ``a`` term
        # via the ``wind_speed * shelter`` multiplication). Values >> 1
        # are physically unusual but possible for funneled flow regimes
        # so we warn rather than reject. CE-QUAL-W2 manual typical
        # range: 0.3 (heavy shelter) to 1.0 (no shelter).
        if not (np.isfinite(wind_shelter) and wind_shelter > 0.0):
            raise ValueError(
                f"wind_shelter must be > 0.0 and finite (set 1.0 for "
                f"no sheltering); got {wind_shelter!r}"
            )
        if wind_shelter > 1.0:
            warnings.warn(
                f"wind_shelter = {wind_shelter} is greater than 1.0. "
                f"CE-QUAL-W2 manual typical range is 0.3-1.0 "
                f"(0.3 heavy shelter, 1.0 no shelter). Values above "
                f"1.0 correspond to wind acceleration (e.g., funneled "
                f"flow) and are physically unusual; verify this is "
                f"intended.",
                UserWarning,
                stacklevel=2,
            )

        # v1's coupling protocol skipped step 1 and started kinetics on step
        # 2. v2 reproduces that behavior so coupled TSM+Riverine runs match
        # v1 outputs exactly. See v3 hotstart contract: when resuming from a
        # saved dataset, ``from_hotstart`` flips this back to False so the
        # next substep is processed normally.
        self.__skip_first_time_step = True

        Process.__init__(self, time_step)

    @ProcessFactory.register("temperature")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "Temperature":
        return Temperature(**config)

    # ---------- v3 hotstart hooks ----------

    def to_hotstart(self) -> dict:
        """Snapshot per-process substep state for a hotstart save.

        Returns the keys this process owns inside the hotstart dataset's
        ``attrs`` mapping. Currently only ``__skip_first_time_step`` is
        process-internal; all other state is in the registry.
        """
        return {
            "temperature.skip_first_time_step": bool(self.__skip_first_time_step),
        }

    def from_hotstart(self, state: dict) -> None:
        """Restore process-internal substep state from a hotstart save.

        After a hotstart, the next substep is processed normally — the
        v1-coupling ``__skip_first_time_step`` semantic only applies to
        a fresh start. If the saved hotstart provides an explicit value,
        we honor it; otherwise we default to ``False`` (don't skip).
        """
        if "temperature.skip_first_time_step" in state:
            self.__skip_first_time_step = bool(
                state["temperature.skip_first_time_step"]
            )
        else:
            self.__skip_first_time_step = False

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Run one TSM substep at ``time``, updating ``water_temperature``."""

        if self.__skip_first_time_step:
            self.__skip_first_time_step = False
            return

        water_temperature = registry.get_at_time("water_temperature", time)
        surface_area = registry.get_at_time("wetted_surface_area", time)
        volume = registry.get_at_time("volume", time)

        cloudiness = registry.get_at_time("cloudiness", time)
        air_temperature = registry.get_at_time("air_temperature", time)
        solar_flux = registry.get_at_time("solar_radiation", time)
        wind_speed = registry.get_at_time("wind_speed", time)
        # Cache the per-cell wind-shelter forcing (if registered) so
        # ``wind_function`` and the equilibrium-T Newton-Raphson loop
        # see the same shelter value across a substep. ``None`` here
        # signals to fall back to the constructor scalar
        # ``self.wind_shelter`` -- this preserves back-compat for
        # applications that don't register the optional forcing.
        if "wind_shelter_coefficient" in registry:
            self._cached_shelter = registry.get_at_time(
                "wind_shelter_coefficient", time
            )
        else:
            self._cached_shelter = None
        atmospheric_pressure = registry.get_at_time("atmospheric_pressure", time)
        atmospheric_vapor_pressure = registry.get_at_time(
            "atmospheric_vapor_pressure", time
        )

        sediment_temperature = registry.get_at_time("sediment_temperature", time)
        sediment_thickness = registry.get_at_time("sediment_thickness", time)

        # F2 fix (audit 2026-05-05): use the with-factors helper so we
        # can apply the SAME depth-ramp factor and rate-cap clip ratio
        # to the sediment-side delta below, preserving the per-cell
        # water-sediment energy pair-cancellation invariant.
        (
            delta_water_temperature,
            ramp_factor,
            cap_clip_ratio,
            components,
        ) = self._temperature_change_with_factors(
            water_temperature=water_temperature,
            surface_area=surface_area,
            volume=volume,
            cloudiness=cloudiness,
            air_temperature=air_temperature,
            solar_flux=solar_flux,
            wind_speed=wind_speed,
            sediment_temperature=sediment_temperature,
            sediment_thickness=sediment_thickness,
            atmospheric_pressure=atmospheric_pressure,
            atmospheric_vapor_pressure=atmospheric_vapor_pressure,
        )

        # m19 fix / C5 (review-findings 2026-05-04; gap-analysis row N8
        # disposition "remove from process bodies"): the per-process
        # ``xr.where(volume > 0, ...)`` guard previously applied here
        # is now removed. ``Model.__apply_wet_mask`` is the single point
        # of dry-cell handling at the orchestration layer; running both
        # is dominated by the Model-level mask (per-process zeros the
        # delta first, then the Model-level mask overwrites the result
        # with NaN on dry cells), so the per-process branch is dead.
        # Trade-off: a run without a configured wet-mask that has
        # ``volume == 0`` cells will produce NaN/inf delta values
        # (division by zero in ``temperature_change``). This is the
        # documented behavior change between Phase 2 and Phase 3 and is
        # acceptable because (a) a run without a wet-mask should have a
        # uniformly wet mesh, and (b) a run with a wet-mask uses
        # ``Model.__apply_wet_mask`` to overwrite dry cells with NaN as
        # the design intent.
        updated_water_temperature = water_temperature + delta_water_temperature

        registry.set_at_time("water_temperature", time, updated_water_temperature)

        # Sediment temperature evolution. The Fortran TSM
        # (``HEC-RAS-WQ/RAS-1D-WQ/Kinetics Libraries/Temperature*/Source files/modTemperature.f90``)
        # gates both ``q_sediment`` and ``dT_sed/dt`` on the same
        # ``use_SedTemp`` flag, so the water-sediment exchange is
        # energy-conservative. The earlier Python ports (v1, v2)
        # dropped the dT_sed/dt update, breaking energy conservation
        # between the water and sediment heat reservoirs. v3 restores
        # the Fortran behavior; ``evolve_sediment_temperature=False``
        # is available for backward-compat against tests that depend
        # on a static sediment forcing.
        if self.use_sediment_temperature and self.evolve_sediment_temperature:
            delta_sediment_temperature = self.sediment_temperature_change(
                water_temperature=water_temperature,
                sediment_temperature=sediment_temperature,
                sediment_thickness=sediment_thickness,
            )
            # F2 fix (audit 2026-05-05): apply the SAME depth-ramp
            # factor and rate-cap clip ratio used on the water-side
            # delta. Without this, the water side absorbs only
            # ``q_sediment * ramp * clip_ratio`` of energy per substep
            # while the sediment side loses an unattenuated
            # ``q_sediment * dt`` worth of energy, breaking the
            # per-cell water-sediment cancellation invariant by a
            # factor of ``(1 - ramp * clip_ratio)``. Symmetric scaling
            # restores per-substep ``dE_water + dE_sediment = 0``
            # (verified in Phase R-3 +
            # ``test_water_sediment_energy_conservation_per_substep``;
            # the guarded-path extension is in
            # ``test_water_sediment_energy_conservation_under_ramp``
            # and ``..._under_cap``). Default values
            # ``q_net_depth_ramp_ref = 0.3 m`` and
            # ``dTdt_max_per_hour = 5 K/hr`` make ``ramp = 1.0`` and
            # ``clip_ratio = 1.0`` for typical deep cells, so the
            # multiplication is a no-op outside the guard regions.
            delta_sediment_temperature = (
                delta_sediment_temperature * ramp_factor * cap_clip_ratio
            )
            # m19 fix / C5 (review-findings 2026-05-04; gap-analysis row
            # N8 disposition "remove from process bodies"): the
            # per-process ``xr.where(volume > 0, ...)`` guard previously
            # applied here is now removed in favor of the orchestration
            # layer's ``Model.__apply_wet_mask``. See the matching
            # comment above the water-temperature branch for the full
            # rationale and trade-off.
            updated_sediment_temperature = (
                sediment_temperature + delta_sediment_temperature
            )
            registry.set_at_time(
                "sediment_temperature", time, updated_sediment_temperature
            )

        # Per-component flux diagnostics (audit 2026-05-05 open
        # question 3). Cache on the process instance for sibling
        # consumers, and write to the registry for any of these names
        # the user has pre-registered. Matching N2.run's "only written
        # when the registry knows the variable" pattern.
        self.q_sensible = components["q_sensible"]
        self.q_latent = components["q_latent"]
        self.q_longwave_up = components["q_longwave_up"]
        self.q_longwave_down = components["q_longwave_down"]
        self.q_solar = components["q_solar"]
        self.q_sediment = components["q_sediment"]
        self.q_net = components["q_net"]
        for diagnostic_name in (
            "q_sensible",
            "q_latent",
            "q_longwave_up",
            "q_longwave_down",
            "q_solar",
            "q_sediment",
            "q_net",
        ):
            if diagnostic_name in registry:
                registry.set_at_time(
                    diagnostic_name, time, components[diagnostic_name]
                )

        # Equilibrium-temperature diagnostic (audit 2026-05-05 open
        # question 2). Computed only when the user has pre-registered
        # ``equilibrium_temperature`` in the registry, so the
        # Newton-Raphson cost stays off the hot path otherwise. Cached
        # on the process when computed.
        if "equilibrium_temperature" in registry:
            self.equilibrium_temperature_c = self.equilibrium_temperature(
                cloudiness=cloudiness,
                air_temperature=air_temperature,
                solar_flux=solar_flux,
                wind_speed=wind_speed,
                atmospheric_pressure=atmospheric_pressure,
                atmospheric_vapor_pressure=atmospheric_vapor_pressure,
                sediment_temperature=sediment_temperature,
                sediment_thickness=sediment_thickness,
            )
            registry.set_at_time(
                "equilibrium_temperature", time, self.equilibrium_temperature_c
            )

    # ---------- Energy-balance fluxes ----------
    #
    # Sign convention (audit 2026-05-05 finding F-sign-convention; v1 +
    # Fortran-A + Fortran-B agreement). All flux methods return
    # **magnitudes** (or signed-by-temperature-gradient values for
    # ``flux_sensible`` and ``flux_sediment``). Signs are applied in
    # ``flux_net`` at composition time::
    #
    #     q_net = (
    #           q_sensible            # +/- by (T_air - T_water)
    #         + q_solar               # always positive (input)
    #         + q_sediment            # +/- by (T_sed - T_water)
    #         + q_atmospheric_LW      # always positive (incoming)
    #         - q_upwelling_LW        # always positive (outgoing)
    #         - q_latent              # always positive (evaporative)
    #     )
    #
    # This matches v1 ``tsm/processes.py:q_net`` and Fortran-A
    # ``modTemperature.f90:257`` exactly. **Do not pre-negate inside
    # any flux method** — it makes the composition asymmetric and
    # breaks the audit invariant. (Pre-2026-05-05, ``flux_upwelling_longwave``
    # and ``flux_latent_heat`` were pre-negated; refactored to
    # magnitudes-only as part of audit finding F-sign-convention.)

    def flux_upwelling_longwave(self, water_temperature: ArrayLike) -> xr.DataArray:
        """Upwelling longwave flux **magnitude** (W/m^2).

        Returns a positive value (Stefan-Boltzmann black-body emission
        from the water surface). Subtract in ``flux_net``: this energy
        leaves the water column.
        """
        return (
            constants.EMISSIVITY_WATER
            * constants.STEFAN_BOLTZMANN
            * conversions.celsius_to_kelvin(water_temperature) ** 4
        )

    def flux_atmospheric_longwave(
        self,
        air_temperature: ArrayLike,
        cloudiness: ArrayLike,
    ) -> ArrayLike:
        """Downwelling atmospheric longwave flux (W/m^2).

        Swinbank (1963) air-emissivity polynomial in absolute air
        temperature combined with a Bolz (1949) cloud-cover correction
        and the Stefan-Boltzmann radiation law. The temperature
        dependence of emissivity is already captured by the polynomial
        in ``T_K``.

        m2 fix (review-findings 2026-05-04): the prior docstring
        attributed the ``9.37e-6 * T_K^2`` form to "Brunt" and the
        cloud correction to "Kiehl", which is incorrect. Brunt's form
        is ``epsilon = a + b * sqrt(e_a)`` (vapor-pressure dependent,
        not ``T_K^2``); the ``T_K^2`` form is Swinbank's, and the
        ``(1 + 0.17 * C^2)`` cloud correction is Bolz's.
        """
        air_temperature_kelvin = conversions.celsius_to_kelvin(air_temperature)
        # Swinbank (1963) air-emissivity polynomial:
        # epsilon_a = 9.37e-6 * T_K^2; Bolz (1949) cloud-cover
        # correction (1 + 0.17 * C^2).
        emissivity_air = 9.37e-6 * air_temperature_kelvin**2
        return (
            emissivity_air
            * (1.0 + 0.17 * cloudiness**2)
            * constants.STEFAN_BOLTZMANN
            * air_temperature_kelvin**4
        )

    def flux_latent_heat(
        self,
        atmospheric_pressure: ArrayLike,
        water_temperature: ArrayLike,
        wind_speed: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
        richardson_function: ArrayLike,
    ) -> xr.DataArray:
        """Latent heat flux **magnitude** (W/m^2).

        Returns a positive value when ``e_sat > e_air`` (the typical
        evaporative regime — water evaporates and the latent heat
        leaves with it). Subtract in ``flux_net``. May be negative if
        ``e_sat < e_air`` (condensation onto the surface, an energy
        gain), in which case subtracting a negative still yields the
        correct sign.
        """
        return (
            0.622
            / atmospheric_pressure
            * self.latent_heat_vaporization(water_temperature)
            * self.water_density(water_temperature)
            * self.wind_function(wind_speed, richardson_function)
            * (
                self.saturation_vapor_pressure(water_temperature)
                - atmospheric_vapor_pressure
            )
        )

    def flux_sensible(
        self,
        water_temperature: ArrayLike,
        air_temperature: ArrayLike,
        wind_speed: ArrayLike,
        richardson_function: ArrayLike,
    ) -> ArrayLike:
        """Sensible heat flux (W/m^2): molecular and turbulent transfer
        between air and water surface."""
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        air_temperature_kelvin = conversions.celsius_to_kelvin(air_temperature)
        return (
            self.air_diffusivity_ratio
            * constants.AIR_SPECIFIC_HEAT
            * self.water_density(water_temperature)
            * self.wind_function(wind_speed, richardson_function)
            * (air_temperature_kelvin - water_temperature_kelvin)
        )

    def flux_sediment(
        self,
        water_temperature: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
    ) -> ArrayLike:
        """Sediment heat flux (W/m^2) under the active-layer model.

        The factor ``/ 0.5`` represents the sediment active-layer
        half-thickness convention: heat is exchanged across the upper
        half of ``sediment_thickness``. Matches v1's
        ``tsm/processes.py:q_sediment``. The ``/ 86400`` converts
        ``sediment_diffusivity`` from its declared units of m^2/day to
        m^2/s so the resulting flux is in W/m^2 = J/(m^2 s). (Audit
        2026-05-05 corrected the prior wording, which mis-stated the
        input units as m^2/s.)
        """
        if not self.use_sediment_temperature:
            # m6 fix (review-findings 2026-05-04): return a same-shape
            # zero-array rather than the Python scalar ``0.0`` so
            # downstream ``flux_net`` arithmetic preserves a consistent
            # dtype and dask-graph structure across enable/disable
            # configurations. ``xr.zeros_like`` is used for ``DataArray``
            # inputs (the production code path); scalar/ndarray inputs
            # (used in v1-port unit tests) fall back to ``np.zeros_like``
            # so this method remains type-stable on both surfaces.
            if isinstance(water_temperature, xr.DataArray):
                return xr.zeros_like(water_temperature)
            return np.zeros_like(water_temperature)
        # M2 fix (review-findings 2026-05-04): degenerate-layer guard.
        # ``sediment_thickness == 0`` (missing data, hotstart artifact)
        # or ``< 0`` (transport-coupling bug) would produce inf/NaN here
        # and poison ``water_temperature`` on every wet cell that depends
        # on this flux. Replace the divisor with 1.0 to keep the
        # expression finite, then return 0.0 for the degenerate cells.
        safe_thickness = xr.where(
            sediment_thickness > 0.0, sediment_thickness, 1.0
        )
        flux = (
            self.sediment_density
            * self.sediment_specific_heat
            * self.sediment_diffusivity
            / 0.5
            / safe_thickness
            * (sediment_temperature - water_temperature)
            / 86400.0
        )
        return xr.where(sediment_thickness > 0.0, flux, 0.0)

    def flux_components(
        self,
        water_temperature: ArrayLike,
        cloudiness: ArrayLike,
        air_temperature: ArrayLike,
        solar_flux: ArrayLike,
        wind_speed: ArrayLike,
        atmospheric_pressure: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
    ) -> dict:
        """Return all per-component fluxes plus the net (W/m^2).

        Audit 2026-05-05 open question 3 resolution. Returns a dict
        with the seven Fortran-A pathway outputs (matching
        ``modTemperature.f90`` ``q_*`` variables, the
        ``TempPathwayOutput`` set):

        * ``q_sensible``: sensible heat flux, signed by
          ``T_air - T_water`` (positive = air heats water).
        * ``q_latent``: latent heat flux returned as a **signed
          subtraction term**, NOT a positive magnitude. Phase H-10
          (2026-05-21) docstring clarification: in the evaporative
          regime (``e_sat > e_air``) the value is positive and
          subtracting it in ``q_net`` cools the water; in the
          condensation regime (``e_sat < e_air``) the value is
          negative and subtracting it (`-(-x) = +x`) correctly adds
          the latent heat of condensation back to the column. Sign
          convention preserved across both regimes; downstream
          consumers doing an offline mass-balance closure should use
          the *subtraction-term* form ``q_net = ... - q_latent``,
          not ``q_net = ... + |q_latent|``. Previously documented as
          "magnitude", which invited bugs in calibration plotters
          that took ``abs(q_latent)``.
        * ``q_longwave_up``: upwelling longwave **magnitude**.
          Subtracted in ``q_net``.
        * ``q_longwave_down``: downwelling atmospheric longwave
          **magnitude**. Added in ``q_net``.
        * ``q_solar``: solar input (passthrough of the registry
          forcing).
        * ``q_sediment``: sediment heat flux, signed by
          ``T_sed - T_water`` (positive = sediment heats water).
        * ``q_net``: composition
          ``q_sensible + q_solar + q_sediment + q_longwave_down -
          q_longwave_up - q_latent``.

        Useful for calibration and validation diagnostics. The values
        are also cached on the process instance after each
        ``Temperature.run`` substep (``self.q_sensible``,
        ``self.q_latent``, etc.) and are written to the registry by
        ``run`` for any of these names that the user has pre-registered.
        """
        mixing_ratio_air = self.mixing_ratio_air(
            atmospheric_vapor_pressure, atmospheric_pressure
        )
        density_air = self.density_air(
            atmospheric_pressure, air_temperature, mixing_ratio_air
        )
        density_air_sat = self.density_air_sat(water_temperature, atmospheric_pressure)

        _, richardson_function = self.richardson_number(
            wind_speed,
            density_air_sat=density_air_sat,
            density_air=density_air,
        )
        sensible = self.flux_sensible(
            water_temperature, air_temperature, wind_speed, richardson_function
        )
        latent = self.flux_latent_heat(
            water_temperature=water_temperature,
            atmospheric_pressure=atmospheric_pressure,
            wind_speed=wind_speed,
            atmospheric_vapor_pressure=atmospheric_vapor_pressure,
            richardson_function=richardson_function,
        )
        sediment = self.flux_sediment(
            water_temperature, sediment_temperature, sediment_thickness
        )
        longwave_down = self.flux_atmospheric_longwave(air_temperature, cloudiness)
        longwave_up = self.flux_upwelling_longwave(water_temperature)
        # Audit 2026-05-05 finding F-sign-convention: signs applied
        # here (composition-time) per the magnitudes-only convention
        # documented at the top of the energy-balance section. Matches
        # v1 ``tsm/processes.py:q_net`` and Fortran-A
        # ``modTemperature.f90:257``.
        net = (
            sensible + solar_flux + sediment + longwave_down - longwave_up - latent
        )
        return {
            "q_sensible": sensible,
            "q_latent": latent,
            "q_longwave_up": longwave_up,
            "q_longwave_down": longwave_down,
            "q_solar": solar_flux,
            "q_sediment": sediment,
            "q_net": net,
        }

    def flux_net(
        self,
        water_temperature: ArrayLike,
        cloudiness: ArrayLike,
        air_temperature: ArrayLike,
        solar_flux: ArrayLike,
        wind_speed: ArrayLike,
        atmospheric_pressure: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
    ) -> ArrayLike:
        """Net heat flux (W/m^2). Backward-compat thin wrapper.

        Returns the ``q_net`` value of :py:meth:`flux_components`. New
        code wanting per-component diagnostics should call
        ``flux_components`` directly.
        """
        components = self.flux_components(
            water_temperature=water_temperature,
            cloudiness=cloudiness,
            air_temperature=air_temperature,
            solar_flux=solar_flux,
            wind_speed=wind_speed,
            atmospheric_pressure=atmospheric_pressure,
            atmospheric_vapor_pressure=atmospheric_vapor_pressure,
            sediment_temperature=sediment_temperature,
            sediment_thickness=sediment_thickness,
        )
        return components["q_net"]

    # ---------- Equilibrium temperature (diagnostic) ----------

    def equilibrium_temperature(
        self,
        cloudiness: ArrayLike,
        air_temperature: ArrayLike,
        solar_flux: ArrayLike,
        wind_speed: ArrayLike,
        atmospheric_pressure: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
        max_iterations: int = 10,
        tolerance_kelvin: float = 0.01,
    ) -> ArrayLike:
        """Equilibrium water temperature (deg C) for the current met
        conditions.

        Newton-Raphson root-finding for ``q_net(T_eq) = 0``: the
        water temperature at which the surface net heat flux vanishes
        under the current meteorological forcing and sediment state.
        Diagnostic only — does not affect the model state.

        Mirrors Fortran-A
        ``modTemperature.f90:209-263``: starts from ``T_eq = T_air``
        and iterates ``T_eq_next = T_eq - q_net / dq_net/dT`` for up
        to ``max_iterations`` iterations or until the per-iteration
        change satisfies ``|T_eq_next - T_eq| < tolerance_kelvin`` for
        every cell. Per-iteration cost is one full flux evaluation
        plus four analytic derivative evaluations
        (``dq_longwave_up/dT``, ``dq_latent/dT``, ``dq_sensible/dT``,
        ``dq_sediment/dT``). At default ``max_iterations = 10`` the
        loop converges to machine precision for any realistic forcing
        from the air-temperature initial guess.

        Args:
            cloudiness, air_temperature, solar_flux, wind_speed,
            atmospheric_pressure, atmospheric_vapor_pressure,
            sediment_temperature, sediment_thickness:
                Same forcing variables consumed by
                :py:meth:`flux_components`. ``solar_flux`` enters as
                a constant in the equilibrium balance (independent of
                ``T_eq``).
            max_iterations: Newton-Raphson iteration cap. Default 10
                matches Fortran-A.
            tolerance_kelvin: Per-cell convergence threshold on the
                iterate change. Default 0.01 K matches Fortran-A.

        Returns:
            Equilibrium temperature in deg C, same shape as the
            forcing inputs.

        Audit 2026-05-05 open question 2 resolution.
        """
        # Initial guess: T_eq = T_air. Copy when array-like so we
        # don't mutate the input registry slice on the first
        # assignment; coerce to a Python float only for true scalar
        # inputs. Earlier the ``else`` branch unconditionally called
        # ``float()``, which crashed on multi-element ``np.ndarray``
        # forcings (Gemini review 2026-05-05, finding 1).
        if isinstance(air_temperature, xr.DataArray):
            teq_c = air_temperature.copy()
        elif isinstance(air_temperature, np.ndarray):
            teq_c = air_temperature.astype(float, copy=True)
        else:
            teq_c = float(air_temperature)

        # Pre-compute air-side quantities that do NOT depend on T_eq.
        # density_air uses air_temperature, so it is constant across
        # iterations. We keep mixing_ratio_air and density_air outside
        # the loop.
        mixing_ratio_air = self.mixing_ratio_air(
            atmospheric_vapor_pressure, atmospheric_pressure
        )
        density_air = self.density_air(
            atmospheric_pressure, air_temperature, mixing_ratio_air
        )
        atmospheric_lw = self.flux_atmospheric_longwave(
            air_temperature, cloudiness
        )

        for _ in range(max_iterations):
            teq_k = conversions.celsius_to_kelvin(teq_c)
            density_water = self.water_density(teq_c)
            esat = self.saturation_vapor_pressure(teq_c)

            # Sat-air density and Richardson are recomputed each
            # iteration because they depend on teq_c.
            density_air_sat = self.density_air_sat(teq_c, atmospheric_pressure)
            _, ri_function = self.richardson_number(
                wind_speed,
                density_air_sat=density_air_sat,
                density_air=density_air,
            )
            wind_fn = self.wind_function(wind_speed, ri_function)

            # Flux components at the current iterate.
            upwelling_lw = (
                constants.EMISSIVITY_WATER * constants.STEFAN_BOLTZMANN * teq_k**4
            )
            sensible = (
                self.air_diffusivity_ratio
                * constants.AIR_SPECIFIC_HEAT
                * density_water
                * wind_fn
                * (
                    conversions.celsius_to_kelvin(air_temperature) - teq_k
                )
            )
            lv = self.latent_heat_vaporization(teq_c)
            latent = (
                0.622
                / atmospheric_pressure
                * lv
                * density_water
                * wind_fn
                * (esat - atmospheric_vapor_pressure)
            )
            sediment = self.flux_sediment(teq_c, sediment_temperature, sediment_thickness)

            qnet = (
                sensible
                + solar_flux
                + sediment
                + atmospheric_lw
                - upwelling_lw
                - latent
            )

            # Analytic derivatives w.r.t. T_eq_K. Mirrors Fortran-A
            # ``modTemperature.f90:225-249``. Approximations matching
            # Fortran: ``Lv`` and ``density_water`` are taken as
            # weakly T-dependent and their T-derivatives are dropped
            # in ``d_latent_dT`` and ``d_sensible_dT``.
            d_upwelling_dT = (
                4.0
                * constants.EMISSIVITY_WATER
                * constants.STEFAN_BOLTZMANN
                * teq_k**3
            )
            d_sensible_dT = (
                -self.air_diffusivity_ratio
                * constants.AIR_SPECIFIC_HEAT
                * density_water
                * wind_fn
            )
            # de_sat/dT_K — derivative of the Brutsaert polynomial.
            d_esat_dT = (
                self.__A1
                + 2.0 * self.__A2 * teq_k
                + 3.0 * self.__A3 * teq_k**2
                + 4.0 * self.__A4 * teq_k**3
                + 5.0 * self.__A5 * teq_k**4
                + 6.0 * self.__A6 * teq_k**5
            )
            # Use ``wind_function`` (which applies the same shelter +
            # log-law transforms as the latent-heat flux itself) so
            # the derivative stays consistent with f(W) inside this
            # Newton-Raphson loop. ``wind_function`` returns
            # ``Ri * (a + b * effective_W^c) / 1e6`` so we get the
            # ri_function multiplication folded in here.
            d_latent_dT = (
                (0.622 / atmospheric_pressure)
                * lv
                * density_water
                * self.wind_function(wind_speed, ri_function)
                * d_esat_dT
            )
            # Sediment derivative gated on use_sediment_temperature so
            # disabled sediment runs see d_sediment/dT = 0 (matches
            # Fortran's ``if (use_SedTemp)`` branch).
            if self.use_sediment_temperature:
                safe_thickness = xr.where(
                    sediment_thickness > 0.0, sediment_thickness, 1.0
                )
                d_sediment_dT_active = -(
                    self.sediment_density
                    * self.sediment_specific_heat
                    * self.sediment_diffusivity
                    / 0.5
                    / safe_thickness
                    / 86400.0
                )
                d_sediment_dT = xr.where(
                    sediment_thickness > 0.0, d_sediment_dT_active, 0.0
                )
            else:
                d_sediment_dT = 0.0

            d_qnet_dT = (
                -d_upwelling_dT
                - d_latent_dT
                + d_sensible_dT
                + d_sediment_dT
            )

            teq_next = teq_c - qnet / d_qnet_dT
            # Vectorized convergence check: stop when every finite
            # cell is within tolerance. Fortran-A's loop uses the
            # same test at the abs-of-difference level. NaN cells
            # (dry cells, missing meteorology) are masked out of the
            # check so they don't block early exit -- ``NaN < tol``
            # evaluates to False, which would otherwise force the
            # loop to run ``max_iterations`` times even when every
            # finite cell has already converged (Gemini review
            # 2026-05-05, finding 3).
            if isinstance(teq_next, xr.DataArray):
                diff = np.abs(teq_next.values - np.asarray(teq_c))
            else:
                diff = np.abs(
                    np.asarray(teq_next) - np.asarray(teq_c)
                )
            valid = ~np.isnan(diff)
            if not np.any(valid):
                converged = True
            else:
                converged = bool(np.max(diff[valid]) < tolerance_kelvin)
            teq_c = teq_next
            if converged:
                break

        return teq_c

    # ---------- Thermodynamic state functions ----------

    def water_specific_heat(self, temperature: ArrayLike) -> ArrayLike:
        """Specific heat of water (J/kg/K) as a function of T (Celsius).

        m4 fix (review-findings 2026-05-04): ``np.select`` evaluates
        every comparison against NaN as False and would silently return
        the ``default=4178.0`` value, masking missing-data defects.
        Wrap the result with ``xr.where`` on ``np.isnan(temperature)``
        so NaN propagates rather than being replaced by a finite value.

        Note (Gemini review 2026-05-05, finding 4): ``np.select`` has
        no dask dispatch via xarray's ufunc registry. If ``temperature``
        is a dask-backed ``xr.DataArray``, this call materializes the
        chunk eagerly. v3's ``Model`` currently runs in-memory per
        chunk so this is moot, but if dask-backed temperature arrays
        ever become the production path, replace ``np.select`` with
        chained ``xr.where`` or ``xr.apply_ufunc(..., dask="allowed")``
        to preserve the computational graph.
        """
        result = np.select(
            condlist=[
                temperature <= 0.0,
                temperature <= 5.0,
                temperature <= 10.0,
                temperature <= 15.0,
                temperature <= 20.0,
                temperature <= 25.0,
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
        return xr.where(np.isnan(temperature), np.nan, result)

    def temperature_change(
        self,
        water_temperature: ArrayLike,
        surface_area: ArrayLike,
        volume: ArrayLike,
        cloudiness: ArrayLike,
        air_temperature: ArrayLike,
        solar_flux: ArrayLike,
        wind_speed: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
        atmospheric_pressure: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
    ) -> ArrayLike:
        """Per-substep change in water temperature (Celsius).

        The base form ``flux_net * surface_area * dt / (V * rho * cp)`` is
        equivalent to ``flux_net * dt / (depth * rho * cp)`` and goes
        numerically stiff at small depth. v3 applies two regularizations
        ported from v1 ``tsm/processes.py:dTdt_water_c``:

        1. **Depth ramp** on ``flux_net``:
           ``ramp = min(1, depth / q_net_depth_ramp_ref)``.
           Set ``q_net_depth_ramp_ref = 0.0`` to disable.
        2. **Rate cap** on the per-substep delta T:
           ``|dT| <= dTdt_max_per_hour * dt_hours``.
           Set ``dTdt_max_per_hour = float("inf")`` to disable.

        With both disabled, the formula reduces to v2's pre-merge form.

        Returns the per-substep delta T (Celsius). For the
        energy-conservative pairing with the sediment-side update
        (audit 2026-05-05 finding F2), see
        ``_temperature_change_with_factors`` and ``Temperature.run``.
        """
        delta, _ramp, _clip_ratio, _components = (
            self._temperature_change_with_factors(
                water_temperature=water_temperature,
                surface_area=surface_area,
                volume=volume,
                cloudiness=cloudiness,
                air_temperature=air_temperature,
                solar_flux=solar_flux,
                wind_speed=wind_speed,
                sediment_temperature=sediment_temperature,
                sediment_thickness=sediment_thickness,
                atmospheric_pressure=atmospheric_pressure,
                atmospheric_vapor_pressure=atmospheric_vapor_pressure,
            )
        )
        return delta

    def _temperature_change_with_factors(
        self,
        water_temperature: ArrayLike,
        surface_area: ArrayLike,
        volume: ArrayLike,
        cloudiness: ArrayLike,
        air_temperature: ArrayLike,
        solar_flux: ArrayLike,
        wind_speed: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
        atmospheric_pressure: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
    ) -> tuple:
        """Internal helper for ``Temperature.run`` (audit 2026-05-05 F2).

        Returns ``(delta_water, ramp, clip_ratio, components)``:

        - ``delta_water`` — the guarded per-substep water-T delta
          (Celsius), identical to ``temperature_change``'s public return.
        - ``ramp`` — the depth-ramp factor applied to ``flux_net``.
          ``min(1, depth / q_net_depth_ramp_ref)`` element-wise, or the
          scalar ``1.0`` when the ramp is disabled
          (``q_net_depth_ramp_ref = 0``).
        - ``clip_ratio`` — per-cell ratio of clipped to unclipped delta
          imposed by the rate cap. ``1.0`` for cells where the cap did
          not fire; ``cap / |delta_unclipped|`` for cells where it did;
          always in ``[0, 1]``.
        - ``components`` — the dict of per-component fluxes from
          :py:meth:`flux_components` (audit 2026-05-05 open question
          3). Returned alongside the delta/factors so ``run`` can
          cache and optionally write them to the registry without
          recomputing.

        The ``ramp`` and ``clip_ratio`` factors are exposed so
        ``Temperature.run`` can apply the same scaling to the
        sediment-side delta. Without that, the per-cell
        water-sediment energy pair-cancellation breaks by a factor of
        ``(1 - ramp * clip_ratio)`` whenever either guard is active —
        a one-way energy sink in shallow cells (audit finding F2,
        2026-05-05).
        """
        components = self.flux_components(
            water_temperature=water_temperature,
            cloudiness=cloudiness,
            air_temperature=air_temperature,
            solar_flux=solar_flux,
            wind_speed=wind_speed,
            sediment_temperature=sediment_temperature,
            sediment_thickness=sediment_thickness,
            atmospheric_pressure=atmospheric_pressure,
            atmospheric_vapor_pressure=atmospheric_vapor_pressure,
        )
        flux_net = components["q_net"]

        # Depth ramp.
        depth = xr.where(surface_area > 0.0, volume / surface_area, 0.0)
        # m9 fix (review-findings 2026-05-04): the prior second clamp
        # ``depth = xr.where(depth > 0.0, depth, 0.0)`` was redundant
        # against the wet/dry guard above; the only inputs it would
        # ever rewrite are negative ``volume`` or NaN ``volume`` --
        # both of which indicate a transport-coupling bug upstream.
        # Let NaN/inf propagate so the defect surfaces rather than
        # silently zeroing the cell. Wet/dry margin is already handled
        # by the ``surface_area > 0`` guard above.
        if self.q_net_depth_ramp_ref > 0.0:
            ramp = np.minimum(1.0, depth / self.q_net_depth_ramp_ref)
        else:
            ramp = 1.0
        flux_ramped = flux_net * ramp

        delta_unclipped = (
            flux_ramped
            * surface_area
            * self.time_step_seconds
            / (
                volume
                * self.water_density(water_temperature)
                * self.water_specific_heat(water_temperature)
            )
        )

        # Rate cap. With ``dTdt_max_per_hour = +inf`` both np.minimum and
        # np.maximum become identity ops, so the cap is a no-op.
        dt_hours = self.time_step_seconds / 3600.0
        cap = self.dTdt_max_per_hour * dt_hours
        delta_clipped = np.maximum(-cap, np.minimum(cap, delta_unclipped))

        # Per-cell rate-cap clip ratio: 1.0 for cells where the cap did
        # not fire, cap/|delta_unclipped| for cells where it did. Used
        # by ``Temperature.run`` to apply the same proportional clip to
        # the sediment-side delta (audit F2). The denominator-safety
        # ``xr.where(abs > 0, abs, 1.0)`` avoids division by zero on
        # cells where the unclipped delta was already 0; in that case
        # the cap doesn't change anything and the ratio degenerates to
        # 1.0 via the outer ``xr.where``.
        abs_unclipped = np.abs(delta_unclipped)
        clip_ratio = xr.where(
            abs_unclipped > cap,
            cap / xr.where(abs_unclipped > 0.0, abs_unclipped, 1.0),
            1.0,
        )

        return delta_clipped, ramp, clip_ratio, components

    def sediment_temperature_change(
        self,
        water_temperature: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
    ) -> ArrayLike:
        """Per-substep change in sediment temperature (Celsius).

        Mirrors the Fortran TSM update at
        ``HEC-RAS-WQ/RAS-1D-WQ/Kinetics Libraries/Temperature*/Source files/modTemperature.f90``::

            if (use_SedTemp) dTsedCdt = alphas(r) / (0.5 * h2(r) * h2(r)) * (TwaterC - TsedC)

        where ``alphas`` is in m^2/day. The relaxation time constant is
        ``tau = 0.5 * h2^2 / alphas``; with the default
        ``alphas = 0.0432 m^2/day`` and ``h2 = 0.1 m``, ``tau ~ 2.78 hours``.

        The water-side flux ``flux_sediment`` and this sediment-side
        update are paired so that the heat exchanged per unit area per
        substep is equal and opposite (energy conservation between the
        water and sediment heat reservoirs).
        """
        # M2 fix (review-findings 2026-05-04): degenerate-layer guard.
        # Mirrors the guard in ``flux_sediment``. Division by
        # ``sediment_thickness ** 2`` is even more sensitive to zero or
        # negative inputs; without the guard, a single bad cell can
        # produce NaN/inf in ``sediment_temperature`` for the rest of the
        # run. Returns 0.0 (no sediment temperature change) for cells
        # whose active layer is degenerate.
        safe_thickness = xr.where(
            sediment_thickness > 0.0, sediment_thickness, 1.0
        )
        delta = (
            self.sediment_diffusivity                       # m^2/day
            / (0.5 * safe_thickness**2)                     # 1/m^2
            * (water_temperature - sediment_temperature)    # Celsius
            * self.time_step_seconds                        # seconds
            / 86400.0                                       # seconds -> days
        )                                                   # = Celsius
        return xr.where(sediment_thickness > 0.0, delta, 0.0)

    def water_density(self, temperature: ArrayLike) -> ArrayLike:
        """Fresh-water density (kg/m^3) as a function of T (Celsius).

        Salt-water correction is out of scope for v3 1.0.0; revisit when
        salinity-coupled chemistry is added.
        """
        return 999.973 * (
            1.0
            - (
                ((temperature - 3.9863) ** 2 * (temperature + 288.9414))
                / (508929.2 * (temperature + 68.12963))
            )
        )

    def latent_heat_vaporization(self, water_temperature: ArrayLike) -> ArrayLike:
        """Latent heat of vaporization (J/kg).

        The polynomial coefficients (2.499999e6 J/kg intercept, -2385.74
        J/kg/K slope) are calibrated for **water temperature in Celsius**
        (Lv ~ 2.50 MJ/kg at 0 C, 2.45 MJ/kg at 20 C). v2's pre-merge form
        applied the polynomial to Kelvin, biasing Lv ~26-27% low across the
        typical surface-water range and underestimating evaporative
        cooling. v3 follows v1: the input is Celsius, so no Kelvin
        conversion is performed before the polynomial.
        """
        return 2499999.0 - 2385.74 * water_temperature

    def saturation_vapor_pressure(self, water_temperature: ArrayLike) -> ArrayLike:
        """Saturation vapor pressure (mb) at water temperature (Celsius).

        Brutsaert (1982) sixth-order polynomial in Kelvin.
        """
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        return self.__A0 + water_temperature_kelvin * (
            self.__A1
            + water_temperature_kelvin
            * (
                self.__A2
                + water_temperature_kelvin
                * (
                    self.__A3
                    + water_temperature_kelvin
                    * (
                        self.__A4
                        + water_temperature_kelvin
                        * (self.__A5 + water_temperature_kelvin * self.__A6)
                    )
                )
            )
        )

    def _compute_effective_wind(self, wind_speed: ArrayLike) -> ArrayLike:
        """Apply wind shelter and log-law height correction to raw wind.

        The composition ``raw_wind * shelter * height_factor`` mirrors
        CE-QUAL-W2's pre-correction at ``w2_4_unix.f90:480``:
        ``WIND2(I) = WIND(JW) * WSC(I) * log(2/Z0) / log(WINDH/Z0)``.

        * ``shelter`` is ``self._cached_shelter`` (per-cell registry
          forcing) when populated by :py:meth:`run`, else
          ``self.wind_shelter`` (constructor scalar).
        * ``height_factor = log(2 / surface_z0) / log(wind_input_height
          / surface_z0)`` when ``self.wind_input_height != 2.0``,
          else 1.0.

        Default constructor values (``wind_input_height=2.0``,
        ``wind_shelter=1.0``, no registry forcing) make this an
        identity transform.
        """
        if self._cached_shelter is not None:
            shelter = self._cached_shelter
        else:
            shelter = self.wind_shelter

        # Phase H-7 (2026-05-21): always apply the log-law factor;
        # at ``wind_input_height == 2.0`` exactly the factor reduces
        # to ``log(2/z0) / log(2/z0) == 1.0`` by construction, so the
        # previous exact-equality branch was dead optimization that
        # introduced a float-drift hazard (a value of 2.0000001 from
        # config-file round-tripping would silently bypass the
        # correction). Always evaluating the log-law produces 1.0 at
        # 2.0 and the correct factor at every other height.
        height_factor = (
            np.log(2.0 / self.surface_z0)
            / np.log(self.wind_input_height / self.surface_z0)
        )

        return wind_speed * shelter * height_factor

    def wind_function(
        self, wind_speed: ArrayLike, richardson_function: ArrayLike
    ) -> ArrayLike:
        """Wind function: stability-corrected coefficient relating wind
        speed to bulk transfer coefficients in the latent and sensible
        heat parameterizations.

        The form is ``f(W) = Ri_function * (a + b * W^c) / 1e6``,
        Edinger, Brady & Geyer (1974). The ``1e6`` divisor places the
        coefficients in convenient O(1) magnitudes; ``a, b, c`` are
        ``self.wind_a, self.wind_b, self.wind_c`` (defaults
        ``0.3, 1.5, 2.0``; see :py:meth:`__init__` docstring for the
        full reference and override path).

        ``wind_speed`` is treated as the raw wind magnitude registered
        by the application. Before evaluating the Edinger formula,
        :py:meth:`_compute_effective_wind` applies the wind-shelter
        coefficient and the log-law height correction, mirroring
        W2's ``WIND2(I)`` pre-correction at
        ``w2_4_unix.f90:480``. With default constructor values the
        composition reduces to the identity transform and the formula
        evaluates as ``Ri * (a + b * wind_speed^c) / 1e6`` exactly.
        """
        effective_wind = self._compute_effective_wind(wind_speed)
        return richardson_function * (
            (self.wind_a / 1_000_000)
            + (self.wind_b / 1_000_000) * (effective_wind**self.wind_c)
        )

    def mixing_ratio_air(
        self, atmospheric_vapor_pressure: ArrayLike, atmospheric_pressure: ArrayLike
    ) -> ArrayLike:
        """Air mixing ratio (unitless).

        Args:
            atmospheric_vapor_pressure: Atmospheric vapor pressure (mb)
            atmospheric_pressure: Atmospheric pressure (mb)

        Edge guard: when ``atmospheric_pressure <= atmospheric_vapor_pressure``
        the denominator is zero or negative; we return 0.0 for those cells.
        The exact-equality case matches upstream v2's posture; the
        ``e_air > P_air`` (negative-denominator) case is the C4 fix
        documented in ``design/clearwater_modules_v3_review_findings.md``.
        Without this extension, a data-entry error, mis-scaled forcing, or
        sensor noise near saturation that produces ``e_air > P_air`` would
        yield a negative mixing ratio, which propagates through
        ``density_air``'s ``(1 + r) / (1 + 1.61 r)`` factor and produces
        sign-flipped or near-singular air densities, poisoning every flux
        that depends on ``density_air``. Implemented via ``xr.where`` so
        the guard works for both scalars and multi-cell DataArray inputs.
        (Upstream v2's scalar ``if`` comparison raises ``ValueError`` for
        arrays of length > 1.)
        """
        # m1 / m7 cleanup (review-findings 2026-05-04): compute a
        # divide-safe denominator once rather than nesting an inverted
        # ``xr.where`` inside the outer guard. Behavior is identical to
        # the prior nested form -- the C4 fix (``denom <= 0`` returns
        # 0.0) is preserved -- but the predicate is no longer repeated.
        denominator = atmospheric_pressure - atmospheric_vapor_pressure
        # C4 fix: guard on ``> 0.0`` (rather than ``!= 0.0``) so the
        # pathological ``e_air > P_air`` case also returns 0.0 rather
        # than a negative mixing ratio. See review-findings C4.
        denominator_safe = xr.where(denominator > 0.0, denominator, 1.0)
        return xr.where(
            denominator > 0.0,
            0.622 * atmospheric_vapor_pressure / denominator_safe,
            0.0,
        )

    def density_air(
        self,
        atmospheric_pressure: ArrayLike,
        air_temperature: ArrayLike,
        mixing_ratio_air: ArrayLike,
    ) -> ArrayLike:
        """Air density (kg/m^3) from ideal-gas law with humidity correction.

        Args:
            atmospheric_pressure: Atmospheric pressure (mb)
            air_temperature: Air temperature (Celsius)
            mixing_ratio_air: Air mixing ratio (unitless)
        """
        air_temperature_kelvin = conversions.celsius_to_kelvin(air_temperature)
        return (
            0.348
            * (atmospheric_pressure / air_temperature_kelvin)
            * (1.0 + mixing_ratio_air)
            / (1.0 + 1.61 * mixing_ratio_air)
        )

    def density_air_sat(
        self, water_temperature: ArrayLike, atmospheric_pressure: ArrayLike
    ) -> ArrayLike:
        """Saturated-air density (kg/m^3) at the water-surface temperature.

        Edge guard (Gemini review 2026-05-05, finding 2): when
        ``atmospheric_pressure <= saturation_vapor_pressure`` the
        denominator of the saturation mixing ratio is zero or
        negative, producing a runaway negative mixing ratio that
        propagates a sign-flipped or singular density through every
        flux that depends on this quantity. This is the same C4 fix
        already applied at :py:meth:`mixing_ratio_air`; restore
        symmetric defense here. Trigger conditions are extreme water
        temperatures (e.g., post-hotstart stabilization with
        unphysically high T_water) or mis-scaled forcing (atm passed
        in atm rather than mb). Returns a 0.0 mixing ratio for
        degenerate cells, so the resulting density is still finite
        (and equals ``0.348 * P_atm / T_K`` -- dry-air density at the
        same pressure and temperature, which is the most defensible
        fallback).
        """
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        saturation_vapor_pressure = self.saturation_vapor_pressure(water_temperature)
        denominator = atmospheric_pressure - saturation_vapor_pressure
        denominator_safe = xr.where(denominator > 0.0, denominator, 1.0)
        mixing_ratio_sat = xr.where(
            denominator > 0.0,
            0.622 * saturation_vapor_pressure / denominator_safe,
            0.0,
        )
        return (
            0.348
            * (atmospheric_pressure / water_temperature_kelvin)
            * (1.0 + mixing_ratio_sat)
            / (1.0 + 1.61 * mixing_ratio_sat)
        )

    def richardson_number(
        self, wind_speed: ArrayLike, density_air_sat: ArrayLike, density_air: ArrayLike
    ) -> tuple[ArrayLike, ArrayLike]:
        """Richardson number and its associated stability function.

        The bulk Richardson number characterizes atmospheric stability
        from the buoyancy-to-shear ratio. Bounds are imposed at
        ``[-1, 2]`` so the piecewise stability function evaluates within
        its calibrated range.

        Regimes:
            -0.01 <= rn < 0.01   -> neutral (function = 1.0)
            rn < -0.01           -> unstable: (1 - 22*rn)^0.80
            rn > 0.01            -> stable:   (1 + 34*rn)^-0.80

        Returns:
            (richardson_number, richardson_function)

        NaN handling (M3, review-findings 2026-05-04): a NaN
        ``wind_speed``, ``density_air``, or ``density_air_sat`` (e.g.
        from missing meteorology forcing) produces a NaN
        ``richardson_number``, which propagates through the stability
        function and the dependent latent and sensible fluxes. This is
        intentional so the defect is visible at the output; silently
        clamping NaN to a finite value would mask bad forcing data.
        """
        # Sign convention (audit 2026-05-05). v3 stores
        # ``constants.GRAVITY = -9.806 m/s^2`` (negative) and uses it
        # directly in the formula below: the formula has no explicit
        # leading minus sign because the sign is carried by the
        # constant. v1 (`tsm/processes.py:150`) follows the same
        # convention. Fortran-A and Fortran-B store ``gravity = +9.806``
        # (positive, SI sign) and apply an explicit ``-gravity`` in the
        # formula. Both produce algebraically identical Richardson
        # numbers. **Do NOT 'normalize' GRAVITY to +9.806 without also
        # reintroducing the explicit ``-`` here**; doing so would
        # silently flip every Richardson regime. The earlier v2 source
        # carried a commented ``-1`` factor with a TODO that was
        # resolved per Jason Rutyna's January 2026 diff investigation
        # (commits ``8218962`` and ``7f4166a`` in the modules repo).
        #
        # M3 fix (review-findings 2026-05-04): suppress the
        # ``RuntimeWarning: divide by zero`` emitted when ``wind_speed``
        # is exactly zero. The resulting -inf is clamped to -1.0 below by
        # ``np.maximum``. ``invalid='ignore'`` also silences the 0/0
        # case which produces NaN. Matches v1's posture.
        with np.errstate(divide="ignore", invalid="ignore"):
            richardson_number: ArrayLike = (
                constants.GRAVITY
                * (density_air - density_air_sat)
                * 2.0
                / (density_air * (wind_speed**2.0))
            )

        # M3 fix: ``xr.where(NaN > 2.0, 2.0, NaN) -> NaN`` correctly
        # propagates NaN, but the equivalent chained ``xr.where`` form
        # below was rewritten to use ``np.minimum`` / ``np.maximum``
        # because both are NaN-aware on both branches and avoid eager
        # evaluation of the regime predicates. The clamp is a no-op for
        # NaN inputs (NaN survives), which is the desired visible-defect
        # behavior documented above.
        richardson_number = np.minimum(
            2.0, np.maximum(-1.0, richardson_number)
        )

        # Stability function. v1 used ``np.select``; v2's chained
        # ``xr.where`` is correct but slower. Performance optimization is
        # tracked as a follow-up.
        richardson_function: ArrayLike = 1.0
        # Neutral, rn < 0.
        richardson_function = xr.where(
            (richardson_number < 0.0) & (richardson_number >= -0.01),
            1.0,
            richardson_function,
        )
        # Unstable.
        richardson_function = xr.where(
            (richardson_number < 0.0) & (richardson_number < -0.01),
            (1.0 - 22.0 * richardson_number) ** 0.80,
            richardson_function,
        )
        # Neutral, rn > 0.
        richardson_function = xr.where(
            (richardson_number >= 0.0) & (richardson_number <= 0.01),
            1.0,
            richardson_function,
        )
        # Stable.
        richardson_function = xr.where(
            (richardson_number >= 0.0) & (richardson_number > 0.01),
            (1.0 + 34.0 * richardson_number) ** (-0.80),
            richardson_function,
        )

        # M3 fix (review-findings 2026-05-04): the chained ``xr.where``
        # comparisons above all evaluate False against NaN, so a NaN
        # ``richardson_number`` would silently leave ``richardson_function``
        # at its initial 1.0. That hides upstream NaN (typically from a
        # missing meteorology forcing) inside a finite stability function
        # and propagates a wrong-but-finite result through the sensible
        # and latent fluxes. Force NaN to survive: where ``richardson_number``
        # is NaN, ``richardson_function`` is NaN too. The visible-defect
        # contract from the M3 review fix.
        richardson_function = xr.where(
            np.isnan(richardson_number),
            np.nan,
            richardson_function,
        )
        return (richardson_number, richardson_function)
