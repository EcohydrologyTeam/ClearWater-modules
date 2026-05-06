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
        wind_c: float = 3.0,
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
                :py:meth:`wind_function`). Defaults are ``a = 0.3``,
                ``b = 1.5``, ``c = 3.0`` — the calibration values
                inherited from v1 ``clearwater_modules.tsm.constants``,
                used across QUAL2K, CE-QUAL-W2, and HEC-RAS-WQ
                derivatives. Pass any subset to override per-instance;
                YAML configs may also override via the ``wind_a /
                wind_b / wind_c`` keys at ``init_from_file`` time.

                Reference: Edinger, J.E., D.K. Brady, and J.C. Geyer
                (1974), *Heat exchange and transport in the
                environment*, Report 14, Cooling Water Discharge
                Research Project (RP-49), Electric Power Research
                Institute, Palo Alto, CA, 125 pp. (Audit 2026-05-05
                open question 4: prior code lineage carried no
                citation; recovered and added here.)
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
        # defaults are ``0.3 / 1.5 / 3.0`` per the constructor docstring;
        # callers may override per-instance or via YAML at config time.
        self.wind_a = wind_a
        self.wind_b = wind_b
        self.wind_c = wind_c
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
        """Net heat flux (W/m^2)."""
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
        atmospheric = self.flux_atmospheric_longwave(air_temperature, cloudiness)
        upwelling = self.flux_upwelling_longwave(water_temperature)
        # Audit 2026-05-05 finding F-sign-convention: signs applied
        # here (composition-time) per the magnitudes-only convention
        # documented at the top of the energy-balance section. Matches
        # v1 ``tsm/processes.py:q_net`` and Fortran-A
        # ``modTemperature.f90:257``.
        return sensible + solar_flux + sediment + atmospheric - upwelling - latent

    # ---------- Thermodynamic state functions ----------

    def water_specific_heat(self, temperature: ArrayLike) -> ArrayLike:
        """Specific heat of water (J/kg/K) as a function of T (Celsius).

        m4 fix (review-findings 2026-05-04): ``np.select`` evaluates
        every comparison against NaN as False and would silently return
        the ``default=4178.0`` value, masking missing-data defects.
        Wrap the result with ``xr.where`` on ``np.isnan(temperature)``
        so NaN propagates rather than being replaced by a finite value.
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
        delta, _ramp, _clip_ratio = self._temperature_change_with_factors(
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

        Returns ``(delta_water, ramp, clip_ratio)``:

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

        The ``ramp`` and ``clip_ratio`` factors are exposed so
        ``Temperature.run`` can apply the same scaling to the
        sediment-side delta. Without that, the per-cell
        water-sediment energy pair-cancellation breaks by a factor of
        ``(1 - ramp * clip_ratio)`` whenever either guard is active —
        a one-way energy sink in shallow cells (audit finding F2,
        2026-05-05).
        """
        flux_net = self.flux_net(
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

        return delta_clipped, ramp, clip_ratio

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
        ``0.3, 1.5, 3.0`` from the v1 calibration; see
        :py:meth:`__init__` docstring for the citation and override
        path).
        """
        return richardson_function * (
            (self.wind_a / 1_000_000)
            + (self.wind_b / 1_000_000) * (wind_speed**self.wind_c)
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
        """Saturated-air density (kg/m^3) at the water-surface temperature."""
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        saturation_vapor_pressure = self.saturation_vapor_pressure(water_temperature)
        mixing_ratio_sat = (
            0.622
            * saturation_vapor_pressure
            / (atmospheric_pressure - saturation_vapor_pressure)
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
