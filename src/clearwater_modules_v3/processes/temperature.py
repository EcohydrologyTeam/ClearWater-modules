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

    def __init__(
        self,
        wind_a: float,
        wind_b: float,
        wind_c: float,
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
            wind_a, wind_b, wind_c: Wind-function parameters.
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
        # Wind-function parameters: physical interpretation under review;
        # legacy values come from the original NSM/TSM implementation.
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

        delta_water_temperature = self.temperature_change(
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

        # Per-process dry-cell guard. Phase 3 makes this redundant by adding
        # registry-level wet-mask handling at the orchestration layer; it is
        # kept here so v3 remains correct on its own without that change.
        delta_water_temperature = xr.where(volume > 0, delta_water_temperature, 0)
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
            # No water in contact -> no heat exchange -> no sediment T change.
            delta_sediment_temperature = xr.where(
                volume > 0, delta_sediment_temperature, 0
            )
            updated_sediment_temperature = (
                sediment_temperature + delta_sediment_temperature
            )
            registry.set_at_time(
                "sediment_temperature", time, updated_sediment_temperature
            )

    # ---------- Energy-balance fluxes ----------

    def flux_upwelling_longwave(self, water_temperature: ArrayLike) -> xr.DataArray:
        """Upwelling longwave flux (W/m^2). Negative because energy leaves
        the water column."""
        return (
            -constants.EMISSIVITY_WATER
            * constants.STEFAN_BOLTZMANN
            * conversions.celsius_to_kelvin(water_temperature) ** 4
        )

    def flux_atmospheric_longwave(
        self,
        air_temperature: ArrayLike,
        cloudiness: ArrayLike,
    ) -> ArrayLike:
        """Downwelling atmospheric longwave flux (W/m^2).

        Brunt-style emissivity (``9.37e-6 * T_K^2``) with a Kiehl
        cloud-cover correction ``(1 + 0.17 * C^2)`` and the standard
        Stefan-Boltzmann radiation. The temperature dependence of
        emissivity is already captured by the polynomial in ``T_K``.
        """
        air_temperature_kelvin = conversions.celsius_to_kelvin(air_temperature)
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
        """Latent heat flux (W/m^2). Negative because evaporative heat
        loss removes energy from the water column."""
        return (
            -0.622
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
        ``tsm/processes.py:q_sediment``. The ``/ 86400`` converts the
        product of diffusivity (m^2/s) and bulk thermal capacity into the
        per-substep flux units expected by the energy balance.
        """
        if not self.use_sediment_temperature:
            return 0.0
        return (
            self.sediment_density
            * self.sediment_specific_heat
            * self.sediment_diffusivity
            / 0.5
            / sediment_thickness
            * (sediment_temperature - water_temperature)
            / 86400.0
        )

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
        return sensible + solar_flux + sediment + atmospheric + upwelling + latent

    # ---------- Thermodynamic state functions ----------

    def water_specific_heat(self, temperature: ArrayLike) -> ArrayLike:
        """Specific heat of water (J/kg/K) as a function of T (Celsius)."""
        return np.select(
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
        depth = xr.where(depth > 0.0, depth, 0.0)
        if self.q_net_depth_ramp_ref > 0.0:
            ramp = np.minimum(1.0, depth / self.q_net_depth_ramp_ref)
        else:
            ramp = 1.0
        flux_ramped = flux_net * ramp

        delta_temperature = (
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
        delta_temperature = np.maximum(-cap, np.minimum(cap, delta_temperature))
        return delta_temperature

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
        return (
            self.sediment_diffusivity                       # m^2/day
            / (0.5 * sediment_thickness**2)                 # 1/m^2
            * (water_temperature - sediment_temperature)    # Celsius
            * self.time_step_seconds                        # seconds
            / 86400.0                                       # seconds -> days
        )                                                   # = Celsius

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
        heat parameterizations."""
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

        Edge guard: when ``atmospheric_pressure == atmospheric_vapor_pressure``
        the denominator is zero; we return 0.0 for those cells. Implemented
        via ``xr.where`` so the guard works for both scalars and multi-cell
        DataArray inputs. (Upstream v2's scalar ``if`` comparison raises
        ``ValueError`` for arrays of length > 1.)
        """
        denominator = atmospheric_pressure - atmospheric_vapor_pressure
        return xr.where(
            denominator == 0.0,
            0.0,
            0.622 * atmospheric_vapor_pressure / xr.where(denominator == 0.0, 1.0, denominator),
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
        """
        # The original v2 source carried a commented "-1" multiplier with
        # a TODO. v1 has no such factor and Jason Rutyna's January 2026
        # diff investigation (see commits 8218962 and 7f4166a in the
        # modules repo) concluded the leading "-1" should not be there.
        # v3 matches v1 exactly: no leading "-1".
        richardson_number: ArrayLike = (
            constants.GRAVITY
            * (density_air - density_air_sat)
            * 2.0
            / (density_air * (wind_speed**2.0))
        )

        richardson_number = xr.where(richardson_number > 2.0, 2.0, richardson_number)
        richardson_number = xr.where(richardson_number < -1.0, -1.0, richardson_number)

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
        return (richardson_number, richardson_function)
