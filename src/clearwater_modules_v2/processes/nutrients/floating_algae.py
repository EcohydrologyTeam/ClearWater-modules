from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from processes.base import Process
from clearwater_data.variables import VariableRegistry
from clearwater_data.custom_types import ArrayLike

import utils.conversions

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from model import Model


class FloatingAlgae(Process):
    variables = ["floating_algae", "solar_radiation"]

    def __init__(
        self,
        time_step_frequency: timedelta = timedelta(minutes=5),
        settling_velocity: float = 0.0,
        repiration_rate: float = 0.0,
        repiration_rate_correction_factor: float = 1.0,
        death_rate: float = 0.0,
        death_rate_correction_factor: float = 1.0,
        growth_rate_option: int = 1,
        growth_rate_max: float = 1.0,
        phosphorus_michaelis_menton_constant: float = 1.0,
        nitrogen_michaelis_menton_constant: float = 1.0,
        light_limitation_option: int = 1,
        light_limitation_constant: float = 1.0,
        light_attenuation_coefficient: float = 1.0,
    ) -> None:
        """
        Initialize the floating algae process.

        Parameters:
            time_step_frequency (timedelta): Time step frequency
            growth_rate_option (int): Growth rate option
                1 = Multiplicative
                2 = Limiting Nutrient
                3 = Harmonic Mean
            settling_velocity (float): Settling velocity of floating algae in units of m/d
            repiration_rate (float): Respiration rate of floating algae in units of ug-Chla/L/d
            repiration_rate_correction_factor (float): Respiration rate correction factor
            death_rate (float): Death rate of floating algae in units of ug-Chla/L/d
            death_rate_correction_factor (float): Death rate correction factor
        """
        self.settling_velocity = settling_velocity
        self.repiration_rate = repiration_rate
        self.repiration_rate_correction_factor = repiration_rate_correction_factor
        self.death_rate = death_rate
        self.death_rate_correction_factor = death_rate_correction_factor
        self.growth_rate_option = growth_rate_option
        self.growth_rate_max = growth_rate_max
        self.phosphorus_michaelis_menton_constant = phosphorus_michaelis_menton_constant
        self.nitrogen_michaelis_menton_constant = nitrogen_michaelis_menton_constant
        self.light_limitation_option = light_limitation_option
        self.light_limitation_constant = light_limitation_constant
        self.light_attenuation_coefficient = light_attenuation_coefficient
        Process.__init__(self, time_step_frequency)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        # check if there is nitrogen process and set flags according
        self.use_nitrate = True
        self.use_ammonium = True

        # check if there is a phosphorus process and set flags according
        # TODO: implement
        self.use_phosphate = True

    def run(self, time_step: datetime, registry: VariableRegistry) -> None:
        """
        Run the floating algae process.
        """

        algae = registry.get_at_time("algae_floating", time_step)
        ammonium = registry.get_at_time("ammonium", time_step)
        nitrate = registry.get_at_time("nitrate", time_step)
        phosphorus_total_inorganic = registry.get_at_time(
            "phosphorus_total_inorganic", time_step
        )
        depth = registry.get_at_time("depth", time_step)
        water_temperature = registry.get_at_time("water_temperature", time_step)
        solar = registry.get_at_time("solar_radiation", time_step)

        # get rate of change
        rate = self.rate(
            algae=algae,
            depth=depth,
            water_temperature=water_temperature,
            phosphorus_total_inorganic=phosphorus_total_inorganic,
            phosphate_fraction_dissolved=0.5,  # TODO: figure out where this value should be coming from
            ammonium=ammonium,
            nitrate=nitrate,
            solar=solar,
        )

        # update algae
        # rate is in ug-Chla/L/d
        # 86400 converts rate to seconds
        algae = 0 + algae * rate * self.time_step_frequency.total_seconds() * 86400

    def rate(
        self,
        algae: ArrayLike,
        depth: ArrayLike,
        water_temperature: ArrayLike,
        phosphorus_total_inorganic: ArrayLike,
        phosphate_fraction_dissolved: ArrayLike,
        ammonium: ArrayLike,
        nitrate: ArrayLike,
        solar: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the rate of change of floating algae.
        """
        # growth limiting factors
        # TODO: confirm arguments are correct
        limit_phosphorus = self.limit_phosphorus(
            concentration=phosphorus_total_inorganic,
            fraction_dissolved=phosphate_fraction_dissolved,
        )
        limit_nitrogen = self.limit_nitrogen(ammonium=ammonium, nitrate=nitrate)
        limit_light = self.limit_light(depth=depth, surface_light_intensity=solar)

        return (
            self.rate_growth(algae, limit_phosphorus, limit_nitrogen, limit_light)
            - self.rate_death(algae, water_temperature)
            - self.rate_respiration(algae, water_temperature)
            - self.rate_settling(algae, depth)
        )

    def rate_growth(
        self,
        algae: ArrayLike,
        limit_phosphorus: ArrayLike,
        limit_nitrogen: ArrayLike,
        limit_light: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the rate of growth of floating algae.
        """

        # Multiplicative
        if self.growth_rate_option == 1:
            rate = (
                self.growth_rate_max * limit_phosphorus * limit_nitrogen * limit_light
            )
        # Limiting Nutrient
        elif self.growth_rate_option == 2:
            rate = xr.where(
                limit_phosphorus > limit_nitrogen,
                self.growth_rate_max * limit_nitrogen * limit_light,
                self.growth_rate_max * limit_phosphorus * limit_light,
            )
        # Harmonic Mean
        elif self.growth_rate_option == 3:
            rate = xr.where(
                limit_nitrogen == 0.0 or limit_phosphorus == 0.0,
                0,
                self.growth_rate_max
                * limit_light
                * 2.0
                / (1.0 / limit_nitrogen + 1.0 / limit_phosphorus),
            )
        else:
            raise ValueError("Invalid growth rate option")

        return rate * algae

    def rate_death(self, algae: ArrayLike, water_temperature: ArrayLike) -> ArrayLike:
        """
        Compute the rate of death of floating algae.
        """
        corrected_death_rate = utils.conversions.arrhenius_correction(
            water_temperature,
            self.death_rate,
            self.death_rate_correction_factor,
        )
        return algae * corrected_death_rate

    def rate_respiration(
        self, algae: ArrayLike, water_temperature: ArrayLike
    ) -> ArrayLike:
        """
        Compute the rate of respiration of floating algae.
        """
        corrected_respiration_rate = utils.conversions.arrhenius_correction(
            water_temperature,
            self.repiration_rate,
            self.repiration_rate_correction_factor,
        )
        return algae * corrected_respiration_rate

    def rate_settling(self, algae: ArrayLike, depth: ArrayLike) -> ArrayLike:
        """
        Compute the rate of settling of floating algae.
        """
        return algae / depth * self.settling_velocity

    def limit_phosphorus(
        self,
        concentration: ArrayLike,
        fraction_dissolved: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the limiting phosphorus for floating algae.
        """

        # if we are not modeling phosphate assume this is not limiting
        if not self.use_phosphate:
            return 1.0

        rate_raw = (
            fraction_dissolved
            * concentration
            / (
                self.phosphorus_michaelis_menton_constant
                + fraction_dissolved * concentration
            )
        )

        # TODO: see if we can combine these conditionals
        # any nan's are not limiting by definition (no p = no limit)
        rate = xr.where(rate_raw == np.nan, 0.0, rate_raw)
        # any rates > 1 are limiting
        rate = xr.where(rate > 1, 1, rate)

        return rate

    def limit_nitrogen(
        self,
        nitrate: ArrayLike,
        ammonium: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the limiting nitrogen for floating algae.
        """
        if not self.use_nitrate and not self.use_ammonium:
            return 1.0

        n_concentration = nitrate if self.use_nitrate else 0.0
        n_concentration += ammonium if self.use_ammonium else 0.0

        rate_raw = n_concentration / (
            self.nitrogen_michaelis_menton_constant + n_concentration
        )

        # TODO: see if we can combine these conditionals
        # any nan's are not limiting by definition (no p = no limit)
        rate = xr.where(rate_raw == np.nan, 0.0, rate_raw)
        # any rates > 1 are limiting
        rate = xr.where(rate > 1, 1, rate)

        return rate

    def limit_light(
        self,
        depth: ArrayLike,
        surface_light_intensity: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the limiting light for floating algae.
        """

        # Half-saturation light limitation
        if self.light_limitation_option == 1:
            raw_rate = (
                (1.0 / (self.light_attenuation_coefficient * depth))
                * np.log(self.light_limitation_constant + surface_light_intensity)
                / (
                    self.light_limitation_constant
                    + surface_light_intensity
                    * np.exp(-(self.light_attenuation_coefficient * depth))
                )
            )
        # Smith Model
        elif self.light_limitation_option == 2:
            raw_rate = xr.where(
                abs(self.light_limitation_constant) < 1e-10,
                1,
                (
                    (1.0 / (self.light_attenuation_coefficient * depth))
                    * np.log(
                        (
                            surface_light_intensity / self.light_limitation_constant
                            + (
                                (
                                    1.0
                                    + (
                                        surface_light_intensity
                                        / self.light_limitation_constant
                                    )
                                    ** 2.0
                                )
                                ** 0.5
                            )
                        )
                        / (
                            surface_light_intensity
                            * np.exp(-self.light_attenuation_coefficient * depth)
                            / self.light_limitation_constant
                            + (
                                (
                                    1.0
                                    + (
                                        surface_light_intensity
                                        * np.exp(
                                            -self.light_attenuation_coefficient * depth
                                        )
                                        / self.light_limitation_constant
                                    )
                                    ** 2.0
                                )
                                ** 0.5
                            )
                        )
                    )
                ),
            )
        # Steele Model
        elif self.light_limitation_option == 3:
            raw_rate = xr.where(
                abs(self.light_limitation_constant) < 1e-10,
                0,
                (
                    (2.718 / (self.light_attenuation_coefficient * depth))
                    * (
                        np.exp(
                            -surface_light_intensity
                            / self.light_limitation_constant
                            * np.exp(-self.light_attenuation_coefficient * depth)
                        )
                        - np.exp(
                            -surface_light_intensity / self.light_limitation_constant
                        )
                    )
                ),
            )
        else:
            raise ValueError("Invalid light limitation option")

        rate = xr.where(raw_rate > 1, 1.0, raw_rate)
        rate = xr.where(rate < 0, 0.0, rate)
        return rate
