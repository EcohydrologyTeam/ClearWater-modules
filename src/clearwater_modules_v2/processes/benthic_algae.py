from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_data.variables import VariableRegistry
from clearwater_data.custom_types import ArrayLike

from .floating_algae import FloatingAlgae

from ..utils.conversions import arrhenius_correction

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from model import Model


class BenthicAlgae(FloatingAlgae):
    variables = ["benthic_algae", "solar_radiation"]

    def __init__(
        self,
        *args,
        density_michaelis_menton_constant: float = 1.0,
        **kwargs,
    ) -> None:
        """
        Initialize the floating algae process.

        Parameters:
            time_step_frequency (timedelta): Time step frequency
            growth_rate_option (int): Growth rate option
                1 = Multiplicative
                2 = Limiting Nutrient
            density_michaelis_menton_constant (float): Michaelis-Menton constant for density
        """
        self.density_michaelis_menton_constant = density_michaelis_menton_constant
        FloatingAlgae.__init__(self, *args, **kwargs)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        # check if there is nitrogen process and set flags according
        self.use_nitrate = True
        self.use_ammonium = True

        # check if there is a phosphorus process and set flags according
        # TODO: implement
        self.use_phosphate = True

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
        limit_phosphorus = self.limit_phosphorus(
            concentration=phosphorus_total_inorganic,
            fraction_dissolved=phosphate_fraction_dissolved,
        )
        limit_nitrogen = self.limit_nitrogen(ammonium=ammonium, nitrate=nitrate)
        limit_light = self.limit_light(
            algae=algae,
            depth=depth,
            surface_light_intensity=solar,
        )
        limit_density = self.limit_density(algae=algae)

        return (
            self.rate_growth(
                algae,
                water_temperature,
                limit_phosphorus,
                limit_nitrogen,
                limit_light,
                limit_density,
            )
            - self.rate_death(algae, water_temperature)
            - self.rate_respiration(algae, water_temperature)
        )

    def rate_growth(
        self,
        algae: ArrayLike,
        water_temperature: ArrayLike,
        limit_phosphorus: ArrayLike,
        limit_nitrogen: ArrayLike,
        limit_light: ArrayLike,
        limit_density: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the rate of growth of floating algae.
        """

        growth_rate = arrhenius_correction(
            water_temperature,
            self.growth_rate_max,
            self.growth_rate_correction,
        )

        # Multiplicative
        if self.growth_rate_option == 1:
            rate = (
                growth_rate
                * limit_phosphorus
                * limit_nitrogen
                * limit_light
                * limit_density
            )
        # Limiting Nutrient
        elif self.growth_rate_option == 2:
            rate = xr.where(
                limit_phosphorus > limit_nitrogen,
                growth_rate * limit_nitrogen * limit_light * limit_density,
                growth_rate * limit_phosphorus * limit_light * limit_density,
            )
        else:
            raise ValueError("Invalid growth rate option")

        return rate * algae

    def limit_light(
        self,
        algae: ArrayLike,
        depth: ArrayLike,
        surface_light_intensity: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the limiting light for floating algae.
        """

        light_at_depth_coefficent = np.exp(-self.light_attenuation_coefficient * depth)

        # Half-saturation light limitation
        if self.light_limitation_option == 1:
            raw_rate = (
                surface_light_intensity
                * light_at_depth_coefficent
                / (
                    self.light_limitation_constant
                    + surface_light_intensity * light_at_depth_coefficent
                )
            )
        # Smith Model
        elif self.light_limitation_option == 2:
            raw_rate = (
                surface_light_intensity
                * light_at_depth_coefficent
                / (
                    (
                        self.light_limitation_constant**2
                        + (surface_light_intensity * light_at_depth_coefficent) ** 2
                    )
                    ** 0.5
                )
            )
        # Steele Model
        elif self.light_limitation_option == 3:
            raw_rate = xr.where(
                abs(self.light_limitation_constant) < 1e-10,
                0,
                (
                    surface_light_intensity
                    * light_at_depth_coefficent
                    / (
                        self.light_limitation_constant
                        * np.exp(
                            1.0
                            - surface_light_intensity
                            * light_at_depth_coefficent
                            / self.light_limitation_constant
                        )
                    )
                ),
            )
        else:
            raise ValueError("Invalid light limitation option")

        # conditions where growth is limited to zero
        # no algae present
        rate = xr.where(algae <= 0, 0, raw_rate)
        # light cannot pententrate deep enough to benethic layer
        rate = xr.where(light_at_depth_coefficent <= 0, 0, rate)
        # there is no light
        rate = xr.where(surface_light_intensity <= 0, 0.0, rate)
        return rate

    def limit_density(self, algae: ArrayLike) -> ArrayLike:
        """
        Compute the limiting density for floating algae.
        """
        limit_raw = 1.0 - (algae / (algae + self.density_michaelis_menton_constant))
        limit = xr.where(limit_raw > 1, 1, limit_raw)
        limit = xr.where(limit_raw == np.nan, 0, limit)

        return limit
