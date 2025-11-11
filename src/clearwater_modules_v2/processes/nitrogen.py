from datetime import timedelta, datetime

import numpy as np
import xarray as xr

from processes.base import Process, ProcessFactory
from clearwater_data.variables import VariableRegistry
from clearwater_data.custom_types import ArrayLike
from clearwater_modules_v2.utils import conversions

arrhenius_correction = conversions.arrhenius_correction

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from model import Model


class Nitrogen(Process):
    variables = ["nitrate", "ammonium"]

    def __init__(
        self,
        time_step: timedelta = timedelta(minutes=5),
        denitrification_rate: ArrayLike = 1.0,
        denitrification_theta: ArrayLike = 1.0,
        nitrification_rate: ArrayLike = 1.0,
        nitrification_theta: ArrayLike = 1.0,
        sediment_denitrification_rate: ArrayLike = 1.0,
        sediment_denitrification_theta: ArrayLike = 1.0,
        sediment_ammonium_release_rate: ArrayLike = 1.0,
        sediment_ammonium_release_theta: ArrayLike = 1.0,
        floating_algae_preference_factor: ArrayLike = 0.5,
        settling_velocity: ArrayLike = 1.0,
        death_rate: ArrayLike = 1.0,
        float_algea_faction_uptake_from_nitrate: ArrayLike = 1.0,
    ) -> None:
        Process.__init__(self, time_step)
        self.denitrification_rate = denitrification_rate
        self.denitrification_theta = denitrification_theta
        self.nitrification_rate = nitrification_rate
        self.nitrification_theta = nitrification_theta
        self.sediment_denitrification_rate = sediment_denitrification_rate
        self.sediment_denitrification_theta = sediment_denitrification_theta
        self.sediment_ammonium_release_rate = sediment_ammonium_release_rate
        self.sediment_ammonium_release_theta = sediment_ammonium_release_theta
        self.floating_algae_preference_factor = floating_algae_preference_factor
        self.settling_velocity = settling_velocity
        # TODO: this should come from floating algae process
        self.death_rate = death_rate
        self.float_algea_faction_uptake_from_nitrate = (
            float_algea_faction_uptake_from_nitrate
        )

    @ProcessFactory.register("nitrogen")
    @staticmethod
    def from_config(config: dict) -> "Nitrogen":
        return Nitrogen(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        self.use_nitrate = True
        self.use_ammonium = True

        if model.has_process("FloatingAlgae"):
            self.use_floating_algae = True
            self.floating_algae_process = model.get_process("FloatingAlgae")
        if model.has_process("BenthicAlgae"):
            self.use_benthic_algae = True
            self.benthic_algae_process = model.get_process("BenthicAlgae")

    def run(self, time_step: datetime, registry: VariableRegistry) -> None:
        # pull data from regsitry
        nitrate = registry.get_at_time("nitrate", time_step)
        ammonium = registry.get_at_time("ammonium", time_step)
        temperature = registry.get_at_time("water_temperature", time_step)
        depth = registry.get_at_time("depth", time_step)
        dissolved_oxygen = registry.get_at_time("oxygen_dissolved", time_step)

        # update ammonium
        ammonium_rate = self.change_ammonium(
            nitrate,
            ammonium,
            temperature,
            depth,
            dissolved_oxygen,
        )
        ammonium = 0 + ammonium * ammonium_rate * self.time_step.total_seconds()
        ammonium = xr.where(ammonium < 0, 0, ammonium)

        # update nitrate
        nitrate_rate = self.change_nitrate(
            nitrate,
            ammonium,
            temperature,
            depth,
            dissolved_oxygen,
        )
        nitrate = 0 + nitrate * nitrate_rate * self.time_step_frequency.total_seconds()
        nitrate = xr.where(nitrate < 0, 0, nitrate)

    def change_ammonium(
        self,
        nitrate: ArrayLike,
        ammonium: ArrayLike,
        temperature: ArrayLike,
        depth: ArrayLike,
        dissolved_oxygen: ArrayLike,
    ) -> None:
        if not self.use_ammonium:
            return 0

        rate = (
            self.ammonium_decay(
                ammonium,
                temperature,
            )
            - self.ammonium_nitrification(
                ammonium,
                temperature,
                # TODO
                # nitrification_inhibition,
            )
            + self.ammonium_from_bed(
                depth=depth,
                temperature=temperature,
            )
            + self.ammonium_floating_respiration()
            - self.ammonium_floating_growth()
            + self.ammonium_benthic_respiration()
            - self.ammonium_benthic_growth()
        )

        # Replace nan's with 0's
        rate = xr.where(rate == np.nan, 0, rate)
        return rate

    def ammonium_floating_respiration(self, rna) -> ArrayLike:
        if not self.use_floating_algae:
            return 0
        return self.floating_algae_process.ammonium_respiration()

    def change_nitrate(
        self,
        nitrate: ArrayLike,
        ammonium: ArrayLike,
        temperature: ArrayLike,
        depth: ArrayLike,
        dissolved_oxygen: ArrayLike,
    ) -> None:
        if not self.use_nitrate:
            return 0

        rate = (
            self.ammonium_nitrification(
                ammonium,
                temperature,
                # TODO
                # nitrification_inhibition,
            )
            - self.nitrate_denitrification(
                dissolved_oxygen,
                # TODO: need argument
                # half_saturation_oxygen,
                nitrate,
                temperature,
            )
            - self.nitrate_bed_denitrification(
                depth,
                nitrate,
                temperature,
            )
            - self.nitrate_uptake_floating_algae(
                nitrate,
                ammonium,
                # TODO: need argument
                # algea_growth_rate,
            )
            - self.nitrate_uptake_benthic_algae(
                nitrate,
                ammonium,
                # TODO: need argument
                # algea_growth_rate,
                depth,
            )
        )

        # Replace nan's with 0's
        rate = xr.where(rate == np.nan, 0, rate)
        return rate

    def ammonium_from_bed(self, depth: ArrayLike, temperature: ArrayLike) -> ArrayLike:
        rate = arrhenius_correction(
            temperature,
            self.sediment_ammonium_release_rate,
            self.sediment_ammonium_release_theta,
        )
        return rate / depth

    def ammonium_uptake_floating_algae(
        self, ammonium: ArrayLike, nitrate: ArrayLike
    ) -> ArrayLike:
        rate = np.nan
        if self.use_ammonium and not self.use_nitrate:
            rate = 1.0
        elif not self.use_ammonium and self.use_nitrate:
            rate = 0.0
        elif not self.use_ammonium and not self.use_nitrate:
            rate = 0.5
        elif self.use_ammonium and self.use_nitrate and self.use_floating_algae:
            rate = (
                self.floating_algae_preference_factor
                * ammonium
                / (
                    self.floating_algae_preference_factor * ammonium
                    + (1.0 - self.floating_algae_preference_factor) * nitrate
                )
            )

        # For cases where NH4 or NO3 are very small, force uptake fractions to ratio
        return xr.where(rate == np.nan, self.floating_algae_preference_factor, rate)

    def ammonium_rate_settling(self, depth: ArrayLike) -> ArrayLike:
        return depth * self.settling_velocity

    def ammonium_rate_death(self, rna: ArrayLike) -> ArrayLike:
        if not self.use_floating_algae:
            return 0.0
        return rna * self.death_rate

    def ammonium_nitrification(
        self,
        ammonium: ArrayLike,
        temperature: ArrayLike,
        nitrification_inhibition: ArrayLike,
    ) -> ArrayLike:
        if not self.use_ammonium:
            return 0.0

        # temperature adjust rate
        rate_corrected = arrhenius_correction(
            temperature, self.nitrification_rate, self.nitrification_theta
        )

        return ammonium * rate_corrected * nitrification_inhibition

    def nitrate_denitrification(
        self,
        dissolved_oxygen: ArrayLike,
        half_saturation_oxygen: ArrayLike,
        nitrate: ArrayLike,
        temperature: ArrayLike,
    ) -> ArrayLike:
        if not self.use_nitrate:
            return 0.0

        # temperature adjust rate
        rate_corrected = arrhenius_correction(
            temperature, self.denitrification_rate, self.denitrification_theta
        )

        rate = (
            nitrate
            * rate_corrected
            * (1.0 - (dissolved_oxygen / (dissolved_oxygen + half_saturation_oxygen)))
        )

        # replace nan's with 0's
        return xr.where(rate == np.nan, 0.0, rate)

    def nitrate_bed_denitrification(
        self,
        depth: ArrayLike,
        nitrate: ArrayLike,
        temperature: ArrayLike,
    ) -> ArrayLike:
        # temperature adjust rate
        rate_corrected = arrhenius_correction(
            temperature,
            self.sediment_denitrification_rate,
            self.sediment_denitrification_theta,
        )

        return nitrate * rate_corrected / depth

    def nitrate_uptake_floating_algae(
        self, nitrate: ArrayLike, ammonium: ArrayLike, algea_growth_rate: ArrayLike
    ) -> ArrayLike:
        if not self.use_floating_algae:
            return 0.0

        return (
            self.floating_algae_nitrogen_weight
            / self.algal_chlorophyll
            * algea_growth_rate
            * self.float_algea_faction_uptake_from_nitrate
        )

    def nitrate_uptake_benthic_algae(
        self,
        nitrate: ArrayLike,
        ammonium: ArrayLike,
        algea_growth_rate: ArrayLike,
        depth: ArrayLike,
    ) -> ArrayLike:
        if not self.use_benthic_algae:
            return 0.0

        return (
            self.benthic_algae_nitrogen_weight
            / self.algal_chlorophyll
            * algea_growth_rate
            * self.benthic_algea_faction_uptake_from_nitrate
            * self.fraction_bottom_area
        )
