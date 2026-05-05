from datetime import timedelta, datetime
import logging

import numpy as np
import xarray as xr

from .base import Process, ProcessFactory
from clearwater_data.variables import VariableRegistry
from clearwater_data.custom_types import ArrayLike
from clearwater_modules_v2.utils import conversions

arrhenius_correction = conversions.arrhenius_correction

logger = logging.getLogger(__name__)


# Defer v3 import to first use to break the v2 <-> v3 circular import
# chain that fires when v2.processes.benthic_algae triggers v3.__init__
# during v2.processes.__init__ enumeration. See floating_algae.py for the
# full discussion.

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from model import Model


class Nitrogen(Process):
    """v2 NSM1 Nitrogen Process (Phase 2.B fixes applied).

    Phase 1.3 (v3 NSM1 design spec, Section 3.4 + Appendix B): adopt the
    ``Process.DEFAULTS`` merge pattern for parameter handling. The v3
    ``NITROGEN_DEFAULTS`` dict is exposed as the class attribute
    ``DEFAULTS`` and merged with a user-supplied ``parameters`` dict at
    construction time; each merged key becomes a ``self.<name>`` attribute
    so subsequent NSM1 Process classes can copy the pattern.

    Phase 2.B (v3 NSM1 design spec, Section 11): the 16 known v2 Nitrogen
    bugs are now fixed (integrator additive Forward Euler in days, NaN
    guards via ``isnull()``, ``half_saturation_oxygen`` wired from
    ``self.KsOxdn``, algal growth/mortality coupling via the rate-variable
    cache on FloatingAlgae/BenthicAlgae, ``set_at_time`` persistence,
    clip-with-log diagnostics). OrgN is added as a third state variable
    (``change_organic_nitrogen``) per design spec Section 5.

    The legacy v2 keyword arguments (``nitrification_rate``,
    ``denitrification_rate``, etc.) are preserved for backward
    compatibility with existing tests and YAML configs. v2 kwargs and v3
    DEFAULTS keys have non-overlapping names, so both naming schemes
    coexist on the instance.
    """

    variables = [
        "nitrate",
        "ammonium",
        "organic_nitrogen",
        "oxygen_dissolved",
        "water_temperature",
        "depth",
    ]

    # Class-level v3 defaults (Section 3.4 of design spec). Lazy-loaded
    # on first instantiation; see the module-level note about the v2
    # <-> v3 circular import chain.
    DEFAULTS: dict[str, float | int | bool] = {}

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
        denitrification_rate: ArrayLike = 1.0,
        denitrification_theta: ArrayLike = 1.0,
        nitrification_rate: ArrayLike = 1.0,
        nitrification_theta: ArrayLike = 1.0,
        sediment_denitrification_rate: ArrayLike = 1.0,
        sediment_denitrification_theta: ArrayLike = 1.0,
        sediment_ammonium_release_rate: ArrayLike = 1.0,
        sediment_ammonium_release_theta: ArrayLike = 1.0,
        ammonium_decay_rate: ArrayLike = 1.0,
        ammonium_decay_theta: ArrayLike = 1.0,
        floating_algae_preference_factor: ArrayLike = 0.5,
        settling_velocity: ArrayLike = 1.0,
        death_rate: ArrayLike = 1.0,
        float_algea_faction_uptake_from_nitrate: ArrayLike = 1.0,
        nitrification_oxygen_inhibition_factor: ArrayLike = 1.0,
    ) -> None:
        Process.__init__(self, time_step)

        # --- Phase 1.3: v3-style parameter merge (DEFAULTS + user overrides) ---
        # Warn (don't error, don't silently ignore) on unknown keys so that
        # YAML typos surface but don't break instantiation.
        # Lazy-load NITROGEN_DEFAULTS to break the v2 <-> v3 circular
        # import chain.
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.nitrogen import DEFAULTS as NITROGEN_DEFAULTS
            type(self).DEFAULTS = NITROGEN_DEFAULTS

        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "Nitrogen: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in NITROGEN_DEFAULTS).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # --- Legacy v2 kwargs (preserved for backward compatibility) ---
        # These names do not collide with NITROGEN_DEFAULTS keys, so both
        # naming schemes coexist on the instance. Phase 2 will rewire the
        # kinetic methods to the v3 names and retire these.
        self.denitrification_rate = denitrification_rate
        self.denitrification_theta = denitrification_theta
        self.nitrification_rate = nitrification_rate
        self.nitrification_theta = nitrification_theta
        self.sediment_denitrification_rate = sediment_denitrification_rate
        self.sediment_denitrification_theta = sediment_denitrification_theta
        self.sediment_ammonium_release_rate = sediment_ammonium_release_rate
        self.sediment_ammonium_release_theta = sediment_ammonium_release_theta
        self.ammonium_decay_rate = ammonium_decay_rate
        self.ammonium_decay_theta = ammonium_decay_theta
        self.floating_algae_preference_factor = floating_algae_preference_factor
        self.settling_velocity = settling_velocity
        # Phase 2.B Bug #12: ``death_rate`` is preserved as a legacy
        # kwarg for back-compat with existing v2 unit tests, but the
        # canonical algal-death routing in v3 reads
        # ``floating_algae_process.algal_death_rate`` /
        # ``algal_orgn_from_mortality_rate`` (Phase 2.A populates these).
        self.death_rate = death_rate
        self.float_algea_faction_uptake_from_nitrate = (
            float_algea_faction_uptake_from_nitrate
        )
        self.nitrification_oxygen_inhibition_factor = (
            nitrification_oxygen_inhibition_factor
        )

        # Phase 2.B: diagnostics handle for clip-with-log. ``init_process``
        # will replace this with the model's run-level Diagnostics if a
        # v3 Model is wired up; tests that drive ``run`` directly without
        # a Model (e.g. Tier 1) keep this local instance.
        from clearwater_modules_v3.utils.numerics import Diagnostics
        self.diagnostics = Diagnostics()

        # Integration Item 1 (v3 NSM1, registry rate-variable convention,
        # spec resolved Q10 / Section 14): step-scoped flux caches read
        # by sibling Processes (N2 denit source, DOX nitrification O2
        # sink, eventually Alkalinity). Initialized to 0 so a downstream
        # Process that runs before Nitrogen.run gets a clean fallback
        # via getattr instead of an AttributeError.
        self.nitrification_flux_rate: ArrayLike = 0.0
        self.denitrification_flux_rate: ArrayLike = 0.0

    @ProcessFactory.register("nitrogen")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "Nitrogen":
        return Nitrogen(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        self.use_nitrate = True
        self.use_ammonium = True
        self.use_floating_algae = model.has_process("FloatingAlgae")
        self.use_benthic_algae = model.has_process("BenthicAlgae")

        if self.use_floating_algae:
            self.floating_algae_process = model.get_process("FloatingAlgae")
        if self.use_benthic_algae:
            self.benthic_algae_process = model.get_process("BenthicAlgae")

        # Phase 2.B: capture run-level Diagnostics from the v3 Model so
        # ``clip_negative_state`` records clip events on the canonical
        # diagnostics. v2's Model has no ``diagnostics`` attribute; in
        # that case the locally-instantiated Diagnostics from __init__
        # is retained.
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Run the Nitrogen process for one time step.

        Phase 2.B fixes:
        * Bug #1 / #2: replaced the multiplicative ``X = 0 + X * rate * dt``
          update with additive Forward Euler ``X_new = X + rate * dt_days``
          (rates are 1/d per the v1 NSM1 convention; ``dt_days`` converts
          the legacy seconds-based ``time_step`` to days). The ``dt``
          attribute name typo (``time_step_frequency``) is also fixed.
        * Bug #16: persist ``ammonium_new`` / ``nitrate_new`` /
          ``organic_nitrogen_new`` via ``registry.set_at_time``.
        * Q7 clip-with-log via ``clip_negative_state`` (with diagnostics).
        * OrgN: third state variable integrated alongside NH4 / NO3.

        Integration Item 1 (v3 NSM1, registry rate-variable convention,
        spec resolved Q10): cache the *step-scoped* nitrification and
        denitrification fluxes onto ``self`` for downstream consumers
        (N2, DOX, eventual Alkalinity). The legacy
        ``self.nitrification_rate`` / ``self.denitrification_rate``
        attribute names remain bound to the kinetic rate *constants*
        (1/d) for back-compat with v2 kwargs; the new fluxes use the
        ``_flux_rate`` suffix to disambiguate, with units mg-N/L/d and
        positive-valued absolute magnitudes.
        """
        from clearwater_modules_v3.utils.numerics import clip_negative_state

        # Pull state from registry.
        nitrate = registry.get_at_time("nitrate", time)
        ammonium = registry.get_at_time("ammonium", time)
        temperature = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)
        oxygen_dissolved = registry.get_at_time("oxygen_dissolved", time)

        # OrgN may not be present in legacy v2 registries; fall back to
        # zeros-like(ammonium) so the integrator and rate calculations
        # still execute (with ``use_OrgN`` gating the source/sink terms).
        if "organic_nitrogen" in registry:
            organic_nitrogen = registry.get_at_time("organic_nitrogen", time)
        else:
            organic_nitrogen = xr.zeros_like(ammonium) if hasattr(ammonium, "dims") else 0.0

        # 1/d -> per-step concentration delta.
        dt_days = self.time_step.total_seconds() / 86400.0

        # --- Step-scoped flux caches (Integration Item 1) ---
        # Compute the nitrification flux (NH4 -> NO3) and denitrification
        # flux (NO3 -> N2) here, before the change-rate decomposition.
        # These are positive-valued absolute magnitudes in mg-N/L/d that
        # downstream Processes (N2 source, DOX O2 sink) read via getattr.
        # NOTE: ``ammonium_nitrification`` and ``nitrate_denitrification``
        # already return non-negative fluxes by construction (kinetic
        # rates >= 0, NH4/NO3 >= 0, inhibition factors in [0, 1]); the
        # signs in ``change_ammonium`` / ``change_nitrate`` come from
        # how they are summed into the change-rates, not from the flux
        # values themselves.
        self.nitrification_flux_rate = self.ammonium_nitrification(
            ammonium,
            temperature,
            oxygen_dissolved,
        )
        self.denitrification_flux_rate = self.nitrate_denitrification(
            oxygen_dissolved,
            self.KsOxdn,
            nitrate,
            temperature,
        )

        # --- Ammonium update ---
        ammonium_rate = self.change_ammonium(
            nitrate,
            ammonium,
            temperature,
            depth,
            oxygen_dissolved,
            organic_nitrogen=organic_nitrogen,
        )
        ammonium_new = ammonium + ammonium_rate * dt_days
        ammonium_new = self._clip(ammonium_new, "ammonium")

        # --- Nitrate update ---
        nitrate_rate = self.change_nitrate(
            nitrate,
            ammonium,
            temperature,
            depth,
            oxygen_dissolved,
        )
        nitrate_new = nitrate + nitrate_rate * dt_days
        nitrate_new = self._clip(nitrate_new, "nitrate")

        # --- Organic Nitrogen update ---
        orgn_rate = self.change_organic_nitrogen(
            organic_nitrogen=organic_nitrogen,
            temperature=temperature,
            depth=depth,
        )
        organic_nitrogen_new = organic_nitrogen + orgn_rate * dt_days
        organic_nitrogen_new = self._clip(organic_nitrogen_new, "organic_nitrogen")

        # --- Persistence (Bug #16) ---
        registry.set_at_time("ammonium", time, ammonium_new)
        registry.set_at_time("nitrate", time, nitrate_new)
        if "organic_nitrogen" in registry:
            registry.set_at_time("organic_nitrogen", time, organic_nitrogen_new)

    def _clip(self, state, name: str):
        """Apply v3 clip-with-log if available; fall back to xr.where."""
        from clearwater_modules_v3.utils.numerics import clip_negative_state
        if isinstance(state, xr.DataArray) and self.diagnostics is not None:
            return clip_negative_state(state, name, self.diagnostics, step=0)
        return xr.where(state < 0, 0, state)

    def change_ammonium(
        self,
        nitrate: ArrayLike,
        ammonium: ArrayLike,
        temperature: ArrayLike,
        depth: ArrayLike,
        oxygen_dissolved: ArrayLike,
        organic_nitrogen: ArrayLike | None = None,
    ) -> None:
        if not self.use_ammonium:
            return 0

        # Phase 2.B: OrgN -> NH4 hydrolysis source (v1 OrgN_NH4_Decay).
        # Adds ``kon_tc * OrgN`` to the NH4 rate when ``use_OrgN`` is on.
        orgn_to_nh4 = self.organic_nitrogen_to_ammonium_hydrolysis(
            organic_nitrogen=organic_nitrogen,
            temperature=temperature,
        )

        rate = (
            self.ammonium_decay_nitrate(
                ammonium,
                temperature,
            )
            - self.ammonium_nitrification(
                ammonium,
                temperature,
                oxygen_dissolved,
            )
            + self.ammonium_from_bed(
                depth=depth,
                temperature=temperature,
            )
            + self.ammonium_floating_respiration()
            - self.ammonium_floating_growth()
            + self.ammonium_benthic_respiration()
            - self.ammonium_benthic_growth()
            + orgn_to_nh4
        )

        # Phase 2.B Bugs #5-#8: replace ``rate == np.nan`` (always False)
        # with a real null check.
        if isinstance(rate, xr.DataArray):
            rate = xr.where(rate.isnull(), 0, rate)
        else:
            rate = np.where(np.isnan(rate), 0, rate)
        return rate

    def ammonium_floating_respiration(self) -> ArrayLike:
        if not self.use_floating_algae:
            return 0
        return self.floating_algae_process.ammonium_respiration()

    def ammonium_benthic_respiration(self) -> ArrayLike:
        if not self.use_benthic_algae:
            return 0
        return self.benthic_algae_process.ammonium_respiration()

    def ammonium_floating_growth(self) -> ArrayLike:
        if not self.use_floating_algae:
            return 0
        return self.floating_algae_process.ammonium_growth()

    def ammonium_benthic_growth(self) -> ArrayLike:
        if not self.use_benthic_algae:
            return 0
        return self.benthic_algae_process.ammonium_growth()

    def change_nitrate(
        self,
        nitrate: ArrayLike,
        ammonium: ArrayLike,
        temperature: ArrayLike,
        depth: ArrayLike,
        oxygen_dissolved: ArrayLike,
    ) -> None:
        if not self.use_nitrate:
            return 0

        # Phase 2.B Bugs #10, #11: read algal growth rates from the
        # FloatingAlgae / BenthicAlgae step-scoped rate cache (Phase 2.A
        # populated these in their ``run`` methods). Falls back to 0 when
        # the algae module is disabled or when ``run`` has not yet been
        # called this step.
        if self.use_floating_algae:
            float_algae_growth = getattr(
                self.floating_algae_process, "algal_growth_rate", 0
            )
        else:
            float_algae_growth = 0
        if self.use_benthic_algae:
            benthic_algae_growth = getattr(
                self.benthic_algae_process, "balgae_growth_rate", 0
            )
        else:
            benthic_algae_growth = 0

        rate = (
            self.ammonium_nitrification(
                ammonium,
                temperature,
                oxygen_dissolved,
            )
            - self.nitrate_denitrification(
                oxygen_dissolved,
                # Phase 2.B Bug #9: wire half-saturation O2 from the v3
                # ``KsOxdn`` parameter (NITROGEN_DEFAULTS).
                self.KsOxdn,
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
                float_algae_growth,
            )
            - self.nitrate_uptake_benthic_algae(
                nitrate,
                ammonium,
                benthic_algae_growth,
                depth,
            )
        )

        # Phase 2.B Bugs #5-#8: real null check (was ``rate == np.nan``).
        if isinstance(rate, xr.DataArray):
            rate = xr.where(rate.isnull(), 0, rate)
        else:
            rate = np.where(np.isnan(rate), 0, rate)
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

        # For cases where NH4 or NO3 are very small, force uptake fractions to ratio.
        # Phase 2.B Bugs #5-#8: real null check (was ``rate == np.nan``).
        if isinstance(rate, xr.DataArray):
            return xr.where(rate.isnull(), self.floating_algae_preference_factor, rate)
        return np.where(np.isnan(rate), self.floating_algae_preference_factor, rate)

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
        oxygen_dissolved: ArrayLike,
    ) -> ArrayLike:
        if not self.use_ammonium:
            return 0.0

        # temperature adjust rate
        rate_corrected = arrhenius_correction(
            temperature, self.nitrification_rate, self.nitrification_theta
        )

        return (
            ammonium * rate_corrected * self.nitrification_inhibition(oxygen_dissolved)
        )

    def ammonium_decay_nitrate(
        self, ammonium: ArrayLike, temperature: ArrayLike
    ) -> ArrayLike:
        if not self.use_ammonium:
            return 0.0

        # temperature adjust rate
        rate_corrected = arrhenius_correction(
            temperature, self.ammonium_decay_rate, self.ammonium_decay_theta
        )

        return ammonium * rate_corrected

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

        # Phase 2.B Bugs #5-#8: real null check (was ``rate == np.nan``).
        if isinstance(rate, xr.DataArray):
            return xr.where(rate.isnull(), 0.0, rate)
        return np.where(np.isnan(rate), 0.0, rate)

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

    def nitrification_inhibition(self, oxygen_dissolved: ArrayLike) -> ArrayLike:
        if not self.use_nitrate:
            return 1.0

        return 1.0 - np.exp(
            -self.nitrification_oxygen_inhibition_factor * oxygen_dissolved
        )

    # ------------------------------------------------------------------
    # Organic Nitrogen (Phase 2.B; v3 NSM1 design spec Section 5)
    # ------------------------------------------------------------------
    # OrgN is a third state variable on the Nitrogen Process. Sources:
    # algal mortality routing from FloatingAlgae/BenthicAlgae (Phase 2.A's
    # ``algal_orgn_from_mortality_rate`` / ``balgae_orgn_from_mortality_rate``).
    # Sinks: hydrolysis to NH4 (``kon_tc * OrgN``) and settling
    # (``vson_tc / depth * OrgN``). v1 source: ``processes.py`` 1173-1420.

    def organic_nitrogen_to_ammonium_hydrolysis(
        self,
        organic_nitrogen: ArrayLike | None,
        temperature: ArrayLike,
    ) -> ArrayLike:
        """OrgN -> NH4 hydrolysis flux (mg-N/L/d). v1 ``OrgN_NH4_Decay``.

        Returns ``kon_tc * OrgN`` when ``use_OrgN`` is on; 0 otherwise.
        """
        if not getattr(self, "use_OrgN", True):
            return 0.0
        if organic_nitrogen is None:
            return 0.0
        kon_tc = arrhenius_correction(temperature, self.kon_20, self.kon_theta)
        return kon_tc * organic_nitrogen

    def organic_nitrogen_settling(
        self,
        organic_nitrogen: ArrayLike,
        temperature: ArrayLike,
        depth: ArrayLike,
    ) -> ArrayLike:
        """OrgN -> bed settling flux (mg-N/L/d). v1 ``OrgN_Settling``.

        ``vson_tc / depth * OrgN``. ``vson_20`` carries m/d units; the
        Arrhenius temperature correction is applied for parity with v1.
        """
        if not getattr(self, "use_OrgN", True):
            return 0.0
        vson_tc = arrhenius_correction(temperature, self.vson_20, self.vson_theta)
        return vson_tc / depth * organic_nitrogen

    def organic_nitrogen_from_floating_algae_mortality(self) -> ArrayLike:
        """Algal-mortality OrgN source (mg-N/L/d). v1 ``ApDeath_OrgN``.

        Routes through Phase 2.A's ``algal_orgn_from_mortality_rate`` cache
        on the FloatingAlgae process.
        """
        if not self.use_floating_algae:
            return 0.0
        return getattr(
            self.floating_algae_process, "algal_orgn_from_mortality_rate", 0.0
        )

    def organic_nitrogen_from_benthic_algae_mortality(self) -> ArrayLike:
        """Benthic-algal-mortality OrgN source (mg-N/L/d). v1 ``AbDeath_OrgN``.

        Routes through Phase 2.A's ``balgae_orgn_from_mortality_rate`` cache
        on the BenthicAlgae process.
        """
        if not self.use_benthic_algae:
            return 0.0
        return getattr(
            self.benthic_algae_process, "balgae_orgn_from_mortality_rate", 0.0
        )

    def change_organic_nitrogen(
        self,
        organic_nitrogen: ArrayLike,
        temperature: ArrayLike,
        depth: ArrayLike,
    ) -> ArrayLike:
        """Net OrgN rate of change (mg-N/L/d). v1 ``dOrgNdt`` (line 1383).

        Sources: algal mortality (floating + benthic).
        Sinks: hydrolysis to NH4, settling to bed.
        Gated by ``use_OrgN`` (default True).
        """
        if not getattr(self, "use_OrgN", True):
            return 0.0

        ap_death_orgn = self.organic_nitrogen_from_floating_algae_mortality()
        ab_death_orgn = self.organic_nitrogen_from_benthic_algae_mortality()
        orgn_to_nh4 = self.organic_nitrogen_to_ammonium_hydrolysis(
            organic_nitrogen=organic_nitrogen,
            temperature=temperature,
        )
        orgn_settling = self.organic_nitrogen_settling(
            organic_nitrogen=organic_nitrogen,
            temperature=temperature,
            depth=depth,
        )

        rate = ap_death_orgn + ab_death_orgn - orgn_to_nh4 - orgn_settling

        # Phase 2.B Bugs #5-#8: real null check.
        if isinstance(rate, xr.DataArray):
            rate = xr.where(rate.isnull(), 0, rate)
        elif isinstance(rate, np.ndarray):
            rate = np.where(np.isnan(rate), 0, rate)
        return rate
