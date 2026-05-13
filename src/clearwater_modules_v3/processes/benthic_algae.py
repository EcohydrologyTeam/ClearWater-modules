"""v2 NSM1 BenthicAlgae Process.

Phase 2.A (v3 NSM1 design spec, Section 11 + Section 6 bug list): apply
the parallel set of bug fixes to ``BenthicAlgae`` that Phase 2.A applies
to ``FloatingAlgae``:

- Integrator fix: replace the broken multiplicative integrator inherited
  from ``FloatingAlgae.run`` with additive Forward Euler.
- Implement ``ammonium_respiration`` and ``ammonium_growth`` per v1
  benthic algae kinetics (NH4_AbRespiration / NH4_AbGrowth) using the
  benthic stoichiometry (BWn) and the depth/Fb area-integration factor.
- NaN guards: replace broken ``rate == np.nan`` with ``isnull`` /
  ``np.isnan``.
- Persistence: persist updated state via ``registry.set_at_time``.

Adopt the v3 ``BALGAE_DEFAULTS`` merge pattern. Algal-mortality routing
methods are added that mirror FloatingAlgae's, using benthic
stoichiometry (BWn, BWp, BWc) and the v1 ``Fw`` (fraction to water
column) and ``Fb`` (fraction of bottom area) factors.
"""

from datetime import datetime, timedelta
import logging

import numpy as np
import xarray as xr

from clearwater_data.variables import VariableRegistry
from clearwater_data.custom_types import ArrayLike

from .floating_algae import FloatingAlgae

from clearwater_modules_v3.utils.conversions import arrhenius_correction
from clearwater_modules_v3.utils.numerics import clip_negative_state


def _sanitize_cache(value):
    """Replace NaN cells with 0 in a downstream-consumer cache.

    Cached rate variables (``balgae_growth_rate``, ``balgae_respiration_rate``)
    are read by sibling processes (DOX, Nitrogen, Phosphorus) that sum many
    sub-fluxes into a per-cell rate. A NaN in one sub-flux poisons the entire
    cell's rate via the final ``where(isnull, 0, rate)`` sanitization, which
    zeroes the cell rather than just the bad term — freezing the cell's
    state at IC indefinitely.

    The semantically correct value when the rate computation produces NaN
    (typically at cells whose inputs were NaN at an earlier dry step, or
    where ``algae == 0`` makes the rate mathematically 0) is 0, not NaN.
    Sanitizing at the cache source ensures all downstream consumers receive
    a finite contribution.
    """
    if isinstance(value, xr.DataArray):
        return xr.where(value.isnull(), 0, value)
    if isinstance(value, np.ndarray):
        return np.where(np.isnan(value), 0, value)
    if isinstance(value, (int, float)):
        if value != value:  # NaN
            return 0
    return value


# Defer v3 imports to first use; see floating_algae.py for the full
# discussion of the v2 <-> v3 circular-import chain.

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


logger = logging.getLogger(__name__)


# Inline fallback defaults for partitioning inputs that live in the
# v3 phosphorus / global_vars / carbon modules. Used only when the user
# does not pass them in the merged parameters dict.
_BENTHIC_FDP_DEFAULTS = {
    "kdpo4": 0.0,
    "Solid": 1.0,
    "use_TIP": True,
    # mortality routing fractions (v1 carbon: f_pocb)
    "f_pocb": 0.5,  # fraction of benthic algal mortality C that goes to POC
}


class BenthicAlgae(FloatingAlgae):
    variables = ["benthic_algae", "solar_radiation"]

    # Class-level v3 defaults for benthic algae parameters. Lazy-loaded
    # on first instantiation to avoid the v2 <-> v3 circular import.
    DEFAULTS: dict[str, float | int | bool] = {}

    # Phase 9.A.1 wiring fix: BenthicAlgae overrides the legacy-kwarg
    # to v3-DEFAULTS mapping inherited from FloatingAlgae so that
    # benthic-specific defaults (e.g. ``KsNb`` instead of ``KsN``,
    # ``mub_max_20`` instead of ``mu_max_20``) flow into the legacy
    # attribute names that the rate methods read.
    _LEGACY_TO_DEFAULTS: dict[str, str] = {
        "growth_rate_max": "mub_max_20",
        "growth_rate_correction": "mub_max_theta",
        "death_rate": "kdb_20",
        "death_rate_correction_factor": "kdb_theta",
        "repiration_rate": "krb_20",
        "repiration_rate_correction_factor": "krb_theta",
        "light_limitation_constant": "KLb",
        "nitrogen_michaelis_menton_constant": "KsNb",
        "phosphorus_michaelis_menton_constant": "KsPb",
        "density_michaelis_menton_constant": "Ksb",
    }

    def __init__(
        self,
        parameters: dict | None = None,
        *args,
        density_michaelis_menton_constant: float | None = None,
        **kwargs,
    ) -> None:
        """Initialize the benthic algae process.

        Args:
            parameters: Optional dict of v3 balgae parameter overrides.
                Merged with the class-level ``DEFAULTS`` (v3
                ``BALGAE_DEFAULTS``). Unknown keys are warned and ignored.
            density_michaelis_menton_constant: Michaelis-Menton constant
                for benthic-density limitation (v1 Ksb). If None, falls
                back to ``Ksb`` (10.0 g-D/m^2) from BALGAE_DEFAULTS.

        Phase 9.A.1 wiring fix: see FloatingAlgae.__init__ for the
        contract. BenthicAlgae overrides ``_LEGACY_TO_DEFAULTS`` so
        the rate methods (most of which are inherited) read benthic-
        specific values (``mub_max_20``, ``kdb_20``, ``KsNb``, etc.)
        rather than the floating-algae shadows.
        """
        # --- Phase 2.A: v3-style parameter merge for BALGAE_DEFAULTS ---
        # Done here rather than relying on FloatingAlgae.__init__ because
        # BenthicAlgae has its own DEFAULTS (BALGAE_DEFAULTS) that
        # override ALGAE_DEFAULTS for benthic-specific parameters.
        # Lazy-load BALGAE_DEFAULTS to break the v2 <-> v3 circular
        # import. Phase 3.5 inter-process coupling: compose ``h2``
        # (active sediment layer thickness) from POM_DEFAULTS so
        # ``balgae_pom_from_mortality_rate`` can be cached without
        # requiring ``h2`` to be passed explicitly. Mirrors the
        # Phase 5.B DOX multi-group composition pattern.
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.balgae import DEFAULTS as BALGAE_DEFAULTS
            from clearwater_modules_v3.parameters.pom import DEFAULTS as POM_DEFAULTS
            composed: dict[str, float | int | bool] = {}
            composed.update(BALGAE_DEFAULTS)
            # h2: active sediment layer thickness (m). Lives in POM_DEFAULTS
            # but is needed by the POM<-BenthicAlgae mortality coupling.
            composed["h2"] = POM_DEFAULTS["h2"]
            # Selector-name remap: BALGAE_DEFAULTS uses
            # ``b_growth_rate_option`` and ``b_light_limitation_option``;
            # the FloatingAlgae rate methods read ``growth_rate_option``
            # and ``light_limitation_option``. Mirror the benthic
            # selectors onto the floating-algae names so the inherited
            # methods read the benthic-specific values.
            if "b_growth_rate_option" in composed:
                composed["growth_rate_option"] = composed["b_growth_rate_option"]
            if "b_light_limitation_option" in composed:
                composed["light_limitation_option"] = composed[
                    "b_light_limitation_option"
                ]
            type(self).DEFAULTS = composed

        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "BenthicAlgae: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in BALGAE_DEFAULTS).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Inline partitioning / routing defaults from non-balgae groups.
        for k, v in _BENTHIC_FDP_DEFAULTS.items():
            if not hasattr(self, k):
                setattr(self, k, user_params.get(k, v))

        # FloatingAlgae.__init__ does NOT take a `parameters` kwarg here;
        # we already applied the BALGAE_DEFAULTS merge above. Pass
        # parameters=None to FloatingAlgae to skip a second merge against
        # ALGAE_DEFAULTS (which would overwrite our BALGAE values for any
        # non-overlapping keys with algae-pelagic defaults). Legacy v2
        # kwargs pass through unchanged. The Phase 9.A.1 wiring fix
        # (resolving the legacy-kwarg shadow with self._LEGACY_TO_DEFAULTS)
        # runs inside FA.__init__ and uses BenthicAlgae's overridden
        # mapping (``mub_max_20``, ``kdb_20``, ``KsNb``, etc.) because
        # ``self._LEGACY_TO_DEFAULTS`` resolves to the leaf class.
        # ``density_michaelis_menton_constant`` is handled by FA.__init__
        # via the inherited legacy-kwarg path; pass it through kwargs.
        if density_michaelis_menton_constant is not None:
            kwargs["density_michaelis_menton_constant"] = (
                density_michaelis_menton_constant
            )
        FloatingAlgae.__init__(self, parameters=None, *args, **kwargs)

        # Step-scoped rate-variable cache for benthic algae, keyed for the
        # benthic-specific consumer side (v1 NH4_AbRespiration etc.).
        self.balgae_growth_rate: ArrayLike = 0.0
        self.balgae_respiration_rate: ArrayLike = 0.0
        self.balgae_death_rate: ArrayLike = 0.0
        self.balgae_nh4_uptake_fraction: ArrayLike = 0.5
        self.balgae_orgn_from_mortality_rate: ArrayLike = 0.0
        self.balgae_orgp_from_mortality_rate: ArrayLike = 0.0
        self.balgae_poc_from_mortality_rate: ArrayLike = 0.0
        self.balgae_doc_from_mortality_rate: ArrayLike = 0.0
        # Phase 3.5 inter-process coupling: POM consumer reads this
        # cache. v1 ``POM_benthic_algae_mortality = Ab * kdb_tc * Fb *
        # (1 - Fw) / h2`` (mg/L/d). Note: this uses ``(1 - Fw)`` (the
        # fraction of mortality NOT released to the water column, which
        # ends up as POM in the active sediment layer); the OrgN / OrgP
        # / POC / DOC routings above use ``Fw`` (the fraction released
        # to the water column) instead.
        self.balgae_pom_from_mortality_rate: ArrayLike = 0.0

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        # Phase 1.D: replace the hardcoded ``use_* = True`` with sibling
        # discovery via ``model.has_process``, matching the FloatingAlgae
        # pattern. Defaults remain ``True`` when ``model`` does not
        # expose ``has_process`` (e.g., legacy test harnesses) for back-compat.
        has_proc = getattr(model, "has_process", None)
        if has_proc is None:
            # Legacy / fixture-only path: keep prior behaviour so isolated
            # unit tests that exercise BenthicAlgae without a real Model
            # continue to work.
            self.use_nitrate = True
            self.use_ammonium = True
            self.use_phosphate = True
        else:
            # When wired through a Model, the per-state-variable flags
            # follow the sibling-Process presence. BenthicAlgae's
            # nutrient-limitation kinetics need NH4 / NO3 / TIP whether
            # or not the matching Process is present, but the *uptake*
            # path (which writes to NH4/NO3/TIP via the rate cache) is
            # only meaningful when a sibling Process consumes those
            # values. Match FloatingAlgae's convention exactly:
            self.use_nitrate = model.has_process("Nitrogen")
            self.use_ammonium = model.has_process("Nitrogen")
            self.use_phosphate = model.has_process("Phosphorus")

        # Capture run-level Diagnostics (Phase 2.0 prerequisite).
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Run the benthic algae process.

        Forward-Euler integrator with v3 clip-with-log diagnostics and
        registry persistence. Mirrors ``FloatingAlgae.run`` but uses the
        benthic algae registry variable name and the additional density
        limitation factor.
        """

        algae = registry.get_at_time("benthic_algae", time)
        ammonium = registry.get_at_time("ammonium", time)
        nitrate = registry.get_at_time("nitrate", time)
        # Phase 8.A: read the canonical v3 inorganic-P state name "tip".
        # The legacy v2 name "phosphorus_total_inorganic" (still used by
        # v2 Riverine and a handful of v2-era YAML configs) is supported
        # as a fallback so older configs continue to run unchanged.
        if "tip" in registry:
            phosphorus_total_inorganic = registry.get_at_time("tip", time)
        else:
            phosphorus_total_inorganic = registry.get_at_time(
                "phosphorus_total_inorganic", time
            )
        depth = registry.get_at_time("depth", time)
        water_temperature = registry.get_at_time("water_temperature", time)
        solar = registry.get_at_time("solar_radiation", time)

        # Use v3 fdp utility for the dissolved P fraction.
        from clearwater_modules_v3.utils.partitioning import fdp as fdp_partition

        phosphate_fraction_dissolved = fdp_partition(
            use_TIP=self.use_TIP,
            Solid=self.Solid,
            kdpo4=self.kdpo4,
        )

        rate = self.rate(
            algae=algae,
            depth=depth,
            water_temperature=water_temperature,
            phosphorus_total_inorganic=phosphorus_total_inorganic,
            phosphate_fraction_dissolved=phosphate_fraction_dissolved,
            ammonium=ammonium,
            nitrate=nitrate,
            solar=solar,
        )

        # Cache depth for v1 NH4_AbRespiration / NH4_AbGrowth (which
        # divide by depth to convert g/m^2/d areal rates into mg-N/L/d
        # volumetric rates).
        self._cached_depth = depth

        # Cache step-scoped routing variables for downstream Processes.
        self._cache_benthic_mortality_rates(algae, water_temperature, depth)

        # Cache the benthic NH4-uptake fraction.
        self.balgae_nh4_uptake_fraction = self._ab_uptake_fr_nh4(
            ammonium=ammonium, nitrate=nitrate
        )

        # Forward Euler in days.
        dt_days = self.time_step.total_seconds() / 86400.0
        algae_new = algae + rate * dt_days

        # Clip-with-log per Q7. Import is module-level (Phase 1.C);
        # step attribution is automatic via ``diagnostics.current_step``
        # (Phase 0.6 Q1).
        algae_new = clip_negative_state(
            algae_new, "benthic_algae", self.diagnostics
        )

        # Persistence (Bug #16 parallel for benthic algae).
        registry.set_at_time("benthic_algae", time, algae_new)

    def _cache_benthic_mortality_rates(
        self,
        algae: ArrayLike,
        water_temperature: ArrayLike,
        depth: ArrayLike,
    ) -> None:
        """Compute and cache benthic-algae mortality routing rates.

        Per v1 (AbDeath_OrgN, AbDeath_OrgP, POC_benthic_algae_mortality,
        DOC_benthic_algae_mortality, POM_benthic_algae_mortality):

        - AbDeath_OrgN = rnb * Fw * Fb * AbDeath / depth
        - AbDeath_OrgP = rpb * Fw * Fb * AbDeath / depth
        - POC_balgae_mortality = (1 / depth) * f_pocb * Fb * Fw * rcb * AbDeath
        - DOC_balgae_mortality = (1 / depth) * (1 - f_pocb) * Fb * Fw * rcb * AbDeath
        - POM_balgae_mortality = Ab * kdb_tc * Fb * (1 - Fw) / h2
                               = AbDeath * Fb * (1 - Fw) / h2

        v1 rnb / rpb / rcb are BWn/BWd, BWp/BWd, BWc/BWd. Here we use
        the v3 BALGAE_DEFAULTS keys directly.

        Note on the POM term: the OrgN / OrgP / POC / DOC routings use
        ``Fw`` (fraction released to water column) divided by ``depth``
        (water column depth). The POM routing uses ``(1 - Fw)`` (the
        complementary fraction NOT released to water; it ends up as
        POM in the active sediment layer) divided by ``h2`` (active
        sediment layer thickness, m). Per the Phase 4 Phosphorus agent
        report, the mortality-routing caches return consumer-ready
        fluxes: the consumer reads them directly without re-multiplying
        by ``Fw``, ``Fb``, ``depth``, or ``h2``.
        """
        ab_death = self.rate_death(algae, water_temperature)

        rnb = self.BWn / self.BWd  # mg-N/mg-D
        rpb = self.BWp / self.BWd  # mg-P/mg-D
        rcb = self.BWc / self.BWd  # mg-C/mg-D

        fw = self.Fw
        fb = self.Fb

        self.balgae_death_rate = ab_death
        self.balgae_orgn_from_mortality_rate = rnb * fw * fb * ab_death / depth
        self.balgae_orgp_from_mortality_rate = rpb * fw * fb * ab_death / depth
        self.balgae_poc_from_mortality_rate = (
            self.f_pocb * fb * fw * rcb * ab_death / depth
        )
        self.balgae_doc_from_mortality_rate = (
            (1.0 - self.f_pocb) * fb * fw * rcb * ab_death / depth
        )
        # POM source from benthic algae mortality (mg/L/d). Uses (1-Fw)
        # and h2 (sediment layer thickness), NOT Fw and depth.
        self.balgae_pom_from_mortality_rate = (
            ab_death * fb * (1.0 - fw) / self.h2
        )

    def _ab_uptake_fr_nh4(
        self, ammonium: ArrayLike, nitrate: ArrayLike
    ) -> ArrayLike:
        """v1 AbUptakeFr_NH4 (lines 1263-1302). Same shape as
        ``ApUptakeFr_NH4`` but uses the benthic preference factor PNb."""
        pnb = getattr(self, "PNb", 0.5)

        if self.use_ammonium and not self.use_nitrate:
            return 1.0
        if not self.use_ammonium and self.use_nitrate:
            return 0.0
        if not self.use_ammonium and not self.use_nitrate:
            return 0.5

        denom = pnb * ammonium + (1.0 - pnb) * nitrate
        result = xr.where(denom > 0, pnb * ammonium / denom, pnb)
        if isinstance(result, xr.DataArray):
            result = xr.where(result.isnull(), pnb, result)
        return result

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
        """Compute the rate of change of benthic algae."""
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
        """Compute the rate of growth of benthic algae."""

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

        result = rate * algae
        # Cache for downstream consumers (Q10 GS-rates pattern).
        # Sanitize NaN at the cache source: a per-cell NaN here propagates
        # via DOX / Nitrogen / Phosphorus rate sums and zeroes the *entire*
        # cell's rate after their final ``sanitize_rate`` step, freezing the
        # state at IC indefinitely. When ``algae == 0`` the contribution is
        # mathematically 0 regardless of upstream computation, so NaN -> 0
        # is the semantically correct value at dry / never-wet cells.
        result = _sanitize_cache(result)
        self.balgae_growth_rate = result
        return result

    def rate_respiration(
        self, algae: ArrayLike, water_temperature: ArrayLike
    ) -> ArrayLike:
        """Compute the rate of respiration of benthic algae.

        Mirrors v1 AbRespiration but uses ``corrected_respiration_rate *
        algae``. We override FloatingAlgae's so we can cache the benthic
        rate variable separately.
        """
        corrected_respiration_rate = arrhenius_correction(
            water_temperature,
            self.repiration_rate,
            self.repiration_rate_correction_factor,
        )
        result = algae * corrected_respiration_rate
        # See ``rate_growth`` — sanitize NaN at the cache source.
        result = _sanitize_cache(result)
        self.balgae_respiration_rate = result
        return result

    def limit_light(
        self,
        algae: ArrayLike,
        depth: ArrayLike,
        surface_light_intensity: ArrayLike,
    ) -> ArrayLike:
        """Compute the limiting light for benthic algae."""

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
        # Phase 9.A.1 audit B6: with x = PAR*KEXT/KLb, the Steele form
        # is ``x * exp(1 - x)``. The previous form used division
        # (``x / exp(1-x)`` = ``x * exp(x-1)``), which is the wrong
        # sign for the exponent argument. v1/Fortran use ``x * exp(1-x)``.
        elif self.light_limitation_option == 3:
            raw_rate = xr.where(
                abs(self.light_limitation_constant) < 1e-10,
                0,
                (
                    surface_light_intensity
                    * light_at_depth_coefficent
                    / self.light_limitation_constant
                    * np.exp(
                        1.0
                        - surface_light_intensity
                        * light_at_depth_coefficent
                        / self.light_limitation_constant
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
        """Compute the limiting density for benthic algae.

        Matches v1 ``FSb(Ab, Ksb) = 1 - Ab / (Ab + Ksb)``.
        """
        limit_raw = 1.0 - (algae / (algae + self.density_michaelis_menton_constant))
        limit = xr.where(limit_raw > 1, 1, limit_raw)
        # Bug fix: x == np.nan is always False per IEEE 754. Use isnull.
        if isinstance(limit, xr.DataArray):
            limit = xr.where(limit.isnull(), 0, limit)
        else:
            limit = xr.where(np.isnan(limit), 0, limit)

        return limit

    # ------------------------------------------------------------------
    # NH4 coupling for benthic algae (v1 NH4_AbRespiration / NH4_AbGrowth)
    # ------------------------------------------------------------------

    def ammonium_respiration(self) -> ArrayLike:
        """v1 NH4_AbRespiration (line 1525): (rnb * AbRespiration * Fb) / depth.

        Returns mg-N/L/d transferred from benthic algal respiration to
        NH4. Uses cached ``balgae_respiration_rate``. We use a stored
        depth value cached during ``run``; if ``run`` has not been
        called yet, returns 0.
        """
        rnb = self.BWn / self.BWd  # mg-N/mg-D
        fb = self.Fb
        depth = getattr(self, "_cached_depth", None)
        if depth is None:
            return 0.0
        return rnb * self.balgae_respiration_rate * fb / depth

    def ammonium_growth(self) -> ArrayLike:
        """v1 NH4_AbGrowth (line 1547): (AbUptakeFr_NH4 * rnb * Fb * AbGrowth) / depth.

        Returns mg-N/L/d removed from NH4 by benthic algal growth.
        Uses the cached ``balgae_growth_rate`` and
        ``balgae_nh4_uptake_fraction`` written by ``run``.
        """
        rnb = self.BWn / self.BWd  # mg-N/mg-D
        fb = self.Fb
        depth = getattr(self, "_cached_depth", None)
        if depth is None:
            return 0.0
        return self.balgae_nh4_uptake_fraction * rnb * fb * self.balgae_growth_rate / depth

    # ------------------------------------------------------------------
    # Algal mortality routing helpers (Q10 GS-rates contract)
    # ------------------------------------------------------------------

    def death_to_orgn(self) -> ArrayLike:
        """v1 AbDeath_OrgN: rnb * Fw * Fb * AbDeath / depth (mg-N/L/d)."""
        return self.balgae_orgn_from_mortality_rate

    def death_to_orgp(self) -> ArrayLike:
        """v1 AbDeath_OrgP: rpb * Fw * Fb * AbDeath / depth (mg-P/L/d)."""
        return self.balgae_orgp_from_mortality_rate

    def death_to_poc(self) -> ArrayLike:
        """v1 POC_benthic_algae_mortality: (1/depth) * f_pocb * Fb * Fw * rcb * AbDeath."""
        return self.balgae_poc_from_mortality_rate

    def death_to_doc(self) -> ArrayLike:
        """v1 DOC_benthic_algae_mortality: (1/depth) * (1-f_pocb) * Fb * Fw * rcb * AbDeath."""
        return self.balgae_doc_from_mortality_rate
