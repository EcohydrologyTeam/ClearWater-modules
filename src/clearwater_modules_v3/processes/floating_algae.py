"""v3 NSM1 FloatingAlgae Process (v3-native; not a v2 re-export).

As-implemented record (Phase 2.A, v3 NSM1 design spec Section 11 +
Section 6 bug list). The following v2-inherited bugs are **fixed and
regression-tested**, and the v3 patterns are adopted:

- Bug #4 (fixed): the broken multiplicative integrator
  ``algae * rate * dt * 86400`` was replaced with additive Forward
  Euler ``algae + rate * dt_days``.
- Bug #13 (fixed): ``ammonium_respiration`` is implemented (was
  returning 0 in v2).
- Bug #14 (fixed): ``ammonium_growth`` is implemented (was returning 0
  in v2).
- Bug #15 (fixed): the hard-coded ``phosphate_fraction_dissolved=0.5``
  was replaced with the v3 ``fdp`` partitioning utility.
- Bug #16: re-add ``set_at_time`` persistence after the integrator step
  (the prior fix was lost when streaming branch was rebased onto
  ``upstream/memory-refactor-pytestUpdate`` via the C8 sync).
- NaN guards: replace broken ``rate == np.nan`` comparisons with
  ``rate.isnull()`` / ``np.isnan`` (IEEE 754 makes ``x == np.nan``
  always False, so the prior guards were no-ops).

Irradiance basis (NSM1-SCI-A3, gold-standard spec B1): the
``solar_radiation`` registry variable is **total broadband shortwave**
(W/m^2). Algal photosynthesis responds to photosynthetically active
radiation only; ``KL`` (``light_limitation_constant``) is a PAR-scale
half-saturation value. ``run`` therefore converts at the process
boundary, ``PAR = solar_radiation * Fr_PAR`` (``Fr_PAR=0.47``), before
``limit_light`` -- mirroring NSM1 v1 ``processes.py:287``. The registry
variable itself is not mutated.

Adopt the Phase 1.3 DEFAULTS-merge pattern established by ``Nitrogen``:
the v3 ``ALGAE_DEFAULTS`` is the class ``DEFAULTS`` and is merged with
a user ``parameters`` dict at construction time. Legacy v2 kwargs are
preserved for backward compatibility.

Also adds algal-mortality routing methods that compute and stash per-cell
mortality rates as instance attributes so downstream Processes
(``Nitrogen``, ``Phosphorus``, ``Carbon``) can read them after
``run`` completes:

- ``algal_orgn_from_mortality_rate``  (mg-N/L/d)
- ``algal_orgp_from_mortality_rate``  (mg-P/L/d)
- ``algal_poc_from_mortality_rate``   (mg-C/L/d)
- ``algal_doc_from_mortality_rate``   (mg-C/L/d)
- ``algal_pom_from_settling_rate``    (mg/L/d) -- consumed by ``POM``;
  v1 ``POM_algal_settling = vsap * Ap * (AWd/AWa) / h2``. Note: this is
  a *settling* flux (Bug-#-fix-symmetry with the mortality flux pattern),
  not a mortality flux. ``h2`` (active sediment layer thickness) is
  composed onto the FloatingAlgae instance from ``parameters.pom`` via
  the DEFAULTS-composition pattern (Phase 5.B DOX precedent).

Per the resolved Q10 GS-rates contract, these are step-scoped (NOT
time-indexed). The full registry-side rate-variable plumbing
(``Registry.set_rate_variable`` / ``clear_rate_variables``) is a Phase
2.A.1 follow-up; for now downstream Nitrogen.run reads them via
``floating_algae_process.<rate_name>`` directly.
"""

from datetime import datetime, timedelta
import logging

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.base import Process, ProcessFactory
from clearwater_data.variables import VariableRegistry
from clearwater_data.custom_types import ArrayLike

from clearwater_modules_v3.utils.conversions import arrhenius_correction
from clearwater_modules_v3.utils.numerics import clip_negative_state


def _sanitize_cache(value):
    """Replace NaN cells with 0 in a downstream-consumer cache.

    Cached rate variables (``algal_growth_rate``, ``algal_respiration_rate``,
    etc.) are read by sibling processes (DOX, Nitrogen, Phosphorus, Carbon)
    that sum many sub-fluxes into a per-cell rate. A NaN in one sub-flux
    poisons the entire cell's rate via the final NaN-guard sanitization,
    which zeroes the cell rather than just the bad term — freezing the
    cell's state at IC indefinitely (observed at newly-wet cells whose
    inputs were NaN at earlier dry steps).

    The semantically correct value when the rate computation produces NaN
    is 0 (mass conservation: ``rate * 0 algae = 0`` regardless of what
    upstream computations produced). Sanitizing at the cache source ensures
    all downstream consumers receive a finite contribution.
    """
    if isinstance(value, xr.DataArray):
        return xr.where(value.isnull(), 0, value)
    if isinstance(value, np.ndarray):
        return np.where(np.isnan(value), 0, value)
    if isinstance(value, (int, float)):
        if value != value:  # NaN
            return 0
    return value

# Defer v3 imports to break a circular-import chain that fires when v2
# is imported FIRST (test path) or via v3.processes:
#   v2.processes.__init__ imports BenthicAlgae early in RUN_ORDER
#   -> v2.benthic_algae imports v3.parameters.balgae
#   -> v3.__init__ runs and triggers v3.processes.benthic_algae
#   -> v3.processes.benthic_algae re-exports v2.benthic_algae (still loading)
# By deferring the v3 imports to first-use inside __init__/run, the v3
# package finishes its own initialization before we reach back to it.

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model
    from clearwater_modules_v3.utils.numerics import Diagnostics


logger = logging.getLogger(__name__)


# Inline fallback defaults for partitioning inputs that live in the
# v3 phosphorus / global_vars modules (not algae). Used only when the
# user does not pass them in the merged parameters dict; the
# Phase 2.A.1 follow-up will pull these via the global parameter
# registry rather than inline.
_FDP_DEFAULTS = {
    "kdpo4": 0.0,   # L/kg; matches v3 phosphorus DEFAULTS (TIP partitioning disabled)
    "Solid": 1.0,   # mg/L; matches v3 global_vars DEFAULTS
    "use_TIP": True,
    # mortality routing fractions (v1 carbon: f_pocp; not in algae DEFAULTS).
    # NSM1-SCI-A2 (gold-standard spec C1, E1 author decision 2026-05-16):
    # 0.5 -> 0.8. Dead algal carbon is predominantly particulate
    # (CE-QUAL-W2 ``APOM`` ~0.8; QUAL2K/Bowie/Chapra; v1 used 0.9). The
    # pre-fix 0.5 mis-routed ~40% of mortality C to DOC, biasing
    # DOC->DIC and DO demand. This is the *operative* value (Carbon
    # consumes the FloatingAlgae-cached rate); keep it consistent with
    # ``parameters/carbon.py`` f_pocp. See parameter_defaults_corrections.md.
    "f_pocp": 0.8,  # fraction of algal mortality C that goes to POC; (1-f_pocp) -> DOC
}


class FloatingAlgae(Process):
    variables = [
        "floating_algae",
        "solar_radiation",
        "depth",
        "water_temperature",
    ]

    # Class-level v3 defaults (Section 3.4 of design spec). Populated
    # lazily on first instantiation from ``v3.parameters.algae`` to
    # avoid the v2 <-> v3 circular import described above. After the
    # first ``FloatingAlgae(...)`` call, this attribute is the same
    # ALGAE_DEFAULTS dict every subsequent instance reads.
    DEFAULTS: dict[str, float | int | bool] = {}

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface FloatingAlgae exposes via the opportunistic-write loop in
    # ``run``. Each name maps to a ``self.<name>`` cache attribute set
    # inside ``_change_with_components`` (or as a side effect of the
    # rate / mortality-routing helpers it calls) and matches the
    # inventory in ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md``
    # §3.
    #
    # All four cache attributes consumed by sibling Processes
    # (DOX / Carbon read ``algal_growth_rate`` and ``algal_respiration_rate``;
    # Nitrogen / Phosphorus read ``algal_nh4_uptake_fraction``;
    # Carbon / Nitrogen / Phosphorus / POM read the
    # ``algal_*_from_mortality_rate`` family) are preserved attribute
    # names — Phase 5 keeps them exact.
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "algal_growth_rate",
        "algal_respiration_rate",
        "algal_death_rate",
        "algal_settling_rate",
        "algal_orgn_from_mortality_rate",
        "algal_orgp_from_mortality_rate",
        "algal_poc_from_mortality_rate",
        "algal_doc_from_mortality_rate",
        "algal_pom_from_settling_rate",
        "algal_nh4_uptake_fraction",
        "algal_light_limitation",
        "algal_nutrient_limitation_n",
        "algal_nutrient_limitation_p",
    )

    # Legacy v2 kwarg name -> v3 DEFAULTS-aligned attribute name.
    # Phase 9.A.1 wiring fix: when a legacy kwarg is omitted, the
    # rate methods read the corresponding DEFAULTS-merged attribute
    # (v1/Fortran-aligned values). When a legacy kwarg is supplied
    # explicitly, both the legacy attribute and the DEFAULTS-aligned
    # attribute are set so older callers continue to work.
    _LEGACY_TO_DEFAULTS: dict[str, str] = {
        "growth_rate_max": "mu_max_20",
        "growth_rate_correction": "mu_max_theta",
        "death_rate": "kdp_20",
        "death_rate_correction_factor": "kdp_theta",
        "repiration_rate": "krp_20",
        "repiration_rate_correction_factor": "krp_theta",
        "settling_velocity": "vsap",
        "light_limitation_constant": "KL",
        "nitrogen_michaelis_menton_constant": "KsN",
        "phosphorus_michaelis_menton_constant": "KsP",
    }

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
        settling_velocity: float | None = None,
        repiration_rate: float | None = None,
        repiration_rate_correction_factor: float | None = None,
        death_rate: float | None = None,
        death_rate_correction_factor: float | None = None,
        growth_rate_option: int | None = None,
        growth_rate_max: float | None = None,
        growth_rate_correction: float | None = None,
        phosphorus_michaelis_menton_constant: float | None = None,
        nitrogen_michaelis_menton_constant: float | None = None,
        density_michaelis_menton_constant: float | None = None,
        light_limitation_option: int | None = None,
        light_limitation_constant: float | None = None,
        light_attenuation_coefficient: float = 1.0,
        ratio_chla_carbon: float = 40.0,
        ratio_chla_nitrogen: float = 7.2,
        ratio_chla_phosphorus: float = 1.0,
    ) -> None:
        """Initialize the floating algae process.

        Args:
            parameters: Optional dict of v3 algae parameter overrides.
                Merged with the class-level ``DEFAULTS`` (v3
                ``ALGAE_DEFAULTS``). Unknown keys are warned and ignored.
            time_step: Substep cadence for this Process.
            settling_velocity: Settling velocity of floating algae (m/d).
                If None, defaults to v3 ``vsap`` (0.15) from ALGAE_DEFAULTS.
            repiration_rate: Respiration rate (1/d). If None, defaults
                to v3 ``krp_20`` (0.2).
            death_rate: Death rate (1/d). If None, defaults to v3
                ``kdp_20`` (0.15).
            growth_rate_option: 1=Multiplicative, 2=Limiting Nutrient,
                3=Harmonic Mean. If None, defaults to v3 ``growth_rate_option``.
            ratio_chla_nitrogen, ratio_chla_phosphorus, ratio_chla_carbon:
                Stoichiometric mass ratios (mg/ug-Chla). NOTE: these
                kwargs are kept for backward compatibility but are NOT
                the per-Chla ratios used by v1's algal-mortality
                routing — those are derived from
                ``self.AWn / self.AWa`` (and AWp/AWc/AWd) per the v3
                ``rna``/``rpa``/``rca``/``rda`` helpers.

        Phase 9.A.1 wiring fix: legacy kwargs that map to v3 DEFAULTS
        keys (see ``_LEGACY_TO_DEFAULTS``) now default to ``None``. When
        a kwarg is None, the DEFAULTS value is used (v1/Fortran-aligned).
        When a kwarg is provided, it overrides both the legacy attribute
        name and the DEFAULTS-aligned name so all rate methods see the
        same value regardless of which name they read.
        """
        # --- Phase 1.3 / 2.A: v3-style parameter merge (DEFAULTS + user overrides) ---
        # Lazy-load ALGAE_DEFAULTS on the first instantiation; see
        # the module-level note about the v2 <-> v3 circular import.
        # Phase 3.5 inter-process coupling: compose ``h2`` (active sediment
        # layer thickness) from POM_DEFAULTS onto the algae DEFAULTS so
        # ``algal_pom_from_settling_rate`` can be cached without requiring
        # the user to pass ``h2`` explicitly. Mirrors the Phase 5.B DOX
        # multi-group composition pattern.
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.algae import DEFAULTS as ALGAE_DEFAULTS
            from clearwater_modules_v3.parameters.pom import DEFAULTS as POM_DEFAULTS
            composed: dict[str, float | int | bool] = {}
            composed.update(ALGAE_DEFAULTS)
            # h2: active sediment layer thickness (m). Lives in POM_DEFAULTS
            # but is needed by the POM<-FloatingAlgae settling coupling.
            composed["h2"] = POM_DEFAULTS["h2"]
            type(self).DEFAULTS = composed

        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "FloatingAlgae: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in ALGAE_DEFAULTS).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Inline partitioning defaults that come from non-algae v3 groups.
        for k, v in _FDP_DEFAULTS.items():
            if not hasattr(self, k):
                setattr(self, k, user_params.get(k, v))

        # --- Legacy v2 kwargs (preserved for backward compatibility) ---
        # Phase 9.A.1 wiring fix: each legacy kwarg either takes the
        # value supplied by the user (and is also mirrored onto the
        # DEFAULTS-aligned attribute), or falls back to the
        # DEFAULTS-merged value (v1/Fortran-aligned). This guarantees
        # that default-instantiated FloatingAlgae() produces the
        # v1/Fortran-aligned rates rather than zero.
        legacy_kwargs = {
            "settling_velocity": settling_velocity,
            "repiration_rate": repiration_rate,
            "repiration_rate_correction_factor": repiration_rate_correction_factor,
            "death_rate": death_rate,
            "death_rate_correction_factor": death_rate_correction_factor,
            "growth_rate_max": growth_rate_max,
            "growth_rate_correction": growth_rate_correction,
            "light_limitation_constant": light_limitation_constant,
            "nitrogen_michaelis_menton_constant": nitrogen_michaelis_menton_constant,
            "phosphorus_michaelis_menton_constant": phosphorus_michaelis_menton_constant,
            "density_michaelis_menton_constant": density_michaelis_menton_constant,
        }
        for legacy_name, defaults_name in self._LEGACY_TO_DEFAULTS.items():
            user_val = legacy_kwargs.get(legacy_name)
            if user_val is None:
                # Fall back to the DEFAULTS-aligned attribute (already
                # set above by the merged-dict loop). Mirror onto the
                # legacy name so downstream code reading either name
                # sees the same value.
                setattr(self, legacy_name, getattr(self, defaults_name))
            else:
                # User supplied an explicit override; sync both names.
                setattr(self, legacy_name, user_val)
                setattr(self, defaults_name, user_val)

        # ``density_michaelis_menton_constant`` is only meaningful for
        # BenthicAlgae (whose ``_LEGACY_TO_DEFAULTS`` includes it). For
        # FloatingAlgae the legacy attribute does not appear in the
        # mapping; if a caller supplied it explicitly, mirror it onto
        # ``self.density_michaelis_menton_constant`` so subclasses that
        # set it directly (legacy behavior) still see it. This is a
        # no-op for default-instantiated FloatingAlgae().
        if (
            density_michaelis_menton_constant is not None
            and "density_michaelis_menton_constant"
            not in self._LEGACY_TO_DEFAULTS
        ):
            self.density_michaelis_menton_constant = (
                density_michaelis_menton_constant
            )

        # Selector kwargs: growth_rate_option / light_limitation_option are
        # already in ALGAE_DEFAULTS so the merge above set them; only
        # override if the caller explicitly passed a value.
        if growth_rate_option is not None:
            self.growth_rate_option = growth_rate_option
        if light_limitation_option is not None:
            self.light_limitation_option = light_limitation_option

        # ``light_attenuation_coefficient`` (lambda) is not in ALGAE_DEFAULTS
        # (Fortran/v1 compute lambda from the POM/Chla sum in modGlobalParam).
        # Keep the v2 scalar default as the wiring-only fallback per audit O4.
        self.light_attenuation_coefficient = light_attenuation_coefficient
        self.ratio_chla_carbon = ratio_chla_carbon
        self.ratio_chla_nitrogen = ratio_chla_nitrogen
        self.ratio_chla_phosphorus = ratio_chla_phosphorus

        # Step-scoped rate-variable cache (resolved Q10 GS-rates pattern,
        # Phase 2.A workaround until full Registry-side rate-variable
        # plumbing lands in 2.A.1). Downstream Nitrogen.run reads these
        # via ``floating_algae_process.<name>``.
        self.algal_growth_rate: ArrayLike = 0.0
        self.algal_respiration_rate: ArrayLike = 0.0
        self.algal_death_rate: ArrayLike = 0.0
        self.algal_settling_rate: ArrayLike = 0.0
        self.algal_nh4_uptake_fraction: ArrayLike = 0.5
        self.algal_orgn_from_mortality_rate: ArrayLike = 0.0
        self.algal_orgp_from_mortality_rate: ArrayLike = 0.0
        self.algal_poc_from_mortality_rate: ArrayLike = 0.0
        self.algal_doc_from_mortality_rate: ArrayLike = 0.0
        # Phase 3.5 inter-process coupling: POM consumer reads this
        # cache to avoid re-computing the v1 settling flux inline.
        # v1 ``POM_algal_settling = vsap * Ap * (AWd/AWa) / h2`` (mg/L/d).
        self.algal_pom_from_settling_rate: ArrayLike = 0.0

        # Diagnostics handle is set in init_process when a v3 Model is
        # available; otherwise a local Diagnostics is used so
        # clip_negative_state has somewhere to record. Lazy-imported
        # for the same circular-import reason as ALGAE_DEFAULTS.
        from clearwater_modules_v3.utils.numerics import Diagnostics
        self.diagnostics = Diagnostics()

        Process.__init__(self, time_step)

    @ProcessFactory.register("floating_algae")
    @staticmethod
    def from_config(
        config: dict, variable_registry: VariableRegistry
    ) -> "FloatingAlgae":
        return FloatingAlgae(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        # The nutrient-availability flags are enabled unconditionally by
        # design (NSM1 floating algae always couples to N and P; the
        # actual limitation is governed by KsN/KsP and the registry
        # concentrations, not by a process-presence gate). This matches
        # the v1 kinetic contract and the parity-test fixtures, which
        # set the same three flags. Not an incomplete stub.
        self.use_nitrate = True
        self.use_ammonium = True
        self.use_phosphate = True

        # Phase 2.A: capture the v3 Model's run-level Diagnostics handle
        # so clip_negative_state can route clip events into the run
        # diagnostics. v2's Model does not have this attribute; in that
        # case the Process keeps its locally-instantiated Diagnostics.
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Run the floating algae process.

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_change_with_components``
        (B); applies Forward Euler with unconditional clip-with-log
        (C, D); persists primary output (E); caches step-scoped rates
        on ``self.<name>`` (F); opportunistically writes diagnostics
        (G).

        """
        # --- State and forcing reads (pattern A) ---
        algae = registry.get_at_time("algae_floating", time)
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
        # ``solar_radiation`` is total broadband shortwave (W/m^2). Algal
        # photosynthesis responds to photosynthetically active radiation
        # only, and ``KL`` (``light_limitation_constant``) is a PAR-scale
        # half-saturation value. Convert at the process boundary, mirroring
        # NSM1 v1 ``processes.py:287`` (``PAR = q_solar * Fr_PAR``).
        # NSM1-SCI-A3 fix (gold-standard spec B1): pre-fix v3 passed total
        # shortwave straight into ``limit_light`` against the PAR-scale
        # ``KL``, under-limiting light and over-predicting algal growth
        # ~30-60% wherever light is the binding constraint (a v1->v3
        # regression -- v1 applied Fr_PAR, v3 had dropped it).
        solar_shortwave = registry.get_at_time("solar_radiation", time)
        solar = solar_shortwave * self.Fr_PAR

        # --- Fused rate composition (pattern B) ---
        rate, components = self._change_with_components(
            algae=algae,
            depth=depth,
            water_temperature=water_temperature,
            phosphorus_total_inorganic=phosphorus_total_inorganic,
            ammonium=ammonium,
            nitrate=nitrate,
            solar=solar,
        )

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        # Names match REGISTRY_DIAGNOSTICS. The growth/respiration/
        # settling/death/mortality-routing/nh4-uptake-fraction caches
        # are *also* set as side effects of the helpers called inside
        # ``_change_with_components``; the setattr loop here is
        # idempotent on those names (same value) and adds the new
        # ``algal_light_limitation`` / ``algal_nutrient_limitation_*``
        # entries that did not have a self-cache before Phase 5.
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])

        # --- Forward Euler in days (pattern C; Bug #4) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        algae_new = algae + rate * dt_days

        # --- Clip-with-log per the resolved Q7 contract (pattern D) ---
        # Clip target is exactly 0 (Monod ratios are well-defined at C=0).
        # Step attribution is automatic via ``diagnostics.current_step``
        # (Phase 0.6 Q1); import is module-level (Phase 1.C).
        algae_new = clip_negative_state(
            algae_new, "algae_floating", self.diagnostics
        )

        # --- Persist primary output (pattern E; Bug #16) ---
        registry.set_at_time("algae_floating", time, algae_new)

        # --- Opportunistic diagnostic registry writes (pattern G) ---
        for name in self.REGISTRY_DIAGNOSTICS:
            if name in registry:
                registry.set_at_time(name, time, components[name])

    # ------------------------------------------------------------------
    # Rate-composition helpers
    # ------------------------------------------------------------------

    def _change_with_components(
        self,
        *,
        algae: ArrayLike,
        depth: ArrayLike,
        water_temperature: ArrayLike,
        phosphorus_total_inorganic: ArrayLike,
        ammonium: ArrayLike,
        nitrate: ArrayLike,
        solar: ArrayLike,
    ) -> tuple[ArrayLike, dict]:
        """Compute ``(rate, components)`` for FloatingAlgae.

        ``rate`` is the net per-day algal rate of change (ug-Chla/L/d).
        ``components`` is the dict[str, ArrayLike] indexed by
        ``REGISTRY_DIAGNOSTICS``.

        Code-motion-only refactor of ``run``'s former inline
        composition (§11.6): operand order, intermediate names, and
        kinetic-helper calls are preserved verbatim. The new
        ``algal_light_limitation`` / ``algal_nutrient_limitation_*``
        diagnostics are computed by separate (pure-function) calls to
        ``limit_light`` / ``limit_nitrogen`` / ``limit_phosphorus``;
        ``rate()`` invokes the same helpers internally and the values
        are bit-identical.

        """
        # Bug #15: compute fdp via the v3 partitioning utility instead of
        # the previous hard-coded 0.5.
        from clearwater_modules_v3.utils.partitioning import fdp as fdp_partition
        phosphate_fraction_dissolved = fdp_partition(
            use_TIP=self.use_TIP,
            Solid=self.Solid,
            kdpo4=self.kdpo4,
        )

        # Net per-day rate. Side effects: rate() invokes rate_growth /
        # rate_respiration / rate_settling, each of which writes
        # ``self.algal_growth_rate`` / ``algal_respiration_rate`` /
        # ``algal_settling_rate`` as cache attributes. _cache_mortality_rates
        # then writes the death / orgn / orgp / poc / doc / pom-settling
        # mortality routing caches.
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

        # Cache mortality routing variables (sets the algal_*_from_*
        # caches as side effects).
        self._cache_mortality_rates(algae, water_temperature)

        # NH4-uptake fraction so Nitrogen.run can split algal N uptake
        # between NH4 and NO3 (v1 ApUptakeFr_NH4). Side-effect parity
        # with the pre-Phase-5 ``run`` body, which assigned this
        # attribute directly. The pattern F setattr loop in ``run``
        # (which iterates REGISTRY_DIAGNOSTICS) is idempotent on this
        # name.
        nh4_uptake_fraction = self._ap_uptake_fr_nh4(
            ammonium=ammonium, nitrate=nitrate
        )
        self.algal_nh4_uptake_fraction = nh4_uptake_fraction

        # Limitation diagnostics (pure-function recomputes; same values
        # as those used inside rate()).
        limit_phosphorus = self.limit_phosphorus(
            concentration=phosphorus_total_inorganic,
            fraction_dissolved=phosphate_fraction_dissolved,
        )
        limit_nitrogen = self.limit_nitrogen(
            ammonium=ammonium, nitrate=nitrate
        )
        limit_light = self.limit_light(
            algae=algae, depth=depth, surface_light_intensity=solar
        )

        components = {
            "algal_growth_rate": self.algal_growth_rate,
            "algal_respiration_rate": self.algal_respiration_rate,
            "algal_death_rate": self.algal_death_rate,
            "algal_settling_rate": self.algal_settling_rate,
            "algal_orgn_from_mortality_rate": self.algal_orgn_from_mortality_rate,
            "algal_orgp_from_mortality_rate": self.algal_orgp_from_mortality_rate,
            "algal_poc_from_mortality_rate": self.algal_poc_from_mortality_rate,
            "algal_doc_from_mortality_rate": self.algal_doc_from_mortality_rate,
            "algal_pom_from_settling_rate": self.algal_pom_from_settling_rate,
            "algal_nh4_uptake_fraction": nh4_uptake_fraction,
            "algal_light_limitation": limit_light,
            "algal_nutrient_limitation_n": limit_nitrogen,
            "algal_nutrient_limitation_p": limit_phosphorus,
        }

        return rate, components
    def _cache_mortality_rates(
        self, algae: ArrayLike, water_temperature: ArrayLike
    ) -> None:
        """Compute and cache the algal mortality routing rates.

        Per v1 (``processes.py`` ApDeath_OrgN, ApDeath_OrgP,
        POC_algal_mortality, DOC_algal_mortality):

        - ApDeath_OrgN = rna * ApDeath
        - ApDeath_OrgP = rpa * ApDeath
        - POC_algal_mortality = f_pocp * kdp_tc * rca * Ap
                              = f_pocp * rca * ApDeath
        - DOC_algal_mortality = (1 - f_pocp) * kdp_tc * rca * Ap
                              = (1 - f_pocp) * rca * ApDeath

        rna/rpa/rca are AWn/AWa, AWp/AWa, AWc/AWa per v1 (lines 308-348).
        ApDeath is computed from rate_death (kdp_tc * Ap).

        Also caches ``algal_pom_from_settling_rate`` (Phase 3.5
        inter-process coupling) per v1
        ``POM_algal_settling = vsap * Ap * (AWd/AWa) / h2`` (mg/L/d).
        Note: this is a *settling* flux, not a mortality flux. It is
        cached here for symmetry with the mortality routing helpers and
        because the POM consumer reads it via the same step-scoped Q10
        GS-rates pattern.
        """
        ap_death = self.rate_death(algae, water_temperature)

        rna = self.AWn / self.AWa  # mg-N/ug-Chla
        rpa = self.AWp / self.AWa  # mg-P/ug-Chla
        rca = self.AWc / self.AWa  # mg-C/ug-Chla

        self.algal_death_rate = ap_death
        self.algal_orgn_from_mortality_rate = rna * ap_death
        self.algal_orgp_from_mortality_rate = rpa * ap_death
        self.algal_poc_from_mortality_rate = self.f_pocp * rca * ap_death
        self.algal_doc_from_mortality_rate = (1.0 - self.f_pocp) * rca * ap_death

        # POM source from settling of algal biomass (mg/L/d).
        # v1 dPOMdt term: vsap * Ap * (AWd/AWa) / h2.
        self.algal_pom_from_settling_rate = (
            self.vsap * algae * (self.AWd / self.AWa) / self.h2
        )

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
        """Compute the rate of change of floating algae (ug-Chla/L/d)."""
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

        return (
            self.rate_growth(
                algae,
                water_temperature,
                limit_phosphorus,
                limit_nitrogen,
                limit_light,
            )
            - self.rate_death(algae, water_temperature)
            - self.rate_respiration(algae, water_temperature)
            - self.rate_settling(algae, depth)
        )

    def rate_growth(
        self,
        algae: ArrayLike,
        water_temperature: ArrayLike,
        limit_phosphorus: ArrayLike,
        limit_nitrogen: ArrayLike,
        limit_light: ArrayLike,
    ) -> ArrayLike:
        """Compute the rate of growth of floating algae."""

        growth_rate = arrhenius_correction(
            water_temperature,
            self.growth_rate_max,
            self.growth_rate_correction,
        )

        # Multiplicative
        if self.growth_rate_option == 1:
            rate = growth_rate * limit_phosphorus * limit_nitrogen * limit_light
        # Limiting Nutrient
        elif self.growth_rate_option == 2:
            rate = xr.where(
                limit_phosphorus > limit_nitrogen,
                growth_rate * limit_nitrogen * limit_light,
                growth_rate * limit_phosphorus * limit_light,
            )
        # Harmonic Mean
        # Phase 9.A.1 audit F14: zero-guard should fire when FN==0 OR
        # FP==0 (avoid division by zero in 1/FN + 1/FP). The previous
        # form fired when FP==1 (i.e. when phosphorus is fully
        # non-limiting), shutting down growth precisely when P is
        # saturating — opposite of intent. v1/Fortran fire on
        # ``(FN == 0) | (FP == 0)``.
        elif self.growth_rate_option == 3:
            rate_raw = (
                growth_rate
                * limit_light
                * 2.0
                / (1.0 / limit_nitrogen + 1.0 / limit_phosphorus)
            )
            rate = xr.where(
                (limit_nitrogen == 0.0) | (limit_phosphorus == 0.0),
                0,
                rate_raw,
            )
        else:
            raise ValueError("Invalid growth rate option")

        result = rate * algae
        # Cache for downstream consumers (Q10 GS-rates pattern).
        # Sanitize NaN at the cache source — see ``_sanitize_cache`` docstring
        # for rationale (NaN here poisons DOX/Nitrogen/Phosphorus rate sums).
        result = _sanitize_cache(result)
        self.algal_growth_rate = result
        return result

    def rate_death(self, algae: ArrayLike, water_temperature: ArrayLike) -> ArrayLike:
        """Compute the rate of death of floating algae."""
        corrected_death_rate = arrhenius_correction(
            water_temperature,
            self.death_rate,
            self.death_rate_correction_factor,
        )
        return _sanitize_cache(algae * corrected_death_rate)

    def rate_respiration(
        self, algae: ArrayLike, water_temperature: ArrayLike
    ) -> ArrayLike:
        """Compute the rate of respiration of floating algae."""
        corrected_respiration_rate = arrhenius_correction(
            water_temperature,
            self.repiration_rate,
            self.repiration_rate_correction_factor,
        )
        result = algae * corrected_respiration_rate
        # See ``rate_growth`` — sanitize NaN at the cache source.
        result = _sanitize_cache(result)
        self.algal_respiration_rate = result
        return result

    def rate_settling(self, algae: ArrayLike, depth: ArrayLike) -> ArrayLike:
        """Compute the rate of settling of floating algae."""
        result = _sanitize_cache(algae / depth * self.settling_velocity)
        self.algal_settling_rate = result
        return result

    def limit_phosphorus(
        self,
        concentration: ArrayLike,
        fraction_dissolved: ArrayLike,
    ) -> ArrayLike:
        """Compute the limiting phosphorus for floating algae."""

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

        # Bug fix: x == np.nan is always False per IEEE 754. Use isnull.
        if isinstance(rate_raw, xr.DataArray):
            rate = xr.where(rate_raw.isnull(), 0.0, rate_raw)
        else:
            rate = xr.where(np.isnan(rate_raw), 0.0, rate_raw)
        # any rates > 1 are limiting
        rate = xr.where(rate > 1, 1, rate)

        return rate

    def limit_nitrogen(
        self,
        nitrate: ArrayLike,
        ammonium: ArrayLike,
    ) -> ArrayLike:
        """Compute the limiting nitrogen for floating algae."""
        if not self.use_nitrate and not self.use_ammonium:
            return 1.0

        n_concentration = nitrate if self.use_nitrate else 0.0
        n_concentration = n_concentration + (ammonium if self.use_ammonium else 0.0)

        rate_raw = n_concentration / (
            self.nitrogen_michaelis_menton_constant + n_concentration
        )

        # Bug fix: x == np.nan is always False per IEEE 754. Use isnull.
        if isinstance(rate_raw, xr.DataArray):
            rate = xr.where(rate_raw.isnull(), 0.0, rate_raw)
        else:
            rate = xr.where(np.isnan(rate_raw), 0.0, rate_raw)
        # any rates > 1 are limiting
        rate = xr.where(rate > 1, 1, rate)

        return rate

    def limit_light(
        self,
        algae: ArrayLike,
        depth: ArrayLike,
        surface_light_intensity: ArrayLike,
    ) -> ArrayLike:
        """Compute the limiting light for floating algae."""

        # Half-saturation light limitation
        # Phase 9.A.1 audit F5: the ``np.log`` argument must be the
        # ratio (KL+PAR) / (KL+PAR*exp(-Ld)). The previous form split
        # numerator and denominator across the *  /  operators so the
        # log received only (KL+PAR), and the divisor multiplied the
        # final result instead of forming the log argument.
        if self.light_limitation_option == 1:
            raw_rate = (
                (1.0 / (self.light_attenuation_coefficient * depth))
                * np.log(
                    (self.light_limitation_constant + surface_light_intensity)
                    / (
                        self.light_limitation_constant
                        + surface_light_intensity
                        * np.exp(-(self.light_attenuation_coefficient * depth))
                    )
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

        # conditions where growth is limited to zero
        # no algae present
        rate = xr.where(algae <= 0, 0, raw_rate)
        # light cannot penetrate deep enough
        rate = xr.where(self.light_attenuation_coefficient * depth <= 0, 0, rate)
        rate = xr.where(raw_rate > 1, 1.0, rate)
        return rate

    # ------------------------------------------------------------------
    # NH4 uptake fractionation (v1 ApUptakeFr_NH4) and ammonium coupling
    # ------------------------------------------------------------------

    def _ap_uptake_fr_nh4(
        self, ammonium: ArrayLike, nitrate: ArrayLike
    ) -> ArrayLike:
        """v1 ApUptakeFr_NH4 (lines 1206-1247).

        Returns the NH4 fraction of total inorganic-N uptake by floating
        algae. Uses the algae preference factor PN
        (= ``floating_algae_preference_factor`` here, mirroring v1 PN) and
        the Monod-style competitive form when both NH4 and NO3 are
        active.
        """
        # PN preference factor: prefer the v3 DEFAULTS-merged name PN if
        # the user supplied it; otherwise fall back to a standard 0.5.
        pn = getattr(self, "PN", 0.5)

        if self.use_ammonium and not self.use_nitrate:
            return xr.ones_like(ammonium) * 1.0 if isinstance(ammonium, xr.DataArray) else 1.0
        if not self.use_ammonium and self.use_nitrate:
            return xr.zeros_like(nitrate) if isinstance(nitrate, xr.DataArray) else 0.0
        if not self.use_ammonium and not self.use_nitrate:
            return 0.5

        # Both NH4 and NO3 active: competitive Monod form.
        denom = pn * ammonium + (1.0 - pn) * nitrate
        # Avoid div/0 at very low N: fall back to PN
        result = xr.where(denom > 0, pn * ammonium / denom, pn)
        if isinstance(result, xr.DataArray):
            result = xr.where(result.isnull(), pn, result)
        return result

    def ammonium_respiration(self) -> ArrayLike:
        """v1 NH4_ApRespiration (line 1486): rna * ApRespiration.

        Returns mg-N/L/d transferred from algal respiration to NH4.
        Uses the cached ``algal_respiration_rate`` written by ``run``;
        falls back to 0 if ``run`` has not yet been called.
        """
        rna = self.AWn / self.AWa  # mg-N/ug-Chla
        return rna * self.algal_respiration_rate

    def ammonium_growth(self) -> ArrayLike:
        """v1 NH4_ApGrowth (line 1504): ApUptakeFr_NH4 * rna * ApGrowth.

        Returns mg-N/L/d removed from NH4 by algal growth.
        Uses the cached ``algal_growth_rate`` and
        ``algal_nh4_uptake_fraction`` written by ``run``.
        """
        rna = self.AWn / self.AWa  # mg-N/ug-Chla
        return self.algal_nh4_uptake_fraction * rna * self.algal_growth_rate

    # ------------------------------------------------------------------
    # Algal mortality routing helpers (Q10 GS-rates contract)
    # ------------------------------------------------------------------

    def death_to_orgn(self) -> ArrayLike:
        """rna * ApDeath -> OrgN (mg-N/L/d). v1 ApDeath_OrgN."""
        return self.algal_orgn_from_mortality_rate

    def death_to_orgp(self) -> ArrayLike:
        """rpa * ApDeath -> OrgP (mg-P/L/d). v1 ApDeath_OrgP."""
        return self.algal_orgp_from_mortality_rate

    def death_to_poc(self) -> ArrayLike:
        """f_pocp * rca * ApDeath -> POC (mg-C/L/d). v1 POC_algal_mortality."""
        return self.algal_poc_from_mortality_rate

    def death_to_doc(self) -> ArrayLike:
        """(1 - f_pocp) * rca * ApDeath -> DOC (mg-C/L/d). v1 DOC_algal_mortality."""
        return self.algal_doc_from_mortality_rate
