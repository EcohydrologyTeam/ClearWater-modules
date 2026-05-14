from datetime import timedelta, datetime
import logging

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.base import Process, ProcessFactory
from clearwater_data.variables import VariableRegistry
from clearwater_data.custom_types import ArrayLike
from clearwater_modules_v3.utils.conversions import arrhenius_correction
from clearwater_modules_v3.utils.numerics import clip_negative_state

logger = logging.getLogger(__name__)


from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


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

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface Nitrogen exposes via the opportunistic-write loop in
    # ``run``. Each name maps to a ``self.<name>`` cache attribute set
    # inside ``_change_with_components`` and matches the inventory in
    # ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3.
    #
    # ``nitrification_flux_rate`` and ``denitrification_flux_rate`` are
    # **preserved** attribute names that DOX, Alkalinity, and N2 already
    # consume via ``getattr(nitrogen_process, ...)``. Renaming them
    # would silently break sibling reads. Phase 4 keeps the names exact.
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "nitrification_flux_rate",
        "denitrification_flux_rate",
        "nh4_from_bed",
        "no3_from_bed_denit",
        "orgn_hydrolysis_rate",
        "orgn_settling_rate",
        "nh4_algal_growth_rate",
        "no3_algal_growth_rate",
        "nh4_algal_resp_rate",
        "nh4_balgae_resp_rate",
    )

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
        denitrification_rate: ArrayLike | None = None,
        denitrification_theta: ArrayLike | None = None,
        nitrification_rate: ArrayLike | None = None,
        nitrification_theta: ArrayLike | None = None,
        sediment_denitrification_rate: ArrayLike | None = None,
        sediment_denitrification_theta: ArrayLike | None = None,
        sediment_ammonium_release_rate: ArrayLike | None = None,
        sediment_ammonium_release_theta: ArrayLike | None = None,
        ammonium_decay_rate: ArrayLike | None = None,
        ammonium_decay_theta: ArrayLike | None = None,
        floating_algae_preference_factor: ArrayLike = 0.5,
        settling_velocity: ArrayLike = 1.0,
        death_rate: ArrayLike = 1.0,
        float_algea_faction_uptake_from_nitrate: ArrayLike | None = None,
        nitrification_oxygen_inhibition_factor: ArrayLike | None = None,
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

        # Phase 9.F.C defensive guard: v3 1.0.0 Nitrogen does NOT gate
        # ammonium_from_bed or nitrate_bed_denitrification by
        # ``use_SedFlux``; those terms are silenced solely by the
        # ``rnh4_20 = 0`` and ``vno3_20 = 0`` defaults. If a user opts
        # into ``use_SedFlux=True`` they are signaling intent to enable
        # the full sediment-flux feature, which requires the NSM2
        # diagenesis path that is not implemented in v3 1.0.0. Refuse
        # explicitly rather than silently producing partial behavior.
        # See corrections doc Section 2.1.
        if user_params.get("use_SedFlux", False):
            raise NotImplementedError(
                "Nitrogen: use_SedFlux=True is not implemented in v3 1.0.0. "
                "The full sediment-flux feature requires the NSM2 diagenesis "
                "path. Set rnh4_20 / vno3_20 directly to specify constant "
                "sediment release rates for site-specific calibration "
                "(without use_SedFlux). See parameter_defaults_corrections.md "
                "Section 2.1."
            )

        unknown_keys = set(user_params) - set(self.DEFAULTS) - {"use_SedFlux"}
        for key in sorted(unknown_keys):
            logger.warning(
                "Nitrogen: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in NITROGEN_DEFAULTS).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        merged.pop("use_SedFlux", None)
        for k, v in merged.items():
            setattr(self, k, v)

        # --- Legacy v2 kwargs (preserved for backward compatibility) ---
        # Phase 9.A.2 (audit findings N1, N2, N4, N10, N11): the kinetic
        # methods now read from the DEFAULTS-key attributes
        # (``self.knit_20``, ``self.kdnit_20``, ``self.rnh4_20``, etc.)
        # rather than the legacy v2 kwargs. Legacy kwargs default to
        # ``None`` and only override the DEFAULTS-key attribute when
        # explicitly provided by the caller. This preserves backwards
        # compatibility with existing tests/YAML configs that pass
        # ``nitrification_rate=0.1`` etc. while making
        # ``Nitrogen()`` (no kwargs) produce v1/Fortran-correct kinetics.
        #
        # Mapping legacy kwarg -> DEFAULTS key (synced both directions):
        #   nitrification_rate                       -> knit_20
        #   nitrification_theta                      -> knit_theta
        #   denitrification_rate                     -> kdnit_20
        #   denitrification_theta                    -> kdnit_theta
        #   sediment_ammonium_release_rate           -> rnh4_20
        #   sediment_ammonium_release_theta          -> rnh4_theta
        #   sediment_denitrification_rate            -> vno3_20
        #   sediment_denitrification_theta           -> vno3_theta
        #   nitrification_oxygen_inhibition_factor   -> KNR
        _legacy_to_defaults = {
            "nitrification_rate": ("knit_20", nitrification_rate),
            "nitrification_theta": ("knit_theta", nitrification_theta),
            "denitrification_rate": ("kdnit_20", denitrification_rate),
            "denitrification_theta": ("kdnit_theta", denitrification_theta),
            "sediment_ammonium_release_rate": ("rnh4_20", sediment_ammonium_release_rate),
            "sediment_ammonium_release_theta": ("rnh4_theta", sediment_ammonium_release_theta),
            "sediment_denitrification_rate": ("vno3_20", sediment_denitrification_rate),
            "sediment_denitrification_theta": ("vno3_theta", sediment_denitrification_theta),
            "nitrification_oxygen_inhibition_factor": ("KNR", nitrification_oxygen_inhibition_factor),
        }
        for legacy_name, (defaults_key, legacy_value) in _legacy_to_defaults.items():
            if legacy_value is not None:
                # User explicitly supplied the legacy kwarg; sync both
                # naming schemes onto the same value.
                setattr(self, legacy_name, legacy_value)
                setattr(self, defaults_key, legacy_value)
            else:
                # Legacy kwarg not supplied; keep DEFAULTS-key value and
                # mirror it onto the legacy attribute for back-compat
                # with any caller still reading the legacy name.
                setattr(self, legacy_name, getattr(self, defaults_key))

        # ammonium_decay_rate / ammonium_decay_theta: phantom NH4 source
        # term (audit finding N2). v1 and Fortran NSM1 have no analogue.
        # The kinetic method ``ammonium_decay_nitrate`` is dropped from
        # the NH4 budget below (audit fix); these attributes are kept on
        # the instance for back-compat with callers that might inspect
        # them, but they no longer feed into ``change_ammonium``.
        self.ammonium_decay_rate = (
            ammonium_decay_rate if ammonium_decay_rate is not None else 0.0
        )
        self.ammonium_decay_theta = (
            ammonium_decay_theta if ammonium_decay_theta is not None else 1.0
        )

        self.floating_algae_preference_factor = floating_algae_preference_factor
        self.settling_velocity = settling_velocity
        # Phase 2.B Bug #12: ``death_rate`` is preserved as a legacy
        # kwarg for back-compat with existing v2 unit tests, but the
        # canonical algal-death routing in v3 reads
        # ``floating_algae_process.algal_death_rate`` /
        # ``algal_orgn_from_mortality_rate`` (Phase 2.A populates these).
        self.death_rate = death_rate

        # Phase 9.A.2 audit finding N12: ``float_algea_faction_uptake_from_nitrate``
        # was a static parameter (default 1.0) that broke algal-N mass
        # balance. The corrected ``nitrate_uptake_floating_algae`` reads
        # the dynamic ``1 - algal_nh4_uptake_fraction`` from FloatingAlgae
        # instead. Retain the attribute for back-compat with tests/YAML
        # that set it, but it is now a deprecated no-op.
        self.float_algea_faction_uptake_from_nitrate = (
            float_algea_faction_uptake_from_nitrate
            if float_algea_faction_uptake_from_nitrate is not None
            else 1.0
        )

        # Phase 8.A: legacy v2 algal-uptake attributes that are read by
        # ``nitrate_uptake_floating_algae`` / ``nitrate_uptake_benthic_algae``
        # / ``ammonium_*`` paths but were previously injected out-of-band by
        # the ``init_from_file`` YAML config path. Without that path,
        # ``Nitrogen()`` instantiates fine but ``Nitrogen.run()`` raised
        # ``AttributeError`` on the first algal-coupled step. Provide
        # sensible defaults sourced from the v3 algae stoichiometric ratios
        # (``ALGAE_DEFAULTS`` / ``BALGAE_DEFAULTS``) so a bare
        # ``Nitrogen()`` is fully runnable. Callers can still override any
        # of these by setting them on the instance after construction
        # (which is what the YAML config path does).
        #
        # NOTE: ``benthic_algea_faction_uptake_from_nitrate`` retains the
        # legacy "algea" typo for back-compat with existing YAML configs
        # and tests; do NOT rename it. A future v3.x can deprecate it.
        from clearwater_modules_v3.parameters.algae import DEFAULTS as ALGAE_DEFAULTS
        from clearwater_modules_v3.parameters.balgae import DEFAULTS as BALGAE_DEFAULTS
        # AWn / BWn / AWa are raw stoichiometric weights "per stoichiometric
        # unit", NOT concentration ratios. The actual mg-N/ug-Chla ratio
        # is rna = AWn/AWa (computed at line 729 below); the actual
        # mg-N/mg-D mass fraction is rnb = BWn/BWd. Treating AWn / BWn
        # directly as mg-N/ug-Chla / mg-N/g-D (as some test helpers and
        # earlier comments did) overstates the algal-N stoichiometry by
        # AWa = 1000x and was the source of a closed-system N
        # conservation bug (Phase 9.G commit ee31218).
        self.floating_algae_nitrogen_weight = ALGAE_DEFAULTS["AWn"]      # mg-N per stoichiometric unit
        self.benthic_algae_nitrogen_weight = BALGAE_DEFAULTS["BWn"]      # mg-N per stoichiometric unit
        self.algal_chlorophyll = ALGAE_DEFAULTS["AWa"]                   # ug-Chla per stoichiometric unit
        self.benthic_algea_faction_uptake_from_nitrate = 0.5             # legacy "algea" typo (preserved)
        self.fraction_bottom_area = 1.0                                  # unitless

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

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_change_with_components``
        (B); applies Forward Euler with unconditional clip-with-log (C,
        D); persists primary outputs (E); caches step-scoped rates on
        ``self.<name>`` (F); opportunistically writes diagnostics (G).


        Phase 2.B fixes (preserved through this refactor):
        * Bug #1 / #2: additive Forward Euler ``X_new = X + rate * dt_days``.
        * Bug #16: persist via ``registry.set_at_time``.
        * Q7 clip-with-log via ``clip_negative_state``.
        * OrgN as third state variable.

        Cached attribute names ``nitrification_flux_rate`` and
        ``denitrification_flux_rate`` are preserved exactly — DOX,
        Alkalinity, and N2 already consume them via ``getattr`` and
        renaming would silently break sibling reads.
        """
        # --- State and forcing reads (pattern A) ---
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

        # --- Fused rate composition (pattern B) ---
        ammonium_rate, nitrate_rate, orgn_rate, components = (
            self._change_with_components(
                nitrate=nitrate,
                ammonium=ammonium,
                organic_nitrogen=organic_nitrogen,
                temperature=temperature,
                depth=depth,
                oxygen_dissolved=oxygen_dissolved,
            )
        )

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        # Names match REGISTRY_DIAGNOSTICS. The two preserved names
        # ``nitrification_flux_rate`` / ``denitrification_flux_rate``
        # carry the same values DOX/Alkalinity/N2 read pre-Phase-4
        # (computed via ``ammonium_nitrification`` and
        # ``nitrate_denitrification`` exactly as before).
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])

        # --- Forward Euler in days (pattern C) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        ammonium_new = ammonium + ammonium_rate * dt_days
        nitrate_new = nitrate + nitrate_rate * dt_days
        organic_nitrogen_new = organic_nitrogen + orgn_rate * dt_days

        # --- Clip-with-log per the resolved Q7 contract (pattern D) ---
        ammonium_new = clip_negative_state(ammonium_new, "ammonium", self.diagnostics)
        nitrate_new = clip_negative_state(nitrate_new, "nitrate", self.diagnostics)
        organic_nitrogen_new = clip_negative_state(
            organic_nitrogen_new, "organic_nitrogen", self.diagnostics
        )

        # --- Persist primary outputs (pattern E; Bug #16) ---
        registry.set_at_time("ammonium", time, ammonium_new)
        registry.set_at_time("nitrate", time, nitrate_new)
        if "organic_nitrogen" in registry:
            registry.set_at_time("organic_nitrogen", time, organic_nitrogen_new)

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
        nitrate: ArrayLike,
        ammonium: ArrayLike,
        organic_nitrogen: ArrayLike,
        temperature: ArrayLike,
        depth: ArrayLike,
        oxygen_dissolved: ArrayLike,
    ) -> tuple[ArrayLike, ArrayLike, ArrayLike, dict]:
        """Compute ``(ammonium_rate, nitrate_rate, orgn_rate, components)``.

        Code-motion-only refactor of ``run``'s former inline composition
        (§11.6): operand order, intermediate names, kinetic-helper calls,
        and per-state rate compositions are preserved verbatim.

        The ``components`` dict is populated from the same intermediates
        the integrator consumes; pure-function sub-flux helpers
        (``ammonium_nitrification``, ``nitrate_denitrification``,
        ``ammonium_from_bed``, ``nitrate_bed_denitrification``,
        ``organic_nitrogen_to_ammonium_hydrolysis``,
        ``organic_nitrogen_settling``, the algal coupling helpers) are
        called once each here for the cache; ``change_ammonium`` /
        ``change_nitrate`` / ``change_organic_nitrogen`` invoke the
        same helpers internally to compose the per-state rates. The
        recomputation is bit-identical and matches the pre-Phase-4
        behaviour exactly.

        """
        # --- Step-scoped flux caches (preserved attribute names) ---
        # ``ammonium_nitrification`` and ``nitrate_denitrification`` are
        # pure functions; the values populated here match what
        # ``change_ammonium`` / ``change_nitrate`` compute internally.
        # Sibling consumers (DOX, Alkalinity, N2) read these attribute
        # names via ``getattr`` and rely on bit-identical values.
        nitrification_flux = self.ammonium_nitrification(
            ammonium,
            temperature,
            oxygen_dissolved,
        )
        denitrification_flux = self.nitrate_denitrification(
            oxygen_dissolved,
            self.KsOxdn,
            nitrate,
            temperature,
        )

        # --- Per-state rate compositions (verbatim from pre-Phase-4) ---
        ammonium_rate = self.change_ammonium(
            nitrate,
            ammonium,
            temperature,
            depth,
            oxygen_dissolved,
            organic_nitrogen=organic_nitrogen,
        )
        nitrate_rate = self.change_nitrate(
            nitrate,
            ammonium,
            temperature,
            depth,
            oxygen_dissolved,
        )
        orgn_rate = self.change_organic_nitrogen(
            organic_nitrogen=organic_nitrogen,
            temperature=temperature,
            depth=depth,
        )

        # --- Sub-fluxes for the components dict ---
        # Each is a pure-function recompute mirroring what the change_*
        # methods invoke internally; same arguments → identical values.
        # When ``use_OrgN`` is False the OrgN sub-fluxes default to 0
        # (matching ``change_organic_nitrogen``'s early-return).
        nh4_from_bed = self.ammonium_from_bed(depth=depth, temperature=temperature)
        no3_from_bed_denit = self.nitrate_bed_denitrification(
            depth, nitrate, temperature
        )
        orgn_hydrolysis = self.organic_nitrogen_to_ammonium_hydrolysis(
            organic_nitrogen=organic_nitrogen,
            temperature=temperature,
        )
        if getattr(self, "use_OrgN", True):
            orgn_settling = self.organic_nitrogen_settling(
                organic_nitrogen=organic_nitrogen,
                temperature=temperature,
                depth=depth,
            )
        else:
            orgn_settling = 0.0

        # Algal NH4 / NO3 coupling diagnostics (sums of floating +
        # benthic). The change_ammonium / change_nitrate methods sum
        # the per-source contributions internally; we mirror that sum
        # here so the components dict carries the consumer-visible
        # totals. Float / benthic per-source values stay accessible via
        # the per-source attribute caches on FloatingAlgae / BenthicAlgae.
        nh4_algal_resp = self.ammonium_floating_respiration()
        nh4_balgae_resp = self.ammonium_benthic_respiration()
        nh4_algal_growth = (
            self.ammonium_floating_growth() + self.ammonium_benthic_growth()
        )

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
        no3_algal_growth = (
            self.nitrate_uptake_floating_algae(nitrate, ammonium, float_algae_growth)
            + self.nitrate_uptake_benthic_algae(
                nitrate, ammonium, benthic_algae_growth, depth
            )
        )

        components = {
            "nitrification_flux_rate": nitrification_flux,
            "denitrification_flux_rate": denitrification_flux,
            "nh4_from_bed": nh4_from_bed,
            "no3_from_bed_denit": no3_from_bed_denit,
            "orgn_hydrolysis_rate": orgn_hydrolysis,
            "orgn_settling_rate": orgn_settling,
            "nh4_algal_growth_rate": nh4_algal_growth,
            "no3_algal_growth_rate": no3_algal_growth,
            "nh4_algal_resp_rate": nh4_algal_resp,
            "nh4_balgae_resp_rate": nh4_balgae_resp,
        }

        return ammonium_rate, nitrate_rate, orgn_rate, components
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

        # Phase 9.A.2 audit finding N2: dropped the phantom
        # ``ammonium_decay_nitrate`` source term. v1 (processes.py:1584)
        # and Fortran NSM1 (modNitrogen.f90:296) have no analogue. The
        # legacy v2 kwarg ``ammonium_decay_rate=1.0/d`` (default) was
        # injecting a first-order NH4 *source* with no matching sink,
        # making default-instantiated NH4 grow exponentially. The
        # ``ammonium_decay_nitrate`` method is retained for back-compat
        # with any caller that invokes it directly, but it is no longer
        # part of the NH4 budget.
        rate = (
            -self.ammonium_nitrification(
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
        # Phase 9.A.2 audit finding N4: read NITROGEN_DEFAULTS attributes
        # (``rnh4_20=0`` v1/Fortran default) instead of the legacy v2 kwarg
        # (was 1.0/d). Legacy kwarg, when supplied, syncs onto ``rnh4_20``
        # in __init__.
        rate = arrhenius_correction(
            temperature,
            self.rnh4_20,
            self.rnh4_theta,
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

        # Phase 9.A.2 audit finding N1: read NITROGEN_DEFAULTS attributes
        # (``knit_20=0.1``, ``knit_theta=1.083`` v1/Fortran defaults).
        # Legacy ``nitrification_rate``/``nitrification_theta`` kwargs,
        # when supplied, sync onto these in __init__.
        rate_corrected = arrhenius_correction(
            temperature, self.knit_20, self.knit_theta
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

        # Phase 9.A.2 audit finding N10: read NITROGEN_DEFAULTS attributes
        # (``kdnit_20=0.002``, ``kdnit_theta=1.08`` v1 defaults). Legacy
        # ``denitrification_rate``/``denitrification_theta`` kwargs, when
        # supplied, sync onto these in __init__.
        rate_corrected = arrhenius_correction(
            temperature, self.kdnit_20, self.kdnit_theta
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
        # Phase 9.A.2 audit finding N11: read NITROGEN_DEFAULTS attributes
        # (``vno3_20=0`` v1/Fortran default). Legacy
        # ``sediment_denitrification_rate``/``sediment_denitrification_theta``
        # kwargs, when supplied, sync onto these in __init__.
        rate_corrected = arrhenius_correction(
            temperature,
            self.vno3_20,
            self.vno3_theta,
        )

        return nitrate * rate_corrected / depth

    def nitrate_uptake_floating_algae(
        self, nitrate: ArrayLike, ammonium: ArrayLike, algea_growth_rate: ArrayLike
    ) -> ArrayLike:
        """v1 NO3_ApGrowth: ApUptakeFr_NO3 * rna * ApGrowth (mg-N/L/d).

        Phase 9.A.2 audit finding N12: previously used the static
        ``float_algea_faction_uptake_from_nitrate`` (default 1.0) for the
        NO3 uptake fraction while the NH4 path read the dynamic
        ``algal_nh4_uptake_fraction`` recomputed each step. The two paths
        therefore did not sum to ``rna * ApGrowth``; algal-N mass balance
        was violated. Per v1 (processes.py:1675) and Fortran
        (modNitrogen.f90:321), the NO3 fraction is exactly
        ``1 - ApUptakeFr_NH4``, computed dynamically from current NH4/NO3.

        We now read ``floating_algae_process.algal_nh4_uptake_fraction``
        (recomputed in FloatingAlgae.run via _ap_uptake_fr_nh4) and use
        ``1 - that`` for the NO3 path. Total NH4 + NO3 algal uptake then
        sums to ``rna * algal_growth_rate`` exactly.
        """
        if not self.use_floating_algae:
            return 0.0

        # Dynamic NO3 uptake fraction = 1 - dynamic NH4 uptake fraction.
        algal_nh4_uptake_fraction = getattr(
            self.floating_algae_process, "algal_nh4_uptake_fraction", 0.5
        )
        algal_no3_uptake_fraction = 1.0 - algal_nh4_uptake_fraction

        # rna = AWn / AWa (mg-N/ug-Chla).
        rna = self.floating_algae_nitrogen_weight / self.algal_chlorophyll
        return rna * algea_growth_rate * algal_no3_uptake_fraction

    def nitrate_uptake_benthic_algae(
        self,
        nitrate: ArrayLike,
        ammonium: ArrayLike,
        algea_growth_rate: ArrayLike,
        depth: ArrayLike,
    ) -> ArrayLike:
        """v1 NO3_AbGrowth: AbUptakeFr_NO3 * rnb * Fb * AbGrowth / depth.

        Phase 9.A.2 audit finding N13: previous implementation had four
        structural defects:
        (1) divided by ``algal_chlorophyll`` (= AWa, the *floating*-algae
            chlorophyll factor of 1000) instead of by ``BWd`` (benthic
            dry-weight). Wrong stoichiometry by orders of magnitude.
        (2) omitted the ``/depth`` divisor (benthic processes are area-
            integrated; v1 and Fortran both convert g/m^2/d to mg-N/L/d
            via /depth).
        (3) multiplied by ``self.fraction_bottom_area`` (default 1.0)
            instead of ``Fb`` (benthic ``fraction of bottom area``,
            default 0.9).
        (4) used static ``benthic_algea_faction_uptake_from_nitrate``
            (default 0.5) instead of dynamic
            ``1 - balgae_nh4_uptake_fraction`` recomputed per step in
            BenthicAlgae.run.

        Corrected form per v1 (processes.py:1697) and Fortran
        (modNitrogen.f90:328):

            NO3_AbGrowth = AbUptakeFr_NO3 * rnb * Fb * AbGrowth / depth

        with ``rnb = BWn / BWd`` (mg-N/mg-D), ``Fb`` from
        BENTHIC_DEFAULTS, ``AbGrowth`` (g/m^2/d) cached on BenthicAlgae,
        and ``AbUptakeFr_NO3 = 1 - balgae_nh4_uptake_fraction`` dynamic.

        Mirrors the v3 Phosphorus benthic-uptake pattern
        (``Phosphorus._tip_uptake_benthic_algae``).
        """
        if not self.use_benthic_algae:
            return 0.0

        # Dynamic NO3 fraction = 1 - dynamic NH4 fraction (recomputed
        # per step in BenthicAlgae.run via _ab_uptake_fr_nh4).
        balgae_nh4_uptake_fraction = getattr(
            self.benthic_algae_process, "balgae_nh4_uptake_fraction", 0.5
        )
        balgae_no3_uptake_fraction = 1.0 - balgae_nh4_uptake_fraction

        # rnb = BWn / BWd (mg-N/mg-D); read directly from BenthicAlgae
        # for parity with the NH4 path (benthic_algae.py:508).
        bwn = getattr(self.benthic_algae_process, "BWn", self.benthic_algae_nitrogen_weight)
        bwd = getattr(self.benthic_algae_process, "BWd", 100.0)
        rnb = bwn / bwd

        # Fb from BenthicAlgae (default 0.9).
        fb = getattr(self.benthic_algae_process, "Fb", 0.9)

        return balgae_no3_uptake_fraction * rnb * fb * algea_growth_rate / depth

    def nitrification_inhibition(self, oxygen_dissolved: ArrayLike) -> ArrayLike:
        if not self.use_nitrate:
            return 1.0

        # Phase 9.A.2 audit finding N1: read NITROGEN_DEFAULTS attribute
        # ``KNR=0.6`` (v1/Fortran default). Legacy
        # ``nitrification_oxygen_inhibition_factor`` kwarg, when supplied,
        # syncs onto ``self.KNR`` in __init__.
        return 1.0 - np.exp(-self.KNR * oxygen_dissolved)

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

        ``vson / depth * OrgN``. ``vson_20`` is the raw settling velocity
        (m/d); per Phase 9.E correction, no Arrhenius temperature
        correction is applied -- this matches both v1
        (``processes.py:1333``) and Fortran (``modNitrogen.f90:233``)
        which both use raw ``vson``. The deliberate Fortran/v1 design
        distinction is that rate constants get Arrhenius corrections but
        settling velocities do not (settling depends on water viscosity
        with weak temperature dependence ~1.009, not on biochemical
        Arrhenius scaling). The ``temperature`` argument is retained for
        API stability with other rate methods on this Process; it is no
        longer consumed.
        """
        if not getattr(self, "use_OrgN", True):
            return 0.0
        return self.vson_20 / depth * organic_nitrogen

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
