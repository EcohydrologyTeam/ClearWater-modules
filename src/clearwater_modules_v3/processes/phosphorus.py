"""v3 NSM1 Phosphorus Process.

Phase 4 of the v3 NSM1 implementation plan
(``design/clearwater_modules_v3_nsm1_design_specification.md`` Section 11
Phase 4, Section 5 Phosphorus design notes). v3-native (NOT a v2
overlay). A single ``Phosphorus`` Process owns both phosphorus state
variables:

* ``tip``                 — Total Inorganic Phosphorus (mg-P/L)
* ``organic_phosphorus``  — Organic Phosphorus (mg-P/L)

Kinetics (mirrors v1 ``processes.py:1833-2168`` and v3 design spec
Section 5):

    dTIP/dt = + kop_tc * OrgP                              (hydrolysis from OrgP)
              - vs / depth * (1 - fdp) * TIP               (settling of particulate fraction)
              + DIPfromBed                                 (sediment release; rpo4_tc / depth)
              - rpa * algal_growth_rate                    (algal P uptake; floating)
              + rpa * algal_respiration_rate               (algal P release; floating)
              - rpb * Fb * algal_growth_rate / depth       (algal P uptake; benthic)
              + rpb * Fb * algal_respiration_rate / depth  (algal P release; benthic)

    dOrgP/dt = + algal_orgp_from_mortality_rate            (floating algae mortality)
               + balgae_orgp_from_mortality_rate           (benthic algae mortality)
               - kop_tc * OrgP                             (hydrolysis to TIP)
               - vsop / depth * OrgP                       (settling)

Where:

* ``kop_tc = arrhenius_correction(kop_20, kop_theta, T)`` — OrgP -> TIP
  hydrolysis rate (1/d).
* ``rpo4_tc = arrhenius_correction(rpo4_20, rpo4_theta, T)`` —
  temperature-corrected sediment P release rate (g-P/m^2/d). The default
  ``rpo4_20=0`` keeps this term silently zero (Phase 0 audit
  ``FIXME(phase1-audit)`` in ``parameters/phosphorus.py``).
* ``fdp = fdp(use_TIP, Solid, kdpo4)`` from
  ``clearwater_modules_v3.utils.partitioning``. With ``kdpo4=0`` (the
  v3 default; FIXME(phase1-audit)), ``fdp = 1.0`` which means TIP is
  fully dissolved and ``(1 - fdp) * TIP = 0``: TIP settling is silently
  disabled at default parameters. Setting ``kdpo4 > 0`` activates the
  particulate fraction.
* ``rpa`` (mg-P/ug-Chla) — floating algae P:Chla stoichiometric ratio.
  Per v1 ``processes.py:308-348`` and v2 floating_algae.py:321,
  ``rpa = AWp / AWa``. Read from ``self.floating_algae_process``.
* ``rpb`` (mg-P/mg-D) — benthic algae P:dry-weight ratio,
  ``rpb = BWp / BWd``. Read from ``self.benthic_algae_process``.

Algal coupling reads (Q10 GS-rates contract; resolved Q10 in spec
Section 14):

* ``floating_algae_process.algal_growth_rate``       — ug-Chla/L/d
* ``floating_algae_process.algal_respiration_rate``  — ug-Chla/L/d
* ``floating_algae_process.algal_orgp_from_mortality_rate``  — mg-P/L/d
* ``benthic_algae_process.balgae_growth_rate``       — g/m^2/d (the
  cached attribute is the per-area benthic biomass rate; the depth
  and Fb area-integration is applied here, mirroring v1
  ``DIP_AbGrowth`` / ``DIP_AbRespiration``).
* ``benthic_algae_process.balgae_respiration_rate``  — g/m^2/d
* ``benthic_algae_process.balgae_orgp_from_mortality_rate``  — mg-P/L/d
  (the cached attribute already includes the ``Fw * Fb / depth`` factor
  per v2 benthic_algae.py:236, so this Phosphorus Process consumes it
  directly without additional area integration).

Forward Euler integrator (resolved Q1 / Q4):

    state_new = state_old + rate * dt_days

with ``dt_days = time_step.total_seconds() / 86400``. Both states are
clipped via ``clip_negative_state`` (resolved Q7 clip-with-log) and
persisted via ``registry.set_at_time``.

Out of scope for Phase 4: ProcessFactory.register / Model integration
(Phase 5.5), DIP partitioning audit (NSM2 territory), full SedFlux
sediment-P budget (NSM2).
"""

from __future__ import annotations

from datetime import datetime, timedelta
import logging
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from clearwater_data.custom_types import ArrayLike
from clearwater_data.variables import VariableRegistry

from clearwater_modules_v3.processes.base import Process, ProcessFactory
from clearwater_modules_v3.utils.conversions import arrhenius_correction
from clearwater_modules_v3.utils.numerics import (
    Diagnostics,
    clip_negative_state,
    sanitize_rate,
)
from clearwater_modules_v3.utils.partitioning import fdp as fdp_partition

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


logger = logging.getLogger(__name__)


# Inline fallback defaults for non-phosphorus parameters consumed by the
# fdp partitioning calculation. The user can override via ``parameters``;
# otherwise these match the v3 ``global_parameters`` / ``global_vars``
# DEFAULTS. Kept inline to avoid coupling Phosphorus to those modules at
# construction time, mirroring the Pathogen pattern.
_PARTITIONING_DEFAULTS: dict[str, float | bool] = {
    "use_TIP": True,    # bool; total inorganic phosphorus module switch
    "Solid": 1.0,       # mg/L; suspended solids concentration
    "use_OrgP": True,   # bool; organic phosphorus module switch
    # Benthic algae geometry (used in DIP_AbGrowth / DIP_AbRespiration
    # depth+Fb area integration). v3 ``balgae`` DEFAULTS values.
    "Fb": 0.9,          # unitless; fraction of bottom area available
}


class Phosphorus(Process):
    # Phase H-9 (2026-05-21): Phosphorus reads step-scoped rate caches
    # from FloatingAlgae (algal_growth_rate, algal_respiration_rate) and
    # BenthicAlgae (balgae_growth_rate, balgae_respiration_rate).
    upstream_processes = ("FloatingAlgae", "BenthicAlgae")

    """v3 NSM1 Phosphorus Process (TIP + OrgP).

    Single Process managing both Phosphorus state variables. Forward
    Euler integrator advances each state by ``rate * dt_days`` per
    substep. Algal coupling is read from FloatingAlgae and BenthicAlgae
    processes through the step-scoped rate cache (Q10 GS-rates
    workaround) when those Processes are present in the Model.

    State variable registry keys (matching the Tier 1 conftest
    fixture and ``total_p`` helper):

    * ``tip`` — Total Inorganic Phosphorus (mg-P/L)
    * ``organic_phosphorus`` — Organic Phosphorus (mg-P/L)

    Phase 5.5 integration note: v2 algae processes register the
    inorganic phosphorus pool under ``phosphorus_total_inorganic``
    (see ``v2/processes/floating_algae.py:244``). The Phase 4
    Phosphorus Process and the Tier 1 fixture both use ``tip``;
    Phase 5 integration should reconcile these two names by
    standardizing on one (the v3 convention is ``tip``).
    """

    variables = [
        "tip",
        "organic_phosphorus",
        "water_temperature",
        "depth",
    ]

    # Class-level v3 defaults. Lazy-loaded on first instantiation to
    # mirror the v2 nitrogen.py:91-93 pattern (the lazy idiom is also
    # used in v3 Pathogen for consistency).
    DEFAULTS: dict[str, float | int | bool] = {}

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface Phosphorus exposes via the opportunistic-write loop in
    # ``run``. Each name maps to a ``self.<name>`` cache attribute set
    # inside ``_change_with_components`` and matches the inventory in
    # ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3.
    #
    # Note: ``orgp_hydrolysis_rate`` is the Appendix A name (matches the
    # NSM1 1.0.0 registry-coupling convention). It aliases the existing
    # ``orgp_to_tip_hydrolysis_rate`` cache attribute that
    # ``test_phosphorus_v1_parity_v3.py`` reads via ``getattr``; both
    # names point at the same value (set as side effects in
    # ``_change_with_components``) for back-compat.
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "orgp_hydrolysis_rate",
        "orgp_settling_rate",
        "tip_settling_rate",
        "dip_from_bed",
        "orgp_algal_mortality_rate",
        "tip_algal_growth_rate",
        "tip_balgae_growth_rate",
    )

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
    ) -> None:
        """Initialize the Phosphorus Process.

        Args:
            parameters: Optional dict of parameter overrides. Merged
                with the class-level ``DEFAULTS`` (v3
                ``PHOSPHORUS_DEFAULTS``). Unknown keys log a warning
                and are ignored. Recognized phosphorus keys include
                ``kop_20``, ``kop_theta``, ``rpo4_20``, ``rpo4_theta``,
                ``kdpo4``, ``vsop``, ``vs``. Recognized partitioning /
                module-switch keys include ``use_TIP``, ``use_OrgP``,
                ``Solid``, ``Fb``.
            time_step: Substep cadence for this Process. Forward Euler
                converts to days as ``time_step.total_seconds() / 86400``.
        """
        Process.__init__(self, time_step)

        # --- v3-style parameter merge (DEFAULTS + user overrides) ---
        # Lazy-load PHOSPHORUS_DEFAULTS (mirror v2 nitrogen.py:91-93).
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.phosphorus import (
                DEFAULTS as PHOSPHORUS_DEFAULTS,
            )
            type(self).DEFAULTS = PHOSPHORUS_DEFAULTS

        user_params = parameters or {}

        # Phase 9.F.C defensive guard: v3 1.0.0 Phosphorus.run() gates
        # ``dip_from_bed`` by ``use_TIP`` only, NOT by ``use_SedFlux``.
        # The bed-flux term is silenced solely by the ``rpo4_20 = 0``
        # default. If a user opts into ``use_SedFlux=True`` they are
        # signaling intent to enable the full sediment-flux feature,
        # which requires the NSM2 diagenesis path that is not
        # implemented in v3 1.0.0. Refuse explicitly rather than
        # silently producing partial behavior. See corrections doc
        # Section 2.1.
        if user_params.get("use_SedFlux", False):
            raise NotImplementedError(
                "Phosphorus: use_SedFlux=True is not implemented in v3 1.0.0. "
                "The full sediment-flux feature requires the NSM2 diagenesis "
                "path. Set rpo4_20 directly to specify a constant sediment "
                "release rate for site-specific calibration (without "
                "use_SedFlux). See parameter_defaults_corrections.md "
                "Section 2.1."
            )

        # Allowed key universe = phosphorus DEFAULTS U partitioning fallbacks.
        allowed_keys = set(self.DEFAULTS) | set(_PARTITIONING_DEFAULTS) | {"use_SedFlux"}
        unknown_keys = set(user_params) - allowed_keys
        for key in sorted(unknown_keys):
            logger.warning(
                "Phosphorus: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in PHOSPHORUS_DEFAULTS or partitioning "
                "fallbacks).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        merged.pop("use_SedFlux", None)
        for k, v in merged.items():
            setattr(self, k, v)

        # Inline partitioning / module-switch defaults from non-phosphorus
        # parameter groups. Allow user override via ``parameters``.
        for k, v in _PARTITIONING_DEFAULTS.items():
            if not hasattr(self, k):
                setattr(self, k, user_params.get(k, v))

        # Diagnostics handle: a v3 Model will replace this with its
        # run-level Diagnostics in ``init_process``. The local default
        # keeps ``run`` callable in unit tests that don't wire a Model.
        self.diagnostics: Diagnostics = Diagnostics()

        # Algae coupling flags (set in ``init_process`` when a Model is
        # available; default to False so Tier 1 standalone runs work).
        self.use_floating_algae: bool = False
        self.use_benthic_algae: bool = False
        self.floating_algae_process = None
        self.benthic_algae_process = None

        # Step-scoped rate-variable cache (Q10 GS-rates contract).
        # Phase 5+ DOX / sediment may consume these; cached for parity
        # with the CBOD pattern.
        self.tip_settling_rate: ArrayLike = 0.0
        self.orgp_settling_rate: ArrayLike = 0.0
        self.orgp_to_tip_hydrolysis_rate: ArrayLike = 0.0

    @ProcessFactory.register("phosphorus")
    @staticmethod
    def from_config(
        config: dict, variable_registry: VariableRegistry
    ) -> "Phosphorus":
        return Phosphorus(**config)

    def init_process(
        self, model: "Model", registry: VariableRegistry
    ) -> None:
        """Wire Phosphorus to the v3 Model.

        Captures algae-process handles so the kinetic terms can read
        ``algal_growth_rate``, ``algal_respiration_rate``, and the
        mortality routing rates from the FloatingAlgae and BenthicAlgae
        Processes (resolved Q10 GS-rates contract). Also captures the
        Model's run-level ``Diagnostics`` handle.
        """
        self.diagnostics = getattr(model, "diagnostics", None) or self.diagnostics
        self.use_floating_algae = model.has_process("FloatingAlgae")
        self.use_benthic_algae = model.has_process("BenthicAlgae")
        if self.use_floating_algae:
            self.floating_algae_process = model.get_process("FloatingAlgae")
        if self.use_benthic_algae:
            self.benthic_algae_process = model.get_process("BenthicAlgae")

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Advance TIP and OrgP by one substep using Forward Euler.

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_change_with_components``
        (B); applies Forward Euler with unconditional clip-with-log
        (C, D); persists primary outputs (E); caches step-scoped rates
        on ``self.<name>`` (F); opportunistically writes diagnostics
        (G).

        """
        # --- State and forcing reads (pattern A) ---
        tip = registry.get_at_time("tip", time)
        orgp = registry.get_at_time("organic_phosphorus", time)
        water_temperature = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)
        # Shared suspended-solids input (clearwater_modules_v3_solid_input_source):
        # prefer the canonical per-cell ``Solid`` registry forcing when a provider
        # registered it; otherwise fall back to the constructor scalar
        # ``self.Solid``. Mirrors Temperature's optional wind_shelter_coefficient
        # read. Runs that register no ``Solid`` are byte-identical to before.
        solid = (
            registry.get_at_time("Solid", time)
            if "Solid" in registry
            else self.Solid
        )

        # --- Fused rate composition (pattern B) ---
        dtip_dt, dorgp_dt, components = self._change_with_components(
            tip=tip,
            orgp=orgp,
            water_temperature=water_temperature,
            depth=depth,
            solid=solid,
        )

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        # The first three names (orgp_settling_rate, tip_settling_rate,
        # orgp_to_tip_hydrolysis_rate / orgp_hydrolysis_rate) were
        # already populated as side effects of _change_with_components;
        # the setattr loop is idempotent on those. The remaining four
        # are new caches with this phase.
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])

        # --- Forward Euler in days (pattern C) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        tip_new = tip + dtip_dt * dt_days
        orgp_new = orgp + dorgp_dt * dt_days

        # --- Clip-with-log per the Q7 contract (pattern D) ---
        # Step attribution is automatic via ``diagnostics.current_step``
        # (Phase 0.6 Q1).
        tip_new = clip_negative_state(tip_new, "tip", self.diagnostics)
        orgp_new = clip_negative_state(
            orgp_new, "organic_phosphorus", self.diagnostics
        )

        # --- Persist primary outputs (pattern E) ---
        registry.set_at_time("tip", time, tip_new)
        registry.set_at_time("organic_phosphorus", time, orgp_new)

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
        tip: ArrayLike,
        orgp: ArrayLike,
        water_temperature: ArrayLike,
        depth: ArrayLike,
        solid: ArrayLike | None = None,
    ) -> tuple[ArrayLike, ArrayLike, dict]:
        """Compute ``(dtip_dt, dorgp_dt, components)`` for Phosphorus.

        Per-day rates (mg-P/L/d). Code-motion-only refactor of ``run``'s
        former inline composition (§11.6): operand order, intermediate
        names, kinetic-helper calls, and the use_TIP / use_OrgP gating
        are preserved verbatim.


        Side effect: sets ``self.orgp_to_tip_hydrolysis_rate``,
        ``self.tip_settling_rate``, ``self.orgp_settling_rate``
        (preserved attribute names that the existing v3 code and
        ``test_phosphorus_v1_parity_v3.py`` depend on). The pattern F
        ``setattr`` loop in ``run`` is idempotent on these names.
        """
        # Temperature-corrected rates.
        kop_tc = arrhenius_correction(
            water_temperature, self.kop_20, self.kop_theta
        )
        rpo4_tc = arrhenius_correction(
            water_temperature, self.rpo4_20, self.rpo4_theta
        )

        # Dissolved fraction. With ``kdpo4=0`` (v3 default), fdp = 1.0
        # which means (1 - fdp) * TIP = 0 and TIP settling is silently
        # zero. Setting ``kdpo4 > 0`` activates the particulate fraction.
        # ``solid`` is resolved registry-first in ``run`` (Solid input source
        # spec); ``None`` only when a direct caller omitted it -> scalar default.
        if solid is None:
            solid = self.Solid
        fdp = fdp_partition(
            use_TIP=self.use_TIP,
            Solid=solid,
            kdpo4=self.kdpo4,
        )

        # OrgP -> TIP hydrolysis (mg-P/L/d). v1 ``OrgP_DIP_decay``.
        if self.use_OrgP:
            orgp_to_tip = kop_tc * orgp
        else:
            orgp_to_tip = 0.0
        self.orgp_to_tip_hydrolysis_rate = orgp_to_tip

        # TIP settling (mg-P/L/d). v1 ``TIP_Settling``.
        if self.use_TIP:
            tip_settling = self.vs / depth * (1.0 - fdp) * tip
        else:
            tip_settling = 0.0
        self.tip_settling_rate = tip_settling

        # OrgP settling (mg-P/L/d). v1 ``OrgP_Settling``.
        if self.use_OrgP:
            orgp_settling = self.vsop / depth * orgp
        else:
            orgp_settling = 0.0
        self.orgp_settling_rate = orgp_settling

        # Sediment P release (mg-P/L/d). v1 ``DIPfromBed``.
        # rpo4_tc has units g-P/m^2/d; dividing by depth (m) gives
        # mg-P/L/d (1 g/m^3 == 1 mg/L). At the v3 default ``rpo4_20=0``
        # this term is silently zero.
        if self.use_TIP:
            dip_from_bed = rpo4_tc / depth
        else:
            dip_from_bed = 0.0

        # Algae coupling. Read step-scoped rate-variable cache from
        # FloatingAlgae / BenthicAlgae per the Q10 GS-rates contract.
        ap_tip_uptake = self._tip_uptake_floating_algae()
        ap_tip_release = self._tip_release_floating_algae_respiration()
        ab_tip_uptake = self._tip_uptake_benthic_algae(depth=depth)
        ab_tip_release = self._tip_release_benthic_algae_respiration(depth=depth)

        ap_orgp_mortality = self._orgp_from_floating_algae_mortality()
        ab_orgp_mortality = self._orgp_from_benthic_algae_mortality()

        # --- TIP rate of change (mg-P/L/d) ---
        # v1 ``dTIPdt`` (line 2058-2091).
        if self.use_TIP:
            dtip_dt = (
                orgp_to_tip
                - tip_settling
                + dip_from_bed
                - ap_tip_uptake
                + ap_tip_release
                - ab_tip_uptake
                + ab_tip_release
            )
        else:
            dtip_dt = xr.zeros_like(tip) if hasattr(tip, "dims") else 0.0

        # --- OrgP rate of change (mg-P/L/d) ---
        # v1 ``dOrgPdt`` (line 1937-1956).
        if self.use_OrgP:
            dorgp_dt = (
                ap_orgp_mortality
                + ab_orgp_mortality
                - orgp_to_tip
                - orgp_settling
            )
        else:
            dorgp_dt = xr.zeros_like(orgp) if hasattr(orgp, "dims") else 0.0

        # NaN guard (mirrors v3 Nitrogen / FloatingAlgae pattern).
        dtip_dt = sanitize_rate(dtip_dt)
        dorgp_dt = sanitize_rate(dorgp_dt)

        # --- Components dict ---
        # ``orgp_hydrolysis_rate`` is the Appendix A name; aliases the
        # existing ``orgp_to_tip_hydrolysis_rate`` attribute set above.
        # ``tip_algal_growth_rate`` and ``tip_balgae_growth_rate``
        # split the algal TIP uptake by source per Appendix A.
        # ``orgp_algal_mortality_rate`` sums floating + benthic
        # mortality contributions to OrgP (mirrors the Nitrogen
        # ``nh4_algal_growth_rate`` total convention).
        components = {
            "orgp_hydrolysis_rate": orgp_to_tip,
            "orgp_settling_rate": orgp_settling,
            "tip_settling_rate": tip_settling,
            "dip_from_bed": dip_from_bed,
            "orgp_algal_mortality_rate": ap_orgp_mortality + ab_orgp_mortality,
            "tip_algal_growth_rate": ap_tip_uptake,
            "tip_balgae_growth_rate": ab_tip_uptake,
        }

        return dtip_dt, dorgp_dt, components
    # ------------------------------------------------------------------
    # Algal-coupling helpers
    # ------------------------------------------------------------------

    def _rpa(self) -> float:
        """Floating-algae P:Chla stoichiometric ratio (mg-P/ug-Chla).

        Per v1 ``processes.py:308-348`` and v2 ``floating_algae.py:321``,
        ``rpa = AWp / AWa``.
        """
        if self.floating_algae_process is None:
            return 0.0
        return self.floating_algae_process.AWp / self.floating_algae_process.AWa

    def _rpb(self) -> float:
        """Benthic-algae P:dry-weight ratio (mg-P/mg-D).

        Per v2 ``benthic_algae.py:228``, ``rpb = BWp / BWd``.
        """
        if self.benthic_algae_process is None:
            return 0.0
        return self.benthic_algae_process.BWp / self.benthic_algae_process.BWd

    def _tip_uptake_floating_algae(self) -> ArrayLike:
        """TIP -> floating-algae uptake (mg-P/L/d).

        v1 ``DIP_ApGrowth = rpa * ApGrowth``.
        """
        if not self.use_floating_algae:
            return 0.0
        algal_growth = getattr(
            self.floating_algae_process, "algal_growth_rate", 0.0
        )
        return self._rpa() * algal_growth

    def _tip_release_floating_algae_respiration(self) -> ArrayLike:
        """Floating-algae respiration -> TIP (mg-P/L/d).

        v1 ``DIP_ApRespiration = rpa * ApRespiration``.
        """
        if not self.use_floating_algae:
            return 0.0
        algal_resp = getattr(
            self.floating_algae_process, "algal_respiration_rate", 0.0
        )
        return self._rpa() * algal_resp

    def _tip_uptake_benthic_algae(self, depth: ArrayLike) -> ArrayLike:
        """TIP -> benthic-algae uptake (mg-P/L/d).

        v1 ``DIP_AbGrowth = rpb * Fb * AbGrowth / depth``.

        ``balgae_growth_rate`` is the cached per-area benthic biomass
        growth (g/m^2/d) from the BenthicAlgae Process; this method
        applies the depth+Fb area-integration factor.
        """
        if not self.use_benthic_algae:
            return 0.0
        balgae_growth = getattr(
            self.benthic_algae_process, "balgae_growth_rate", 0.0
        )
        return self._rpb() * self.Fb * balgae_growth / depth

    def _tip_release_benthic_algae_respiration(
        self, depth: ArrayLike
    ) -> ArrayLike:
        """Benthic-algae respiration -> TIP (mg-P/L/d).

        v1 ``DIP_AbRespiration = rpb * Fb * AbRespiration / depth``.
        """
        if not self.use_benthic_algae:
            return 0.0
        balgae_resp = getattr(
            self.benthic_algae_process, "balgae_respiration_rate", 0.0
        )
        return self._rpb() * self.Fb * balgae_resp / depth

    def _orgp_from_floating_algae_mortality(self) -> ArrayLike:
        """Floating-algae mortality -> OrgP (mg-P/L/d).

        Reads the cached ``algal_orgp_from_mortality_rate`` populated by
        FloatingAlgae.run (Phase 2.A; v2 floating_algae.py:326).
        """
        if not self.use_floating_algae:
            return 0.0
        return getattr(
            self.floating_algae_process,
            "algal_orgp_from_mortality_rate",
            0.0,
        )

    def _orgp_from_benthic_algae_mortality(self) -> ArrayLike:
        """Benthic-algae mortality -> OrgP (mg-P/L/d).

        Reads the cached ``balgae_orgp_from_mortality_rate`` populated by
        BenthicAlgae.run (v2 benthic_algae.py:236). The cached value
        already includes the ``Fw * Fb / depth`` factor, so this
        Phosphorus Process consumes it directly.
        """
        if not self.use_benthic_algae:
            return 0.0
        return getattr(
            self.benthic_algae_process,
            "balgae_orgp_from_mortality_rate",
            0.0,
        )

