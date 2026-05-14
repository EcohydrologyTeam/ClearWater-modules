"""v3 NSM1 Particulate Organic Matter (POM) Process.

Phase 3.2 of the v3 NSM1 implementation plan (design spec Section 11
Phase 3, Section 4 component inventory). v3-native (no v2 overlay
exists; v2 has no POM Process).

Kinetics follow v1 ``clearwater_modules/nsm1/processes.py`` lines
2185-2329 (``kpom_tc``, ``POM_algal_settling``, ``POM_dissolution``,
``POM_POC_settling``, ``POM_benthic_algae_mortality``, ``POM_burial``,
``dPOMdt``, ``POM``):

    dPOM/dt =   POM_algal_settling             (floating algae source)
              - POM_dissolution                (POM -> DOC sink)
              + POM_POC_settling               (POC settling source)
              + POM_benthic_algae_mortality    (benthic algae source)
              - POM_burial                     (sediment loss)

Where (per v1):

    POM_algal_settling          = vsap * Ap * rda / h2          ; rda = AWd/AWa
    POM_dissolution             = kpom_tc * POM
    POM_POC_settling            = vsoc * POC / h2 / fcom
    POM_benthic_algae_mortality = Ab * kdb_tc * Fb * (1 - Fw) / h2
    POM_burial                  = vb * POM / h2

POM is normalized by ``h2`` (active sediment layer thickness, m), not by
water depth. The Phase 3.2 task brief mentioned a ``vsom_20`` settling
velocity divided by ``depth`` -- v1 has no such parameter. The v1 burial
term ``vb / h2 * POM`` is the closest analogue and is used here. See the
"v1/v2 ambiguities" note in the Phase 3.2 report.

Conceptual note (Phase 9.F.C documentation fix, corrections doc
Section 2.5): NSM1's POM state variable represents the **bed-sediment
POM compartment** (Fortran ``POM2`` -- the "2" suffix denotes Di Toro's
"layer 2", i.e. the lower anaerobic sediment layer), NOT the
water-column POM. v1 and v3 dropped the ``2`` subscript when porting
from Fortran, but the conceptual identity is preserved via the
``h2`` divisor. ``h2 = 0.1`` m matches the Di Toro (2001) / QUAL2K
v2.11 §5.6 convention for the lower anaerobic layer thickness ``H_2``
(approx 10 cm), and ``h2`` plays the role of converting areal
water-column fluxes (m * mg/L/d) into bed volumetric concentration
changes (mg/L/d). The water-column algal/POC/benthic-mortality source
terms in the equations above carry the ``/ h2`` factor for exactly
this dimensional reason. Implementing the full two-layer Di Toro
diagenesis (separate ``H_1`` aerobic layer, full nutrient flux model)
is the future NSM2 sediment-diagenesis scope; v3 NSM1 1.0.0 carries
only the ``H_2`` layer with first-order burial/dissolution kinetics,
matching v1 and Fortran exactly.

Forward-Euler integrator pattern (Phase 2.A / 2.B): rates are 1/d,
``dt_days = time_step.total_seconds() / 86400``, ``state_new = state +
rate * dt_days``, then ``clip_negative_state`` with diagnostics, then
``set_at_time``.

Phase 3.5 inter-process coupling: the Phase 2.A FloatingAlgae /
BenthicAlgae Processes now cache:

* ``algal_pom_from_settling_rate``    -- v1 ``vsap * Ap * rda / h2``
* ``balgae_pom_from_mortality_rate``  -- v1 ``kdb_tc * Ab * Fb * (1 - Fw)
                                          / h2``

so POM.run reads them via ``getattr(...)`` (consumer-ready fluxes,
already in mg/L/d). POM also caches its own ``pom_hydrolysis_rate``
for the Phase 5.A Carbon consumer to read as a DOC source term. The
cache is in **mg-C/L/d (water-column volumetric)**, computed as
``fcom * kpom_tc * pom * h2 / depth``: the raw POM-mass dissolution
rate (``kpom_tc * pom``, mg-D/L_sed/d) is multiplied by ``fcom`` (mg-C
per mg-D, the carbon mass fraction of organic dry weight) and by
``h2 / depth`` (sediment-volume to water-column-volume conversion).
This matches the symmetric ``poc_hydrolysis_rate`` cache on Carbon
(mg-C/L/d) and ensures the POM dissolution pathway is mass-conserving
in carbon terms. Without these factors, every mg of dry POM dissolved
would produce 1 mg of DOC at default depth=1 m / h2=0.1 m / fcom=0.4
-- a 25x overcount that fails closed-system C conservation.
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

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from model import Model


logger = logging.getLogger(__name__)


# Inline fallback defaults for parameters that live in v3 ``global_vars``
# / ``global_parameters`` rather than in ``pom`` itself. Used only when
# the user does not pass them in the ``parameters`` dict; Phase 3.5
# integration may wire these via the model-level parameter registry.
_POM_GLOBAL_DEFAULTS = {
    "vsoc": 0.01,    # m/d; POC settling velocity (v3 global_vars)
    "fcom": 0.4,     # unitless; fraction of sediment as combustible OM (v3 global_vars)
    "vb": 6.85e-6,   # m/d; sediment burial velocity (= 0.0025 m/yr = 0.25 cm/yr; matches v3 global_vars after Phase 9.F.A correction)
    "use_Algae": True,
    "use_Balgae": True,
    "use_POC": True,
}


class POM(Process):
    """v3 NSM1 Particulate Organic Matter Process.

    State variable: ``pom`` (mg/L).

    Sources / sinks (per v1 NSM1):
      * Floating algae settling -> POM (source)
      * POC settling            -> POM (source)
      * Benthic algae mortality -> POM (source, fraction routed to water)
      * POM dissolution         -> DOC (sink)
      * POM burial              -> sediment (sink)
    """

    variables = ["pom", "water_temperature", "depth"]

    # Class-level v3 defaults (Section 3.4 of design spec). Lazy-loaded
    # on first instantiation to be consistent with the lazy pattern used
    # by Nitrogen / FloatingAlgae (originally introduced to break the
    # v2 <-> v3 circular import). POM is v3-native, but the lazy pattern
    # keeps consistency with the rest of v3.
    DEFAULTS: dict[str, float | int | bool] = {}

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface POM exposes via the opportunistic-write loop in ``run``.
    # Each name maps to a ``self.<name>`` cache attribute set inside
    # ``_change_with_components`` and matches the inventory in
    # ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3.
    #
    # Note: ``pom_doc_source_rate`` (consumer-ready unit-converted DOC
    # source rate that Carbon reads via getattr) is NOT in
    # REGISTRY_DIAGNOSTICS — it is set as a side effect of the helper
    # for sibling consumption but not exposed to the registry; the
    # raw POM dissolution rate is exposed instead under
    # ``pom_hydrolysis_rate``.
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "pom_hydrolysis_rate",
        "pom_settling_rate",
        "pom_algal_mortality_rate",
        "pom_balgae_mortality_rate",
    )

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
    ) -> None:
        """Initialize the POM process.

        Args:
            parameters: Optional dict of v3 POM parameter overrides.
                Merged with the class-level ``DEFAULTS`` (v3
                ``POM_DEFAULTS``). Unknown keys (relative to the union of
                POM_DEFAULTS and the inline coupling defaults) are warned
                and ignored.
            time_step: Substep cadence for this Process.
        """
        Process.__init__(self, time_step)

        # --- v3-style parameter merge (DEFAULTS + user overrides) ---
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.pom import DEFAULTS as POM_DEFAULTS
            type(self).DEFAULTS = POM_DEFAULTS

        user_params = parameters or {}
        # Allowed user keys: POM_DEFAULTS keys plus the inline coupling
        # defaults (vsoc, fcom, vb, use_*) that POM reads from the
        # global parameter groups.
        known_keys = set(self.DEFAULTS) | set(_POM_GLOBAL_DEFAULTS)
        unknown_keys = set(user_params) - known_keys
        for key in sorted(unknown_keys):
            logger.warning(
                "POM: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in POM_DEFAULTS or POM coupling defaults).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Inline coupling defaults that come from non-POM v3 groups.
        for k, v in _POM_GLOBAL_DEFAULTS.items():
            if not hasattr(self, k):
                setattr(self, k, user_params.get(k, v))

        # Diagnostics handle: defaults to a local Diagnostics so
        # ``clip_negative_state`` has somewhere to record when the
        # Process is driven directly from a unit test (no Model).
        from clearwater_modules_v3.utils.numerics import Diagnostics
        self.diagnostics = Diagnostics()

        # Coupling flags filled in by ``init_process``; default False so
        # ``run`` works in unit-test mode without a Model.
        self.use_floating_algae: bool = False
        self.use_benthic_algae: bool = False
        self.floating_algae_process = None
        self.benthic_algae_process = None

        # Phase 3.5 inter-process coupling: step-scoped rate cache
        # (Q10 GS-rates contract). Carbon consumes this as a DOC source
        # term (POM hydrolysis -> DOC). Units: **mg-C/L/d** (water-column
        # volumetric), with ``fcom`` and ``h2/depth`` already applied.
        # Initialized to 0 so a Carbon consumer that reads this attribute
        # before POM.run is ever called gets a defined value.
        self.pom_doc_source_rate: ArrayLike = 0.0

    @ProcessFactory.register("pom")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "POM":
        return POM(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """Capture cross-process refs and the run-level Diagnostics handle."""
        self.use_floating_algae = model.has_process("FloatingAlgae")
        self.use_benthic_algae = model.has_process("BenthicAlgae")

        if self.use_floating_algae:
            self.floating_algae_process = model.get_process("FloatingAlgae")
        if self.use_benthic_algae:
            self.benthic_algae_process = model.get_process("BenthicAlgae")

        # Capture the v3 Model's run-level Diagnostics handle so
        # ``clip_negative_state`` records clip events on the canonical
        # diagnostics. v2's Model has no ``diagnostics`` attribute; in
        # that case the locally-instantiated Diagnostics is retained.
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Run the POM process for one time step.

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_change_with_components``
        (B); applies Forward Euler with unconditional clip-with-log
        (C, D); persists primary output (E); caches step-scoped rates
        on ``self.<name>`` (F); opportunistically writes diagnostics
        (G).

        See ``_change_legacy_inline`` for the pre-Phase-7 inline
        composition (retained through Phase 10 for the helper-vs-inline
        parity test under §11.3). Phase 7 also resolves the Phase 1
        deferred item: the ``pom_doc_source_rate`` cache (consumed by
        Carbon via getattr) now lives in ``_change_with_components``
        rather than inside the public ``rate()`` method, so the cache
        does not depend on a side effect of an externally-callable
        helper.
        """
        # --- State and forcing reads (pattern A) ---
        pom = registry.get_at_time("pom", time)
        water_temperature = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)

        # POC source (only present if POC Process / state is registered).
        # In closed-system Tier 1 mode without a Carbon Process the POC
        # state may still be present in the conftest fixture; reading it
        # is a no-op as long as ``use_POC`` is set appropriately by the
        # Phase 3.5 wiring. To keep Phase 3.2 mass-conserving in the
        # closed-system test (vsoc=0), we read POC defensively.
        if "poc" in registry:
            poc = registry.get_at_time("poc", time)
        else:
            poc = xr.zeros_like(pom)

        # --- Fused rate composition (pattern B) ---
        rate, components = self._change_with_components(
            pom=pom,
            water_temperature=water_temperature,
            poc=poc,
            depth=depth,
        )

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])

        # --- Forward Euler in days (pattern C) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        pom_new = pom + rate * dt_days

        # --- Clip-with-log per Q7 (pattern D) ---
        # Step attribution is automatic via ``diagnostics.current_step``
        # (Phase 0.6 Q1).
        pom_new = clip_negative_state(pom_new, "pom", self.diagnostics)

        # --- Persist primary output (pattern E) ---
        registry.set_at_time("pom", time, pom_new)

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
        pom: ArrayLike,
        water_temperature: ArrayLike,
        poc: ArrayLike,
        depth: ArrayLike,
    ) -> tuple[ArrayLike, dict]:
        """Compute ``(rate, components)`` for POM.

        ``rate`` is the net per-day POM rate of change (mg/L/d).
        ``components`` is the dict[str, ArrayLike] indexed by
        ``REGISTRY_DIAGNOSTICS``.

        Code-motion-only refactor of ``run``'s former inline
        composition (§11.6): operand order, intermediate names, and
        sub-flux helpers are preserved verbatim from the body of
        ``rate()`` (the pre-Phase-7 entry point).

        **Phase 7 cache relocation**: the ``self.pom_doc_source_rate``
        cache that Carbon reads via getattr is now set here as a side
        effect (not inside the public ``rate()`` method). The ``rate()``
        method retains its prior behaviour for back-compat with
        external tests that call ``pom.rate(...)`` directly, but the
        canonical write path is via this helper.

        The companion shadow ``_change_legacy_inline`` returns just
        the net rate and is used by
        ``tests/v3/nsm1/test_pom_helper_vs_inline.py`` to verify this
        helper produces a bit-identical net rate through Phase 10.
        """
        # Temperature-corrected dissolution rate.
        kpom_tc = arrhenius_correction(
            water_temperature, self.kpom_20, self.kpom_theta
        )

        # Sink: dissolution. POM mass loss is in mg-D/L_sed/d (POM state
        # is dry-mass-per-sediment-volume). Used as a sink in dPOM/dt.
        rate_dissolution = kpom_tc * pom

        # Phase 5.A consumer-ready DOC source rate (mg-C/L_water/d).
        # The raw POM dissolution rate (mg-D/L_sed/d) becomes a DOC
        # source by (1) converting dry mass to carbon mass via
        # ``fcom`` (mg-C/mg-D) and (2) converting sediment volumetric
        # to water-column volumetric via ``h2 / depth``. Matches the
        # convention of ``poc_hydrolysis_rate`` on Carbon (mg-C/L/d
        # water-column). Required for closed-system C conservation.
        # Phase 7: cache moved here from ``rate()`` per spec §6 Phase 7.
        self.pom_doc_source_rate = (
            self.fcom * rate_dissolution * self.h2 / depth
        )

        # Sink: burial.
        rate_burial = self.vb * pom / self.h2

        # Source: POC settling. Disabled when ``use_POC`` is False or
        # when ``vsoc == 0`` (closed-system mode).
        if self.use_POC:
            rate_poc_settling = self.vsoc * poc / self.h2 / self.fcom
        else:
            rate_poc_settling = xr.zeros_like(pom)

        # Source: floating-algae settling. Phase 3.5 inter-process
        # coupling: the FloatingAlgae Process caches
        # ``algal_pom_from_settling_rate`` (mg/L/d) on each ``run``
        # invocation. Read it via ``getattr`` so a missing cache
        # (FloatingAlgae not yet run, or absent) degrades to zero.
        if self.use_floating_algae and self.use_Algae and self.floating_algae_process is not None:
            rate_algal_settling = getattr(
                self.floating_algae_process,
                "algal_pom_from_settling_rate",
                0,
            )
        else:
            rate_algal_settling = xr.zeros_like(pom)

        # Source: benthic-algae mortality. The BenthicAlgae Process
        # caches ``balgae_pom_from_mortality_rate`` (mg/L/d). Same
        # getattr-with-zero fallback.
        if self.use_benthic_algae and self.use_Balgae and self.benthic_algae_process is not None:
            rate_benthic_mortality = getattr(
                self.benthic_algae_process,
                "balgae_pom_from_mortality_rate",
                0,
            )
        else:
            rate_benthic_mortality = xr.zeros_like(pom)

        rate = (
            rate_algal_settling
            - rate_dissolution
            + rate_poc_settling
            + rate_benthic_mortality
            - rate_burial
        )

        components = {
            "pom_hydrolysis_rate": rate_dissolution,
            "pom_settling_rate": rate_burial,
            "pom_algal_mortality_rate": rate_algal_settling,
            "pom_balgae_mortality_rate": rate_benthic_mortality,
        }

        return rate, components

    def _change_legacy_inline(
        self,
        *,
        pom: ArrayLike,
        water_temperature: ArrayLike,
        poc: ArrayLike,
        depth: ArrayLike,
    ) -> ArrayLike:
        """Pre-Phase-7 inline rate composition. **Verbatim copy** of
        the body that ``run`` invoked via ``self.rate(...)`` before
        Phase 7 of the pattern-alignment spec landed.

        The pre-Phase-7 ``run`` called ``self.rate(pom=..., ...,
        time=..., registry=...)`` which in turn read ``depth`` from
        the registry inside the helper. This shadow takes ``depth``
        as a direct argument (consistent with the new helper
        signature) but performs the same arithmetic on the same
        intermediates, so the per-substep net rate is bit-identical
        to the pre-Phase-7 path.

        Retained through Phase 10 for the helper-vs-inline parity
        test per §11.3. Deleted in Phase 10 alongside its parity test
        once the final end-to-end baseline parity passes.

        Returns just the net per-day rate (mg/L/d).
        """
        kpom_tc = arrhenius_correction(
            water_temperature, self.kpom_20, self.kpom_theta
        )
        rate_dissolution = kpom_tc * pom

        # Pre-Phase-7 ``rate()`` set this side effect; preserved here so
        # the shadow leaves the same self state.
        self.pom_doc_source_rate = (
            self.fcom * rate_dissolution * self.h2 / depth
        )

        rate_burial = self.vb * pom / self.h2

        if self.use_POC:
            rate_poc_settling = self.vsoc * poc / self.h2 / self.fcom
        else:
            rate_poc_settling = xr.zeros_like(pom)

        if self.use_floating_algae and self.use_Algae and self.floating_algae_process is not None:
            rate_algal_settling = getattr(
                self.floating_algae_process,
                "algal_pom_from_settling_rate",
                0,
            )
        else:
            rate_algal_settling = xr.zeros_like(pom)

        if self.use_benthic_algae and self.use_Balgae and self.benthic_algae_process is not None:
            rate_benthic_mortality = getattr(
                self.benthic_algae_process,
                "balgae_pom_from_mortality_rate",
                0,
            )
        else:
            rate_benthic_mortality = xr.zeros_like(pom)

        return (
            rate_algal_settling
            - rate_dissolution
            + rate_poc_settling
            + rate_benthic_mortality
            - rate_burial
        )

    # ------------------------------------------------------------------
    # Kinetic helpers
    # ------------------------------------------------------------------

    def rate(
        self,
        pom: ArrayLike,
        water_temperature: ArrayLike,
        poc: ArrayLike,
        time: datetime,
        registry: VariableRegistry,
    ) -> ArrayLike:
        """Compute the rate of change of POM (mg/L/d).

        Returns the sum of source/sink terms per v1 ``dPOMdt``.
        """
        # Temperature-corrected dissolution rate.
        kpom_tc = arrhenius_correction(
            water_temperature, self.kpom_20, self.kpom_theta
        )

        # Sink: dissolution. POM mass loss is in mg-D/L_sed/d (POM state
        # is dry-mass-per-sediment-volume). Used as a sink in dPOM/dt.
        rate_dissolution = kpom_tc * pom

        # Phase 5.A consumer-ready DOC source rate (mg-C/L_water/d).
        # The raw POM dissolution rate (mg-D/L_sed/d) becomes a DOC
        # source by (1) converting dry mass to carbon mass via
        # ``fcom`` (mg-C/mg-D) and (2) converting sediment volumetric
        # to water-column volumetric via ``h2 / depth``. Matches the
        # convention of ``poc_hydrolysis_rate`` on Carbon (mg-C/L/d
        # water-column). Required for closed-system C conservation.
        depth = registry.get_at_time("depth", time)
        self.pom_doc_source_rate = (
            self.fcom * rate_dissolution * self.h2 / depth
        )

        # Sink: burial.
        rate_burial = self.vb * pom / self.h2

        # Source: POC settling. Disabled when ``use_POC`` is False or
        # when ``vsoc == 0`` (closed-system mode).
        if self.use_POC:
            rate_poc_settling = self.vsoc * poc / self.h2 / self.fcom
        else:
            rate_poc_settling = xr.zeros_like(pom)

        # Source: floating-algae settling. Phase 3.5 inter-process
        # coupling: the Phase 2.A FloatingAlgae Process caches
        # ``algal_pom_from_settling_rate`` (v1 ``vsap * Ap * (AWd/AWa) /
        # h2``, mg/L/d) on each ``run`` invocation. Read it via
        # ``getattr`` so a missing cache (FloatingAlgae not yet run, or
        # absent) degrades to zero.
        if self.use_floating_algae and self.use_Algae and self.floating_algae_process is not None:
            rate_algal_settling = getattr(
                self.floating_algae_process,
                "algal_pom_from_settling_rate",
                0,
            )
        else:
            rate_algal_settling = xr.zeros_like(pom)

        # Source: benthic-algae mortality. Phase 3.5 inter-process
        # coupling: the Phase 2.A BenthicAlgae Process caches
        # ``balgae_pom_from_mortality_rate`` (v1 ``Ab * kdb_tc * Fb *
        # (1 - Fw) / h2``, mg/L/d). Same getattr-with-zero fallback.
        if self.use_benthic_algae and self.use_Balgae and self.benthic_algae_process is not None:
            rate_benthic_mortality = getattr(
                self.benthic_algae_process,
                "balgae_pom_from_mortality_rate",
                0,
            )
        else:
            rate_benthic_mortality = xr.zeros_like(pom)

        return (
            rate_algal_settling
            - rate_dissolution
            + rate_poc_settling
            + rate_benthic_mortality
            - rate_burial
        )
