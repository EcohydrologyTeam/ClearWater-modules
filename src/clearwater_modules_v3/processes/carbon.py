"""v3 NSM1 ``Carbon`` Process (POC, DOC, DIC).

Phase 5.A of the v3 NSM1 implementation plan
(``design/clearwater_modules_v3_nsm1_design_specification.md`` Section 11
Phase 5, Section 5 Carbon design notes). v3-native (NOT a v2 overlay; v2
has no Carbon Process).

A single Carbon Process owns three state variables (mg/L):

* ``poc``  -- particulate organic carbon
* ``doc``  -- dissolved organic carbon
* ``dic``  -- dissolved inorganic carbon

Kinetics (mirrors v1 ``processes.py:2439-2870`` and Fortran
``modCarbon.f90``):

    dPOC/dt =   algal_poc_from_mortality_rate                # FloatingAlgae
              + balgae_poc_from_mortality_rate               # BenthicAlgae
              - kpoc_tc * POC                                # POC -> DOC hydrolysis
                                                             # (no DOX-Monod, per
                                                             # Fortran modCarbon.f90:170
                                                             # and v1 processes.py:2455)
              - vsoc / depth * POC                           # settling

    dDOC/dt = + kpoc_tc * POC                                # POC hydrolysis source
              + algal_doc_from_mortality_rate                # FloatingAlgae
              + balgae_doc_from_mortality_rate               # BenthicAlgae
              + getattr(pom_process, "pom_doc_source_rate",
                        0)                                   # POM dissolution -> DOC (if POM); mg-C/L_water/d
              - DOX_attenuation * kdoc_tc * DOC              # DOC -> DIC oxidation

    dDIC/dt = + DOX_attenuation * kdoc_tc * DOC              # DOC oxidation -> DIC (mg-C/L/d)
              + 0.923 * ka_tc * (KH * pCO2 / 1e6 * 12000
                                  - FCO2 * DIC)              # CO2 reaeration (atm); Henry's law converted from mol-C/L to mg-C/L
              + rca * (algal_respiration_rate
                       - algal_growth_rate)                  # FloatingAlgae C
              + rcb * Fb / depth *
                       (balgae_respiration_rate
                        - balgae_growth_rate)                # BenthicAlgae C
              + cbod_oxidation_rate / roc                    # CBOD -> DIC oxidation
              + JDIC / depth                                 # sediment release

Where:

* ``DOX_attenuation = DOX / (KsOxmc + DOX)`` (Monod attenuation by
  dissolved oxygen). Reads DOX from the Phase 5.B ``DOX`` Process via
  ``getattr(self.dox_process, "dox", _DOX_FALLBACK_MG_PER_L)``; the
  fallback ``8.0 mg/L`` is a typical surface-water saturation value used
  while DOX is being built concurrently.
* ``kpoc_tc``, ``kdoc_tc`` -- Arrhenius-corrected rate constants from
  ``CARBON_DEFAULTS`` (``kpoc_20``, ``kdoc_20``, ``kpoc_theta``,
  ``kdoc_theta``).
* ``vsoc`` -- POC settling velocity (m/d), from ``GLOBAL_VAR_DEFAULTS``.
* ``ka_tc`` -- effective reaeration coefficient from
  ``utils.reaeration.ka_tc`` (Phase 3.4 N2 pattern). Reused for CO2 since
  the v1 model assumes CO2 and O2 share the gas-transfer coefficient
  (a customary simplification for natural waters).
* ``KH`` -- Henry's-law constant for CO2 (mol/L/atm), temperature-
  corrected via the v1 empirical formula
  ``10**(2385.73 / Tk + 0.0152642 * Tk - 14.0184)``.
* ``pCO2`` -- atmospheric CO2 partial pressure (ppm; default 383 per the
  Phase 0.2 audit).
* ``FCO2`` -- fraction of DIC as free CO2 (depends on pH; per the
  resolved Q14 simple-tracer assumption, v3 1.0.0 does not solve for pH
  and uses a constant ``FCO2`` from the parameter file). Documented as a
  simple-tracer placeholder until the carbonate solver lands (v3 1.x).
* ``rca = AWc / AWa`` -- floating algal C:Chl-a stoichiometric ratio
  (mg-C/ug-Chla; per v1 ``rca`` helper at ``processes.py:337-348``).
  ``rcb = BWc / BWd`` -- benthic algal C:dry-weight ratio (mg-C/mg-D;
  per v1 ``rcb`` helper at ``processes.py:776-786``). v3 (Phase 9.B
  audit) computes these once at the top of ``run`` rather than passing
  the raw weights ``AWc`` / ``BWc``; the prior v3 implementation used
  ``AWc`` / ``BWc`` directly which yielded a 1000x / 100x scaling error
  in DIC algal coupling. ``Fb`` -- bottom-area fraction.
* ``JDIC`` -- DIC sediment release flux (g/m^2/d). v1 derives this from
  SOD (``SOD_tc / roc``); a Phase 5.B sediment integration may rewire
  this. For Phase 5.A standalone the term is identically zero unless the
  user passes ``JDIC > 0`` in the parameters dict OR ``use_SedFlux`` is
  ``True`` AND the SOD-derived path is enabled.

Q10 GS-rates contract: the per-step DOC->DIC oxidation flux
``self.doc_dic_oxidation_rate`` (mg-C/L/d) is cached as an instance
attribute after ``run`` completes. Phase 5.B DOX consumes it as an O2
sink term (multiplied by ``roc = 32/12``).

Forward Euler integrator pattern: rates are 1/d, ``dt_days =
time_step.total_seconds() / 86400``, ``state_new = state +
rate * dt_days``, then ``clip_negative_state`` with diagnostics, then
``set_at_time``.

DEFAULTS sources: the v3 ``CARBON_DEFAULTS`` dict (``parameters.carbon``)
holds carbon-specific kinetics; reaeration and sediment-flux parameters
live in ``parameters.dox`` and ``parameters.global_vars`` respectively.
The merged DEFAULTS are applied to the instance as ``self.<name>``
attributes, mirroring the Phase 3.4 N2 composed-defaults pattern.

Per design spec Section 11, this Process is **not** wired into the
package ``__init__`` here -- registration is Phase 5.5 integration.
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
from clearwater_modules_v3.utils.reaeration import kah_20, kaw_20, ka_tc


if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


logger = logging.getLogger(__name__)


# Stub DOX concentration used when the registry/DOX Process is not yet
# wired in (Phase 5.A standalone; Phase 5.B DOX is being built
# concurrently). 8.0 mg/L is a typical surface-water saturation value.
_DOX_FALLBACK_MG_PER_L: float = 8.0


def henrys_k_co2(t_water_c: ArrayLike) -> ArrayLike:
    """Temperature-dependent Henry's-law constant for CO2 (mol/L/atm).

    v1 source: ``clearwater_modules/nsm1/processes.py:Henrys_k`` (2687-2695).

    Note: this is the v1 empirical formula, retained here because no
    closed-form Henry's-law constant for CO2 is exposed in
    ``utils.reaeration``. Tk is water temperature in Kelvin.
    """
    t_water_k = t_water_c + 273.15
    return 10.0 ** (
        2385.73 / t_water_k + 0.0152642 * t_water_k - 14.0184
    )


class Carbon(Process):
    """v3 NSM1 Carbon Process (POC, DOC, DIC).

    Forward Euler integrator advances POC, DOC, DIC by ``dt_days =
    time_step / 86400`` each substep. Negative cells are clipped to zero
    via ``clip_negative_state`` with diagnostics on
    ``self.diagnostics``.

    Step-scoped rate caches consumed by sibling Processes:

    * ``self.doc_dic_oxidation_rate`` -- mg-C/L/d, consumed by Phase 5.B
      DOX (sink, multiplied by ``roc``).
    * ``self.poc_hydrolysis_rate``    -- mg-C/L/d, exposed for diagnostic
      / future POM coupling symmetry.
    """

    variables = [
        "poc",
        "doc",
        "dic",
        "water_temperature",
        "depth",
    ]

    # Class-level v3 defaults. Lazy-loaded on first instantiation; the
    # composed dict pulls from ``parameters.carbon`` (kinetics),
    # ``parameters.dox`` (reaeration menu), ``parameters.global_vars``
    # (POC settling velocity, hydraulic forcings), and
    # ``parameters.algae`` / ``parameters.balgae`` (stoichiometric
    # ratios used in the algal-coupling terms).
    DEFAULTS: dict[str, float | int | bool] = {}

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface Carbon exposes via the opportunistic-write loop in
    # ``run``. Each name maps to a ``self.<name>`` cache attribute set
    # inside ``_change_with_components`` and matches the inventory in
    # ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3.
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "poc_hydrolysis_rate",
        "doc_dic_oxidation_rate",
        "dic_atm_exchange_rate",
        "dic_sed_release_rate",
        "carbon_algal_resp_rate",
        "carbon_balgae_resp_rate",
        "carbon_algal_photo_rate",
        "carbon_balgae_photo_rate",
        "carbon_cbod_oxidation_rate",
    )

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
    ) -> None:
        """Initialize the Carbon process.

        Args:
            parameters: Optional dict of v3 Carbon parameter overrides.
                Merged with the class-level composed ``DEFAULTS``.
                Unknown keys are warned and ignored.
            time_step: Substep cadence for this Process.
        """
        Process.__init__(self, time_step)

        # Lazy-load composed DEFAULTS (parameters.carbon is augmented by
        # selected entries from dox / global_vars / algae / balgae /
        # global_parameters so that Carbon can compute the CO2
        # reaeration term and the algal-coupling terms without requiring
        # those to be passed in by the user). Phase 5.5 integration may
        # consolidate via the model-level parameter registry.
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.carbon import (
                DEFAULTS as CARBON_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.dox import (
                DEFAULTS as DOX_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.global_vars import (
                DEFAULTS as GLOBAL_VAR_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.global_parameters import (
                DEFAULTS as GLOBAL_PARAM_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.algae import (
                DEFAULTS as ALGAE_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.balgae import (
                DEFAULTS as BALGAE_DEFAULTS,
            )

            composed: dict[str, float | int | bool] = {}
            # Carbon-specific kinetics (kpoc_20, kdoc_20, kpoc_theta,
            # kdoc_theta, KsOxmc, pCO2, FCO2, roc, f_pocp, f_pocb).
            composed.update(CARBON_DEFAULTS)
            # Reaeration menu (used for the CO2 atmospheric exchange term).
            for k in (
                "kah_20_user",
                "kaw_20_user",
                "kah_theta",
                "kaw_theta",
                "hydraulic_reaeration_option",
                "wind_reaeration_option",
            ):
                composed[k] = DOX_DEFAULTS[k]
            # POC settling velocity + hydraulic forcings.
            for k in (
                "vsoc",
                "velocity",
                "flow",
                "topwidth",
                "slope",
                "shear_velocity",
                "wind_speed",
            ):
                composed[k] = GLOBAL_VAR_DEFAULTS[k]
            # Feature flags.
            for k in (
                "use_DOX",
                "use_Algae",
                "use_Balgae",
                "use_POC",
                "use_DOC",
                "use_DIC",
                "use_POM",
                "use_SedFlux",
            ):
                composed[k] = GLOBAL_PARAM_DEFAULTS.get(k, True)
            # Algal stoichiometric ratios for DIC photosynthesis /
            # respiration coupling. v3 derives ``rca = AWc / AWa`` and
            # ``rcb = BWc / BWd`` at run time per the Phase 9.B audit
            # (Fortran ``modAlgae.f90`` / ``modBenthicAlgae.f90`` use the
            # same identities); the raw weights ``AWc`` / ``BWc`` are
            # composed here so the user can override either the raw
            # weights or the derived ratio knobs.
            composed["AWc"] = ALGAE_DEFAULTS["AWc"]   # mg-C raw weight
            composed["AWa"] = ALGAE_DEFAULTS["AWa"]   # ug-Chla per algal unit
            composed["BWc"] = BALGAE_DEFAULTS["BWc"]  # mg-C / g-D raw weight
            composed["BWd"] = BALGAE_DEFAULTS["BWd"]  # mg-D / g-D dry-weight
            # Sediment release flux placeholder (g-C/m^2/d). Defaults to
            # zero; user can opt in via ``parameters={"JDIC": ...}``.
            composed["JDIC"] = 0.0
            # Benthic-area fraction (used in the BenthicAlgae DIC
            # coupling term). v2 BenthicAlgae stores this as ``Fb``.
            composed["Fb"] = float(BALGAE_DEFAULTS.get("Fb", 1.0))
            type(self).DEFAULTS = composed

        # Merge user overrides over the composed defaults.
        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "Carbon: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in CARBON_DEFAULTS / dox / global_vars / "
                "global_parameters / algae / balgae).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Step-scoped rate caches consumed by sibling Processes
        # (Q10 GS-rates contract). doc_dic_oxidation_rate is the sink
        # term Phase 5.B DOX adds to its O2 integrator after multiplying
        # by ``roc = 32/12 mg-O2/mg-C``.
        self.doc_dic_oxidation_rate: ArrayLike = 0.0
        self.poc_hydrolysis_rate: ArrayLike = 0.0

        # Diagnostics handle: a fresh v3 ``Diagnostics`` until
        # ``init_process`` overrides with the v3 Model's run-level
        # handle.
        self.diagnostics: Diagnostics = Diagnostics()

        # Coupling flags. Defaulted False so ``run`` works without an
        # explicit ``init_process`` call (e.g. Tier 1 unit tests).
        self.use_floating_algae: bool = False
        self.use_benthic_algae: bool = False
        self.use_pom: bool = False
        self.use_dox: bool = False
        self.use_cbod: bool = False
        self.floating_algae_process = None
        self.benthic_algae_process = None
        self.pom_process = None
        self.dox_process = None
        self.cbod_process = None

    @ProcessFactory.register("carbon")
    @staticmethod
    def from_config(
        config: dict, variable_registry: VariableRegistry
    ) -> "Carbon":
        return Carbon(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """Capture cross-process refs and the run-level Diagnostics handle."""
        # Capture v3 Model's run-level Diagnostics if present.
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

        # Optional couplings. Each Process is read with
        # ``hasattr/has_process`` so absent sibling Processes degrade to
        # zero source/sink contributions in ``run``.
        if hasattr(model, "has_process"):
            self.use_floating_algae = model.has_process("FloatingAlgae")
            self.use_benthic_algae = model.has_process("BenthicAlgae")
            self.use_pom = model.has_process("POM")
            self.use_dox = model.has_process("DOX")
            self.use_cbod = model.has_process("CBOD")
            if self.use_floating_algae:
                self.floating_algae_process = model.get_process("FloatingAlgae")
            if self.use_benthic_algae:
                self.benthic_algae_process = model.get_process("BenthicAlgae")
            if self.use_pom:
                self.pom_process = model.get_process("POM")
            if self.use_dox:
                self.dox_process = model.get_process("DOX")
            if self.use_cbod:
                self.cbod_process = model.get_process("CBOD")

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Advance POC, DOC, DIC by one substep using Forward Euler.

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_change_with_components``
        (B); applies Forward Euler with unconditional clip-with-log (C,
        D); persists primary outputs (E); caches step-scoped rates on
        ``self.<name>`` (F); opportunistically writes diagnostics (G).

        (The Phase-10 ``_change_legacy_inline`` shadow method and its
        ``test_carbon_helper_vs_inline.py`` parity test were removed once
        the helper refactor was settled; the coupled-demo bit-identical
        baseline parity test is the standing regression guard.)
        """
        # --- State reads (pattern A) ---
        poc = registry.get_at_time("poc", time)
        doc = registry.get_at_time("doc", time)
        dic = registry.get_at_time("dic", time)
        t_water_c = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)

        # DOX coupling: prefer the registry (Phase 5.B DOX persists its
        # state under ``oxygen_dissolved``), then probe a hypothetical
        # ``self.dox`` cache (not present on Phase 5.B DOX as of writing,
        # but probed for forward compatibility), then fall back to the
        # 8.0 mg/L stub. The stub serves Phase 5.A standalone tests
        # where neither DOX nor the registry has the variable.
        if self.use_dox and self.dox_process is not None:
            cached_dox = getattr(self.dox_process, "dox", None)
            if cached_dox is not None:
                dox = cached_dox
            else:
                dox = _dox_from_registry(registry, time, poc)
        else:
            dox = _dox_from_registry(registry, time, poc)

        # --- Fused rate composition (pattern B) ---
        d_poc, d_doc, d_dic, components = self._change_with_components(
            poc=poc,
            doc=doc,
            dic=dic,
            t_water_c=t_water_c,
            depth=depth,
            dox=dox,
        )

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        # Names match the REGISTRY_DIAGNOSTICS tuple and the spec §4
        # inventory. Preserved attribute names: ``doc_dic_oxidation_rate``
        # and ``poc_hydrolysis_rate`` were already published by Carbon
        # for DOX consumption; do not rename.
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])

        # --- Forward Euler in days (pattern C) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        poc_new = poc + d_poc * dt_days
        doc_new = doc + d_doc * dt_days
        dic_new = dic + d_dic * dt_days

        # --- Clip-with-log per the resolved Q7 contract (pattern D) ---
        poc_new = clip_negative_state(poc_new, "poc", self.diagnostics)
        doc_new = clip_negative_state(doc_new, "doc", self.diagnostics)
        dic_new = clip_negative_state(dic_new, "dic", self.diagnostics)

        # --- Persist primary outputs (pattern E) ---
        registry.set_at_time("poc", time, poc_new)
        registry.set_at_time("doc", time, doc_new)
        registry.set_at_time("dic", time, dic_new)

        # --- Opportunistic diagnostic registry writes (pattern G) ---
        # Each Appendix-A name is written ONLY if the user has
        # pre-registered it. Zero cost when not subscribed; lets
        # calibration / validation runs subscribe to any subset.
        for name in self.REGISTRY_DIAGNOSTICS:
            if name in registry:
                registry.set_at_time(name, time, components[name])

    # ------------------------------------------------------------------
    # Rate-composition helpers
    # ------------------------------------------------------------------

    def _change_with_components(
        self,
        *,
        poc: ArrayLike,
        doc: ArrayLike,
        dic: ArrayLike,
        t_water_c: ArrayLike,
        depth: ArrayLike,
        dox: ArrayLike,
    ) -> tuple[ArrayLike, ArrayLike, ArrayLike, dict]:
        """Compute ``(d_poc, d_doc, d_dic, components)``.

        Code-motion-only refactor of ``run``'s former inline composition
        (§11.6): operand order, intermediate names, and arithmetic are
        preserved verbatim from the pre-refactor body. The ``components``
        dict is populated from the same intermediates the integrator
        consumes (no recomputation).

        Phase 2 pattern-alignment spec §6 deliverable. The companion
        ``_change_legacy_inline`` shadow and its
        ``test_carbon_helper_vs_inline.py`` parity test were removed once
        the helper refactor was settled and proven bit-identical; the
        standing regression guard is the coupled-demo bit-identical
        baseline parity test (``test_coupled_demo_parity.py``).
        """
        # --- Temperature-corrected rate constants (Arrhenius / van't Hoff) ---
        kpoc_tc_value = arrhenius_correction(
            t_water_c, self.kpoc_20, self.kpoc_theta
        )
        kdoc_tc_value = arrhenius_correction(
            t_water_c, self.kdoc_20, self.kdoc_theta
        )

        # --- DOX-Monod attenuation (DOC -> DIC oxidation only) ---
        # Per Fortran ``modCarbon.f90:198`` and v1 ``processes.py:2629``,
        # DOC oxidation uses ``DOX / (KsOxmc + DOX)`` Monod attenuation.
        # POC hydrolysis is *not* attenuated by DOX in either reference;
        # earlier v3 applied the same factor to POC hydrolysis but Phase
        # 9.B audit removed that factor for parity.
        dox_attenuation = dox / (self.KsOxmc + dox)

        # --- Stoichiometric C-to-Chla / C-to-D ratios derived from raw
        # weights (Phase 9.B audit fix; Fortran ``modAlgae.f90`` and v1
        # ``processes.py`` derive these the same way). Cached once so the
        # arithmetic below stays in mg-C/L/d cleanly. ---
        rca = self.AWc / self.AWa     # mg-C / ug-Chla
        rcb = self.BWc / self.BWd     # mg-C / mg-D

        # --- POC kinetic terms (mg-C/L/d) ---
        # POC hydrolysis: pure first-order (no DOX-Monod), per Fortran
        # and v1.
        poc_hydrolysis = kpoc_tc_value * poc
        poc_settling = self.vsoc / depth * poc
        poc_algal_mortality = self._poc_algal_mortality(poc)
        poc_balgae_mortality = self._poc_balgae_mortality(depth, poc)

        d_poc = (
            poc_algal_mortality
            + poc_balgae_mortality
            - poc_hydrolysis
            - poc_settling
        )

        # --- DOC kinetic terms (mg-C/L/d) ---
        doc_oxidation = kdoc_tc_value * doc * dox_attenuation
        doc_algal_mortality = self._doc_algal_mortality(doc)
        doc_balgae_mortality = self._doc_balgae_mortality(depth, doc)
        # POM dissolution -> DOC source (mg-C/L_water/d). The cache
        # ``pom_doc_source_rate`` on POM Process is in water-column
        # volumetric mg-C/L/d, with ``fcom`` (mg-C/mg-D) and ``h2/depth``
        # (sediment->water-column volume) already applied -- consumer-
        # ready for direct addition to dDOC/dt. Without those factors
        # the v3 implementation produces a closed-system C conservation
        # leak proportional to ``(1 - fcom * h2 / depth)`` per dissolved
        # mg-D of POM (about 25x overcount at default depth=1 m / h2=0.1
        # m / fcom=0.4). Defaults to 0 when POM is absent.
        if self.use_pom and self.pom_process is not None:
            pom_doc_source = getattr(
                self.pom_process, "pom_doc_source_rate", 0
            )
        else:
            pom_doc_source = 0

        d_doc = (
            poc_hydrolysis
            + doc_algal_mortality
            + doc_balgae_mortality
            + pom_doc_source
            - doc_oxidation
        )

        # --- DIC kinetic terms (mg-C/L/d) ---
        # Phase 9.E DIC unit reconciliation: the v3 (and v1) DIC budget
        # previously inherited Fortran ``modCarbon.f90:268``'s
        # ``/ 12000.0`` divisions on every explicit-formula term, which
        # produces a rate in mol-C/L/d. v1 and v3 store DIC as mg-C/L
        # (matches Fortran's labeling at ``modMain.f90:301``: "mg-C/L"),
        # so the mol-C/L/d rates were implicitly being added to a
        # mg-C/L state -- a 12000x scaling error that effectively froze
        # DIC dynamics. Phase 9.E removes the ``/ 12000.0`` from every
        # explicit-formula DIC source/sink and converts the Henry's-law
        # atmospheric-equilibrium term from mol-C/L to mg-C/L by
        # multiplying by ``MG_C_PER_MOL_C = 12 g/mol * 1000 mg/g =
        # 12000``. After the fix every dDIC/dt term is in mg-C/L/d,
        # consistent with the mg-C/L state. This is a v3 correction
        # over Fortran/v1; see ``parameter_defaults_corrections.md``
        # Section 1.11 (Phase 9.E DIC unit reconciliation).
        MG_C_PER_MOL_C = 12000.0  # 12 g-C/mol * 1000 mg-C/g = mg-C per mol-C

        # CO2 atmospheric reaeration: 0.923 * ka_tc * ([CO2*]_eq - [CO2]).
        # KH is mol/L/atm; pCO2 is ppm (= 1e-6 atm). Henry's-law product
        # ``KH * pCO2 / 1e6`` is in mol-C/L; multiply by MG_C_PER_MOL_C
        # to convert to mg-C/L for unit consistency with ``FCO2 * DIC``
        # (both terms now mg-C/L; resulting rate mg-C/L/d).
        ka_tc_value = self._ka_tc(t_water_c, depth)
        kh_co2 = henrys_k_co2(t_water_c)
        co2_reaeration = (
            0.923 * ka_tc_value
            * (kh_co2 * self.pCO2 / 1.0e6 * MG_C_PER_MOL_C - self.FCO2 * dic)
        )

        # Algal photosynthesis / respiration -> DIC source/sink.
        # Floating algae: rates are stored as ug-Chla/L/d on the
        # FloatingAlgae Process; ``rca = AWc / AWa`` is mg-C per
        # ug-Chla. Product is mg-C/L/d directly (no further unit
        # conversion needed; Phase 9.E removed the ``/ 12000`` from
        # this term).
        algae_growth = self._floating_algae_growth_rate()
        algae_respiration = self._floating_algae_respiration_rate()
        dic_algal_resp = algae_respiration * rca
        dic_algal_photo = algae_growth * rca

        # Benthic algae: rates are g-D/m^2/d on the BenthicAlgae
        # Process; ``rcb = BWc / BWd`` is mg-C per mg-D = g-C/g-D.
        # ``rcb * balgae * Fb / depth``: g-C/m^2/d / m = g-C/m^3/d =
        # mg-C/L/d directly. Phase 9.E removed the ``/ 12000``.
        balgae_growth = self._benthic_algae_growth_rate()
        balgae_respiration = self._benthic_algae_respiration_rate()
        dic_balgae_resp = balgae_respiration * rcb * self.Fb / depth
        dic_balgae_photo = balgae_growth * rcb * self.Fb / depth

        # Sediment release: ``JDIC`` is g-C/m^2/d. ``JDIC / depth``:
        # g-C/m^3/d = mg-C/L/d directly (1 g/m^3 = 1 mg/L). Phase 9.E
        # removed the ``/ 12000`` from this term.
        if self.use_SedFlux:
            dic_sed_release = self.JDIC / depth
        else:
            dic_sed_release = 0.0

        # CBOD oxidation -> DIC source. Per Fortran
        # ``modCarbon.f90:262-266`` and v1 ``processes.py:2793-2814``.
        # ``cbod_oxidation_rate`` is mg-O2/L/d; dividing by ``roc =
        # 32/12 g-O2/g-C`` gives mg-C/L/d directly. Phase 9.E removed
        # the ``/ 12000``. (Phase 9.B audit C3 added the missing source
        # term itself; Phase 9.E corrects its unit scaling.)
        if self.use_cbod and self.cbod_process is not None:
            cbod_ox_rate = getattr(
                self.cbod_process, "cbod_oxidation_rate", 0
            )
            if cbod_ox_rate is None:
                cbod_ox_rate = 0
            dic_cbod_oxidation = cbod_ox_rate / self.roc
        else:
            dic_cbod_oxidation = 0.0

        d_dic = (
            doc_oxidation
            + co2_reaeration
            + dic_algal_resp
            - dic_algal_photo
            + dic_balgae_resp
            - dic_balgae_photo
            + dic_cbod_oxidation
            + dic_sed_release
        )

        # --- NaN guards on the rates (mirrors Nitrogen) ---
        d_poc = sanitize_rate(d_poc)
        d_doc = sanitize_rate(d_doc)
        d_dic = sanitize_rate(d_dic)

        # --- Components dict (pattern G + F single source of truth) ---
        # Each name matches REGISTRY_DIAGNOSTICS and the spec §4
        # inventory. Preserved attribute names ``doc_dic_oxidation_rate``
        # and ``poc_hydrolysis_rate`` carry the same sanitize_rate
        # treatment as in the pre-refactor code so DOX (sibling consumer)
        # reads identical values.
        #
        # Sanitize NaN at the cache source: a NaN here propagates via
        # DOX's rate sum and zeroes the entire cell's DOX rate, freezing
        # the cell at IC. When DOC/POC are 0 the contribution is 0
        # regardless.
        components = {
            "poc_hydrolysis_rate": sanitize_rate(poc_hydrolysis),
            "doc_dic_oxidation_rate": sanitize_rate(doc_oxidation),
            "dic_atm_exchange_rate": sanitize_rate(co2_reaeration),
            "dic_sed_release_rate": sanitize_rate(dic_sed_release),
            "carbon_algal_resp_rate": sanitize_rate(dic_algal_resp),
            "carbon_balgae_resp_rate": sanitize_rate(dic_balgae_resp),
            "carbon_algal_photo_rate": sanitize_rate(dic_algal_photo),
            "carbon_balgae_photo_rate": sanitize_rate(dic_balgae_photo),
            "carbon_cbod_oxidation_rate": sanitize_rate(dic_cbod_oxidation),
        }

        return d_poc, d_doc, d_dic, components
    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _ka_tc(self, t_water_c: ArrayLike, depth: ArrayLike) -> ArrayLike:
        """Effective reaeration coefficient (1/d), temperature-corrected.

        Uses the same hydraulic + wind reaeration menu as DOX/N2. The
        ``kah_20`` and ``kaw_20`` utilities wrap ``np.select`` outputs in
        bare ``xr.DataArray`` constructors that produce a ``dim_0``
        anonymous axis rather than preserving the ``depth`` dim. To keep
        the downstream rate arithmetic in the proper per-cell shape, we
        re-wrap the combined ``ka_tc`` result with ``depth``'s dims if
        the value count matches; otherwise we fall through to scalar
        broadcasting.
        """
        kah_20_value = kah_20(
            kah_20_user=self.kah_20_user,
            hydraulic_reaeration_option=self.hydraulic_reaeration_option,
            velocity=self.velocity,
            depth=depth,
            flow=self.flow,
            topwidth=self.topwidth,
            slope=self.slope,
            shear_velocity=self.shear_velocity,
        )
        kaw_20_value = kaw_20(
            kaw_20_user=self.kaw_20_user,
            wind_speed=self.wind_speed,
            wind_reaeration_option=self.wind_reaeration_option,
        )
        return ka_tc(
            kah_20=kah_20_value,
            kaw_20=kaw_20_value,
            kah_theta=self.kah_theta,
            kaw_theta=self.kaw_theta,
            T_water_C=t_water_c,
            depth=depth,
        )

    # --- Floating algae coupling (rates cached on the FloatingAlgae Process) ---

    def _floating_algae_growth_rate(self) -> ArrayLike:
        if not (self.use_floating_algae and self.use_Algae):
            return 0.0
        if self.floating_algae_process is None:
            return 0.0
        return getattr(
            self.floating_algae_process, "algal_growth_rate", 0.0
        )

    def _floating_algae_respiration_rate(self) -> ArrayLike:
        if not (self.use_floating_algae and self.use_Algae):
            return 0.0
        if self.floating_algae_process is None:
            return 0.0
        return getattr(
            self.floating_algae_process, "algal_respiration_rate", 0.0
        )

    def _poc_algal_mortality(self, poc: ArrayLike) -> ArrayLike:
        """Floating algal mortality routed to POC (mg-C/L/d).

        Reads the cache ``algal_poc_from_mortality_rate`` populated by
        the Phase 2.A FloatingAlgae Process (``f_pocp * rca * ap_death``).
        """
        if not (self.use_floating_algae and self.use_Algae):
            return _zeros_like(poc)
        if self.floating_algae_process is None:
            return _zeros_like(poc)
        return getattr(
            self.floating_algae_process,
            "algal_poc_from_mortality_rate",
            _zeros_like(poc),
        )

    def _doc_algal_mortality(self, doc: ArrayLike) -> ArrayLike:
        """Floating algal mortality routed to DOC (mg-C/L/d).

        Reads the cache ``algal_doc_from_mortality_rate`` populated by
        the Phase 2.A FloatingAlgae Process
        (``(1 - f_pocp) * rca * ap_death``).
        """
        if not (self.use_floating_algae and self.use_Algae):
            return _zeros_like(doc)
        if self.floating_algae_process is None:
            return _zeros_like(doc)
        return getattr(
            self.floating_algae_process,
            "algal_doc_from_mortality_rate",
            _zeros_like(doc),
        )

    # --- Benthic algae coupling (rates cached on the BenthicAlgae Process) ---

    def _benthic_algae_growth_rate(self) -> ArrayLike:
        if not (self.use_benthic_algae and self.use_Balgae):
            return 0.0
        if self.benthic_algae_process is None:
            return 0.0
        return getattr(
            self.benthic_algae_process, "balgae_growth_rate", 0.0
        )

    def _benthic_algae_respiration_rate(self) -> ArrayLike:
        if not (self.use_benthic_algae and self.use_Balgae):
            return 0.0
        if self.benthic_algae_process is None:
            return 0.0
        return getattr(
            self.benthic_algae_process, "balgae_respiration_rate", 0.0
        )

    def _poc_balgae_mortality(
        self, depth: ArrayLike, poc: ArrayLike
    ) -> ArrayLike:
        """Benthic algal mortality routed to POC (mg-C/L/d).

        Reads the cache ``balgae_poc_from_mortality_rate`` populated by
        the Phase 2.A BenthicAlgae Process
        (``f_pocb * rcb * (kdb_tc * Ab * Fb * Fw) / depth``).
        """
        if not (self.use_benthic_algae and self.use_Balgae):
            return _zeros_like(poc)
        if self.benthic_algae_process is None:
            return _zeros_like(poc)
        return getattr(
            self.benthic_algae_process,
            "balgae_poc_from_mortality_rate",
            _zeros_like(poc),
        )

    def _doc_balgae_mortality(
        self, depth: ArrayLike, doc: ArrayLike
    ) -> ArrayLike:
        """Benthic algal mortality routed to DOC (mg-C/L/d).

        Reads the cache ``balgae_doc_from_mortality_rate`` populated by
        the Phase 2.A BenthicAlgae Process.
        """
        if not (self.use_benthic_algae and self.use_Balgae):
            return _zeros_like(doc)
        if self.benthic_algae_process is None:
            return _zeros_like(doc)
        return getattr(
            self.benthic_algae_process,
            "balgae_doc_from_mortality_rate",
            _zeros_like(doc),
        )


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _dox_from_registry(
    registry: VariableRegistry,
    time: datetime,
    template: ArrayLike,
) -> ArrayLike:
    """Read DOX from the registry; fall back to a constant stub.

    Phase 5.A standalone fallback used when the DOX Process is absent
    AND ``oxygen_dissolved`` is not in the registry. 8.0 mg/L is a
    typical surface-water saturation value.
    """
    if "oxygen_dissolved" in registry:
        return registry.get_at_time("oxygen_dissolved", time)
    if isinstance(template, xr.DataArray):
        return xr.full_like(template, _DOX_FALLBACK_MG_PER_L)
    return _DOX_FALLBACK_MG_PER_L


def _zeros_like(template: ArrayLike) -> ArrayLike:
    """Return zeros with the same dims as ``template`` (or scalar 0.0)."""
    if isinstance(template, xr.DataArray):
        return xr.zeros_like(template)
    return 0.0


