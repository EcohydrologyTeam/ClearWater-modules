"""v3 NSM1 ``DOX`` (dissolved oxygen) Process.

Phase 5.B (v3 NSM1 design spec, Section 11 Phase 5, Section 5 DOX design
notes). v3-native (NOT a v2 overlay). DOX is the largest single
constituent in NSM1 — it is a sink (and partial source) for nearly every
other Process. This Process is therefore terminal in the rate-DAG: it
*reads* step-scoped rate caches from algae / Nitrogen / Carbon / CBOD
upstream Processes and integrates the DOX state forward; it does NOT
publish any rate cache that downstream Processes consume (other than the
registry-state ``oxygen_dissolved`` itself, which several upstream
Processes Monod-attenuate against).

Kinetics (mirrors v1 ``processes.py`` lines 2876-3135):

    dDOX/dt = + ka_tc * (O2sat - DOX)                 (atmospheric reaeration)
              + ApGrowth * rca * roc * (138/106 - 32/106 * ApUptakeFr_NH4)
                                                       (floating algae photosynthesis)
              - ApRespiration * rca * roc              (floating algae respiration)
              + AbGrowth * rcb * roc * Fb / depth * (138/106 - 32/106 * AbUptakeFr_NH4)
                                                       (benthic algae photosynthesis,
                                                        depth-integrated)
              - AbRespiration * rcb * roc * Fb / depth (benthic algae respiration)
              - ron * nitrification_flux_rate          (nitrification O2 sink;
                                                        flux read from
                                                        Nitrogen Process per
                                                        Integration Item 1)
              - roc * doc_dic_oxidation_rate           (DOC -> DIC oxidation O2 sink)
              - cbod_oxidation_rate                    (CBOD oxidation: 1 mg-O2 per
                                                        mg-CBOD by definition)
              - SOD_tc / depth                         (sediment oxygen demand,
                                                        area-integrated)

Where:

* ``ka_tc`` -- ``utils.reaeration.ka_tc`` (combined hydraulic + wind
  reaeration, temperature-corrected). **NSM1-DOX-F2 (spec C4):**
  ``hydraulic_reaeration_option == 1`` with ``kah_20_user == 0.0``
  yields zero hydraulic reaeration (silently zero atmospheric
  reaeration if wind is also off). v3 preserves the unfloored value by
  default for v1/Fortran parity but warns once and offers an opt-in
  CE-QUAL-W2-``MINKL``-style floor ``min_reaeration_ka`` (default
  ``0.0`` = OFF; default behaviour unchanged).
* ``O2sat`` -- APHA / QUAL2E formulation; v1 lines 2901-2923. Requires
  ``T_water_K`` (Kelvin), ``pressure_mb`` (mb), and the empirical
  ``DOs_atm_alpha`` correction. The full APHA form is implemented (no
  simpler approximation is used). **Fresh-water assumption (NSM1-DOX-F1,
  spec C3):** the APHA *salinity* correction is intentionally omitted
  (matches v1; exact for fresh water). DOsat is overstated for
  brackish/marine water (~18% at 35 ppt); the salinity-corrected form
  is a documented deferral (audit C6 / NSM2). ``run`` warns once if a
  nonzero ``salinity`` is in the registry. See ``dox_sat_apha``.
* ``ron = 4.57`` -- v1 default ratio of O2 mass per N nitrified
  (= 2.0 * 32.0 / 14.0). Lives in ``DOX_DEFAULTS``.
* ``roc = 32/12`` -- O2 mass per C oxidized. Lives in
  ``CARBON_DEFAULTS``; DOX composes it on the instance.
* ``rca = AWc / AWa`` -- algal C:Chla stoichiometric ratio
  (mg-C/ug-Chla). Lives in ``ALGAE_DEFAULTS`` as the raw weights;
  derived at run time per the Phase 9.B audit (Fortran ``modAlgae.f90``
  and v1 ``processes.py:337-348`` derive the same way).
* ``rcb = BWc / BWd`` -- benthic algae C:dry-weight ratio
  (mg-C/mg-D). Lives in ``BALGAE_DEFAULTS`` as raw weights; derived
  at run time per the same audit (Fortran ``modBenthicAlgae.f90``
  and v1 ``processes.py:776-786``).
* ``138/106 - 32/106 * X`` -- v1 Redfield-derived photosynthesis
  stoichiometric factor that fractionates O2 production by NH4 vs NO3
  uptake. With ``ApUptakeFr_NH4 == 1`` (all NH4) the factor is 1.0;
  with ``ApUptakeFr_NH4 == 0`` (all NO3) the factor is 138/106 ≈ 1.30,
  reflecting the extra O2 produced when reducing NO3 back to organic-N.
* ``KNR`` -- O2 inhibition factor for nitrification (v1 ``KNR``,
  default 0.6 mg-O2/L). Lives in ``NITROGEN_DEFAULTS``; DOX composes it.
* ``knit_tc`` -- temperature-corrected nitrification rate constant
  (Arrhenius). Lives in ``NITROGEN_DEFAULTS`` as
  ``knit_20`` / ``knit_theta``; DOX composes them.
* ``SOD_tc`` -- ``utils.sediment.SOD_tc`` (Arrhenius-corrected SOD).
  Phase 1.1 made the utility pure Arrhenius (no DOX-Monod
  attenuation); Phase 9.B audit re-applies the Fortran
  ``modGlobalParam.f90:254`` form ``SOD_tc *= DOX / (DOX + KsSOD)``
  here in the DOX consumer. Under hypoxia this drops SOD toward zero
  (the sediment cannot deplete oxygen that is not there); without the
  correction, v3 over-consumed O2 at low DOX.

Coupling pattern (mirrors Phase 2.B Nitrogen for the multi-source-sink
integrator). Every coupled Process is consumed via ``getattr(...)``
with a default of ``0`` so a missing sibling degrades gracefully:

* FloatingAlgae:  ``algal_growth_rate``, ``algal_respiration_rate``,
                  ``algal_nh4_uptake_fraction``
* BenthicAlgae:   ``balgae_growth_rate``, ``balgae_respiration_rate``,
                  ``balgae_nh4_uptake_fraction``
* Nitrogen:       ``nitrification_flux_rate`` -- step-scoped
                  nitrification flux (mg-N/L/d, positive magnitude)
                  cached by Nitrogen.run per Integration Item 1
                  (registry rate-variable convention, spec resolved
                  Q10). DOX multiplies by ``ron`` to get the O2 sink.
* Carbon:         ``doc_dic_oxidation_rate`` (Phase 5.A sibling; may not
                  yet exist — ``getattr(..., 0)`` falls back gracefully).
* CBOD:           ``cbod_oxidation_rate`` (Phase 3.3, exists).

DEFAULTS sources (v3 ``parameters.dox`` provides the DOX-specific
kinetics: ``ron``, ``KsSOD``, ``SOD_20``, ``SOD_theta``, ``kaw_20_user``,
``kah_20_user``, ``kaw_theta``, ``kah_theta``,
``hydraulic_reaeration_option``, ``wind_reaeration_option``). DOX
additionally composes onto the instance:

* ``parameters.global_parameters``: ``pressure_mb`` (atmospheric pressure),
  ``use_NH4``, ``use_DOC``, ``use_Algae``, ``use_Balgae``
* ``parameters.global_vars``: hydraulic forcings (``velocity``, ``flow``,
  ``topwidth``, ``slope``, ``shear_velocity``, ``wind_speed``)
* ``parameters.algae``: ``AWc``, ``AWa`` (rca = AWc / AWa)
* ``parameters.balgae``: ``BWc``, ``BWd``, ``Fb``, ``Fw``
  (rcb = BWc / BWd)
* ``parameters.nitrogen``: ``KNR``, ``knit_20``, ``knit_theta``
* ``parameters.carbon``:   ``roc``

Per design spec Section 11, this Process is **not** wired into the
package ``__init__`` here — registration is Phase 5.5.
"""

from __future__ import annotations

from datetime import datetime, timedelta
import logging
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from clearwater_data.variables import VariableRegistry
from clearwater_data.custom_types import ArrayLike

from clearwater_modules_v3.processes.base import Process, ProcessFactory
from clearwater_modules_v3.utils.numerics import (
    Diagnostics,
    clip_negative_state,
    sanitize_rate,
)
from clearwater_modules_v3.utils.reaeration import kah_20, kaw_20, ka_tc
from clearwater_modules_v3.utils.sediment import SOD_tc as sod_tc_util

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# DO saturation primitives (stateless; v1 NSM1 lines 2876-2923, APHA form)
# ---------------------------------------------------------------------------


def _kelvin(t_c: ArrayLike) -> ArrayLike:
    """deg C -> K."""
    return t_c + 273.15


def _pwv_atm(t_water_k: ArrayLike) -> ArrayLike:
    """Partial pressure of water vapor (atm), v1 ``pwv`` (lines 2878-2886)."""
    return np.exp(11.8571 - 3840.70 / t_water_k - 216961.0 / t_water_k**2)


def _do_atm_alpha(t_water_c: ArrayLike) -> ArrayLike:
    """DO-saturation atmospheric correction coefficient.

    v1 ``DOs_atm_alpha`` (lines 2890-2898). Empirical fit in degrees C
    (note: v1 docstring incorrectly labels the input as Kelvin; the
    formula and the calling site both use Celsius).
    """
    return 0.000975 - 1.426e-5 * t_water_c + 6.436e-8 * t_water_c**2


def dox_sat_apha(
    t_water_c: ArrayLike,
    pressure_mb: ArrayLike,
) -> ArrayLike:
    """Dissolved oxygen saturation concentration (mg-O2/L), APHA form.

    Mirrors v1 ``DOX_sat`` (lines 2901-2923). The uncorrected (1 atm,
    fresh water) saturation is the Benson-Krause form; the pressure
    correction multiplies by ``P_atm * (1 - p_wv/P_atm) * (1 - alpha *
    P_atm) / ((1 - p_wv) * (1 - alpha))``.

    **Freshwater assumption (NSM1-DOX-F1, gold-standard spec C3).**
    This is the **fresh-water** APHA saturation. The APHA salinity
    correction term is intentionally **not** applied (matches v1; exact
    for fresh water — multiplicative factor 1.0). For brackish/marine
    water DOsat is overstated (~18% at 35 ppt salinity). CE-QUAL-W2
    applies the APHA-exact salinity correction; the salinity-corrected
    form is a documented deferral (audit C6 / NSM2 scope). ``DOX.run``
    emits a one-time warning if a nonzero ``salinity`` is present in the
    registry so brackish input cannot pass silently. No salinity input
    is taken here by design.

    Args:
        t_water_c | deg C | water temperature.
        pressure_mb | mb | atmospheric pressure (use 1013.25 for sea level).

    Returns:
        mg-O2/L | dissolved oxygen saturation concentration.
    """
    t_k = _kelvin(t_water_c)
    pressure_atm = pressure_mb * 0.000986923
    pwv_atm = _pwv_atm(t_k)
    alpha = _do_atm_alpha(t_water_c)

    # Uncorrected saturation (Benson-Krause).
    log_uncorrected = (
        -139.34410
        + 1.575701e5 / t_k
        - 6.642308e7 / (t_k**2)
        + 1.243800e10 / (t_k**3)
        - 8.621949e11 / (t_k**4)
    )
    dox_sat_uncorrected = np.exp(log_uncorrected)

    # Pressure correction.
    return (
        dox_sat_uncorrected
        * pressure_atm
        * (1.0 - pwv_atm / pressure_atm)
        * (1.0 - alpha * pressure_atm)
        / ((1.0 - pwv_atm) * (1.0 - alpha))
    )


# ---------------------------------------------------------------------------
# DOX Process
# ---------------------------------------------------------------------------


class DOX(Process):
    # Phase H-9 (2026-05-21): DOX reads step-scoped rate caches from
    # Nitrogen (nitrification_flux_rate), FloatingAlgae (algal_growth_rate,
    # algal_respiration_rate, algal_nh4_uptake_fraction), BenthicAlgae
    # (balgae_growth_rate, balgae_respiration_rate, balgae_nh4_uptake_fraction),
    # Carbon (doc_dic_oxidation_rate), and CBOD (cbod_oxidation_rate).
    # Each must run BEFORE DOX in the same substep. Model.__init_model
    # validates the registered process order against this list.
    upstream_processes = (
        "Nitrogen", "FloatingAlgae", "BenthicAlgae", "Carbon", "CBOD",
    )

    """v3 NSM1 dissolved oxygen Process.

    State variable: ``oxygen_dissolved`` (mg-O2/L).

    Sources / sinks (mg-O2/L/d, integrated by Forward Euler):

    * Atmospheric reaeration:                   ``ka_tc * (O2sat - DOX)``
    * Algal photosynthesis (floating + benthic): O2 produced, fractionated
      by the NH4 vs NO3 uptake split via the Redfield 138/106 factor.
    * Algal respiration (floating + benthic):    O2 consumed.
    * Nitrification (NH4 -> NO3):                ``ron * knit_tc * NH4 *
      (1 - exp(-KNR * DOX))``
    * DOC oxidation (DOC -> DIC):                ``roc *
      doc_dic_oxidation_rate`` from Carbon (Phase 5.A sibling).
    * CBOD oxidation:                            ``cbod_oxidation_rate``
      from CBOD (Phase 3.3 sibling, 1 mg-O2 per mg-CBOD by definition).
    * Sediment oxygen demand (SOD):              ``SOD_tc / depth``.

    Caches step-scoped quantities for diagnostics and Phase 5.5
    integration tests:

    * ``self.dox_sat``                  -- saturation concentration (mg/L)
    * ``self.atm_reaeration_rate``      -- atmospheric reaeration flux (mg/L/d)
    * ``self.dox_nitrification_rate``   -- O2 sink from nitrification (mg/L/d)
    * ``self.dox_sod_rate``             -- O2 sink from SOD (mg/L/d)
    """

    variables = [
        "oxygen_dissolved",
        "water_temperature",
        "depth",
    ]

    # Class-level v3 defaults. Lazy-loaded on first instantiation;
    # composed from the relevant v3 parameter groups (see module
    # docstring).
    DEFAULTS: dict[str, float | int | bool] = {}

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface DOX exposes via the opportunistic-write loop in ``run``.
    # Each name maps to a ``self.<name>`` cache attribute set inside
    # ``_change_with_components`` and matches the inventory in
    # ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3.
    #
    # Note: ``sod_rate`` and ``dox_sod_rate`` are aliases for the same
    # volumetric SOD sink (mg-O2/L/d). ``dox_sod_rate`` is the v3 cache
    # attribute DOX has always written; ``sod_rate`` is the Appendix A
    # name inherited from the 1.0.0 "sediment-globals" naming. Both are
    # exposed so legacy consumers (and the Appendix A catalog) line up.
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "dox_sat",
        "atm_reaeration_rate",
        "dox_nitrification_rate",
        "dox_sod_rate",
        "sod_rate",
        "dox_doc_oxidation_rate",
        "dox_cbod_oxidation_rate",
        "dox_algal_photo_rate",
        "dox_algal_resp_rate",
        "dox_balgae_photo_rate",
        "dox_balgae_resp_rate",
    )

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
    ) -> None:
        Process.__init__(self, time_step)

        # Lazy-load DEFAULTS by composing from the relevant v3 parameter
        # groups. Mirrors the N2 Process pattern (this is the second
        # composing-Process in v3 NSM1 after N2).
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.dox import (
                DEFAULTS as DOX_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.global_parameters import (
                DEFAULTS as GLOBAL_PARAM_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.global_vars import (
                DEFAULTS as GLOBAL_VAR_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.algae import (
                DEFAULTS as ALGAE_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.balgae import (
                DEFAULTS as BALGAE_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.nitrogen import (
                DEFAULTS as NITROGEN_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.carbon import (
                DEFAULTS as CARBON_DEFAULTS,
            )

            composed: dict[str, float | int | bool] = {}
            # DOX-specific kinetics (DOX is the canonical home of these).
            composed.update(DOX_DEFAULTS)
            # Atmospheric pressure + use_* feature flags.
            composed["pressure_mb"] = GLOBAL_PARAM_DEFAULTS["pressure_mb"]
            for k in (
                "use_NH4",
                "use_DOC",
                "use_Algae",
                "use_Balgae",
                "use_DOX",
            ):
                composed[k] = GLOBAL_PARAM_DEFAULTS.get(k, True)
            # Hydraulic forcings (toy values; should be overridden per
            # cell in production).
            for k in (
                "velocity",
                "flow",
                "topwidth",
                "slope",
                "shear_velocity",
                "wind_speed",
            ):
                composed[k] = GLOBAL_VAR_DEFAULTS[k]
            # Algal stoichiometry. v3 derives ``rca = AWc / AWa`` and
            # ``rcb = BWc / BWd`` at run time (Phase 9.B audit fix; the
            # raw weights ``AWc``/``BWc`` are 40 mg-C while ``rca`` is
            # 0.04 mg-C/ug-Chla and ``rcb`` is 0.4 mg-C/mg-D — using the
            # raw weights as ratios under-stated O2 photosynthesis /
            # respiration coupling by 100-1000x).
            composed["AWc"] = ALGAE_DEFAULTS["AWc"]   # mg-C raw weight
            composed["AWa"] = ALGAE_DEFAULTS["AWa"]   # ug-Chla algal unit
            # Benthic algal stoichiometry and the bottom-area / wet-
            # fraction split.
            composed["BWc"] = BALGAE_DEFAULTS["BWc"]  # mg-C / g-D raw
            composed["BWd"] = BALGAE_DEFAULTS["BWd"]  # mg-D / g-D raw
            composed["Fb"] = BALGAE_DEFAULTS.get("Fb", 0.9)
            composed["Fw"] = BALGAE_DEFAULTS.get("Fw", 0.9)
            # Nitrogen kinetics for the local nitrification flux fallback.
            composed["KNR"] = NITROGEN_DEFAULTS["KNR"]
            composed["knit_20"] = NITROGEN_DEFAULTS["knit_20"]
            composed["knit_theta"] = NITROGEN_DEFAULTS["knit_theta"]
            # Carbon stoichiometry.
            composed["roc"] = CARBON_DEFAULTS["roc"]

            type(self).DEFAULTS = composed

        # Merge user overrides over the composed defaults.
        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "DOX: unknown parameter %r in 'parameters' dict; ignoring "
                "(not in DOX_DEFAULTS / global_parameters / global_vars / "
                "algae / balgae / nitrogen / carbon).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Step-scoped caches for diagnostics / future Phase 5.5 work.
        self.dox_sat: ArrayLike = 0.0
        self.atm_reaeration_rate: ArrayLike = 0.0
        self.dox_nitrification_rate: ArrayLike = 0.0
        self.dox_sod_rate: ArrayLike = 0.0
        # The total dDOX/dt rate at the current step (mg/L/d). Useful
        # for downstream debugging and for Phase 5.5 semi-implicit
        # opt-in (which would split this into source / sink halves).
        self.dox_rate: ArrayLike = 0.0

        # NSM1-DOX-F1 (gold-standard spec C3): one-time warn-latch for
        # the freshwater DO-saturation assumption (the APHA salinity
        # correction is not applied; see ``dox_sat_apha`` docstring and
        # audit C6). Tripped if a nonzero ``salinity`` is in the registry.
        self._salinity_freshwater_warned: bool = False

        # NSM1-DOX-F2 (gold-standard spec C4): one-time warn-latch for
        # the silent-zero atmospheric-reaeration path
        # (hydraulic_reaeration_option == 1 with kah_20_user == 0.0 and
        # no wind path). CE-QUAL-W2 enforces a MINKL minimum-reaeration
        # floor on every branch; v3 keeps the unfloored value by default
        # for v1 parity but warns and offers the opt-in
        # ``min_reaeration_ka`` floor.
        self._reaeration_silent_zero_warned: bool = False

        # Diagnostics handle (replaced by the model's run-level handle in
        # ``init_process`` if available).
        self.diagnostics = Diagnostics()

        # Coupling flags. Defaulted here so ``run`` works without an
        # explicit ``init_process`` call (e.g. Tier 1 unit tests).
        self.use_floating_algae: bool = False
        self.use_benthic_algae: bool = False
        self.use_nitrogen: bool = False
        self.use_carbon: bool = False
        self.use_cbod: bool = False
        self.floating_algae_process = None
        self.benthic_algae_process = None
        self.nitrogen_process = None
        self.carbon_process = None
        self.cbod_process = None

    @ProcessFactory.register("dox")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "DOX":
        return DOX(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """Wire up sibling-process couplings and capture run-level diagnostics."""
        # Capture run-level Diagnostics if the v3 Model exposes one.
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

        # Optional couplings — every sibling Process is detected via
        # ``has_process``. v2's Model also exposes this method so the
        # same code path serves both.
        has_proc = getattr(model, "has_process", None)
        if has_proc is None:
            return

        if model.has_process("FloatingAlgae"):
            self.use_floating_algae = True
            self.floating_algae_process = model.get_process("FloatingAlgae")
        if model.has_process("BenthicAlgae"):
            self.use_benthic_algae = True
            self.benthic_algae_process = model.get_process("BenthicAlgae")
        if model.has_process("Nitrogen"):
            self.use_nitrogen = True
            self.nitrogen_process = model.get_process("Nitrogen")
        if model.has_process("Carbon"):
            self.use_carbon = True
            self.carbon_process = model.get_process("Carbon")
        if model.has_process("CBOD"):
            self.use_cbod = True
            self.cbod_process = model.get_process("CBOD")

        # Phase H-8 (2026-05-21): adopt the Temperature module's
        # ``wind_input_height`` if Temperature is in the model, so the
        # internal log-law correction inside ``kaw_20`` agrees with
        # Temperature's wind-function. Without this, a user passing
        # ``Temperature(wind_input_height=10.0)`` for KSLE ASOS wind
        # would have Temperature apply the log-law correction
        # downstream to 2 m, but DOX's reaeration would treat the
        # same 10-m wind as if it were at 2 m, applying an additional
        # ``(10/2)**0.143 == 1.35`` factor on top -- silently
        # inflating DOX gas exchange by ~35%. With this lookup, both
        # modules see the same convention.
        self.wind_input_height = 2.0  # legacy default
        if model.has_process("Temperature"):
            temperature_process = model.get_process("Temperature")
            self.wind_input_height = float(
                getattr(temperature_process, "wind_input_height", 2.0)
            )

    # ------------------------------------------------------------------
    # Per-source helpers (each returns mg-O2/L/d; missing siblings -> 0)
    # ------------------------------------------------------------------

    def _atm_reaeration_flux(
        self,
        dox: ArrayLike,
        dox_sat: ArrayLike,
        ka_tc_value: ArrayLike,
    ) -> ArrayLike:
        """Atmospheric O2 reaeration flux (mg/L/d). v1 ``Atm_O2_reaeration``."""
        return ka_tc_value * (dox_sat - dox)

    def _floating_algae_growth_flux(self) -> ArrayLike:
        """Floating-algae photosynthesis O2 source (mg/L/d). v1 ``DOX_ApGrowth``.

        ``ApGrowth * rca * roc * (138/106 - 32/106 * ApUptakeFr_NH4)``.

        Reads the cached ``algal_growth_rate`` (ug-Chla/L/d) and
        ``algal_nh4_uptake_fraction`` written by FloatingAlgae.run.
        Returns 0 when FloatingAlgae is not wired up or ``use_Algae`` is
        False.
        """
        if not (self.use_floating_algae and self.use_Algae):
            return 0
        if self.floating_algae_process is None:
            return 0
        ap_growth = getattr(
            self.floating_algae_process, "algal_growth_rate", 0
        )
        ap_uptake_fr_nh4 = getattr(
            self.floating_algae_process, "algal_nh4_uptake_fraction", 0.5
        )
        # rca: algal C:Chla stoichiometric ratio = AWc / AWa
        # (mg-C/ug-Chla). algal_growth_rate is in ug-Chla/L/d, so rca *
        # algal_growth_rate is in mg-C/L/d. Phase 9.B audit fix: prior
        # v3 used ``self.AWc`` directly (40 mg-C) instead of the derived
        # ratio (0.04 mg-C/ug-Chla); Fortran ``modDOX.f90:135`` and v1
        # ``processes.py:2942-2959`` derive the ratio.
        rca = self.AWc / self.AWa
        return (
            ap_growth
            * rca
            * self.roc
            * (138.0 / 106.0 - 32.0 / 106.0 * ap_uptake_fr_nh4)
        )

    def _floating_algae_respiration_flux(self) -> ArrayLike:
        """Floating-algae respiration O2 sink (mg/L/d). v1 ``DOX_ApRespiration``.

        ``ApRespiration * rca * roc``. Reads cached
        ``algal_respiration_rate`` (ug-Chla/L/d).
        """
        if not (self.use_floating_algae and self.use_Algae):
            return 0
        if self.floating_algae_process is None:
            return 0
        ap_resp = getattr(
            self.floating_algae_process, "algal_respiration_rate", 0
        )
        # rca = AWc / AWa per Phase 9.B audit (Fortran modDOX.f90:136 +
        # v1 processes.py:2962-2977). Prior v3 used ``self.AWc`` directly.
        rca = self.AWc / self.AWa
        return ap_resp * rca * self.roc

    def _benthic_algae_growth_flux(self, depth: ArrayLike) -> ArrayLike:
        """Benthic-algae photosynthesis O2 source (mg/L/d). v1 ``DOX_AbGrowth``.

        ``(138/106 - 32/106 * AbUptakeFr_NH4) * roc * rcb * AbGrowth * Fb / depth``.

        Reads cached ``balgae_growth_rate`` and ``balgae_nh4_uptake_fraction``.
        """
        if not (self.use_benthic_algae and self.use_Balgae):
            return 0
        if self.benthic_algae_process is None:
            return 0
        ab_growth = getattr(
            self.benthic_algae_process, "balgae_growth_rate", 0
        )
        ab_uptake_fr_nh4 = getattr(
            self.benthic_algae_process, "balgae_nh4_uptake_fraction", 0.5
        )
        # rcb: benthic algae C-to-dry-weight stoichiometric ratio =
        # BWc / BWd (mg-C/mg-D). Phase 9.B audit fix: prior v3 used
        # ``self.BWc`` directly (40 mg-C/g-D raw weight) instead of the
        # derived ratio (0.4 mg-C/mg-D); Fortran
        # ``modDOX.f90:143`` and v1 ``processes.py:3032-3054`` derive
        # the ratio.
        rcb = self.BWc / self.BWd
        return (
            (138.0 / 106.0 - 32.0 / 106.0 * ab_uptake_fr_nh4)
            * self.roc
            * rcb
            * ab_growth
            * self.Fb
            / depth
        )

    def _benthic_algae_respiration_flux(self, depth: ArrayLike) -> ArrayLike:
        """Benthic-algae respiration O2 sink (mg/L/d). v1 ``DOX_AbRespiration``.

        ``roc * rcb * AbRespiration * Fb / depth``.
        """
        if not (self.use_benthic_algae and self.use_Balgae):
            return 0
        if self.benthic_algae_process is None:
            return 0
        ab_resp = getattr(
            self.benthic_algae_process, "balgae_respiration_rate", 0
        )
        # rcb = BWc / BWd per Phase 9.B audit (Fortran modDOX.f90:144 +
        # v1 processes.py:3057-3078). Prior v3 used ``self.BWc`` directly.
        rcb = self.BWc / self.BWd
        return self.roc * rcb * ab_resp * self.Fb / depth

    def _nitrification_flux(
        self,
        ammonium: ArrayLike,
        t_water_c: ArrayLike,
        dox: ArrayLike,
    ) -> ArrayLike:
        """Nitrification O2 sink (mg/L/d). v1 ``DOX_Nitrification`` (line 2997).

        ``ron * nitrification_flux_rate``.

        Integration Item 1 (registry rate-variable convention, spec
        resolved Q10): when the Nitrogen Process is wired up, read the
        step-scoped ``nitrification_flux_rate`` (mg-N/L/d, positive
        magnitude) cached by ``Nitrogen.run`` and multiply by ``ron``
        (O2 mass per N nitrified, default 4.57). When Nitrogen is not
        wired the sink is identically zero, even if NH4 is in the
        registry.

        Returns 0 when ``use_NH4`` is False or Nitrogen is not wired
        (regardless of whether NH4 happens to be present).
        """
        if not self.use_NH4:
            return 0
        if not (self.use_nitrogen and self.nitrogen_process is not None):
            return 0
        nitrification_flux = getattr(
            self.nitrogen_process, "nitrification_flux_rate", 0
        )
        if nitrification_flux is None:
            nitrification_flux = 0
        # ``ron`` from DOX_DEFAULTS = 4.5714... (= 2 * 32 / 14).
        return self.ron * nitrification_flux

    def _doc_oxidation_flux(self) -> ArrayLike:
        """DOC -> DIC oxidation O2 sink (mg/L/d). v1 ``DOX_DOC_oxidation``.

        ``roc * doc_dic_oxidation_rate`` from the Carbon process. The v1
        function (line 3013) multiplies by ``roc`` (carbon-to-O2 ratio).

        Returns 0 when Carbon is not wired up or ``use_DOC`` is False.
        Phase 5.A may not yet be complete — the ``getattr`` default of
        0 covers the case where ``Carbon.run`` has not yet populated
        ``doc_dic_oxidation_rate`` for this step.
        """
        if not (self.use_carbon and self.use_DOC):
            return 0
        if self.carbon_process is None:
            return 0
        doc_dic = getattr(self.carbon_process, "doc_dic_oxidation_rate", 0)
        return self.roc * doc_dic

    def _cbod_oxidation_flux(self) -> ArrayLike:
        """CBOD oxidation O2 sink (mg/L/d). v1 ``DOX_CBOD_oxidation``.

        Returns ``cbod_oxidation_rate`` directly because by definition
        1 mg-CBOD == 1 mg-O2 (CBOD is an oxygen-demand quantity).
        """
        if not self.use_cbod:
            return 0
        if self.cbod_process is None:
            return 0
        return getattr(self.cbod_process, "cbod_oxidation_rate", 0)

    def _sod_flux(
        self,
        t_water_c: ArrayLike,
        depth: ArrayLike,
        dox: ArrayLike,
    ) -> ArrayLike:
        """Sediment oxygen demand O2 sink (mg/L/d). v1 ``DOX_SOD``.

        ``SOD_tc * DOX / (DOX + KsSOD) / depth``, where ``SOD_tc`` is
        Arrhenius-corrected ``SOD_20`` (g-O2/m^2/d). Dividing by
        ``depth`` (m) yields mg-O2/L/d (since 1 g/m^2/d / 1 m ==
        1 g/m^3/d == 1 mg/L/d).

        Phase 9.B audit fix (C2): the DOX-Monod attenuation is applied
        here in the DOX consumer, not inside the
        ``utils.sediment.SOD_tc`` primitive (Phase 1.1 made that
        primitive pure Arrhenius for architectural reasons). Matches
        Fortran ``modGlobalParam.f90:254``
        (``SOD_tc *= DOX / (DOX + KsSod)`` when ``use_DOX``) and v1
        ``shared.processes.SOD_tc:200``. Under hypoxia (DOX -> 0) the
        attenuation drives the SOD sink to zero, reflecting that the
        sediment cannot deplete oxygen that is not present.
        """
        sod = sod_tc_util(self.SOD_20, self.SOD_theta, t_water_c)
        if self.use_DOX:
            sod = sod * dox / (dox + self.KsSOD)
        return sod / depth

    # ------------------------------------------------------------------
    # Forward-Euler integrator
    # ------------------------------------------------------------------

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Integrate the DOX state by one Forward Euler substep.

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_change_with_components``
        (B); applies Forward Euler with unconditional clip-with-log (C,
        D); persists primary outputs (E); caches step-scoped rates on
        ``self.<name>`` (F); opportunistically writes diagnostics (G).


        Reads (at ``t = time``):
        * ``oxygen_dissolved`` (mg-O2/L)
        * ``water_temperature`` (deg C)
        * ``depth`` (m)
        * ``ammonium`` (mg-N/L) — optional; falls back to 0 when absent.
        * ``atmospheric_pressure`` (mb) — optional; falls back to
          ``self.pressure_mb``.

        Writes (at ``t = time``):
        * ``oxygen_dissolved`` (mg-O2/L) — in-place state update.
        * Each Appendix A diagnostic in ``REGISTRY_DIAGNOSTICS`` — only
          if the user has pre-registered it (pattern G).
        """
        # --- State and forcing reads (pattern A) ---
        dox = registry.get_at_time("oxygen_dissolved", time)
        t_water_c = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)

        # Optional inputs.
        if "ammonium" in registry:
            ammonium = registry.get_at_time("ammonium", time)
        else:
            # Zero-NH4 fallback: nitrification flux is 0 regardless of DOX.
            ammonium = xr.zeros_like(dox) if hasattr(dox, "dims") else 0.0
        if "atmospheric_pressure" in registry:
            pressure_mb = registry.get_at_time("atmospheric_pressure", time)
        else:
            pressure_mb = self.pressure_mb

        # Phase G-3 (2026-05-21): read hydraulic and wind forcings from
        # the registry when present, falling back to the constructor-
        # time scalars (which are toy defaults from
        # ``parameters/global_vars.py``). Prior to this, reaeration
        # used the frozen constructor scalars even when the registry
        # carried time-varying ``wind_speed`` / ``velocity`` / ``flow``
        # / ``topwidth`` / ``slope`` / ``shear_velocity`` --
        # silently producing wrong gas-exchange rates in any coupled
        # run with realistic met forcing. The Temperature module
        # already reads ``wind_speed`` from the registry, so
        # pre-G-3 the heat budget responded to wind while DOX did
        # not -- a confusing internal asymmetry that this fix closes.
        wind_speed = (
            registry.get_at_time("wind_speed", time)
            if "wind_speed" in registry else self.wind_speed
        )
        # Wind shelter (wind_shelter_coefficient): the same optional per-cell
        # forcing the TSM wind function consumes. Threaded into kaw_20 so a
        # sheltered cell gets reduced wind for gas transfer as well as heat
        # exchange (CE-QUAL-W2 applies WSC once to drive both). Absent -> 1.0
        # (no shelter), preserving prior numerical output.
        wind_shelter = (
            registry.get_at_time("wind_shelter_coefficient", time)
            if "wind_shelter_coefficient" in registry else 1.0
        )
        velocity = (
            registry.get_at_time("velocity", time)
            if "velocity" in registry else self.velocity
        )
        flow = (
            registry.get_at_time("flow", time)
            if "flow" in registry else self.flow
        )
        topwidth = (
            registry.get_at_time("topwidth", time)
            if "topwidth" in registry else self.topwidth
        )
        slope = (
            registry.get_at_time("slope", time)
            if "slope" in registry else self.slope
        )
        shear_velocity = (
            registry.get_at_time("shear_velocity", time)
            if "shear_velocity" in registry else self.shear_velocity
        )

        # --- NSM1-DOX-F1 (spec C3): freshwater DO-saturation guard ---
        # ``dox_sat_apha`` computes the *freshwater* APHA saturation (the
        # salinity correction is intentionally omitted; matches v1 and is
        # exact for fresh water, ~factor 1.0). If the registry carries a
        # nonzero salinity, the freshwater assumption is violated and
        # DOsat is overstated (~18% at 35 ppt). Warn once rather than let
        # brackish input pass silently. The salinity-corrected APHA form
        # is a documented deferral (audit C6 / NSM2 scope); no numeric
        # change here for fresh water.
        if not self._salinity_freshwater_warned and "salinity" in registry:
            _sal = registry.get_at_time("salinity", time)
            if bool(np.any(np.asarray(_sal) > 0.0)):
                logger.warning(
                    "DOX: nonzero 'salinity' present in the registry, but "
                    "the APHA DO-saturation (dox_sat_apha) applies the "
                    "FRESHWATER form only (no salinity correction; "
                    "NSM1-DOX-F1 / audit C6, deferred). DOsat is "
                    "overstated for brackish/marine water (~18%% at 35 "
                    "ppt). This warning is emitted once per process."
                )
                self._salinity_freshwater_warned = True

        # --- Fused rate composition (pattern B) ---
        delta_dox, rate, components = self._change_with_components(
            dox=dox,
            t_water_c=t_water_c,
            depth=depth,
            ammonium=ammonium,
            pressure_mb=pressure_mb,
            wind_speed=wind_speed,
            velocity=velocity,
            flow=flow,
            topwidth=topwidth,
            slope=slope,
            shear_velocity=shear_velocity,
            wind_shelter=wind_shelter,
        )

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        # Preserved names: ``dox_sat``, ``atm_reaeration_rate``,
        # ``dox_nitrification_rate``, ``dox_sod_rate`` were already
        # published by DOX for the Phase 5.5 / Tier 1 / Alkalinity
        # diagnostics. The remaining seven are new with this phase.
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])
        # Net-rate cache retained for downstream debugging (not exposed
        # via REGISTRY_DIAGNOSTICS — it is the integrator's own argument,
        # not an Appendix A diagnostic).
        self.dox_rate = rate

        # --- Forward Euler integration (pattern C) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        dox_new = dox + delta_dox

        # Clip-with-log (pattern D). Phase 0.6 Q1+Q2: ``clip_negative_state``
        # accepts non-DataArray and None-diagnostics inputs; step
        # attribution is automatic via ``diagnostics.current_step``.
        dox_new = clip_negative_state(dox_new, "oxygen_dissolved", self.diagnostics)

        # --- Persist primary state (pattern E) ---
        registry.set_at_time("oxygen_dissolved", time, dox_new)

        # --- Opportunistic diagnostic registry writes (pattern G) ---
        # Each Appendix A name is written ONLY if the user has
        # pre-registered it. Zero cost when not subscribed.
        for name in self.REGISTRY_DIAGNOSTICS:
            if name in registry:
                registry.set_at_time(name, time, components[name])

    # ------------------------------------------------------------------
    # Rate-composition helpers
    # ------------------------------------------------------------------

    def _change_with_components(
        self,
        *,
        dox: ArrayLike,
        t_water_c: ArrayLike,
        depth: ArrayLike,
        ammonium: ArrayLike,
        pressure_mb: ArrayLike,
        wind_speed: ArrayLike | None = None,
        velocity: ArrayLike | None = None,
        flow: ArrayLike | None = None,
        topwidth: ArrayLike | None = None,
        slope: ArrayLike | None = None,
        shear_velocity: ArrayLike | None = None,
        wind_shelter: ArrayLike | float = 1.0,
    ) -> tuple[ArrayLike, ArrayLike, dict]:
        """Compute ``(delta_dox, rate, components)``.

        Code-motion-only refactor of ``run``'s former inline composition
        (§11.6): operand order, intermediate names, and arithmetic are
        preserved verbatim from the pre-refactor body. The ``components``
        dict is populated from the same intermediates the integrator
        consumes (no recomputation).

        ``delta_dox`` is the Forward Euler increment ``rate * dt_days``
        (in mg-O2/L per substep). ``rate`` is the net rate in mg-O2/L/d.
        Both are returned so ``run`` can apply Forward Euler without
        re-multiplying by ``dt_days`` and the caller can cache ``rate``
        on ``self.dox_rate`` for downstream debugging.

        Phase G-3 (2026-05-21): the hydraulic and wind kwargs accept
        registry-driven time-varying forcings from ``run()``. When
        passed as ``None`` (default, for backward compat with callers
        that have not yet been updated), each falls back to the
        constructor-time scalar on ``self``.
        """
        # Phase G-3 kwarg-to-self fallback.
        if wind_speed is None:
            wind_speed = self.wind_speed
        if velocity is None:
            velocity = self.velocity
        if flow is None:
            flow = self.flow
        if topwidth is None:
            topwidth = self.topwidth
        if slope is None:
            slope = self.slope
        if shear_velocity is None:
            shear_velocity = self.shear_velocity
        # --- O2 saturation (APHA / Benson-Krause) ---
        dox_sat = dox_sat_apha(t_water_c, pressure_mb)

        # --- Effective reaeration coefficient ka_tc (1/d) ---
        # Mirrors the N2 fast path: when both menu options are
        # user-defined (==1) and both user values are zero, ``ka_tc``
        # is identically zero and we short-circuit to avoid the
        # numpy.select shape quirk in ``utils.reaeration.kah_20``.
        is_user_hydraulic_zero = (
            self.hydraulic_reaeration_option == 1 and self.kah_20_user == 0.0
        )
        is_user_wind_zero = (
            self.wind_reaeration_option == 1 and self.kaw_20_user == 0.0
        )
        if is_user_hydraulic_zero and is_user_wind_zero:
            ka_tc_value = 0.0
        else:
            kah_20_value = kah_20(
                kah_20_user=self.kah_20_user,
                hydraulic_reaeration_option=self.hydraulic_reaeration_option,
                velocity=velocity,
                depth=depth,
                flow=flow,
                topwidth=topwidth,
                slope=slope,
                shear_velocity=shear_velocity,
            )
            kaw_20_value = kaw_20(
                kaw_20_user=self.kaw_20_user,
                wind_speed=wind_speed,
                wind_reaeration_option=self.wind_reaeration_option,
                wind_input_height=getattr(self, "wind_input_height", 2.0),
                wind_shelter=wind_shelter,
            )
            ka_tc_value = ka_tc(
                kah_20=kah_20_value,
                kaw_20=kaw_20_value,
                kah_theta=self.kah_theta,
                kaw_theta=self.kaw_theta,
                T_water_C=t_water_c,
                depth=depth,
            )

        # --- NSM1-DOX-F2 (spec C4): silent-zero reaeration guard + floor ---
        # ``hydraulic_reaeration_option == 1`` with ``kah_20_user == 0.0``
        # yields zero hydraulic reaeration; if the wind path is also off,
        # atmospheric reaeration is silently zero. CE-QUAL-W2 enforces a
        # ``MINKL`` minimum-reaeration floor on every branch. v3 keeps
        # the unfloored value by default for v1/Fortran parity but
        # (a) warns once, and (b) offers the opt-in ``min_reaeration_ka``
        # floor (default 0.0 = OFF, so default behaviour is unchanged).
        if (
            not self._reaeration_silent_zero_warned
            and is_user_hydraulic_zero
        ):
            logger.warning(
                "DOX: hydraulic_reaeration_option == 1 with "
                "kah_20_user == 0.0 gives zero hydraulic reaeration; "
                "atmospheric reaeration may be silently zero. v3 "
                "preserves v1/Fortran parity (no implicit floor); "
                "CE-QUAL-W2 enforces a MINKL floor. Set "
                "'min_reaeration_ka' > 0 to opt into a minimum-"
                "reaeration floor. NSM1-DOX-F2. Emitted once per process."
            )
            self._reaeration_silent_zero_warned = True
        if self.min_reaeration_ka > 0.0:
            ka_tc_value = np.maximum(ka_tc_value, self.min_reaeration_ka)

        # --- Compute per-source / per-sink fluxes (mg/L/d) ---
        # Each sub-flux is individually sanitized so that a NaN at a single
        # cell in one term cannot poison the entire ``rate`` array via the
        # subsequent sum. Without this, a stale NaN in (e.g.) the cached
        # ``balgae_growth_rate`` at a newly-wet cell — where benthic algae
        # density is 0 and the rate should be 0 — produces ``0 + ... + NaN
        # = NaN`` for the whole cell, which ``sanitize_rate`` then drops to
        # 0, leaving DOX integration as ``dox + 0 = dox`` and freezing the
        # cell at IC indefinitely. Sanitizing per sub-flux preserves the
        # contributions of the well-behaved terms.
        atm_reaer  = sanitize_rate(self._atm_reaeration_flux(dox, dox_sat, ka_tc_value))
        algal_grow = sanitize_rate(self._floating_algae_growth_flux())
        algal_resp = sanitize_rate(self._floating_algae_respiration_flux())
        balgae_grow = sanitize_rate(self._benthic_algae_growth_flux(depth))
        balgae_resp = sanitize_rate(self._benthic_algae_respiration_flux(depth))
        nitr_sink  = sanitize_rate(self._nitrification_flux(ammonium, t_water_c, dox))
        doc_sink   = sanitize_rate(self._doc_oxidation_flux())
        cbod_sink  = sanitize_rate(self._cbod_oxidation_flux())
        sod_sink   = sanitize_rate(self._sod_flux(t_water_c, depth, dox))

        # --- Net rate (mg/L/d). Mirrors v1 ``dDOXdt`` (line 3119) ---
        rate = (
            atm_reaer
            + algal_grow
            - algal_resp
            + balgae_grow
            - balgae_resp
            - nitr_sink
            - doc_sink
            - cbod_sink
            - sod_sink
        )

        # NaN/inf guard (final defense-in-depth, in addition to per-sub-flux
        # sanitization above). Catches any residual NaN from the sum.
        rate = sanitize_rate(rate)

        # --- Forward Euler delta (caller applies as ``dox + delta_dox``) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        delta_dox = rate * dt_days

        # --- Components dict (single source of truth for pattern F+G) ---
        # ``sod_rate`` is exposed as an alias for ``dox_sod_rate`` (the
        # volumetric SOD sink in mg-O2/L/d). Both names map to the
        # same sanitized value.
        components = {
            "dox_sat": dox_sat,
            "atm_reaeration_rate": atm_reaer,
            "dox_nitrification_rate": nitr_sink,
            "dox_sod_rate": sod_sink,
            "sod_rate": sod_sink,
            "dox_doc_oxidation_rate": doc_sink,
            "dox_cbod_oxidation_rate": cbod_sink,
            "dox_algal_photo_rate": algal_grow,
            "dox_algal_resp_rate": algal_resp,
            "dox_balgae_photo_rate": balgae_grow,
            "dox_balgae_resp_rate": balgae_resp,
        }

        return delta_dox, rate, components
