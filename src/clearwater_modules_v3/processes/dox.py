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
  reaeration, temperature-corrected)
* ``O2sat`` -- APHA / QUAL2E formulation; v1 lines 2901-2923. Requires
  ``T_water_K`` (Kelvin), ``pressure_mb`` (mb), and the empirical
  ``DOs_atm_alpha`` correction. The full APHA form is implemented (no
  simpler approximation is used).
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

from clearwater_modules_v2.processes.base import Process, ProcessFactory
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

        Reads (at ``t = time``):
        * ``oxygen_dissolved`` (mg-O2/L)
        * ``water_temperature`` (deg C)
        * ``depth`` (m)
        * ``ammonium`` (mg-N/L) — optional; falls back to 0 when absent.
        * ``atmospheric_pressure`` (mb) — optional; falls back to
          ``self.pressure_mb``.

        Writes (at ``t = time``):
        * ``oxygen_dissolved`` (mg-O2/L) — in-place state update.

        Caches step-scoped quantities on ``self`` for diagnostics:
        * ``dox_sat``, ``atm_reaeration_rate``, ``dox_nitrification_rate``,
          ``dox_sod_rate``, ``dox_rate``.
        """
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
            ka_tc_value = ka_tc(
                kah_20=kah_20_value,
                kaw_20=kaw_20_value,
                kah_theta=self.kah_theta,
                kaw_theta=self.kaw_theta,
                T_water_C=t_water_c,
                depth=depth,
            )

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

        # --- Forward Euler integration ---
        dt_days = self.time_step.total_seconds() / 86400.0
        dox_new = dox + rate * dt_days

        # Clip-with-log per the Q7 contract.
        if isinstance(dox_new, xr.DataArray) and self.diagnostics is not None:
            dox_new = clip_negative_state(
                dox_new, "oxygen_dissolved", self.diagnostics, step=0
            )
        else:
            dox_new = xr.where(dox_new < 0, 0, dox_new)

        # Persist updated state.
        registry.set_at_time("oxygen_dissolved", time, dox_new)

        # Cache step-scoped quantities for diagnostics / Phase 5.5.
        self.dox_sat = dox_sat
        self.atm_reaeration_rate = atm_reaer
        self.dox_nitrification_rate = nitr_sink
        self.dox_sod_rate = sod_sink
        self.dox_rate = rate
