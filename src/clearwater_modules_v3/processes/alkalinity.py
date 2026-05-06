"""v3 NSM1 ``Alkalinity`` Process.

Phase 6 (v3 NSM1 design spec, Section 11 Phase 6, Section 5 Alkalinity
design notes, Section 14 resolved Q). v3-native (NOT a v2 overlay; v2
has no Alkalinity Process).

Per Section 14 resolved Q ("Alkalinity simple-tracer"), v3 1.0.0 treats
Alkalinity as a *simple tracer* with source/sink coupling to Nitrogen
(nitrification consumes alk, denitrification produces alk) and the algae
Processes (algal growth/respiration coupling, fractionated by NH4 vs NO3
uptake). There is **no** carbonate equilibrium / pH solver in v3 1.0.0;
that is NSM2 territory in v3 1.1+.

Kinetics (mirrors v1 ``processes.py:3246-3447``):

    dAlk/dt = + r_alkden * (1 - DOX/(DOX+KsOxdn)) * kdnit_tc * NO3 * 50000   (denit source)
              - r_alkn  * (1 - exp(-KNR*DOX)) * knit_tc * NH4 * 50000        (nitrif sink)
              - (r_alkaa * fNH4 - r_alkan * (1 - fNH4)) * ApGrowth * rca * 50000
                                                                              (algal growth coupling
                                                                               -- net sink under NH4
                                                                               uptake, source under
                                                                               NO3 uptake)
              + r_alkaa * ApRespiration * rca * 50000                         (algal respiration source)
              - (r_alkba * fbNH4 - r_alkbn * (1 - fbNH4)) * AbGrowth * rcb * Fb / depth * 50000
                                                                              (benthic algae growth)
              + r_alkba * AbRespiration * rcb * Fb / depth * 50000            (benthic algae respiration)

In v3 the production-side terms are read from the upstream Processes' cached
*step-scoped flux rates* rather than recomputed locally:

* Nitrogen: ``nitrification_flux_rate`` (mg-N/L/d), ``denitrification_flux_rate``
  (mg-N/L/d). These are absolute-value magnitudes — see Phase 2.B Item 1.
  Multiplying by ``r_alkn * 50000`` and ``r_alkden * 50000`` yields the
  alkalinity change rate in mg-CaCO3/L/d. The DOX-Monod attenuation that
  v1 applied locally inside ``Alk_nitrification`` / ``Alk_denitrification``
  is already baked into the upstream flux (Nitrogen.run computes
  ``ammonium_nitrification`` / ``nitrate_denitrification`` with the same
  Monod term), so we do not re-apply it here.
* FloatingAlgae: ``algal_growth_rate`` (ug-Chla/L/d),
  ``algal_respiration_rate`` (ug-Chla/L/d), ``algal_nh4_uptake_fraction``.
* BenthicAlgae: ``balgae_growth_rate`` (g-D/m^2/d),
  ``balgae_respiration_rate`` (g-D/m^2/d), ``balgae_nh4_uptake_fraction``.

Stoichiometric ratios (from ``ALKALINITY_DEFAULTS``):

* ``r_alkaa`` = 14/106/12/1000 eq/mg-C — algal photosynthesis, NH4 path.
* ``r_alkan`` = 18/106/12/1000 eq/mg-C — algal photosynthesis, NO3 path.
* ``r_alkn``  = 2/14/1000      eq/mg-N — nitrification (Alk sink).
* ``r_alkden`` = 4/14/1000     eq/mg-N — denitrification (Alk source).
* ``r_alkba`` = 14/106/12/1000 eq/mg-C — benthic algae photosynthesis, NH4 path.
* ``r_alkbn`` = 18/106/12/1000 eq/mg-C — benthic algae photosynthesis, NO3 path.

The literal factor of ``50000`` converts eq/L to mg-CaCO3/L (CaCO3
equivalent weight = 50 g/eq = 50000 mg/eq).

DEFAULTS sources (composed from four parameter groups):

* ``parameters.alkalinity``: 6 stoichiometric ratios (above).
* ``parameters.global_parameters``: feature flags (``use_NH4``, ``use_NO3``,
  ``use_DOX``, ``use_Algae``, ``use_Balgae``).
* ``parameters.algae``: ``AWc`` (== rca, mg-C/ug-Chla).
* ``parameters.balgae``: ``BWc`` (== rcb, mg-C/g-D), ``Fb`` (bottom-area
  fraction).

Coupling pattern: every cross-Process value is read via
``getattr(..., default=0)`` so missing siblings degrade gracefully.

Per design spec Section 11, this Process **is** wired into the package
``__init__`` here (registration under Phase 6 acceptance criteria).
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

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


logger = logging.getLogger(__name__)


# Conversion: eq/L -> mg-CaCO3/L (CaCO3 equivalent weight = 50 g/eq).
EQ_TO_MG_CACO3: float = 50000.0


# ---------------------------------------------------------------------------
# Alkalinity Process
# ---------------------------------------------------------------------------


class Alkalinity(Process):
    """v3 NSM1 alkalinity Process (simple-tracer model).

    State variable: ``alkalinity`` (mg-CaCO3/L).

    Sources / sinks (mg-CaCO3/L/d, integrated by Forward Euler):

    * Denitrification (NO3 -> N2):           Alk source.
    * Nitrification (NH4 -> NO3):            Alk sink.
    * Algal photosynthesis:                  Alk sink under NH4 uptake,
      Alk source under NO3 uptake; net depends on
      ``algal_nh4_uptake_fraction``.
    * Algal respiration:                     Alk source (DIC release).
    * Benthic-algae photosynthesis & respiration: same form as floating
      algae, depth-integrated and bottom-area-fractionated.

    Caches step-scoped quantities for diagnostics:

    * ``self.alk_nitrification_rate`` -- Alk sink from nitrification (mg-CaCO3/L/d)
    * ``self.alk_denitrification_rate`` -- Alk source from denitrification (mg-CaCO3/L/d)
    * ``self.alk_algal_growth_rate`` -- Alk source/sink from algal growth (mg-CaCO3/L/d)
    * ``self.alk_algal_respiration_rate`` -- Alk source from algal respiration (mg-CaCO3/L/d)
    * ``self.alk_benthic_algae_growth_rate`` -- benthic algae growth (mg-CaCO3/L/d)
    * ``self.alk_benthic_algae_respiration_rate`` -- benthic algae respiration (mg-CaCO3/L/d)
    * ``self.alk_rate`` -- net dAlk/dt at the current step (mg-CaCO3/L/d)
    """

    variables = [
        "alkalinity",
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
        # groups. Mirrors the Phase 5.B DOX multi-group composition.
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.alkalinity import (
                DEFAULTS as ALKALINITY_DEFAULTS,
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
            # Alkalinity-specific stoichiometric ratios.
            composed.update(ALKALINITY_DEFAULTS)
            # Feature flags.
            for k in (
                "use_NH4",
                "use_NO3",
                "use_DOX",
                "use_Algae",
                "use_Balgae",
            ):
                composed[k] = GLOBAL_PARAM_DEFAULTS.get(k, True)
            # Algal stoichiometry (rca = AWc, mg-C/ug-Chla).
            composed["AWc"] = ALGAE_DEFAULTS["AWc"]
            # Benthic algae stoichiometry (rcb = BWc, mg-C/g-D) and
            # bottom-area fraction.
            composed["BWc"] = BALGAE_DEFAULTS["BWc"]
            composed["Fb"] = BALGAE_DEFAULTS.get("Fb", 0.9)

            type(self).DEFAULTS = composed

        # Merge user overrides over the composed defaults.
        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "Alkalinity: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in ALKALINITY_DEFAULTS / global_parameters / "
                "algae / balgae).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Step-scoped caches for diagnostics.
        self.alk_nitrification_rate: ArrayLike = 0.0
        self.alk_denitrification_rate: ArrayLike = 0.0
        self.alk_algal_growth_rate: ArrayLike = 0.0
        self.alk_algal_respiration_rate: ArrayLike = 0.0
        self.alk_benthic_algae_growth_rate: ArrayLike = 0.0
        self.alk_benthic_algae_respiration_rate: ArrayLike = 0.0
        self.alk_rate: ArrayLike = 0.0

        # Diagnostics handle (replaced by the model's run-level handle in
        # ``init_process`` if available).
        self.diagnostics = Diagnostics()

        # Coupling flags. Defaulted here so ``run`` works without an
        # explicit ``init_process`` call (e.g. Tier 1 unit tests).
        self.use_floating_algae: bool = False
        self.use_benthic_algae: bool = False
        self.use_nitrogen: bool = False
        self.floating_algae_process = None
        self.benthic_algae_process = None
        self.nitrogen_process = None

    @ProcessFactory.register("alkalinity")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "Alkalinity":
        return Alkalinity(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """Wire up sibling-process couplings and capture run-level diagnostics."""
        # Capture run-level Diagnostics if the v3 Model exposes one.
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

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

    # ------------------------------------------------------------------
    # Per-source helpers (each returns mg-CaCO3/L/d; missing siblings -> 0)
    # ------------------------------------------------------------------

    def _nitrification_alk_sink(self) -> ArrayLike:
        """Alk sink from nitrification (mg-CaCO3/L/d). v1 ``Alk_nitrification``.

        ``r_alkn * nitrification_flux_rate * 50000``.

        The v1 routine applies a ``(1 - exp(-KNR*DOX))`` Monod-style
        oxygen attenuation locally; in v3 that attenuation is already
        baked into ``Nitrogen.ammonium_nitrification`` (which the
        Nitrogen Process calls in ``run`` to populate
        ``nitrification_flux_rate``), so we do NOT re-apply it.

        Returns 0 when ``use_NH4`` is False, or when Nitrogen is not
        wired up.
        """
        if not self.use_NH4:
            return 0
        if not (self.use_nitrogen and self.nitrogen_process is not None):
            return 0
        nitr_flux = getattr(
            self.nitrogen_process, "nitrification_flux_rate", 0
        )
        if nitr_flux is None:
            nitr_flux = 0
        return self.r_alkn * nitr_flux * EQ_TO_MG_CACO3

    def _denitrification_alk_source(self) -> ArrayLike:
        """Alk source from denitrification (mg-CaCO3/L/d). v1 ``Alk_denitrification``.

        ``r_alkden * denitrification_flux_rate * 50000``.

        v1 applies a ``(1 - DOX/(DOX+KsOxdn))`` oxygen-inhibition Monod
        factor locally; in v3 that factor is already baked into
        ``Nitrogen.nitrate_denitrification`` (called in Nitrogen.run to
        populate ``denitrification_flux_rate``), so we do NOT re-apply
        it here.

        Returns 0 when ``use_NO3`` is False, or when Nitrogen is not
        wired up.
        """
        if not self.use_NO3:
            return 0
        if not (self.use_nitrogen and self.nitrogen_process is not None):
            return 0
        denit_flux = getattr(
            self.nitrogen_process, "denitrification_flux_rate", 0
        )
        if denit_flux is None:
            denit_flux = 0
        return self.r_alkden * denit_flux * EQ_TO_MG_CACO3

    def _floating_algae_growth_alk_flux(self) -> ArrayLike:
        """Floating-algae growth coupling (mg-CaCO3/L/d). v1 ``Alk_algal_growth``.

        Returns the *net Alk SINK* contribution from algal growth (i.e.
        the term that is *subtracted* from ``dAlk/dt``):

            (r_alkaa * fNH4 - r_alkan * (1 - fNH4)) * ApGrowth * rca * 50000

        With all NH4 uptake (fNH4 == 1) the term is + r_alkaa * ... ; this
        is *subtracted* from dAlk/dt -> net Alk sink (correct: NH4 uptake
        consumes alk).

        With all NO3 uptake (fNH4 == 0) the term is - r_alkan * ... ;
        subtracting yields + r_alkan * ... in dAlk/dt -> net Alk source
        (correct: NO3 uptake produces alk).

        Reads cached ``algal_growth_rate`` (ug-Chla/L/d) and
        ``algal_nh4_uptake_fraction`` from FloatingAlgae.

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
        # rca: algal C:Chla. ApGrowth (ug-Chla/L/d) * rca (mg-C/ug-Chla)
        # = mg-C/L/d. Multiplied by stoich (eq/mg-C) and 50000 (mg-CaCO3/eq)
        # yields mg-CaCO3/L/d.
        rca = self.AWc
        return (
            (self.r_alkaa * ap_uptake_fr_nh4
             - self.r_alkan * (1.0 - ap_uptake_fr_nh4))
            * ap_growth
            * rca
            * EQ_TO_MG_CACO3
        )

    def _floating_algae_respiration_alk_source(self) -> ArrayLike:
        """Floating-algae respiration Alk source (mg-CaCO3/L/d). v1 ``Alk_algal_respiration``.

        ``ApRespiration * r_alkaa * rca * 50000`` -- algal respiration
        always produces alk (DIC release).

        Reads cached ``algal_respiration_rate`` (ug-Chla/L/d).
        """
        if not (self.use_floating_algae and self.use_Algae):
            return 0
        if self.floating_algae_process is None:
            return 0
        ap_resp = getattr(
            self.floating_algae_process, "algal_respiration_rate", 0
        )
        return ap_resp * self.r_alkaa * self.AWc * EQ_TO_MG_CACO3

    def _benthic_algae_growth_alk_flux(self, depth: ArrayLike) -> ArrayLike:
        """Benthic-algae growth coupling (mg-CaCO3/L/d). v1 ``Alk_benthic_algae_growth``.

        Returns the *net Alk SINK* contribution (sign convention as for
        floating algae):

            (1/depth) * (r_alkba * fbNH4 - r_alkbn * (1-fbNH4))
            * AbGrowth * Fb * rcb * 50000

        AbGrowth is in g-D/m^2/d; rcb (mg-C/g-D) converts to mg-C/m^2/d;
        dividing by depth (m) gives mg-C/m^3/d == mg-C/L/d (per
        DOX coupling docstring).
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
        rcb = self.BWc
        return (
            (1.0 / depth)
            * (self.r_alkba * ab_uptake_fr_nh4
               - self.r_alkbn * (1.0 - ab_uptake_fr_nh4))
            * ab_growth
            * self.Fb
            * rcb
            * EQ_TO_MG_CACO3
        )

    def _benthic_algae_respiration_alk_source(
        self, depth: ArrayLike
    ) -> ArrayLike:
        """Benthic-algae respiration Alk source (mg-CaCO3/L/d).
        v1 ``Alk_benthic_algae_respiration``.

        ``(1/depth) * r_alkba * AbRespiration * rcb * Fb * 50000``.
        """
        if not (self.use_benthic_algae and self.use_Balgae):
            return 0
        if self.benthic_algae_process is None:
            return 0
        ab_resp = getattr(
            self.benthic_algae_process, "balgae_respiration_rate", 0
        )
        return (
            (1.0 / depth)
            * self.r_alkba
            * ab_resp
            * self.BWc
            * self.Fb
            * EQ_TO_MG_CACO3
        )

    # ------------------------------------------------------------------
    # Forward-Euler integrator
    # ------------------------------------------------------------------

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Integrate the alkalinity state by one Forward Euler substep.

        Reads (at ``t = time``):
        * ``alkalinity`` (mg-CaCO3/L)
        * ``water_temperature`` (deg C) -- not used directly here but
          required to be present in the registry for symmetry with
          other Processes.
        * ``depth`` (m) -- used to depth-integrate benthic algae
          contributions.

        Writes (at ``t = time``):
        * ``alkalinity`` (mg-CaCO3/L) -- in-place state update.

        Caches step-scoped quantities on ``self`` for diagnostics.
        """
        alk = registry.get_at_time("alkalinity", time)
        # water_temperature is read for symmetry; the simple-tracer
        # model has no direct temperature dependence in v3 1.0.0
        # (all temperature-correction lives upstream in Nitrogen / algae).
        _ = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)

        # --- Compute per-source / per-sink fluxes (mg-CaCO3/L/d) ---
        nitr_sink = self._nitrification_alk_sink()
        denit_source = self._denitrification_alk_source()
        algal_growth_sink = self._floating_algae_growth_alk_flux()
        algal_resp_source = self._floating_algae_respiration_alk_source()
        balgae_growth_sink = self._benthic_algae_growth_alk_flux(depth)
        balgae_resp_source = self._benthic_algae_respiration_alk_source(depth)

        # --- Net rate (mg-CaCO3/L/d). Mirrors v1 ``dAlkdt`` (line 3431):
        #   denit - nitrif - algal_growth + algal_resp
        #          - benthic_growth + benthic_resp
        # ``algal_growth`` and ``benthic_growth`` here carry the v1 sign
        # convention (sink-positive); subtracting them lets the NH4-vs-NO3
        # split flip the sign correctly.
        rate = (
            denit_source
            - nitr_sink
            - algal_growth_sink
            + algal_resp_source
            - balgae_growth_sink
            + balgae_resp_source
        )

        # NaN/inf guard (defense-in-depth; primary dry-cell defense
        # is the orchestration-layer wet-mask in Model). Catches
        # ``inf`` from depth divisions and ``NaN`` from missing
        # forcings.
        rate = sanitize_rate(rate)

        # --- Forward Euler integration ---
        dt_days = self.time_step.total_seconds() / 86400.0
        alk_new = alk + rate * dt_days

        # Clip-with-log per the Q7 contract.
        if isinstance(alk_new, xr.DataArray) and self.diagnostics is not None:
            alk_new = clip_negative_state(
                alk_new, "alkalinity", self.diagnostics, step=0
            )
        else:
            alk_new = xr.where(alk_new < 0, 0, alk_new)

        # Persist updated state.
        registry.set_at_time("alkalinity", time, alk_new)

        # Cache step-scoped quantities for diagnostics.
        self.alk_nitrification_rate = nitr_sink
        self.alk_denitrification_rate = denit_source
        self.alk_algal_growth_rate = algal_growth_sink
        self.alk_algal_respiration_rate = algal_resp_source
        self.alk_benthic_algae_growth_rate = balgae_growth_sink
        self.alk_benthic_algae_respiration_rate = balgae_resp_source
        self.alk_rate = rate
