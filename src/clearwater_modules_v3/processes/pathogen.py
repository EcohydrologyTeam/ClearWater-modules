"""v3 NSM1 Pathogen (PX) Process.

Phase 3.1 (v3 NSM1 design spec, Section 11 + Section 5 Pathogen design
notes): direct port of v1's pathogen kinetics as a v3-native Process
(no v2 overlay exists). Pathogen is an independent constituent — it has
no coupling to other Phase 3 Processes — so the implementation only
reads its own state plus the standard environmental forcings (water
temperature, depth, solar radiation) and any optional light-extinction
contributions (Solid / POC / Ap) the registry may carry.

Kinetics (v1 ``processes.py`` lines 3141-3227):

    dPX/dt = - kdx_tc * PX                                       (natural decay)
             - apx * I0 * (1 - exp(-KEXT * depth)) / (KEXT * depth) * PX
                                                                 (light-induced decay)
             - vx / depth * PX                                   (settling)

Where:

* ``kdx_tc = arrhenius_correction(kdx_20, kdx_theta, T)`` — temperature-
  corrected natural decay rate (1/d).
* ``apx`` — pathogen sunlight-inactivation efficiency
  ((W/m^2)^-1 d^-1). v3 default ``apx = 0.017`` is the Auer & Niehaus
  (1993) / Chapra (1997) canonical value (Phase 9.F.B; was a 1.0
  placeholder in v1/v3 pre-9.F.B).
* ``I0`` — surface broadband solar irradiance ``q_solar`` (W/m^2).
  Phase 9.F.B reverted the earlier Phase 3.1 substitution
  ``I0 = q_solar * Fr_PAR`` because pathogen inactivation is largely
  UVA/UVB-mediated rather than PAR-mediated; the canonical
  Auer/Niehaus/Chapra/QUAL2K formulation operates on total broadband
  ``q_solar``.
* ``KEXT`` — Beer-Lambert composite extinction coefficient (1/m),
  computed via ``utils.light.L``. Inputs (``Solid``, ``POC``, ``Ap``)
  are read from the registry when present; missing inputs default to
  zero with a one-time warning so the Pathogen Process can run
  stand-alone for the Tier 1 conservation harness.
* ``vx`` — pathogen net settling velocity (m/d). v3 default
  ``vx = 1.38`` is the Auer & Niehaus (1993) / Chapra (1997) canonical
  value (Phase 9.F.B; was a 1.0 placeholder in v1/v3 pre-9.F.B).

State variable name: ``pathogen`` (matches Tier 1 conftest fixture and
constituent diff Section 4 / Appendix A).

Per the resolved Q7 clip-with-log contract, the post-Forward-Euler
state is passed through ``clip_negative_state`` and any clip events
are recorded on the run-level ``Diagnostics`` handle (resolved Q10:
state at ``time=t_current`` is the Jacobi state).
"""

from __future__ import annotations

from datetime import datetime, timedelta
import logging
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from clearwater_data.custom_types import ArrayLike
from clearwater_data.variables import VariableRegistry

from clearwater_modules_v3.utils.conversions import arrhenius_correction
from clearwater_modules_v3.processes.base import Process, ProcessFactory
from clearwater_modules_v3.utils.light import L
from clearwater_modules_v3.utils.numerics import (
    Diagnostics,
    clip_negative_state,
    sanitize_rate,
)

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


logger = logging.getLogger(__name__)


# Inline fallback defaults for non-pathogen parameters that participate in
# the light-extinction calculation but live in other parameter groups
# (global_vars / global_parameters / algae). These are read from the merged
# DEFAULTS dict if the user supplied them, otherwise from this fallback.
# The default values match v3 ``global_vars.DEFAULTS`` /
# ``global_parameters.DEFAULTS`` / ``algae.DEFAULTS``.
_LIGHT_DEFAULTS: dict[str, float | bool] = {
    "lambda0": 0.02,        # 1/m; background extinction
    "lambda1": 0.0088,      # 1/m / (ug-Chla/L); linear chlorophyll self-shading
    "lambda2": 0.054,       # unitless; non-linear chlorophyll self-shading
    "lambdas": 0.052,       # L/mg/m; ISS extinction parameter
    "lambdam": 0.174,       # L/mg/m; POM extinction parameter (matches v3 global_vars and Fortran/QUAL2K; corrected from 0.0174 in Phase 9.C)
    "fcom": 0.4,            # unitless; carbon-to-organic-matter mass ratio
    "Fr_PAR": 0.47,         # unitless; PAR fraction of total solar radiation
    "use_Algae": True,      # bool; floating-algae module switch
    "use_POC": True,        # bool; POC module switch
}


class Pathogen(Process):
    """v3 NSM1 Pathogen (PX) Process.

    Independent constituent: reads only its own state plus environmental
    forcings; writes only ``pathogen`` back to the registry. Forward-Euler
    integrator on a per-cell rate of change in 1/d, mirroring the Phase
    2.A FloatingAlgae pattern.
    """

    variables = [
        "pathogen",
        "water_temperature",
        "depth",
        "solar_radiation",
    ]

    # Class-level v3 defaults. Lazy-loaded on first instantiation to
    # mirror the v2 nitrogen.py pattern (the lazy-load is a defensive
    # match for the v2 <-> v3 circular import; in v3-native modules the
    # circular import is not active, but the same idiom keeps the file
    # consistent with the established style in v2 process classes).
    DEFAULTS: dict[str, float | int | bool] = {}

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface Pathogen exposes via the opportunistic-write loop in
    # ``run``. Each name maps to a ``self.<name>`` cache attribute set
    # inside ``_rate_with_components`` and matches the inventory in
    # ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3.
    #
    # Per pattern-alignment spec §10 Q5: Pathogen is a single-state
    # rate-form integrator, so the canonical helper is named
    # ``_rate_with_components`` (returns ``(rate, components)``) rather
    # than ``_change_with_components`` (which returns
    # ``(delta_state(s)..., components)``). The naming reflects the
    # integrator's argument shape.
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "pathogen_natural_death_rate",
        "pathogen_light_death_rate",
        "pathogen_settling_rate",
    )

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
    ) -> None:
        """Initialize the Pathogen Process.

        Args:
            parameters: Optional dict of parameter overrides. Merged
                with the class-level ``DEFAULTS`` (v3
                ``PATHOGEN_DEFAULTS``). Unknown keys log a warning and
                are ignored.
            time_step: Substep cadence for this Process. Forward Euler
                converts to days as ``time_step.total_seconds() / 86400``.
        """
        Process.__init__(self, time_step)

        # --- v3-style parameter merge (DEFAULTS + user overrides) ---
        # Lazy-load PATHOGEN_DEFAULTS (mirror v2 nitrogen.py:91-93 pattern).
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.pathogen import (
                DEFAULTS as PATHOGEN_DEFAULTS,
            )
            type(self).DEFAULTS = PATHOGEN_DEFAULTS

        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS) - set(_LIGHT_DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "Pathogen: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in PATHOGEN_DEFAULTS).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Inline light-extinction defaults from non-pathogen groups. The
        # user can override any of these via ``parameters``; otherwise we
        # use the values from v3 global_vars/global_parameters/algae
        # DEFAULTS (kept inline to avoid coupling Pathogen to those
        # modules at construction time).
        for k, v in _LIGHT_DEFAULTS.items():
            if not hasattr(self, k):
                setattr(self, k, user_params.get(k, v))

        # Diagnostics handle: a v3 ``Model`` will replace this with its
        # run-level Diagnostics in ``init_process``. The local default
        # keeps ``run`` callable in unit tests that don't wire a Model.
        self.diagnostics: Diagnostics = Diagnostics()

        # One-time-warning latch for missing optional registry inputs.
        self._warned_missing_solid = False
        self._warned_missing_poc = False
        self._warned_missing_ap = False

    @ProcessFactory.register("pathogen")
    @staticmethod
    def from_config(
        config: dict, variable_registry: VariableRegistry
    ) -> "Pathogen":
        return Pathogen(**config)

    def init_process(
        self, model: "Model", registry: VariableRegistry
    ) -> None:
        """Wire Pathogen to the v3 Model.

        Captures the Model's run-level ``Diagnostics`` handle so that
        ``clip_negative_state`` records clip events on the canonical
        per-run diagnostics container.
        """
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Advance the pathogen state by one substep via Forward Euler.

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_rate_with_components``
        (B — rate-form per spec §10 Q5); applies Forward Euler with
        unconditional clip-with-log (C, D); persists primary output
        (E); caches step-scoped rates on ``self.<name>`` (F);
        opportunistically writes diagnostics (G).

        """
        # --- State and forcing reads (pattern A) ---
        px = registry.get_at_time("pathogen", time)
        water_temperature = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)
        q_solar = registry.get_at_time("solar_radiation", time)

        # Optional light-extinction inputs. Default to zero (with a
        # one-time warning per missing input) so the Pathogen Process
        # can run stand-alone for the Tier 1 closed-system harness,
        # where these reservoirs may not be populated.
        solid = self._get_optional(registry, "Solid", "_warned_missing_solid", time)
        poc = self._get_optional(registry, "poc", "_warned_missing_poc", time)
        # Canonical floating-algae state is ``algae_floating`` (registered
        # by FloatingAlgae and mapped from the riverine "Ap" constituent).
        # Prefer it; fall back to the legacy ``ap`` name only during
        # migration (pathogen_algae_name design spec, 2026-05-30).
        ap = self._get_optional(
            registry, "algae_floating", "_warned_missing_ap", time, fallback="ap"
        )

        # --- Fused rate composition (pattern B; rate-form per §10 Q5) ---
        rate, components = self._rate_with_components(
            px=px,
            water_temperature=water_temperature,
            depth=depth,
            q_solar=q_solar,
            solid=solid,
            poc=poc,
            ap=ap,
        )

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])

        # --- Forward Euler (pattern C) ---
        # Rates are 1/d; convert dt from seconds to days.
        dt_days = self.time_step.total_seconds() / 86400.0
        px_new = px + rate * dt_days

        # --- Clip-with-log per Q7 (pattern D) ---
        # Step attribution is automatic via ``diagnostics.current_step``
        # (Phase 0.6 Q1).
        px_new = clip_negative_state(px_new, "pathogen", self.diagnostics)

        # --- Persist primary output (pattern E) ---
        registry.set_at_time("pathogen", time, px_new)

        # --- Opportunistic diagnostic registry writes (pattern G) ---
        for name in self.REGISTRY_DIAGNOSTICS:
            if name in registry:
                registry.set_at_time(name, time, components[name])

    # ------------------------------------------------------------------
    # Rate-composition helpers (rate-form per spec §10 Q5)
    # ------------------------------------------------------------------

    def _rate_with_components(
        self,
        *,
        px: ArrayLike,
        water_temperature: ArrayLike,
        depth: ArrayLike,
        q_solar: ArrayLike,
        solid: ArrayLike,
        poc: ArrayLike,
        ap: ArrayLike,
    ) -> tuple[ArrayLike, dict]:
        """Compute ``(rate, components)`` for Pathogen.

        ``rate`` is the net per-day PX rate of change (cfu/100mL/d or
        whatever PX-unit is supplied). ``components`` carries the three
        positive-magnitude sub-rates: natural death, light-induced
        death, and settling. The integrator applies them as sinks via
        ``rate = -(natural + light + settling)``.

        Code-motion-only refactor (§11.6): operand order and
        sub-flux helpers preserved verbatim from the pre-Phase-8
        ``rate()`` body. ``rate()`` is retained for back-compat with
        any external caller; it now delegates to this helper.

        The companion shadow ``_rate_legacy_inline`` returns just
        ``rate`` and is used by
        ``tests/v3/nsm1/test_pathogen_helper_vs_inline.py``.
        """
        natural_death = self._rate_natural_decay(px, water_temperature)
        light_death = self._rate_light_decay(px, depth, q_solar, solid, poc, ap)
        settling = self._rate_settling(px, depth)

        rate = -(natural_death + light_death + settling)

        # NaN/inf guard at the net rate (defense-in-depth; primary
        # dry-cell defense is the orchestration-layer wet-mask in
        # Model).
        rate = sanitize_rate(rate)

        components = {
            "pathogen_natural_death_rate": natural_death,
            "pathogen_light_death_rate": light_death,
            "pathogen_settling_rate": settling,
        }

        return rate, components

    # ------------------------------------------------------------------
    # Kinetic helpers
    # ------------------------------------------------------------------

    def rate(
        self,
        px: ArrayLike,
        water_temperature: ArrayLike,
        depth: ArrayLike,
        q_solar: ArrayLike,
        solid: ArrayLike,
        poc: ArrayLike,
        ap: ArrayLike,
    ) -> ArrayLike:
        """Compute the per-cell rate of change of PX (1/d * PX).

        Returns the *absolute* rate of change ``dPX/dt`` (cfu/100mL/d
        or whatever unit PX is supplied in), not the lumped 1/d
        coefficient. v1 source: ``dPXdt`` (line 3209-3224).
        """
        return -(
            self._rate_natural_decay(px, water_temperature)
            + self._rate_light_decay(px, depth, q_solar, solid, poc, ap)
            + self._rate_settling(px, depth)
        )

    def _rate_natural_decay(
        self, px: ArrayLike, water_temperature: ArrayLike
    ) -> ArrayLike:
        """Natural pathogen decay (cfu/100mL/d).

        v1 ``kdx_tc`` (line 3141) and ``PathogenDeath`` (line 3158).
        """
        kdx_tc = arrhenius_correction(
            water_temperature, self.kdx_20, self.kdx_theta
        )
        return kdx_tc * px

    def _rate_light_decay(
        self,
        px: ArrayLike,
        depth: ArrayLike,
        q_solar: ArrayLike,
        solid: ArrayLike,
        poc: ArrayLike,
        ap: ArrayLike,
    ) -> ArrayLike:
        """Light-induced pathogen decay (cfu/100mL/d).

        v1 ``PathogenDecay`` (line 3172):

            apx * q_solar / (KEXT * depth)
                * (1 - exp(-KEXT * depth)) * PX

        Units note: ``q_solar`` is total incident solar radiation in
        W/m^2 (consistent with ``utils.light.PAR`` and
        ``parameters/global_vars.py``). v1's docstring mislabeled the
        parameter as 1/d but the value (500) and consumption pattern
        (Beer-Lambert scaling) are unambiguously W/m^2; the v3 port
        treats it as W/m^2 throughout (Phase 9.F docstring fix; see
        ``parameter_defaults_corrections.md`` Section 2.7).

        Phase 9.F.B reverted the Phase 3.1 substitution
        ``I0 = q_solar * Fr_PAR`` and now passes total broadband
        ``q_solar`` directly into the depth-averaged Beer-Lambert
        integral, matching v1, Fortran, and the canonical
        Auer & Niehaus (1993) / Chapra (1997) / QUAL2K formulation.
        Pathogen inactivation is largely UVA/UVB-mediated, not
        PAR-mediated, so the PAR substitution was a v3-introduced
        deviation rather than a port improvement (the v3 default
        ``apx = 0.017`` is the Auer/Niehaus canonical value tied to
        broadband irradiance, so reverting the substitution restores
        the canonical literature calibration). See
        ``parameter_defaults_corrections.md`` Section 1.15.
        """
        i0 = q_solar
        kext = L(
            lambda0=self.lambda0,
            lambda1=self.lambda1,
            lambda2=self.lambda2,
            lambdas=self.lambdas,
            lambdam=self.lambdam,
            Solid=solid,
            POC=poc,
            fcom=self.fcom,
            Ap=ap,
            use_Algae=self.use_Algae,
            use_POC=self.use_POC,
        )
        # Guard against (KEXT*depth -> 0) by short-circuiting the
        # depth-averaged light-availability factor to 1 in that limit
        # (per L'Hopital: lim_{x->0} (1 - e^{-x}) / x = 1). For typical
        # production cells with KEXT > 0 and depth > 0 this branch is
        # not exercised; it exists so a degenerate registry (depth==0
        # or KEXT==0) does not produce NaN.
        kd = kext * depth
        light_avail = xr.where(
            kd > 0.0,
            (1.0 - np.exp(-kd)) / xr.where(kd > 0.0, kd, 1.0),
            1.0,
        )
        return self.apx * i0 * light_avail * px

    def _rate_settling(self, px: ArrayLike, depth: ArrayLike) -> ArrayLike:
        """Pathogen settling (cfu/100mL/d).

        v1 ``PathogenSettling`` (line 3193):

            vx / depth * PX
        """
        return self.vx / depth * px

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get_optional(
        self,
        registry: VariableRegistry,
        name: str,
        warn_attr: str,
        time: datetime,
        fallback: str | None = None,
    ) -> ArrayLike:
        """Read an optional registry variable; default to 0 with a
        one-time warning when missing.

        Used for ``Solid`` / ``poc`` / ``algae_floating`` —
        light-extinction inputs that may not be populated when Pathogen
        runs stand-alone.

        When ``fallback`` is given, the canonical ``name`` is preferred
        and ``fallback`` is consulted only if ``name`` is absent. This
        mirrors the ``tip`` / ``phosphorus_total_inorganic`` precedence
        pattern used in the algae processes for migration back-compat
        (here: canonical ``algae_floating`` preferred, legacy ``ap``
        fallback). The zero-with-warning default fires only when neither
        name is present.
        """
        if name in registry:
            return registry.get_at_time(name, time)
        if fallback is not None and fallback in registry:
            return registry.get_at_time(fallback, time)
        if not getattr(self, warn_attr):
            logger.warning(
                "Pathogen: optional registry variable %r not present; "
                "treating as 0 for the light-extinction calculation.",
                name,
            )
            setattr(self, warn_attr, True)
        return 0.0
