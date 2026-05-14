"""v3 NSM1 ``N2`` (dissolved nitrogen gas) Process.

Phase 3.4 (v3 NSM1 design spec, Section 11 Phase 3, Section 5 N2/TDG
design notes). This is a v3-native Process (not a v2 overlay):

- N2 saturation via Henry's law:
  ``N2sat = 2.8e4 * KHN2_tc * 0.79 * (pressure_atm - p_wv)``
  (mg-N/L), with ``KHN2_tc`` from a temperature-corrected Henry's-law
  coefficient (Van't Hoff form, NIST reference 0.00065 mol/(kg*bar) at
  298.15 K, dH/R = 1300 K) and ``p_wv`` from the empirical
  ``exp(11.8571 - 3840.7/T - 216961/T^2)`` water-vapor formula.
- Atmospheric exchange (1/d): ``1.034 * ka_tc * (N2sat - N2)``, where
  ``ka_tc`` is the temperature-corrected combined hydraulic + wind
  reaeration coefficient from ``utils.reaeration``.
- Denitrification source: read from the Nitrogen process when present
  (``getattr(nitrogen_process, "denitrification_flux_rate", 0)``). N2
  is the product of NO3 -> N2 denitrification. Per Integration Item 1
  (registry rate-variable convention, spec resolved Q10), Nitrogen
  caches the *step-scoped* denitrification flux (mg-N/L/d, positive
  magnitude) on this attribute after each ``run``.
- Forward Euler integrator: ``N2_new = N2 + (atm_exchange + denit_source) * dt_days``
- Derived TDG diagnostic: ``TDG = N2 / N2sat`` (Phase 3 simple form;
  the O2-weighted form ``0.79 * N2/N2sat + 0.21 * DOX/DOX_sat`` requires
  the DOX Process and is deferred to Phase 5).

DEFAULTS sources (v3 ``N2_DEFAULTS`` from ``parameters.n2`` is empty per
the Phase 1.2 audit). N2 pulls parameters from three other groups:

- ``parameters.global_parameters``: ``pressure_mb`` (atmospheric pressure)
- ``parameters.dox``: reaeration menu (``kah_20_user``, ``kaw_20_user``,
  ``kah_theta``, ``kaw_theta``, ``hydraulic_reaeration_option``,
  ``wind_reaeration_option``)
- ``parameters.global_vars``: hydraulic forcings (``velocity``, ``flow``,
  ``topwidth``, ``slope``, ``shear_velocity``, ``wind_speed``)

The merged DEFAULTS are applied to the instance as ``self.<name>``
attributes, mirroring the Phase 2.A FloatingAlgae and Phase 2.B Nitrogen
patterns. Per design spec Section 11, this Process is **not** wired into
the package ``__init__`` here -- registration is Phase 3.5.
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

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


logger = logging.getLogger(__name__)


# Conversion: 1 mb -> atm (1 atm == 1013.25 mb), used inside N2sat.
MB_TO_ATM: float = 1.0 / 1013.25


# ---------------------------------------------------------------------------
# Henry's law / saturation primitives (stateless; v1 NSM1 lines 3452-3540)
# ---------------------------------------------------------------------------


def _kelvin(t_c: ArrayLike) -> ArrayLike:
    """deg C -> K."""
    return t_c + 273.15


def khn2_tc(t_water_k: ArrayLike) -> ArrayLike:
    """Henry's-law constant for N2 (mol/L/atm), Van't Hoff temperature correction.

    Reference: NIST. ``KH(298.15K) = 0.00065 mol/(kg*bar)`` and the
    enthalpy-of-dissolution constant ``dH/R = 1300 K``.

    v1 source: ``clearwater_modules/nsm1/processes.py:KHN2_tc`` (3452-3467).
    """
    return 0.00065 * np.exp(1300.0 * (1.0 / t_water_k - 1.0 / 298.15))


def pwv(t_water_k: ArrayLike) -> ArrayLike:
    """Partial pressure of water vapor (atm), empirical fit.

    v1 source: ``clearwater_modules/nsm1/processes.py:pwv`` (2878-2886).
    """
    return np.exp(11.8571 - 3840.70 / t_water_k - 216961.0 / t_water_k**2)


def n2sat_henry(
    khn2_tc_value: ArrayLike,
    pressure_mb: ArrayLike,
    pwv_atm: ArrayLike,
) -> ArrayLike:
    """N2 saturation concentration (mg-N/L) from Henry's law.

    ``N2sat = 2.8e4 * KHN2_tc * 0.79 * (P_atm - p_wv)``

    v1 source: ``clearwater_modules/nsm1/processes.py:N2sat`` (3470-3487).
    """
    pressure_atm = pressure_mb * MB_TO_ATM
    raw = 2.8e4 * khn2_tc_value * 0.79 * (pressure_atm - pwv_atm)
    if isinstance(raw, xr.DataArray):
        return xr.where(raw < 0.0, 1e-6, raw)
    return np.where(raw < 0.0, 1e-6, raw)


# ---------------------------------------------------------------------------
# N2 Process
# ---------------------------------------------------------------------------


class N2(Process):
    """v3 N2 (dissolved nitrogen gas) Process.

    State variable: ``n2`` (mg-N/L).
    Derived variable: ``total_dissolved_gas`` (TDG, fraction of saturation).

    Sources / sinks (1/d, integrated by Forward Euler):
    * Atmospheric exchange: ``1.034 * ka_tc * (N2sat - N2)``
    * Denitrification: read from ``nitrogen_process.denitrification_rate``
      when the Nitrogen process is wired up; 0 otherwise.

    Caches computed step-scoped quantities for downstream consumers (TDG-
    aware Processes in v3 1.x):
    * ``self.n2_sat``                -- saturation concentration (mg-N/L)
    * ``self.n2_atm_exchange_rate``  -- atm exchange flux (mg-N/L/d)
    * ``self.tdg``                   -- total dissolved gas (fraction)
    """

    variables = [
        "n2",
        "total_dissolved_gas",
        "water_temperature",
        "depth",
        "atmospheric_pressure",
    ]

    # Class-level v3 defaults. Lazy-loaded on first instantiation;
    # populated from a curated subset of global_parameters / dox /
    # global_vars (see module docstring).
    DEFAULTS: dict[str, float | int | bool] = {}

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface N2 exposes via the opportunistic-write loop in ``run``.
    # Each name maps to a ``self.<name>`` cache attribute set inside
    # ``_change_with_components`` and matches the inventory in
    # ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3.
    #
    # N2 already had ``total_dissolved_gas`` opportunistically exposed
    # pre-Phase 8 (the sole pre-existing example of pattern G in v3
    # NSM1 1.0.0); Phase 8 *extends* the loop to cover the full
    # Appendix A set rather than replacing it.
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "n2_atm_exchange_rate",
        "n2_sat",
        "total_dissolved_gas",
        "n2_denit_source_rate",
    )

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
    ) -> None:
        Process.__init__(self, time_step)

        # Lazy-load DEFAULTS by composing from the relevant v3 parameter
        # groups. ``parameters.n2`` is empty per the Phase 1.2 audit;
        # N2 pulls from global_parameters (pressure_mb), dox (reaeration
        # menu), and global_vars (hydraulic forcings).
        if not type(self).DEFAULTS:
            from clearwater_modules_v3.parameters.n2 import (
                DEFAULTS as N2_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.global_parameters import (
                DEFAULTS as GLOBAL_PARAM_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.dox import (
                DEFAULTS as DOX_DEFAULTS,
            )
            from clearwater_modules_v3.parameters.global_vars import (
                DEFAULTS as GLOBAL_VAR_DEFAULTS,
            )

            composed: dict[str, float | int | bool] = {}
            # Atmospheric pressure (mb) lives in global_parameters.
            composed["pressure_mb"] = GLOBAL_PARAM_DEFAULTS["pressure_mb"]
            composed["use_DOX"] = GLOBAL_PARAM_DEFAULTS.get("use_DOX", True)
            # Reaeration menu lives in dox.
            for k in (
                "kah_20_user",
                "kaw_20_user",
                "kah_theta",
                "kaw_theta",
                "hydraulic_reaeration_option",
                "wind_reaeration_option",
            ):
                composed[k] = DOX_DEFAULTS[k]
            # Hydraulic forcings live in global_vars (toy values; should
            # be overridden per cell in production).
            for k in (
                "velocity",
                "flow",
                "topwidth",
                "slope",
                "shear_velocity",
                "wind_speed",
            ):
                composed[k] = GLOBAL_VAR_DEFAULTS[k]

            # Anything in N2_DEFAULTS itself wins over the composed fallbacks.
            composed.update(N2_DEFAULTS)
            type(self).DEFAULTS = composed

        # Merge user overrides over the composed defaults.
        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "N2: unknown parameter %r in 'parameters' dict; ignoring "
                "(not in N2_DEFAULTS / global_parameters / dox / global_vars).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Step-scoped rate / derived caches for downstream consumers.
        self.n2_sat: ArrayLike = 0.0
        self.n2_atm_exchange_rate: ArrayLike = 0.0
        self.tdg: ArrayLike = 0.0

        # Diagnostics handle (replaced by the model's run-level handle in
        # ``init_process`` if available).
        self.diagnostics = Diagnostics()

        # Coupling flags. Defaulted here so ``run`` works without an
        # explicit ``init_process`` call (e.g. Tier 1 unit tests).
        self.use_nitrogen: bool = False
        self.nitrogen_process = None

    @ProcessFactory.register("n2")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "N2":
        return N2(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """Wire up Nitrogen coupling and capture the run-level Diagnostics."""
        # Capture run-level Diagnostics if the v3 Model exposes one.
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

        # Optional coupling to Nitrogen for the denitrification source.
        if hasattr(model, "has_process") and model.has_process("Nitrogen"):
            self.use_nitrogen = True
            self.nitrogen_process = model.get_process("Nitrogen")
        else:
            self.use_nitrogen = False
            self.nitrogen_process = None

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Integrate the N2 state by one Forward Euler step.

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_change_with_components``
        (B); applies Forward Euler with unconditional clip-with-log
        (C, D); persists primary output (E); caches step-scoped rates
        on ``self.<name>`` (F); opportunistically writes diagnostics
        (G).

        ``total_dissolved_gas`` opportunistic write (the sole v3 1.0.0
        example of pattern G) to cover the full Appendix A set —
        ``n2_atm_exchange_rate``, ``n2_sat``, and the new
        ``n2_denit_source_rate``.

        ``self.tdg`` attribute name is **preserved** for back-compat
        with the existing v1-parity and Tier 1 tests; ``total_dissolved_gas``
        (Appendix A name) is set as an alias side effect via the
        pattern F setattr loop.
        """
        # --- State and forcing reads (pattern A) ---
        n2_state = registry.get_at_time("n2", time)
        t_water_c = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)

        # Atmospheric pressure: prefer registry-provided per-cell value
        # (e.g. weather forcing) over the scalar parameter.
        if "atmospheric_pressure" in registry:
            pressure_mb = registry.get_at_time("atmospheric_pressure", time)
        else:
            pressure_mb = self.pressure_mb

        # --- Fused rate composition (pattern B) ---
        rate, components = self._change_with_components(
            n2_state=n2_state,
            t_water_c=t_water_c,
            depth=depth,
            pressure_mb=pressure_mb,
        )

        # --- Forward Euler integration (pattern C) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        n2_new = n2_state + rate * dt_days

        # --- Clip-with-log per Q7 (pattern D) ---
        n2_new = clip_negative_state(n2_new, "n2", self.diagnostics)

        # --- Persist primary output (pattern E) ---
        registry.set_at_time("n2", time, n2_new)

        # --- Derived TDG (depends on n2_new; computed after the
        # integrator step). Phase 3 simple form: ``N2_new / N2_sat``.
        # Phase 5 will swap for the O2-weighted form once DOX is
        # available.
        with np.errstate(divide="ignore", invalid="ignore"):
            tdg = n2_new / components["n2_sat"]
        tdg = sanitize_rate(tdg)
        # ``self.tdg`` is the back-compat attribute name (v1-parity +
        # Tier 1 tests); ``total_dissolved_gas`` is the Appendix A
        # name. Both point at the same value.
        self.tdg = tdg
        components["total_dissolved_gas"] = tdg

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        # Sets self.n2_atm_exchange_rate, self.n2_sat,
        # self.total_dissolved_gas, self.n2_denit_source_rate. The
        # tdg alias above is idempotent on self.total_dissolved_gas.
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])

        # --- Opportunistic diagnostic registry writes (pattern G) ---
        # Extends the pre-Phase-8 ``total_dissolved_gas`` write to the
        # full Appendix A set.
        for name in self.REGISTRY_DIAGNOSTICS:
            if name in registry:
                registry.set_at_time(name, time, components[name])

    # ------------------------------------------------------------------
    # Rate-composition helper
    # ------------------------------------------------------------------

    def _change_with_components(
        self,
        *,
        n2_state: ArrayLike,
        t_water_c: ArrayLike,
        depth: ArrayLike,
        pressure_mb: ArrayLike,
    ) -> tuple[ArrayLike, dict]:
        """Compute ``(rate, components)`` for N2 (without TDG; TDG is
        computed post-integrator-step in ``run`` because it depends on
        ``n2_new``).

        ``rate`` is the net per-day N2 rate of change (mg-N/L/d).
        ``components`` carries ``n2_atm_exchange_rate``, ``n2_sat``,
        ``n2_denit_source_rate`` (and ``total_dissolved_gas`` is added
        by ``run`` after the integrator step).

        Code-motion-only refactor of ``run``'s former inline
        composition (§11.6): operand order, intermediate names, the
        ``is_user_*_zero`` fast-path short-circuit, and the sanitize
        guard preserved verbatim.
        """
        # --- Henry's-law saturation ---
        t_water_k = _kelvin(t_water_c)
        khn2 = khn2_tc(t_water_k)
        pwv_atm = pwv(t_water_k)
        n2_sat = n2sat_henry(khn2, pressure_mb, pwv_atm)

        # --- Effective reaeration coefficient ka_tc (1/d) ---
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

        # --- Atmospheric exchange flux (mg-N/L/d) ---
        atm_exchange = 1.034 * ka_tc_value * (n2_sat - n2_state)

        # --- Denitrification source (mg-N/L/d) ---
        if self.use_nitrogen and self.nitrogen_process is not None:
            denit_source = getattr(
                self.nitrogen_process, "denitrification_flux_rate", 0
            )
            if denit_source is None:
                denit_source = 0
        else:
            denit_source = 0

        # --- Net rate (mg-N/L/d) ---
        rate = atm_exchange + denit_source

        # NaN/inf guard.
        rate = sanitize_rate(rate)

        components = {
            "n2_atm_exchange_rate": atm_exchange,
            "n2_sat": n2_sat,
            "n2_denit_source_rate": denit_source,
        }

        return rate, components
