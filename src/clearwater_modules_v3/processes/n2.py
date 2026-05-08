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

        Reads:
        * ``n2`` (mg-N/L) at ``t = time``
        * ``water_temperature`` (deg C) at ``t = time``
        * ``depth`` (m) at ``t = time``
        * ``atmospheric_pressure`` (mb) at ``t = time`` -- optional;
          falls back to ``self.pressure_mb`` when absent.

        Writes:
        * ``n2`` (mg-N/L) at ``t = time`` (in-place state update)
        * ``total_dissolved_gas`` (fraction) at ``t = time`` -- optional;
          only written when the registry knows the variable.

        Caches (for downstream Processes):
        * ``self.n2_sat``, ``self.n2_atm_exchange_rate``, ``self.tdg``
        """
        n2_state = registry.get_at_time("n2", time)
        t_water_c = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)

        # Atmospheric pressure: prefer registry-provided per-cell value
        # (e.g. weather forcing) over the scalar parameter.
        if "atmospheric_pressure" in registry:
            pressure_mb = registry.get_at_time("atmospheric_pressure", time)
        else:
            pressure_mb = self.pressure_mb

        # --- Henry's-law saturation ---
        t_water_k = _kelvin(t_water_c)
        khn2 = khn2_tc(t_water_k)
        pwv_atm = pwv(t_water_k)
        n2_sat = n2sat_henry(khn2, pressure_mb, pwv_atm)

        # --- Effective reaeration coefficient ka_tc (1/d) ---
        # Fast path: when both menu options are user-defined (==1) and
        # both user values are zero, ``ka_tc`` is identically zero and
        # we can short-circuit. This also avoids a known numpy.select
        # broadcast quirk in ``utils.reaeration.kah_20`` where mixing
        # scalar option-checks with DataArray-shaped depth conditions
        # introduces a spurious second dimension on the output even
        # when the ``option == 1`` branch is selected (causing shape
        # (cell, dim_0) instead of (cell,)).
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
        # Integration Item 1: read the step-scoped denitrification flux
        # cached on the Nitrogen Process after its ``run``. Units are
        # mg-N/L/d, positive-valued (absolute magnitude of NO3 -> N2);
        # ``denitrification_flux_rate`` is distinct from the legacy
        # ``denitrification_rate`` kinetic rate-constant attribute.
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

        # NaN/inf guard (defense-in-depth; primary dry-cell defense
        # is the orchestration-layer wet-mask in Model). Catches
        # ``inf`` from ``x / depth`` at ``depth == 0`` and ``NaN``
        # from missing forcings.
        rate = sanitize_rate(rate)

        # --- Forward Euler integration ---
        dt_days = self.time_step.total_seconds() / 86400.0
        n2_new = n2_state + rate * dt_days

        # Clip-with-log per the Q7 contract.
        if isinstance(n2_new, xr.DataArray) and self.diagnostics is not None:
            n2_new = clip_negative_state(n2_new, "n2", self.diagnostics, step=0)
        else:
            n2_new = xr.where(n2_new < 0, 0, n2_new)

        # Persist updated state.
        registry.set_at_time("n2", time, n2_new)

        # --- Derived TDG (fraction of saturation; Phase 3 simple form) ---
        # Phase 5 will swap this for the O2-weighted form
        # ``0.79 * N2/N2sat + 0.21 * DOX/DOX_sat`` once the DOX Process
        # is available.
        with np.errstate(divide="ignore", invalid="ignore"):
            tdg = n2_new / n2_sat
        tdg = sanitize_rate(tdg)

        if "total_dissolved_gas" in registry:
            registry.set_at_time("total_dissolved_gas", time, tdg)

        # Cache step-scoped quantities for downstream consumers.
        self.n2_sat = n2_sat
        self.n2_atm_exchange_rate = atm_exchange
        self.tdg = tdg
