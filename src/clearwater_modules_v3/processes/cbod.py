"""v3 NSM1 CBOD Process.

Phase 3.3 of the v3 NSM1 implementation plan
(``design/clearwater_modules_v3_nsm1_design_specification.md`` Section 11
Phase 3, Section 4 component inventory). v3-native (NOT a v2 overlay).

CBOD (Carbonaceous Biological Oxygen Demand) is unique among NSM1 v3
constituents in that v1 supports a multi-group representation: each
labile group (e.g., CBOD1, CBOD2, CBOD3) carries its own oxidation rate
constant. This Phase 3.3 implementation supports a single CBOD group by
default (the typical default in v1 usage and what the Tier 1 fixture
provides). The multi-group extension path is documented inline at the
top of ``run`` so a Phase 4+ multi-group port can re-use the same
DEFAULTS-merge / Forward-Euler scaffolding without restructuring the
class.

Kinetics (mirrors v1 ``processes.py:2334-2422``):

    dCBOD/dt = -kbod_tc * DOX / (KsOxbod + DOX) * CBOD     # oxidation
               -ksbod_tc * CBOD                             # settling (1/d)

where:

* ``kbod_tc = arrhenius_correction(kbod_20, kbod_theta, T)``
* ``ksbod_tc = arrhenius_correction(ksbod_20, ksbod_theta, T)``
* ``KsOxbod`` is the DOX half-saturation for CBOD oxidation
* ``DOX`` is read from the registry; if absent (e.g., Phase 3 standalone
  CBOD test before the Phase 5 DOX Process lands), a stub value of
  ``8.0 mg/L`` is used and a warning is logged once per ``run`` call.

NSM1-SCI-CB1 (gold-standard spec C2; research doc
``docs/clearwater_modules_v3_nsm1_research_2_3_ksbod.md``;
``parameter_defaults_corrections.md`` §2.3): CBOD settling is a
**first-order rate constant (1/d at 20 °C)**, NOT a settling velocity.
Fortran NSM1 ``modCBOD.f90:114`` and QUAL2E apply ``ksbod_tc * CBOD``
with **no depth division**; pre-fix v3 divided by depth (treating it as
m/d), which silently diverged from the Fortran 1/d convention by a
factor of ``1/depth`` for any nonzero ``ksbod_20``. The Arrhenius
coefficient is the **settling** value ``ksbod_theta = 1.024``
(Bowie 1985 / QUAL2E), not the oxidation ``1.047``.

Note (Phase 0 audit): the v3 ``ksbod_20=0`` default means CBOD does not
settle by default. The settling code path is wired in for completeness
but is identically zero unless the user passes ``ksbod_20 > 0`` — so the
SCI-CB1 form/θ correction is dormant at the shipped default and does not
perturb the coupled-demo trajectory.

Q10 GS-rates contract: the per-step oxidation rate
``self.cbod_oxidation_rate`` is cached as an instance attribute after
``run`` completes. Phase 5 DOX will read it as an upstream rate
variable (sink term in the DOX integrator).
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
from clearwater_modules_v3.parameters.cbod import DEFAULTS as CBOD_DEFAULTS
from clearwater_modules_v3.utils.conversions import arrhenius_correction
from clearwater_modules_v3.utils.numerics import (
    Diagnostics,
    clip_negative_state,
    sanitize_rate,
)


if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


logger = logging.getLogger(__name__)


# Stub DOX concentration used when the registry does not yet contain a
# ``oxygen_dissolved`` variable (Phase 3 standalone CBOD before Phase 5
# DOX Process). 8.0 mg/L is a typical surface-water saturation value.
_DOX_STUB_MG_PER_L: float = 8.0


class CBOD(Process):
    """v3 NSM1 single-group CBOD Process.

    Forward Euler integrator advances CBOD by ``dt_days = time_step / 86400``
    each substep:

        cbod_new = cbod + (-oxidation_rate - settling_rate) * dt_days

    Negative cells are clipped to zero via ``clip_negative_state`` with
    diagnostics on ``self.diagnostics`` (a v3 ``Diagnostics`` instance,
    or the v3 Model's run-level diagnostics if ``init_process`` is
    called with a v3 Model).

    Multi-group extension path
    --------------------------
    A future Phase 4+ multi-group CBOD Process can keep this class
    intact and add a sibling ``CBODMultiGroup`` that owns a list of
    ``CBOD`` instances (one per group). Alternatively, the registry can
    store ``cbod_1``, ``cbod_2``, ... and this class is parameterized
    by ``self.group_index`` to read/write the appropriate variable. The
    current single-group implementation reads ``"cbod"`` for parity
    with the Tier 1 fixture.
    """

    variables = [
        "cbod",
        "water_temperature",
        "depth",
        "oxygen_dissolved",
    ]

    # Class-level v3 defaults (design spec Section 3.4). Loaded from
    # the v3 ``parameters.cbod`` module.
    DEFAULTS: dict[str, float | int | bool] = CBOD_DEFAULTS

    # Pattern-alignment spec §4 / Appendix A diff: the registry-diagnostics
    # surface CBOD exposes via the opportunistic-write loop in ``run``.
    # Each name maps to a ``self.<name>`` cache attribute set inside
    # ``_change_with_components`` and matches the inventory in
    # ``design/clearwater_modules_v3_nsm1_appendix_a_diff.md`` §3.
    #
    # ``cbod_oxidation_rate`` is the preserved attribute name DOX and
    # Carbon read via getattr (sibling-consumer contract).
    REGISTRY_DIAGNOSTICS: tuple[str, ...] = (
        "cbod_oxidation_rate",
        "cbod_settling_rate",
    )

    def __init__(
        self,
        parameters: dict | None = None,
        time_step: timedelta = timedelta(minutes=5),
    ) -> None:
        """Initialize the CBOD process.

        Args:
            parameters: Optional dict of v3 CBOD parameter overrides.
                Merged with the class-level ``DEFAULTS``
                (``clearwater_modules_v3.parameters.cbod.DEFAULTS``).
                Unknown keys are warned and ignored. Recognized keys
                include ``KsOxbod``, ``kbod_20``, ``ksbod_20``,
                ``kbod_theta``, ``ksbod_theta``.
            time_step: Substep cadence for this Process.
        """
        user_params = parameters or {}
        unknown_keys = set(user_params) - set(self.DEFAULTS)
        for key in sorted(unknown_keys):
            logger.warning(
                "CBOD: unknown parameter %r in 'parameters' dict; "
                "ignoring (not in CBOD_DEFAULTS).",
                key,
            )
        merged = {**self.DEFAULTS, **user_params}
        for k, v in merged.items():
            setattr(self, k, v)

        # Step-scoped rate-variable cache per the resolved Q10 GS-rates
        # contract (design spec Section 3.5, Appendix A). Phase 5 DOX
        # will consume ``self.cbod_oxidation_rate`` as a sink term in
        # its integrator (mg-O2/L/d, equivalent to the CBOD oxidation
        # rate because 1 mg-CBOD == 1 mg-O2 by definition).
        self.cbod_oxidation_rate: ArrayLike = 0.0
        self.cbod_settling_rate: ArrayLike = 0.0

        # Whether to scale oxidation by DOX/(KsOxbod+DOX). Mirrors v1's
        # ``use_DOX`` flag (``processes.py:2386``). Defaults to True so
        # the Tier 1 closed-system test exercises the full kinetics; a
        # standalone CBOD test that does not have DOX in the registry
        # falls back to the stub value (see ``run``) rather than
        # disabling the attenuation.
        self.use_DOX: bool = True

        # Diagnostics handle: a fresh v3 ``Diagnostics`` until
        # ``init_process`` overrides with the v3 Model's run-level
        # handle. v2 Models do not provide ``model.diagnostics``.
        self.diagnostics: Diagnostics = Diagnostics()

        Process.__init__(self, time_step)

    @ProcessFactory.register("cbod")
    @staticmethod
    def from_config(
        config: dict, variable_registry: VariableRegistry
    ) -> "CBOD":
        return CBOD(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """Capture the v3 Model's run-level Diagnostics handle if present.

        Phase 5 DOX coupling note: CBOD reads DOX from the registry in
        ``run`` (path (a) in the Phase 3.3 design discussion). When the
        DOX Process is not yet wired in (Phase 3 standalone), the
        registry will not contain ``oxygen_dissolved`` and ``run``
        falls back to the stub value with a logged warning.
        """
        model_diagnostics = getattr(model, "diagnostics", None)
        if model_diagnostics is not None:
            self.diagnostics = model_diagnostics

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """Advance CBOD by one substep using Forward Euler.

        Pattern-alignment spec §3 patterns A–J: reads forcings at top
        (A); delegates rate composition to ``_change_with_components``
        (B); applies Forward Euler with unconditional clip-with-log
        (C, D); persists primary output (E); caches step-scoped rates
        on ``self.<name>`` (F); opportunistically writes diagnostics
        (G).


        Multi-group note: a future extension would loop over
        ``cbod_<group_index>`` keys here, computing per-group
        oxidation/settling rates and summing across groups for
        downstream DOX consumption. The single-group implementation
        hard-codes the ``"cbod"`` key.
        """
        # --- State and forcing reads (pattern A) ---
        cbod = registry.get_at_time("cbod", time)
        water_temperature = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)

        # DOX coupling: read from the registry if present; fall back to
        # the stub value otherwise (Phase 5 has not landed yet OR the
        # test is running CBOD standalone). Path (a) per the Phase 3.3
        # discussion. Phase 1.E: switched from ``try/except KeyError`` to
        # the ``if name in registry`` idiom.
        if "oxygen_dissolved" in registry:
            dox = registry.get_at_time("oxygen_dissolved", time)
        else:
            logger.warning(
                "CBOD.run: 'oxygen_dissolved' not found in registry; "
                "falling back to stub value %.2f mg/L. This is expected "
                "for Phase 3 standalone CBOD tests; Phase 5 DOX wiring "
                "will populate the registry.",
                _DOX_STUB_MG_PER_L,
            )
            dox = xr.full_like(cbod, _DOX_STUB_MG_PER_L)

        # --- Fused rate composition (pattern B) ---
        rate, components = self._change_with_components(
            cbod=cbod,
            water_temperature=water_temperature,
            depth=depth,
            dox=dox,
        )

        # --- Cache step-scoped rates on ``self.<name>`` (pattern F) ---
        # ``cbod_oxidation_rate`` is consumed by DOX and Carbon via
        # getattr; preserved attribute name.
        for name in self.REGISTRY_DIAGNOSTICS:
            setattr(self, name, components[name])

        # --- Forward Euler in days (pattern C) ---
        dt_days = self.time_step.total_seconds() / 86400.0
        cbod_new = cbod + rate * dt_days

        # --- Clip-with-log per Q7 (pattern D) ---
        cbod_new = clip_negative_state(cbod_new, "cbod", self.diagnostics)

        # --- Persist primary output (pattern E) ---
        registry.set_at_time("cbod", time, cbod_new)

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
        cbod: ArrayLike,
        water_temperature: ArrayLike,
        depth: ArrayLike,
        dox: ArrayLike,
    ) -> tuple[ArrayLike, dict]:
        """Compute ``(rate, components)`` for CBOD.

        ``rate`` is the net per-day CBOD rate of change (mg-O2/L/d) —
        the sum of the negated oxidation and settling sinks. The
        positive-magnitude oxidation and settling terms are exposed in
        the components dict.

        Code-motion-only refactor of ``run``'s former inline
        composition (§11.6): operand order, intermediate names, the
        ``use_DOX`` Monod gating, and the per-sub-flux ``sanitize_rate``
        calls are preserved verbatim.

        The companion shadow ``_change_legacy_inline`` returns just the
        net rate and is used by
        ``tests/v3/nsm1/test_cbod_helper_vs_inline.py``.
        """
        # Temperature-corrected rate constants (Arrhenius / van't Hoff).
        kbod_tc = arrhenius_correction(
            water_temperature, self.kbod_20, self.kbod_theta
        )
        ksbod_tc = arrhenius_correction(
            water_temperature, self.ksbod_20, self.ksbod_theta
        )

        # Oxidation rate (mg-O2/L/d). Monod attenuation by DOX when
        # use_DOX is True; otherwise first-order in CBOD only.
        if self.use_DOX:
            oxidation_rate = (
                kbod_tc * dox / (self.KsOxbod + dox) * cbod
            )
        else:
            oxidation_rate = kbod_tc * cbod

        # Settling rate (mg-O2/L/d). NSM1-SCI-CB1 (spec C2): ``ksbod`` is
        # a first-order rate (1/d at 20 °C), so the form is
        # ``ksbod_tc * cbod`` with NO depth division -- matching Fortran
        # ``modCBOD.f90:114`` and QUAL2E. Pre-fix v3 used
        # ``ksbod_tc / depth * cbod`` (a velocity, m/d), diverging by
        # 1/depth for any nonzero ``ksbod_20``. ``depth`` is still read
        # from the registry (variables contract / symmetry) but is no
        # longer used by this term. With the shipped ``ksbod_20=0`` this
        # is identically zero regardless of form.
        settling_rate = ksbod_tc * cbod

        # Sanitize per-sub-flux at the cache source: a NaN/inf here
        # propagates via DOX's rate sum and zeroes the entire cell's
        # DOX rate via the downstream ``sanitize_rate``, freezing the
        # cell at IC. Phase 1.E adopted ``sanitize_rate`` for parity
        # with the other v3 NSM1 Processes.
        oxidation_rate = sanitize_rate(oxidation_rate)
        settling_rate = sanitize_rate(settling_rate)

        rate = -oxidation_rate - settling_rate

        components = {
            "cbod_oxidation_rate": oxidation_rate,
            "cbod_settling_rate": settling_rate,
        }

        return rate, components
