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
               -ksbod_tc / depth * CBOD                      # settling

where:

* ``kbod_tc = arrhenius_correction(kbod_20, kbod_theta, T)``
* ``ksbod_tc = arrhenius_correction(ksbod_20, ksbod_theta, T)``
* ``KsOxbod`` is the DOX half-saturation for CBOD oxidation
* ``DOX`` is read from the registry; if absent (e.g., Phase 3 standalone
  CBOD test before the Phase 5 DOX Process lands), a stub value of
  ``8.0 mg/L`` is used and a warning is logged once per ``run`` call.

Note (Phase 0 audit): the v3 ``ksbod_20=0`` default means CBOD does not
settle by default. The settling code path is wired in for completeness
but is identically zero unless the user passes ``ksbod_20 > 0``.

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

from clearwater_modules_v2.processes.base import Process, ProcessFactory
from clearwater_modules_v3.parameters.cbod import DEFAULTS as CBOD_DEFAULTS
from clearwater_modules_v3.utils.conversions import arrhenius_correction
from clearwater_modules_v3.utils.numerics import (
    Diagnostics,
    clip_negative_state,
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

        Reads state at ``t = time`` and writes the updated state back at
        the same key (in-place semantic, matching the v2 NSM1 Process
        convention used by the Tier 1 ``InMemoryRegistry``). The v3
        production registry uses time-indexed reads/writes; both
        interfaces share the ``get_at_time`` / ``set_at_time`` API so
        the same code path serves both.

        Multi-group note: a Phase 4+ extension would loop over
        ``cbod_<group_index>`` keys here, computing per-group
        oxidation/settling rates and summing
        ``self.cbod_oxidation_rate`` across groups for downstream DOX
        consumption. The single-group implementation hard-codes the
        ``"cbod"`` key.
        """
        cbod = registry.get_at_time("cbod", time)
        water_temperature = registry.get_at_time("water_temperature", time)
        depth = registry.get_at_time("depth", time)

        # DOX coupling: try to read from the registry; fall back to the
        # stub value if absent (Phase 5 has not landed yet OR the test
        # is running CBOD standalone). Path (a) per the Phase 3.3
        # discussion.
        try:
            dox = registry.get_at_time("oxygen_dissolved", time)
        except KeyError:
            logger.warning(
                "CBOD.run: 'oxygen_dissolved' not found in registry; "
                "falling back to stub value %.2f mg/L. This is expected "
                "for Phase 3 standalone CBOD tests; Phase 5 DOX wiring "
                "will populate the registry.",
                _DOX_STUB_MG_PER_L,
            )
            dox = xr.full_like(cbod, _DOX_STUB_MG_PER_L)

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

        # Settling rate (mg-O2/L/d). v1 ``CBOD_sedimentation`` is
        # ``CBOD * ksbod_tc`` (units m/d * 1/m -> 1/d) but the divide
        # by depth makes it m/d * (1/depth) * CBOD which is the v3
        # convention (see Phase 3.3 spec). With the Phase 0 audit
        # default ``ksbod_20=0``, this term is identically zero.
        settling_rate = ksbod_tc / depth * cbod

        # Cache step-scoped rate variables for Phase 5 DOX consumption
        # (Q10 GS-rates contract). cbod_oxidation_rate is the sink term
        # DOX adds to its integrator; cbod_settling_rate is exposed for
        # diagnostic / future sediment coupling.
        self.cbod_oxidation_rate = oxidation_rate
        self.cbod_settling_rate = settling_rate

        # Forward Euler in days. time_step is stored as a timedelta;
        # convert to days the same way the other v3 Processes do.
        dt_days = self.time_step.total_seconds() / 86400.0
        cbod_new = cbod + (-oxidation_rate - settling_rate) * dt_days

        # Clip-with-log per the resolved Q7 contract. Tier 1 closed-
        # system tests assert clip_events is empty under physically
        # reasonable inputs.
        cbod_new = clip_negative_state(
            cbod_new, "cbod", self.diagnostics, step=0
        )

        registry.set_at_time("cbod", time, cbod_new)
