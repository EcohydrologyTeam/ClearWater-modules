"""Contract guards for NSM1 process ``variables`` / ``output_variables``.

These pin the two declaration bugs the ``init_from_file`` provider-coverage
check (``model.py`` Step 8) exposed, and the wet-mask output scope
(``model.py.__apply_wet_mask``):

1. FloatingAlgae declared the orphan ``"floating_algae"`` in ``variables``
   while it actually reads/writes the canonical ``"algae_floating"`` — so
   the provider check required a name nothing supplies.
2. N2 listed the write-only derived diagnostic ``"total_dissolved_gas"``
   among its required inputs.
3. ``output_variables`` (the wet-mask scope) must be the written state,
   not the input forcings, and BenthicAlgae — which subclasses
   FloatingAlgae — must override it (its state is ``benthic_algae``, not
   ``algae_floating``).

The end-to-end provider-coverage behavior is covered by
``tests/v3/test_riverine_init_from_file_integration_v3.py`` (riverine env).
"""
from __future__ import annotations

from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae
from clearwater_modules_v3.processes.n2 import N2


def test_floating_algae_declares_canonical_state_name():
    """variables names the canonical state ``algae_floating`` (read/written
    at runtime), not the orphan ``floating_algae``."""
    assert "algae_floating" in FloatingAlgae.variables
    assert "floating_algae" not in FloatingAlgae.variables
    assert FloatingAlgae.output_variables == ["algae_floating"]


def test_n2_excludes_derived_diagnostic_from_inputs():
    """total_dissolved_gas is a write-only derived diagnostic — it belongs
    in REGISTRY_DIAGNOSTICS, not in variables (required inputs)."""
    assert "total_dissolved_gas" not in N2.variables
    assert "total_dissolved_gas" in N2.REGISTRY_DIAGNOSTICS
    assert "n2" in N2.variables
    assert N2.output_variables == ["n2"]


def test_benthic_output_variables_overrides_floating():
    """BenthicAlgae subclasses FloatingAlgae; its written state is
    ``benthic_algae``, so it must NOT inherit FloatingAlgae's
    ``output_variables`` (which would NaN-mask the wrong state on dry cells)."""
    assert BenthicAlgae.output_variables == ["benthic_algae"]
    assert BenthicAlgae.output_variables != FloatingAlgae.output_variables
    assert "benthic_algae" in BenthicAlgae.variables


def test_output_variables_are_states_not_forcings():
    """Each output is also declared as a variable (the process reads its own
    state), and input forcings are never in output_variables (so the wet-mask
    does not NaN-mask forcings)."""
    for cls in (FloatingAlgae, BenthicAlgae, N2):
        for out in cls.output_variables:
            assert out in cls.variables, (
                f"{cls.__name__}: output {out!r} not declared in variables"
            )
    assert "solar_radiation" not in FloatingAlgae.output_variables
    assert "solar_radiation" not in BenthicAlgae.output_variables
    for forcing in ("water_temperature", "depth", "atmospheric_pressure"):
        assert forcing not in N2.output_variables
