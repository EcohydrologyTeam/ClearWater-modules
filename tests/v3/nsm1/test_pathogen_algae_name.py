"""Pathogen reads the canonical floating-algae name (algae_floating).

Guards the fix from ``design/clearwater_modules_v3_pathogen_algae_name.md``:
the Pathogen process used to read the phytoplankton state for its
light-extinction / shading term under the name ``ap``, which no process
registers (the canonical name is ``algae_floating``, registered by
FloatingAlgae and mapped from the riverine ``Ap`` constituent). In a
coupled run this silently used zero algal shading and emitted a
fallback warning every step.

The fix prefers ``algae_floating`` and falls back to the legacy ``ap``
only when the canonical name is absent (the ``tip`` /
``phosphorus_total_inorganic`` precedence pattern). These tests pin:

- canonical ``algae_floating`` drives the shading term (no warning;
  result differs from the zero-algae case);
- the legacy ``ap`` fallback is value-equivalent during migration;
- ``algae_floating`` is preferred when both are present;
- when neither name is present, the input is treated as zero with a
  single one-time warning.

Light-extinction physics check (utils.light.L): more algae -> larger
KEXT -> smaller (1 - exp(-KEXT*d)) / (KEXT*d) factor -> SMALLER
light-death magnitude. So a bloom reduces the (positive-magnitude)
``pathogen_light_death_rate`` relative to zero algae.
"""
from __future__ import annotations

import logging
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.examples import InMemoryRegistry
from clearwater_modules_v3.processes.pathogen import Pathogen


N = 5
START = datetime(2026, 1, 1, 0, 0, 0)
PATHOGEN_LOGGER = "clearwater_modules_v3.processes.pathogen"


def _da(value) -> xr.DataArray:
    return xr.DataArray(np.full(N, value, dtype=float), dims=["nface"])


def _registry(*, algae_floating=None, ap=None, solar=300.0) -> InMemoryRegistry:
    """Build a registry with the required pathogen forcings plus the
    requested floating-algae name(s). ``Solid`` and ``poc`` are registered
    as zero so the only shading input that varies is the algae, and so
    they emit no missing-input warnings."""
    reg = InMemoryRegistry()
    reg.register("pathogen", _da(1.0e4))
    reg.register("water_temperature", _da(20.0))
    reg.register("depth", _da(1.0))
    reg.register("solar_radiation", _da(solar))
    reg.register("Solid", _da(0.0))
    reg.register("poc", _da(0.0))
    if algae_floating is not None:
        reg.register("algae_floating", _da(algae_floating))
    if ap is not None:
        reg.register("ap", _da(ap))
    return reg


def _light_rate(registry) -> np.ndarray:
    """Run a fresh Pathogen one step against ``registry`` and return the
    cached positive-magnitude light-death rate (pattern-F diagnostic)."""
    p = Pathogen()
    p.run(START, registry)
    return np.asarray(p.pathogen_light_death_rate)


def test_algae_floating_drives_shading_no_warning(caplog):
    """Canonical ``algae_floating`` is read and shades the light term; no
    missing-input warning; result differs from (is below) zero-algae."""
    with caplog.at_level(logging.WARNING, logger=PATHOGEN_LOGGER):
        rate_bloom = _light_rate(_registry(algae_floating=50.0))
    assert not any(
        "not present" in r.getMessage() for r in caplog.records
    ), "canonical algae_floating present -> must not emit a fallback warning"

    rate_zero = _light_rate(_registry(algae_floating=0.0))
    assert np.all(rate_bloom > 0.0)
    # Shading reduces the light-death magnitude (see module docstring).
    assert np.all(rate_bloom < rate_zero)


def test_legacy_ap_fallback_matches_canonical():
    """During migration, reading the legacy ``ap`` (canonical absent) is
    value-equivalent to reading ``algae_floating`` with the same value."""
    rate_canonical = _light_rate(_registry(algae_floating=50.0))
    rate_legacy = _light_rate(_registry(ap=50.0))
    np.testing.assert_array_equal(rate_canonical, rate_legacy)


def test_canonical_preferred_when_both_present():
    """When both names are present, ``algae_floating`` wins over ``ap``."""
    rate_both = _light_rate(_registry(algae_floating=50.0, ap=999.0))
    rate_canonical_only = _light_rate(_registry(algae_floating=50.0))
    np.testing.assert_array_equal(rate_both, rate_canonical_only)

    # And it must NOT match the result of using the legacy value.
    rate_legacy_value = _light_rate(_registry(ap=999.0))
    assert not np.allclose(rate_both, rate_legacy_value)


def test_missing_both_treated_as_zero():
    """With neither name present, the algae input is treated as zero —
    same light term as an explicit zero-algae registry."""
    rate_missing = _light_rate(_registry())  # neither algae_floating nor ap
    rate_zero = _light_rate(_registry(algae_floating=0.0))
    np.testing.assert_array_equal(rate_missing, rate_zero)


def test_missing_both_warns_once(caplog):
    """The missing-input warning names the canonical variable and fires
    only once (the one-time latch), even across multiple substeps."""
    reg = _registry()  # neither algae_floating nor ap
    p = Pathogen()
    with caplog.at_level(logging.WARNING, logger=PATHOGEN_LOGGER):
        p.run(START, reg)
        p.run(START, reg)
    algae_warnings = [
        r.getMessage()
        for r in caplog.records
        if "algae_floating" in r.getMessage() and "not present" in r.getMessage()
    ]
    assert len(algae_warnings) == 1
