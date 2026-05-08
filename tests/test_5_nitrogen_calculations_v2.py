"""Parity tests: v2 Nitrogen sub-rate methods vs v1 nsm1.processes helpers.

Each test constructs a v2 Nitrogen instance, calls one of its rate
sub-methods directly, and compares to the equivalent v1 helper-function
output computed with the same inputs.

This file is the most important of the three v1-parity tests because the
Nitrogen process had 3 of the 4 LimnoTech bugs (missing set_at_time for
NH4 and NO3, time_step_frequency typo, broken NaN check). These tests
lock in the post-fix sub-rate behavior.

Scope: only the v1-equivalent sub-rate computations are exercised
(ammonium_nitrification, ammonium_from_bed, ammonium_decay_nitrate,
nitrate_denitrification, nitrate_bed_denitrification, nitrification_inhibition).
The full ``run()`` time-stepping is intentionally NOT exercised here because
the update equation has unrelated open issues that need design input.

Synthetic mesh: 5-cell numpy arrays, single time step.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.nitrogen import Nitrogen


@pytest.fixture(scope="function")
def nh4_5cell():
    # Function-scoped: some v2 sub-rate methods alias inputs (e.g., FloatingAlgae.limit_nitrogen
    # uses an in-place += that mutates its input). Function scope keeps tests independent.
    return xr.DataArray(np.array([0.05, 0.10, 0.15, 0.20, 0.30]))


@pytest.fixture(scope="function")
def no3_5cell():
    return xr.DataArray(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))


@pytest.fixture(scope="function")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="function")
def depth_5cell():
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))


@pytest.fixture(scope="function")
def dox_5cell():
    return xr.DataArray(np.array([4.0, 6.0, 7.0, 8.0, 10.0]))


@pytest.fixture(scope="function")
def n_instance():
    """Nitrogen instance with v1-aligned defaults."""
    inst = Nitrogen(
        time_step=timedelta(minutes=5),
        denitrification_rate=0.002,
        denitrification_theta=1.045,
        nitrification_rate=0.1,
        nitrification_theta=1.083,
        sediment_denitrification_rate=0.1,
        sediment_denitrification_theta=1.045,
        sediment_ammonium_release_rate=0.1,
        sediment_ammonium_release_theta=1.074,
        ammonium_decay_rate=0.0,
        ammonium_decay_theta=1.0,
        nitrification_oxygen_inhibition_factor=0.6,
    )
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False
    return inst


def test_nitrification_inhibition_matches_v1_NitrificationInhibition(n_instance, dox_5cell):
    """v2 nitrification_inhibition(DOX) == v1 NitrificationInhibition(use_DOX, KNR, DOX)."""
    v2_val = n_instance.nitrification_inhibition(dox_5cell)
    v1_val = v1.NitrificationInhibition(use_DOX=True, KNR=0.6, DOX=dox_5cell)

    np.testing.assert_allclose(v2_val.values, np.asarray(v1_val), rtol=1e-6)


def test_ammonium_nitrification_matches_v1_NH4_Nitrification(
    n_instance, nh4_5cell, water_temp_5cell, dox_5cell
):
    """v2 ammonium_nitrification == v1 NH4_Nitrification(NitrificationInhibition, NH4, knit_tc, use_NH4)."""
    v2_rate = n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)

    knit_tc_val = v1.knit_tc(water_temp_5cell, 0.1, 1.083)
    inhib = v1.NitrificationInhibition(use_DOX=True, KNR=0.6, DOX=dox_5cell)
    v1_rate = v1.NH4_Nitrification(
        NitrificationInhibition=inhib,
        NH4=nh4_5cell,
        knit_tc=knit_tc_val,
        use_NH4=True,
    )

    np.testing.assert_allclose(v2_rate.values, np.asarray(v1_rate), rtol=1e-6)


def test_ammonium_from_bed_matches_v1_NH4fromBed(n_instance, depth_5cell, water_temp_5cell):
    """v2 ammonium_from_bed == v1 NH4fromBed(depth, rnh4_tc)."""
    v2_rate = n_instance.ammonium_from_bed(depth=depth_5cell, temperature=water_temp_5cell)

    rnh4_tc_val = v1.rnh4_tc(water_temp_5cell, 0.1, 1.074)
    v1_rate = v1.NH4fromBed(depth=depth_5cell, rnh4_tc=rnh4_tc_val)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_nitrate_bed_denitrification_matches_v1_NO3_BedDenit(
    n_instance, depth_5cell, no3_5cell, water_temp_5cell
):
    """v2 nitrate_bed_denitrification == v1 NO3_BedDenit(depth, vno3_tc, NO3)."""
    v2_rate = n_instance.nitrate_bed_denitrification(
        depth=depth_5cell, nitrate=no3_5cell, temperature=water_temp_5cell
    )

    vno3_tc_val = v1.vno3_tc(water_temp_5cell, 0.1, 1.045)
    v1_rate = v1.NO3_BedDenit(depth=depth_5cell, vno3_tc=vno3_tc_val, NO3=no3_5cell)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_nitrate_denitrification_water_column_matches_v1_NO3_Denit(
    n_instance, no3_5cell, water_temp_5cell, dox_5cell
):
    """v2 nitrate_denitrification == v1 NO3_Denit(use_DOX, DOX, KsOxdn, kdnit_tc, NO3)
    (when DOX values are well-defined and yield no NaN denominator).
    """
    half_sat_o2 = 1.0
    v2_rate = n_instance.nitrate_denitrification(
        dissolved_oxygen=dox_5cell,
        half_saturation_oxygen=half_sat_o2,
        nitrate=no3_5cell,
        temperature=water_temp_5cell,
    )

    kdnit_tc_val = v1.kdnit_tc(water_temp_5cell, 0.002, 1.045)
    v1_rate = v1.NO3_Denit(
        use_DOX=True,
        DOX=dox_5cell,
        KsOxdn=half_sat_o2,
        kdnit_tc=kdnit_tc_val,
        NO3=no3_5cell,
    )

    np.testing.assert_allclose(np.asarray(v2_rate), np.asarray(v1_rate), rtol=1e-6)


def test_change_ammonium_no_algae_drops_to_minus_nitrification_plus_bed(
    n_instance, nh4_5cell, no3_5cell, water_temp_5cell, depth_5cell, dox_5cell
):
    """When use_floating_algae=use_benthic_algae=False, change_ammonium reduces to:
    -ammonium_nitrification + ammonium_from_bed.

    Phase 9.A.2 audit finding N2: the phantom ``ammonium_decay_nitrate`` term
    (no v1 or Fortran NSM1 analogue) was dropped from change_ammonium. Pre-fix,
    this test asserted the rate equaled ``decay - nitrification + bed``; post-fix
    it asserts ``-nitrification + bed``. The legacy decay method is retained
    on the instance for back-compat but is no longer part of the NH4 budget.

    Verifies the post-fix isnull NaN replacement (bug #4) does not perturb finite inputs.
    """
    n_instance.use_floating_algae = False
    n_instance.use_benthic_algae = False

    v2_rate = n_instance.change_ammonium(
        nitrate=no3_5cell,
        ammonium=nh4_5cell,
        temperature=water_temp_5cell,
        depth=depth_5cell,
        oxygen_dissolved=dox_5cell,
    )

    expected = (
        - n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
        + n_instance.ammonium_from_bed(depth=depth_5cell, temperature=water_temp_5cell)
    )

    np.testing.assert_allclose(v2_rate.values, expected.values, rtol=1e-6)
    # And no NaNs survive the isnull replacement
    assert not np.any(np.isnan(v2_rate.values))


# ---------------------------------------------------------------------------
# Phase 9.A.2 regression tests: default-instantiated Nitrogen() with v3 wiring.
# ---------------------------------------------------------------------------
#
# These tests pin the post-Phase-9.A.2 contract: a bare ``Nitrogen()`` (no
# kwargs) reads NITROGEN_DEFAULTS for nitrification/denitrification/sediment
# rate constants and produces v1/Fortran-correct kinetics. Pre-fix, the
# legacy v2 kwargs all defaulted to ``1.0`` (5x to 500x larger than NSM1
# defaults) and the kinetic methods read from the legacy attributes; this
# made default-instantiated Nitrogen unsafe at any non-zero state.


def test_default_nitrogen_uses_v3_defaults_for_nitrification(
    nh4_5cell, water_temp_5cell, dox_5cell
):
    """Phase 9.A.2 audit finding N1: default Nitrogen() reads ``knit_20=0.1`` /
    ``knit_theta=1.083`` / ``KNR=0.6`` from NITROGEN_DEFAULTS, not the legacy
    kwarg (which defaulted to 1.0 / 1.0 / 1.0).
    """
    inst = Nitrogen()
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False

    assert inst.knit_20 == 0.1
    assert inst.knit_theta == 1.083
    assert inst.KNR == 0.6

    # Inhibition uses KNR (0.6), not the legacy default (1.0).
    inhib = inst.nitrification_inhibition(dox_5cell)
    expected_inhib = 1.0 - np.exp(-0.6 * dox_5cell.values)
    np.testing.assert_allclose(np.asarray(inhib), expected_inhib, rtol=1e-12)

    # Nitrification flux uses knit_20=0.1 (not 1.0).
    rate = inst.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
    knit_tc = 0.1 * 1.083 ** (water_temp_5cell.values - 20.0)
    expected = nh4_5cell.values * knit_tc * expected_inhib
    np.testing.assert_allclose(np.asarray(rate), expected, rtol=1e-6)


def test_default_nitrogen_uses_v3_defaults_for_denitrification(
    no3_5cell, water_temp_5cell, dox_5cell
):
    """Phase 9.A.2 audit finding N10: default Nitrogen() reads ``kdnit_20=0.002``
    / ``kdnit_theta=1.045`` from NITROGEN_DEFAULTS (not legacy 1.0/1.0).

    Phase 9.E correction: ``kdnit_theta`` was 1.08 in v1/v3 but the
    Fortran-aligned canonical value is 1.045 (modNitrogen.f90:95). The
    1.08 was transposed with ``vno3_theta`` during v1's port from
    Fortran. Updated by Phase 9.E with regression coverage in
    ``test_phase9e_nitrogen_theta_corrections``.
    """
    inst = Nitrogen()
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False

    assert inst.kdnit_20 == 0.002
    assert inst.kdnit_theta == 1.045

    rate = inst.nitrate_denitrification(
        dissolved_oxygen=dox_5cell,
        half_saturation_oxygen=inst.KsOxdn,
        nitrate=no3_5cell,
        temperature=water_temp_5cell,
    )
    kdnit_tc = 0.002 * 1.045 ** (water_temp_5cell.values - 20.0)
    expected = (
        no3_5cell.values
        * kdnit_tc
        * (1.0 - dox_5cell.values / (dox_5cell.values + inst.KsOxdn))
    )
    np.testing.assert_allclose(np.asarray(rate), expected, rtol=1e-6)


def test_default_nitrogen_sediment_rates_zero_at_v3_defaults(
    depth_5cell, water_temp_5cell, no3_5cell
):
    """Phase 9.A.2 audit findings N4, N11: default Nitrogen() reads
    ``rnh4_20=0`` / ``vno3_20=0`` (v1/Fortran defaults), so sediment NH4
    release and sediment NO3 denitrification are silently zero by default
    (matches v1/Fortran). Pre-fix the legacy default 1.0/d injected a
    1/depth source/sink at every step.
    """
    inst = Nitrogen()
    inst.use_nitrate = True
    inst.use_ammonium = True

    assert inst.rnh4_20 == 0.0
    assert inst.vno3_20 == 0.0

    nh4_from_bed = inst.ammonium_from_bed(depth=depth_5cell, temperature=water_temp_5cell)
    np.testing.assert_allclose(np.asarray(nh4_from_bed), 0.0)

    no3_bed_denit = inst.nitrate_bed_denitrification(
        depth=depth_5cell, nitrate=no3_5cell, temperature=water_temp_5cell
    )
    np.testing.assert_allclose(np.asarray(no3_bed_denit), 0.0)


def test_phantom_ammonium_decay_term_dropped_from_change_ammonium(
    nh4_5cell, no3_5cell, water_temp_5cell, depth_5cell, dox_5cell
):
    """Phase 9.A.2 audit finding N2: default-instantiated Nitrogen() does NOT
    inject the phantom ``ammonium_decay_nitrate`` source into the NH4 budget.

    Pre-fix, the term ``ammonium_decay_rate=1.0/d * NH4`` was added as a
    *positive source* with no v1/Fortran analogue, causing NH4 to grow
    without bound at default kwargs. Post-fix, the term is removed from
    ``change_ammonium`` regardless of the legacy ``ammonium_decay_rate``
    value.
    """
    # Build a default Nitrogen but intentionally set ammonium_decay_rate=99
    # to confirm the decay term is fully dropped (not just zeroed).
    inst = Nitrogen(ammonium_decay_rate=99.0)
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False

    rate = inst.change_ammonium(
        nitrate=no3_5cell,
        ammonium=nh4_5cell,
        temperature=water_temp_5cell,
        depth=depth_5cell,
        oxygen_dissolved=dox_5cell,
    )

    # Expected (post-fix): -nitrification + bed (+ orgn hydrolysis if present).
    # With default NITROGEN_DEFAULTS (rnh4_20=0), bed term = 0; OrgN absent.
    expected = -inst.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
    np.testing.assert_allclose(np.asarray(rate), np.asarray(expected), rtol=1e-6)


def test_legacy_kwargs_still_override_defaults():
    """Phase 9.A.2 wiring contract: when a legacy v2 kwarg is explicitly
    supplied, it overrides the corresponding NITROGEN_DEFAULTS value AND
    syncs onto both naming schemes (legacy attribute and DEFAULTS-key
    attribute end up with the same value).
    """
    inst = Nitrogen(
        nitrification_rate=0.5,
        nitrification_theta=1.1,
        denitrification_rate=0.05,
        sediment_ammonium_release_rate=0.2,
        sediment_denitrification_rate=0.3,
        nitrification_oxygen_inhibition_factor=0.7,
    )

    # Both attribute names share the user value.
    assert inst.nitrification_rate == 0.5
    assert inst.knit_20 == 0.5
    assert inst.nitrification_theta == 1.1
    assert inst.knit_theta == 1.1
    assert inst.denitrification_rate == 0.05
    assert inst.kdnit_20 == 0.05
    assert inst.sediment_ammonium_release_rate == 0.2
    assert inst.rnh4_20 == 0.2
    assert inst.sediment_denitrification_rate == 0.3
    assert inst.vno3_20 == 0.3
    assert inst.nitrification_oxygen_inhibition_factor == 0.7
    assert inst.KNR == 0.7


def test_nitrate_uptake_floating_algae_uses_dynamic_split(
    no3_5cell, nh4_5cell
):
    """Phase 9.A.2 audit finding N12: NO3 algal uptake uses dynamic
    ``1 - algal_nh4_uptake_fraction`` from FloatingAlgae, NOT the static
    ``float_algea_faction_uptake_from_nitrate=1.0``. NH4 + NO3 paths must
    sum to ``rna * algal_growth_rate``.
    """
    inst = Nitrogen()
    inst.use_floating_algae = True
    inst.use_benthic_algae = False
    inst.use_nitrate = True
    inst.use_ammonium = True

    # Stand-in FloatingAlgae with the two cached attributes the post-fix
    # code reads.
    class FakeFloatingAlgae:
        AWn = 7.2
        AWa = 1000.0
        algal_growth_rate = xr.DataArray(np.array([0.5, 0.6, 0.7, 0.8, 1.0]))
        algal_nh4_uptake_fraction = xr.DataArray(
            np.array([0.2, 0.3, 0.5, 0.7, 0.9])
        )

        def ammonium_growth(self):
            rna = self.AWn / self.AWa
            return self.algal_nh4_uptake_fraction * rna * self.algal_growth_rate

    fake_falgae = FakeFloatingAlgae()
    inst.floating_algae_process = fake_falgae

    nh4_uptake = inst.ammonium_floating_growth()
    no3_uptake = inst.nitrate_uptake_floating_algae(
        nitrate=no3_5cell,
        ammonium=nh4_5cell,
        algea_growth_rate=fake_falgae.algal_growth_rate,
    )

    rna = fake_falgae.AWn / fake_falgae.AWa
    expected_total = rna * fake_falgae.algal_growth_rate
    actual_total = nh4_uptake + no3_uptake

    # Mass-balance invariant: NH4_uptake + NO3_uptake == rna * AlgalGrowth.
    np.testing.assert_allclose(
        np.asarray(actual_total), np.asarray(expected_total), rtol=1e-12
    )


def test_nitrate_uptake_benthic_algae_uses_dynamic_split_and_correct_units(
    no3_5cell, nh4_5cell, depth_5cell
):
    """Phase 9.A.2 audit finding N13: rebuilt ``nitrate_uptake_benthic_algae``
    uses ``rnb = BWn/BWd`` (NOT BWn/AWa), divides by ``depth`` (was missing),
    multiplies by ``Fb`` (NOT fraction_bottom_area), and uses dynamic
    ``1 - balgae_nh4_uptake_fraction`` (NOT static 0.5).

    Mirrors the v3 Phosphorus benthic-uptake pattern.
    """
    inst = Nitrogen()
    inst.use_floating_algae = False
    inst.use_benthic_algae = True
    inst.use_nitrate = True
    inst.use_ammonium = True

    class FakeBenthicAlgae:
        BWn = 7.2
        BWd = 100.0
        Fb = 0.9
        balgae_growth_rate = xr.DataArray(np.array([0.5, 0.6, 0.7, 0.8, 1.0]))
        balgae_nh4_uptake_fraction = xr.DataArray(
            np.array([0.2, 0.3, 0.5, 0.7, 0.9])
        )

    fake_balgae = FakeBenthicAlgae()
    inst.benthic_algae_process = fake_balgae

    no3_uptake = inst.nitrate_uptake_benthic_algae(
        nitrate=no3_5cell,
        ammonium=nh4_5cell,
        algea_growth_rate=fake_balgae.balgae_growth_rate,
        depth=depth_5cell,
    )

    rnb = fake_balgae.BWn / fake_balgae.BWd  # 0.072
    expected = (
        (1.0 - fake_balgae.balgae_nh4_uptake_fraction)
        * rnb
        * fake_balgae.Fb
        * fake_balgae.balgae_growth_rate
        / depth_5cell
    )
    np.testing.assert_allclose(np.asarray(no3_uptake), np.asarray(expected), rtol=1e-12)


def test_change_nitrate_no_algae_drops_to_nitrification_minus_denit_minus_beddenit(
    n_instance, nh4_5cell, no3_5cell, water_temp_5cell, depth_5cell, dox_5cell
):
    """When use_floating_algae=use_benthic_algae=False, change_nitrate reduces to:
    ammonium_nitrification - nitrate_denitrification(half_sat=KsOxdn) - nitrate_bed_denitrification.

    Phase 2.B Bug #9 fix: ``change_nitrate`` now wires ``half_saturation_oxygen``
    from ``self.KsOxdn`` (NITROGEN_DEFAULTS, default 0.1) rather than the
    legacy hard-coded literal ``1``. The test mirrors that wiring so the
    expected matches the fixed behavior.

    Verifies the post-fix isnull NaN replacement and that the time_step_frequency
    typo (bug #3) is no longer reachable from change_nitrate.
    """
    n_instance.use_floating_algae = False
    n_instance.use_benthic_algae = False

    v2_rate = n_instance.change_nitrate(
        nitrate=no3_5cell,
        ammonium=nh4_5cell,
        temperature=water_temp_5cell,
        depth=depth_5cell,
        oxygen_dissolved=dox_5cell,
    )

    expected = (
        n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
        - n_instance.nitrate_denitrification(
            dox_5cell, n_instance.KsOxdn, no3_5cell, water_temp_5cell
        )
        - n_instance.nitrate_bed_denitrification(depth_5cell, no3_5cell, water_temp_5cell)
    )

    np.testing.assert_allclose(np.asarray(v2_rate), np.asarray(expected), rtol=1e-6)
    assert not np.any(np.isnan(np.asarray(v2_rate)))


# ---------------------------------------------------------------------------
# Phase 9.E regression: nitrogen Arrhenius theta transposition fix
# ---------------------------------------------------------------------------
# The four nitrogen Arrhenius theta values were transposed in pairs during
# v1's port from Fortran. v3 inherited the transposition until Phase 9.E.
# Evidence summary (full discussion in parameters/nitrogen.py module
# docstring and parameter_defaults_corrections.md Section 1.10):
#
#   Parameter      v1/v3 (pre-9.E)    v3 (Phase 9.E)    Fortran modNitrogen.f90
#   kon_theta      1.074              1.047             1.047 (line 89)
#   rnh4_theta     1.047              1.074             1.074 (line 82)
#   kdnit_theta    1.08               1.045             1.045 (line 95)
#   vno3_theta     1.045              1.08              1.08  (line 100)
#
# Three independent lines of evidence:
#   (1) Direct Fortran source confirms the canonical values.
#   (2) Phosphorus parallel: kop_theta=1.047 / rpo4_theta=1.074 agree
#       across v1/v3/Fortran; the nitrogen pair should mirror this and
#       does in Fortran but did not in v1/v3 pre-9.E.
#   (3) Literature convention (Chapra 1997, QUAL2K manual, EPA Bowie 1985):
#       organic-matter hydrolysis uses theta=1.047 (universal NSM1 default,
#       matches mu_max_theta, kdp_theta, krp_theta, kpoc_theta, kdoc_theta,
#       kop_theta, kpom_theta, kbod_theta); sediment-water exchange
#       velocities use ~1.074-1.08; water-column denitrification ~1.045.

def test_phase9e_kon_theta_matches_fortran():
    """Phase 9.E: kon_theta = 1.047 (matches Fortran modNitrogen.f90:89,
    matches kop_theta phosphorus-parallel, matches universal NSM1
    organic-matter Arrhenius convention)."""
    inst = Nitrogen()
    assert inst.kon_theta == 1.047


def test_phase9e_rnh4_theta_matches_fortran():
    """Phase 9.E: rnh4_theta = 1.074 (matches Fortran modNitrogen.f90:82,
    matches rpo4_theta phosphorus-parallel for sediment release)."""
    inst = Nitrogen()
    assert inst.rnh4_theta == 1.074


def test_phase9e_kdnit_theta_matches_fortran():
    """Phase 9.E: kdnit_theta = 1.045 (matches Fortran modNitrogen.f90:95
    and Chapra 1997 water-column denitrification convention)."""
    inst = Nitrogen()
    assert inst.kdnit_theta == 1.045


def test_phase9e_vno3_theta_matches_fortran():
    """Phase 9.E: vno3_theta = 1.08 (matches Fortran modNitrogen.f90:100
    sediment-denitrification settling-velocity convention)."""
    inst = Nitrogen()
    assert inst.vno3_theta == 1.08


def test_phase9e_nitrogen_theta_pairs_consistent_with_phosphorus():
    """Phase 9.E: the kon/rnh4 nitrogen pair mirrors the kop/rpo4
    phosphorus pair. Pre-9.E the pairs were transposed for nitrogen; this
    test pins the corrected parallel-process consistency."""
    from clearwater_modules_v3.parameters.phosphorus import (
        DEFAULTS as PHOSPHORUS_DEFAULTS,
    )
    from clearwater_modules_v3.parameters.nitrogen import (
        DEFAULTS as NITROGEN_DEFAULTS,
    )
    # Organic-matter hydrolysis (kon_theta vs kop_theta): must match.
    assert NITROGEN_DEFAULTS["kon_theta"] == PHOSPHORUS_DEFAULTS["kop_theta"]
    # Sediment release (rnh4_theta vs rpo4_theta): must match.
    assert NITROGEN_DEFAULTS["rnh4_theta"] == PHOSPHORUS_DEFAULTS["rpo4_theta"]


# ---------------------------------------------------------------------------
# Phase 9.E follow-up: vson_theta removed; OrgN settling matches v1/Fortran
# ---------------------------------------------------------------------------
# v3 had previously applied an Arrhenius temperature correction to OrgN
# settling (vson_tc = arrhenius_correction(T, vson_20, vson_theta), with
# vson_theta=1.024) added by Phase 1.2 by analogy with rate-constant theta
# values. Phase 2.B's docstring claimed "parity with v1" but the parity
# claim was false: both v1 (processes.py:1333) and Fortran
# (modNitrogen.f90:233) use raw `vson` without Arrhenius correction.
# Fortran's deliberate type distinction confirms the convention: rate
# constants are TempCorrectionStruct (get Arrhenius); settling velocities
# are plain real (no Arrhenius). Phase 9.E removed vson_theta and changed
# the Process to use raw vson_20. See parameter_defaults_corrections.md
# Section 1.12.

def test_phase9e_orgn_settling_matches_v1_no_arrhenius(
    n_instance, water_temp_5cell, depth_5cell
):
    """Phase 9.E: v3 OrgN settling must match v1 OrgN_Settling exactly,
    with no Arrhenius temperature correction. Verifies parity at multiple
    temperatures (Arrhenius would have produced a temperature-varying
    deviation; raw vson produces an exact match)."""
    from clearwater_modules.nsm1 import processes as v1

    # Set vson_20 to a non-trivial value so the settling term is non-zero.
    n_instance.vson_20 = 0.05  # m/d
    n_instance.use_OrgN = True
    organic_nitrogen = xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))

    # v3 settling at multiple water temperatures: should be invariant
    # because the formula is vson_20 / depth * OrgN (no temperature
    # term).
    rate_at_5C = n_instance.organic_nitrogen_settling(
        organic_nitrogen=organic_nitrogen,
        temperature=xr.DataArray(np.full(5, 5.0)),
        depth=depth_5cell,
    )
    rate_at_30C = n_instance.organic_nitrogen_settling(
        organic_nitrogen=organic_nitrogen,
        temperature=xr.DataArray(np.full(5, 30.0)),
        depth=depth_5cell,
    )
    np.testing.assert_allclose(
        np.asarray(rate_at_5C),
        np.asarray(rate_at_30C),
        rtol=1e-12,
        err_msg=(
            "Phase 9.E: OrgN settling must be temperature-invariant "
            "(matches v1 raw-vson form, not the v3-pre-9.E "
            "Arrhenius-corrected form)."
        ),
    )

    # And confirm v3 == v1 OrgN_Settling exactly at any temperature.
    v1_rate = v1.OrgN_Settling(
        vson=n_instance.vson_20, depth=depth_5cell, OrgN=organic_nitrogen
    )
    np.testing.assert_allclose(
        np.asarray(rate_at_5C),
        np.asarray(v1_rate),
        rtol=1e-12,
        err_msg=(
            "Phase 9.E: v3 OrgN settling must match v1 "
            "OrgN_Settling(vson, depth, OrgN) exactly."
        ),
    )


def test_phase9e_vson_theta_removed_from_defaults():
    """Phase 9.E: vson_theta removed from NITROGEN_DEFAULTS (was 1.024,
    a v3-only addition with no v1/Fortran counterpart). Pin so any future
    re-addition requires explicit reconciliation against the Fortran/v1
    convention that settling velocities don't get Arrhenius corrections."""
    from clearwater_modules_v3.parameters.nitrogen import (
        DEFAULTS as NITROGEN_DEFAULTS,
    )

    assert "vson_theta" not in NITROGEN_DEFAULTS, (
        f"Phase 9.E removed vson_theta from NITROGEN_DEFAULTS; if a future "
        f"phase re-adds it, document the rationale and update this test. "
        f"Current keys: {sorted(NITROGEN_DEFAULTS.keys())}"
    )
    # And confirm vson_20 is still present (the actual settling velocity).
    assert "vson_20" in NITROGEN_DEFAULTS
