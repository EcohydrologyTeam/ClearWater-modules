"""Phase 9.A.1 wiring + formula-fix verification tests.

Tier-1 anchored tests for FloatingAlgae and BenthicAlgae:

- Default-instantiated rates must match v1's reference rate functions
  evaluated at the v3 ALGAE_DEFAULTS / BALGAE_DEFAULTS values. The tests
  parameterize v1's reference functions with the current default
  constants (read from ALGAE_DEFAULTS) so updates to the defaults toward
  literature consensus do not require synchronized hand-edits to the
  hardcoded values in this file. The original purpose of these tests is
  to verify *wiring* — i.e., that the default-instantiated Process reads
  from the DEFAULTS dict rather than being shadowed to 0/1 by legacy
  kwargs (audit F1-F4 / B1-B3, B7-B9). The reference function comes
  from v1 because v1's per-rate kinetic formulas are the canonical
  reference implementation; the constants come from ALGAE_DEFAULTS
  because that is the v3 contract.

- FloatingAlgae limit_light option 1 (half-saturation) must equal v1 FL
  (audit F5: parenthesization fix).

- FloatingAlgae rate_growth option 3 (harmonic mean) must zero out only
  when FN==0 OR FP==0 (audit F14: zero-guard fix).

- BenthicAlgae limit_light option 3 (Steele) must equal v1 FLb option 3
  (audit B6: exponent sign fix).

- BenthicAlgae default-instantiated limit_nitrogen / limit_phosphorus /
  limit_density consume KsNb / KsPb / Ksb (not the inherited KsN / KsP /
  Ksb=1 from FloatingAlgae) per audit B7-B9.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae
from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae
from clearwater_modules_v3.parameters.algae import DEFAULTS as ALGAE_DEFAULTS


# ----- Fixtures shared across audit tests --------------------------------

@pytest.fixture
def algae_5cell():
    return xr.DataArray(np.array([10.0, 20.0, 30.0, 40.0, 50.0]))


@pytest.fixture
def benthic_5cell():
    return xr.DataArray(np.array([5.0, 10.0, 15.0, 20.0, 25.0]))


@pytest.fixture
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture
def depth_5cell():
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))


@pytest.fixture
def par_5cell():
    return xr.DataArray(np.array([20.0, 50.0, 100.0, 200.0, 400.0]))


# ----- F1-F4 wiring: default-instantiated FloatingAlgae rates match v1 ----

class TestFloatingAlgaeDefaultsAreLiteratureAligned:
    """Audit F1-F4: default-instantiated FloatingAlgae() must read the v3
    ALGAE_DEFAULTS values rather than the legacy-kwarg defaults (which
    previously shadowed DEFAULTS to 0 / 1).

    The reference rates are computed via v1's per-rate kinetic functions
    (the canonical reference implementation) parameterized with the
    current default constants pulled from ALGAE_DEFAULTS. The defaults
    themselves are central within published literature ranges (Bowie
    1985 EPA/600/3-85/040; QUAL2K v2.11; CE-QUAL-W2 v4.5) — see
    ``parameters/algae.py`` for the source citations. Reading from
    ALGAE_DEFAULTS keeps these tests in sync if the defaults move within
    the literature range in a future update."""

    def test_default_respiration_matches_v1(self, algae_5cell, water_temp_5cell):
        fa = FloatingAlgae()
        v2_rate = fa.rate_respiration(algae_5cell, water_temp_5cell)
        krp_tc = v1.krp_tc(
            water_temp_5cell, ALGAE_DEFAULTS["krp_20"], ALGAE_DEFAULTS["krp_theta"]
        )
        v1_rate = v1.ApRespiration(krp_tc, algae_5cell)
        np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)

    def test_default_death_matches_v1(self, algae_5cell, water_temp_5cell):
        fa = FloatingAlgae()
        v2_rate = fa.rate_death(algae_5cell, water_temp_5cell)
        kdp_tc = v1.kdp_tc(
            water_temp_5cell, ALGAE_DEFAULTS["kdp_20"], ALGAE_DEFAULTS["kdp_theta"]
        )
        v1_rate = v1.ApDeath(kdp_tc, algae_5cell)
        np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)

    def test_default_settling_matches_v1(self, algae_5cell, depth_5cell):
        fa = FloatingAlgae()
        v2_rate = fa.rate_settling(algae_5cell, depth_5cell)
        v1_rate = v1.ApSettling(ALGAE_DEFAULTS["vsap"], algae_5cell, depth_5cell)
        np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)

    def test_default_KsN_KsP_match_v3_defaults(self):
        fa = FloatingAlgae()
        # ALGAE_DEFAULTS: KsN=0.04, KsP=0.0012
        assert fa.KsN == 0.04
        assert fa.KsP == 0.0012
        assert fa.nitrogen_michaelis_menton_constant == 0.04
        assert fa.phosphorus_michaelis_menton_constant == 0.0012

    def test_default_growth_rate_uses_default_constants(
        self, algae_5cell, water_temp_5cell
    ):
        fa = FloatingAlgae()
        # All-1 limit factors so the result == mu_max_tc * Ap.
        ones = xr.ones_like(algae_5cell)
        v2_rate = fa.rate_growth(
            algae_5cell, water_temp_5cell, ones, ones, ones
        )
        mu_max_tc_val = v1.mu_max_tc(
            water_temp_5cell,
            ALGAE_DEFAULTS["mu_max_20"],
            ALGAE_DEFAULTS["mu_max_theta"],
        )
        v1_rate = v1.ApGrowth(mu_max_tc_val, algae_5cell)
        np.testing.assert_allclose(v2_rate.values, np.asarray(v1_rate), rtol=1e-6)


# ----- B1-B3, B7-B9 wiring: default-instantiated BenthicAlgae rates -------

class TestBenthicAlgaeDefaultsAreV1Aligned:
    """Audit B1-B3, B7-B9: default-instantiated BenthicAlgae() must read
    the v3 BALGAE_DEFAULTS values (mub_max_20=0.4/mub_max_theta=1.047,
    kdb_20=0.3/kdb_theta=1.047, krb_20=0.2/krb_theta=1.06, KsNb=0.25,
    KsPb=0.125, Ksb=10.0) rather than inheriting FloatingAlgae's
    pelagic-algae kwargs."""

    def test_default_respiration_matches_v1(self, benthic_5cell, water_temp_5cell):
        ba = BenthicAlgae()
        v2_rate = ba.rate_respiration(benthic_5cell, water_temp_5cell)
        # v1 krb_tc(krb_20=0.2, TwaterC, krb_theta=1.06).
        krb_tc = v1.krb_tc(0.2, water_temp_5cell, 1.06)
        v1_rate = v1.AbRespiration(krb_tc, benthic_5cell)
        np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)

    def test_default_death_matches_v1(self, benthic_5cell, water_temp_5cell):
        ba = BenthicAlgae()
        v2_rate = ba.rate_death(benthic_5cell, water_temp_5cell)
        kdb_tc = v1.kdb_tc(0.3, water_temp_5cell, 1.047)
        v1_rate = v1.AbDeath(kdb_tc, benthic_5cell)
        np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)

    def test_default_KsNb_KsPb_Ksb_match_v3_defaults(self):
        ba = BenthicAlgae()
        # BALGAE_DEFAULTS: KsNb=0.25, KsPb=0.125, Ksb=10.0
        assert ba.KsNb == 0.25
        assert ba.KsPb == 0.125
        assert ba.Ksb == 10.0
        # And the inherited legacy attribute names mirror the benthic
        # values (not the floating defaults).
        assert ba.nitrogen_michaelis_menton_constant == 0.25
        assert ba.phosphorus_michaelis_menton_constant == 0.125
        assert ba.density_michaelis_menton_constant == 10.0

    def test_default_limit_density_matches_v1_FSb(self, benthic_5cell):
        ba = BenthicAlgae()
        v2_limit = ba.limit_density(algae=benthic_5cell)
        v1_limit = v1.FSb(Ab=benthic_5cell, Ksb=10.0)
        np.testing.assert_allclose(v2_limit.values, np.asarray(v1_limit), rtol=1e-6)

    def test_default_limit_nitrogen_uses_KsNb(self):
        ba = BenthicAlgae()
        ba.use_nitrate = True
        ba.use_ammonium = True
        nh4 = xr.DataArray(np.array([0.05, 0.10, 0.15]))
        no3 = xr.DataArray(np.array([0.10, 0.20, 0.30]))
        v2 = ba.limit_nitrogen(nitrate=no3.copy(deep=True), ammonium=nh4.copy(deep=True))
        v1_ = v1.FNb(use_NH4=True, use_NO3=True, NH4=nh4, NO3=no3, KsNb=0.25)
        np.testing.assert_allclose(np.asarray(v2), np.asarray(v1_), rtol=1e-6)

    def test_default_limit_phosphorus_uses_KsPb(self):
        ba = BenthicAlgae()
        ba.use_phosphate = True
        tip = xr.DataArray(np.array([0.05, 0.10, 0.15]))
        fdp = 0.5
        v2 = ba.limit_phosphorus(concentration=tip, fraction_dissolved=fdp)
        v1_ = v1.FPb(fdp=fdp, TIP=tip, use_TIP=True, KsPb=0.125)
        np.testing.assert_allclose(np.asarray(v2), np.asarray(v1_), rtol=1e-6)

    def test_default_growth_rate_uses_v1_aligned_constants(
        self, benthic_5cell, water_temp_5cell
    ):
        ba = BenthicAlgae()
        ones = xr.ones_like(benthic_5cell)
        v2_rate = ba.rate_growth(
            benthic_5cell, water_temp_5cell, ones, ones, ones, ones
        )
        # v1 mub_max_tc(0.4, TwaterC, 1.047).
        mub_max_tc_val = v1.mub_max_tc(0.4, water_temp_5cell, 1.047)
        v1_rate = v1.AbGrowth(mub_max_tc_val, benthic_5cell)
        np.testing.assert_allclose(v2_rate.values, np.asarray(v1_rate), rtol=1e-6)


# ----- F5: FloatingAlgae limit_light option 1 parenthesization fix --------

def test_floating_algae_limit_light_option1_matches_v1_FL(
    algae_5cell, depth_5cell, par_5cell
):
    """F5: ``(1/(L*d)) * log( (KL+PAR) / (KL+PAR*exp(-Ld)) )`` --
    v3 must group log-numerator over log-denominator inside the log,
    not split across the * and / operators."""
    fa = FloatingAlgae()
    fa.use_nitrate = True
    fa.use_ammonium = True
    fa.use_phosphate = True
    fa.light_limitation_option = 1
    # Default values: KL=10.0 (DEFAULTS), light_attenuation_coefficient=1.0.
    v2_limit = fa.limit_light(
        algae=algae_5cell, depth=depth_5cell, surface_light_intensity=par_5cell
    )
    # v1 FL with the same parameters and the same inputs.
    v1_limit = v1.FL(
        Ap=algae_5cell, depth=depth_5cell, PAR=par_5cell,
        light_limitation_option=1, KL=10.0, L=1.0,
    )
    np.testing.assert_allclose(np.asarray(v2_limit), np.asarray(v1_limit), rtol=1e-6)


# ----- F14: FloatingAlgae rate_growth option 3 harmonic-mean guard --------

class TestFloatingAlgaeHarmonicGrowthGuard:
    """F14: harmonic-mean zero-guard must fire on FN==0 OR FP==0, not
    on the previous (incorrect) ``FP==1`` condition."""

    def test_harmonic_mean_FP_one_does_not_zero_growth(
        self, algae_5cell, water_temp_5cell
    ):
        """When FP==1 (P fully non-limiting), growth should NOT zero out.
        The pre-fix code shut growth down at FP==1, the most common
        steady-state P-replete condition."""
        fa = FloatingAlgae(growth_rate_option=3)
        ones = xr.ones_like(algae_5cell)
        half = xr.full_like(algae_5cell, 0.5)
        # FN=0.5, FP=1.0, FL=1.0 -> mu = mu_max_tc * 1 * 2 / (1/0.5 + 1/1) = mu_max_tc * 2/3.
        v2_rate = fa.rate_growth(algae_5cell, water_temp_5cell, ones, half, ones)
        mu_max_tc_val = v1.mu_max_tc(
            water_temp_5cell,
            ALGAE_DEFAULTS["mu_max_20"],
            ALGAE_DEFAULTS["mu_max_theta"],
        )
        v1_mu = v1.mu(
            mu_max_tc=mu_max_tc_val, growth_rate_option=3,
            FL=ones, FP=ones, FN=half,
        )
        v1_rate = v1.ApGrowth(v1_mu, algae_5cell)
        np.testing.assert_allclose(v2_rate.values, np.asarray(v1_rate), rtol=1e-6)
        # And ensure result is non-zero.
        assert (v2_rate.values > 0).all()

    def test_harmonic_mean_FP_zero_zeros_growth(
        self, algae_5cell, water_temp_5cell
    ):
        """When FP==0, growth must zero out (avoids 1/FP division)."""
        fa = FloatingAlgae(growth_rate_option=3)
        ones = xr.ones_like(algae_5cell)
        zeros = xr.zeros_like(algae_5cell)
        # FP==0 case.
        v2_rate = fa.rate_growth(algae_5cell, water_temp_5cell, zeros, ones, ones)
        np.testing.assert_array_equal(v2_rate.values, np.zeros_like(v2_rate.values))


# ----- B6: BenthicAlgae limit_light option 3 Steele exponent fix ----------

def test_benthic_algae_limit_light_option3_steele_matches_v1(
    benthic_5cell, depth_5cell, par_5cell
):
    """B6: Steele form is ``x * exp(1 - x)`` with ``x = PAR*KEXT/KLb``.
    Pre-fix used division (``x / exp(1-x)`` = ``x * exp(x-1)``)."""
    ba = BenthicAlgae(light_limitation_option=3)
    ba.use_nitrate = True
    ba.use_ammonium = True
    ba.use_phosphate = True
    v2_limit = ba.limit_light(
        algae=benthic_5cell, depth=depth_5cell, surface_light_intensity=par_5cell
    )
    # v1 FLb with KLb=10 (BALGAE_DEFAULTS), L=1.0 (default).
    v1_limit = v1.FLb(
        Ab=benthic_5cell, depth=depth_5cell, PAR=par_5cell,
        b_light_limitation_option=3, KLb=10.0, L=1.0,
    )
    np.testing.assert_allclose(np.asarray(v2_limit), np.asarray(v1_limit), rtol=1e-6)
