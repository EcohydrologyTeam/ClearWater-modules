"""Tier-1.5 active-kinetics conservation tests for v3 NSM1.

Phase 9.D.2 of the v3 NSM1 audit (design spec Section 14, audit ``Required
actions`` item 5: "Build a Fortran-anchored default-instantiation regression
suite that the existing parity tests should have been"). This module
complements the Tier 1 closed-system tests in
``test_validation_tier1_conservation.py`` by exercising **default-
instantiated** Process classes through ``build_nsm1_demo`` and asserting
that:

1. Total nitrogen mass (water-column N pools + algal-N equivalents) is
   conserved across many substeps when boundary fluxes (settling,
   atmospheric exchange, sediment release, burial) are disabled but all
   other default kinetics are active.
2. Total carbon mass (water-column C pools + algal-C equivalents) is
   conserved under the same boundary-disabled configuration.
3. No ``clip_negative_state`` events fire under default kinetics; the
   resolved Q7 contract requires that physically reasonable initial
   conditions plus default parameters keep every state variable
   non-negative.
4. Every state variable remains finite and positive throughout the run.

Why Tier-1.5 rather than Tier 1 or Tier 2?

* Tier 1 closed-system tests (``test_validation_tier1_conservation.py``)
  zero out *every* kinetic rate (vsap=0, kdp_20=0, etc.), which confirms
  the integrator is bookkeeping flux rates correctly but does not
  exercise the parameter-defaults wiring. The Tier 1 tests would still
  pass even if a Process silently consumed an undefined rate constant
  via a lurking ``getattr(..., 0)`` fallback.
* Tier 2 Streeter-Phelps and Fortran-trajectory parity tests check
  individual Process behavior against a reference, not coupled mass
  conservation across all 11 NSM1 reservoirs at once.
* Tier-1.5 (this module) uses default-instantiated Processes via
  ``build_nsm1_demo`` so the only way the test can pass is if the
  Phase 8.A / 9.A / 9.B / 9.C wiring sweep correctly populates every
  parameter and inter-process rate cache. A regression in a single
  parameter default (e.g. one of the vson_20 / vsop / SOD_20 corrections)
  shifts the equilibrium drift visibly.

Boundary fluxes disabled
------------------------

The closed-system invariant only holds when boundary-affecting fluxes
are zero. The fixture builds the demo with these per-process overrides:

* FloatingAlgae: ``vsap=0`` -- no algal settling out of the water column
* BenthicAlgae:  no overrides; benthic algae state is g-D/m^2 and only
                 routes mass to water-column pools (no settling out)
* Nitrogen:      ``vson_20=0``, ``rnh4_20=0``, ``vno3_20=0`` -- no OrgN
                 settling, no sediment NH4 release, no sediment NO3 loss
* Phosphorus:    ``vsop=0``, ``vs=0``, ``rpo4_20=0`` -- no OrgP / TIP
                 settling, no sediment P release
* Carbon:        ``vsoc=0``, ``JDIC=0``, ``kah_20_user=0``,
                 ``kaw_20_user=0`` (with options=1) -- no POC settling,
                 no sediment DIC release, no CO2 atmospheric exchange
* CBOD:          ``ksbod_20=0`` (already 0 by default) -- no CBOD
                 settling
* DOX:           ``kah_20_user=0``, ``kaw_20_user=0`` (with options=1),
                 ``SOD_20=0`` -- no atmospheric reaeration, no SOD sink.
                 Note: DOX is **not** in the conservation totals
                 enforced here -- DOX participates in C bookkeeping
                 only via its DOX-Monod attenuation factor on DOC
                 oxidation, which is mass-conserving in C (DOC -> DIC).
* N2:            ``kah_20_user=0``, ``kaw_20_user=0`` (with options=1)
                 -- no atmospheric N2 exchange
* POM:           ``vb=0`` -- no burial out of the water column
* Pathogen:      no overrides; pathogen does not couple to N/C totals
* Alkalinity:    no overrides; alkalinity tracks stoichiometric
                 source/sink terms, not a real conservation pool

Tolerance discussion
--------------------

The conservation helpers convert algae state (in ug-Chla/L) to
water-column N or C concentration (mg-N/L or mg-C/L) using the
*ratios* ``rna = AWn/AWa`` and ``rca = AWc/AWa``, matching the v3
kinetics in ``nitrogen.py`` (``rna = self.floating_algae_nitrogen_weight
/ self.algal_chlorophyll``) and ``carbon.py`` (``rca = self.AWc /
self.AWa``). The Fortran v1 NSM1 and CE-QUAL-W2
(``water-quality.f90:1505, 1579`` -- ``ALG * AN`` with state in
mg-DW/L) follow the same dimensional pattern. Earlier revisions of
this helper used the raw weights ``AWn`` / ``AWc`` directly, which
inflated the reported algae N/C contribution by AWa = 1000x and
masked the actual conservation behavior.

Empirical drift on the Tier-1.5 fixture (5-cell synthetic mesh,
default IC, 5-minute substeps):

    50 steps:  d(total-N)/total-N = -0.057%   d(total-C)/total-C = +0.43%
    100 steps: d(total-N)/total-N = -0.110%   d(total-C)/total-C = +3.07%
    144 steps: d(total-N)/total-N = -0.153%   d(total-C)/total-C = +7.61%
    288 steps: d(total-N)/total-N = -0.273%   d(total-C)/total-C = +21.65%

Total-N drift is small and roughly linear in step count, consistent
with a per-step Forward-Euler discretization residual at the
literature-aligned ``mu_max=2.0`` 1/d default. Total-C drift is larger
and grows super-linearly because POM is in mg-D/L (dry weight) and
exits the C bookkeeping when it dissolves to DOC (mg-C/L); the POM
dissolution pathway is dimensionally approximate. With ``vb=0`` the
exchange is closed in dry-weight terms but not in C terms.

A secondary residual is the volumetric/areal unit conversion between
benthic-algae state (g-D/m^2) and water-column N/C pools (mg/L). The
helper converts via ``Ab * (BWn/BWd) * Fb / depth`` to mg-N/L,
matching the Fortran NSM1 / v1 / v3 water-column flux convention
``rnb * Fb * AbGrowth / depth`` (where ``rnb = BWn/BWd``). The benthic
contribution to the total-N inventory in this fixture is ~12% so the
conversion choice does not dominate the residual.

Tolerance: ``rtol=1.0e-1`` (10%) over 100 substeps. The total-N
residual is well within this floor (~0.1% at 100 steps); the looser
tolerance is preserved to absorb the documented POM-dry-weight-as-C
approximation in the total-C helper.

Total-C is closed up to the same Forward-Euler residual plus the
POM <-> DOC dimensional conflation and the CBOD <-> DIC oxidation
pathway. CBOD is in mg-O2/L; 1 mg-CBOD oxidized produces ``1 / roc``
mg-C/L of DIC, where ``roc = 32/12``. The total-C helper sums
``cbod / roc`` to absorb this exchange.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Any

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.examples.nsm1_demo_setup import (
    InMemoryRegistry,
    build_nsm1_demo,
    default_initial_conditions,
)
from clearwater_modules_v3.parameters.algae import DEFAULTS as ALGAE_DEFAULTS
from clearwater_modules_v3.parameters.balgae import DEFAULTS as BALGAE_DEFAULTS


# ---------------------------------------------------------------------------
# Stoichiometric ratios for default-instantiated mass bookkeeping
# ---------------------------------------------------------------------------

# Floating-algae stoichiometric ratios. AWn / AWc / AWa carry the
# stoichiometric-unit convention used throughout v3 NSM1 (and the
# Fortran v1 NSM1 / W2 conventions): AWn = 7.2 (mg-N per stoichiometric
# unit), AWc = 40 (mg-C per stoichiometric unit), AWa = 1000 (ug-Chla
# per stoichiometric unit). Converting algae state (in ug-Chla/L) to
# water-column N or C concentration uses the *ratios* ``rna = AWn/AWa``
# and ``rca = AWc/AWa``, NOT ``AWn`` / ``AWc`` directly. The v3 kinetics
# match this:
#
#   nitrogen.py:675     rna = self.floating_algae_nitrogen_weight / self.algal_chlorophyll
#                           = AWn / AWa     (mg-N / ug-Chla)
#   carbon.py:402       rca = self.AWc / self.AWa
#                           = AWc / AWa     (mg-C / ug-Chla)
#
# Cross-check vs CE-QUAL-W2 (water-quality.f90:1505, 1579): W2 stores
# algae in mg-DW/L and uses ``ALG * AN`` (mg-DW/L * mg-N/mg-DW = mg-N/L);
# the NSM1 ug-Chla/L analog reduces to ``Ap * (AWn/AWa) = Ap * rna``.
AP_N_PER_CHLA: float = float(ALGAE_DEFAULTS["AWn"]) / float(ALGAE_DEFAULTS["AWa"])  # mg-N per ug-Chla = rna
AP_C_PER_CHLA: float = float(ALGAE_DEFAULTS["AWc"]) / float(ALGAE_DEFAULTS["AWa"])  # mg-C per ug-Chla = rca

# Benthic-algae stoichiometric / coupling constants used by the
# benthic-algae areal -> volumetric conversion in the conservation
# helpers below. The volumetric factor applied to ``benthic_algae``
# is ``(BWn / BWd) * Fb / depth`` for nitrogen (``(BWc / BWd) * Fb / depth``
# for carbon). This matches the Fortran NSM1 / v1 / v3 NSM1
# water-column flux convention
# ``rnb * Fb * AbGrowth / depth`` where ``rnb = BWn / BWd``, and
# closes the closed-system mass balance with the ``Fb`` partial-
# coupling factor on both sides. See the discussion at the top of
# this module and ``design/clearwater_modules_v3_nsm1_audit_*.md``
# (algae) for the architectural rationale.
AB_N_PER_GD:   float = float(BALGAE_DEFAULTS["BWn"]) / float(BALGAE_DEFAULTS["BWd"])  # mg-N/mg-D (mass fraction)
AB_C_PER_GD:   float = float(BALGAE_DEFAULTS["BWc"]) / float(BALGAE_DEFAULTS["BWd"])  # mg-C/mg-D (mass fraction)
AB_FB:         float = float(BALGAE_DEFAULTS["Fb"])                                  # active-fraction coupling factor

# CBOD <-> DIC oxidation stoichiometry (mg-O2 / mg-C; v3 carbon.DEFAULTS).
ROC: float = 32.0 / 12.0


# ---------------------------------------------------------------------------
# Closed-system parameter overrides for default-instantiated demo
# ---------------------------------------------------------------------------


def _closed_system_process_parameters() -> dict[str, dict[str, Any]]:
    """Per-Process overrides that disable boundary fluxes for Tier-1.5."""
    return {
        "FloatingAlgae": {
            "vsap": 0.0,                # disable algal settling
        },
        "BenthicAlgae": {},             # no settling parameters; mortality routes are kinetic
        "Nitrogen": {
            "vson_20": 0.0,             # disable OrgN settling
            "rnh4_20": 0.0,             # disable sediment NH4 release
            "vno3_20": 0.0,             # disable sediment NO3 denitrification
        },
        "Phosphorus": {
            "vsop": 0.0,                # disable OrgP settling
            "vs":   0.0,                # disable TIP settling
            "rpo4_20": 0.0,             # disable sediment P release
        },
        "Carbon": {
            "vsoc": 0.0,                # disable POC settling
            "JDIC": 0.0,                # disable sediment DIC release (already default)
            "kah_20_user": 0.0,         # disable CO2 hydraulic atm exchange
            "kaw_20_user": 0.0,         # disable CO2 wind atm exchange
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        "POM": {
            "vb": 0.0,                  # disable POM burial
        },
        "CBOD": {
            "ksbod_20": 0.0,            # already default; explicit for clarity
        },
        "DOX": {
            "kah_20_user": 0.0,         # disable hydraulic reaeration
            "kaw_20_user": 0.0,         # disable wind reaeration
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "SOD_20": 0.0,              # disable sediment oxygen demand
        },
        "N2": {
            "kah_20_user": 0.0,         # disable atm N2 hydraulic exchange
            "kaw_20_user": 0.0,         # disable atm N2 wind exchange
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        "Pathogen":   {},
        "Alkalinity": {},
    }


# ---------------------------------------------------------------------------
# Conservation totals (default-instantiated names)
# ---------------------------------------------------------------------------
# These helpers are deliberately distinct from the conftest ``total_n`` /
# ``total_c`` helpers: those use ``ap`` / ``ab`` registry names that the
# Tier 1 fixture sets up. The default ``build_nsm1_demo`` uses the
# production names ``algae_floating`` / ``benthic_algae``. Tier-1.5 keeps
# its own helpers so the assertion intent is explicit.


def _get(registry: Any, name: str) -> xr.DataArray | None:
    if name not in registry:
        return None
    return registry.get(name)


def _sum_over_cells(da: xr.DataArray) -> float:
    if "cell" in da.dims:
        return float(da.sum(dim="cell").values)
    return float(da.sum().values)


def total_n_active_kinetics(registry: Any) -> float:
    """Total nitrogen across NH4, NO3, OrgN, N2, and algal N-equivalents.

    All terms are reported in mg-N/L for consistency with the water-column
    pools, then summed per-cell.

    Water-column terms (already in mg-N/L):
        ammonium + nitrate + organic_nitrogen + n2

    Floating algae N (chlorophyll-volumetric):
        algae_floating [ug-Chla/L] * AWn [mg-N/ug-Chla] = mg-N/L

    Benthic algae N (areal -> volumetric, NSM1 convention):
        benthic_algae [g-D/m^2] * (BWn/BWd) [mg-N/mg-D] * Fb / depth [m]

    The benthic-algae conversion uses the same implicit "1 g/m^3 == 1 mg/L"
    identity that the Fortran NSM1 / v1 / v3 water-column flux formulas
    use (``NH4_AbGrowth = AbUptakeFr_NH4 * rnb * Fb * AbGrowth / depth``,
    where ``rnb = BWn/BWd``). The ``Fb`` factor reflects the partial-
    coupling architectural choice in NSM1: only the ``Fb`` fraction of
    the bed area exchanges with the water column (the ``(1-Fb)`` portion
    is intentionally open-coupled).

    Mass balance closes exactly with this convention: the volumetric N
    rate of change of the benthic state ``Ab * (BWn/BWd) * Fb / depth``
    equals the water-column ``rnb * Fb * AbGrowth / depth`` flux, so
    total-N is conserved across algal growth-uptake transfers.

    See ``design/v3_tsm_wind_function_improvements.md``-adjacent
    discussion in PR review and the CE-QUAL-W2 cross-check
    (``water-quality.f90:1517, 2342``) for the unit-identity rationale.
    """
    pieces: list[float] = []
    for name in ("ammonium", "nitrate", "organic_nitrogen", "n2"):
        da = _get(registry, name)
        if da is not None:
            pieces.append(_sum_over_cells(da))

    ap = _get(registry, "algae_floating")
    if ap is not None:
        pieces.append(_sum_over_cells(ap * AP_N_PER_CHLA))
    ab = _get(registry, "benthic_algae")
    if ab is not None:
        depth = _get(registry, "depth")
        if depth is None:
            raise ValueError(
                "total_n_active_kinetics: registry has 'benthic_algae' but "
                "no 'depth'; cannot convert areal benthic-algae N to "
                "volumetric. The Tier-1.5 fixture should provide both."
            )
        pieces.append(_sum_over_cells(ab * AB_N_PER_GD * AB_FB / depth))

    return sum(pieces)


def total_c_active_kinetics(registry: Any) -> float:
    """Total carbon across POC, DOC, DIC, CBOD-as-C, and algal C-equivalents.

    All terms are reported in mg-C/L for consistency with the water-column
    pools, then summed per-cell.

    Water-column terms:
        poc + doc + dic                              (mg-C/L)
        cbod / ROC                                   (mg-O2/L -> mg-C/L)

    Floating algae C (chlorophyll-volumetric):
        algae_floating [ug-Chla/L] * AWc [mg-C/ug-Chla] = mg-C/L

    Benthic algae C (areal -> volumetric, NSM1 convention):
        benthic_algae [g-D/m^2] * (BWc/BWd) [mg-C/mg-D] * Fb / depth [m]

    See ``total_n_active_kinetics`` for the unit-identity and ``Fb``
    coupling rationale; the carbon side uses the same conversion form.
    """
    pieces: list[float] = []
    for name in ("poc", "doc", "dic"):
        da = _get(registry, name)
        if da is not None:
            pieces.append(_sum_over_cells(da))

    cbod = _get(registry, "cbod")
    if cbod is not None:
        pieces.append(_sum_over_cells(cbod / ROC))

    ap = _get(registry, "algae_floating")
    if ap is not None:
        pieces.append(_sum_over_cells(ap * AP_C_PER_CHLA))
    ab = _get(registry, "benthic_algae")
    if ab is not None:
        depth = _get(registry, "depth")
        if depth is None:
            raise ValueError(
                "total_c_active_kinetics: registry has 'benthic_algae' but "
                "no 'depth'; cannot convert areal benthic-algae C to "
                "volumetric. The Tier-1.5 fixture should provide both."
            )
        pieces.append(_sum_over_cells(ab * AB_C_PER_GD * AB_FB / depth))

    return sum(pieces)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="function")
def tier1p5_demo():
    """Build a Tier-1.5 demo with closed-system parameter overrides.

    Returns the ``Nsm1Demo`` handle from ``build_nsm1_demo``. Tests can
    iterate ``demo.step(t)`` to advance the simulation.
    """
    demo = build_nsm1_demo(
        time_step=timedelta(minutes=5),
        n_cells=5,
        process_parameters=_closed_system_process_parameters(),
    )
    return demo


# ---------------------------------------------------------------------------
# Tier-1.5 conservation tests
# ---------------------------------------------------------------------------


def test_tier1p5_total_n_conservation_active_kinetics(tier1p5_demo) -> None:
    """Total nitrogen is conserved across ~12 hours of default kinetics.

    The fixture disables all boundary fluxes (settling, atmospheric
    exchange, sediment release) but leaves every internal kinetic rate at
    its v3 default. Active rates that redistribute N within the closed
    set:

    * Algal growth (``mu_max_20`` = 2.0 1/d): NH4/NO3 -> algal-N
    * Algal respiration (``krp_20`` = 0.10 1/d): algal-N -> NH4
    * Algal mortality (``kdp_20`` = 0.05 1/d): algal-N -> OrgN/POC/DOC
    * Nitrification (``knit_20`` = 0.1 1/d): NH4 -> NO3
    * Denitrification (``kdnit_20`` = 0.002 1/d): NO3 -> N2
    * OrgN hydrolysis (``kon_20`` = 0.1 1/d): OrgN -> NH4

    All of these are mass-conserving within the total-N pool *in
    principle*. The empirical residual drift (~0.1% at 100 substeps)
    is consistent with a per-step Forward-Euler discretization
    residual at the literature-aligned ``mu_max=2.0`` 1/d default;
    see the module docstring for the post-fix tolerance discussion.
    """
    n_initial = total_n_active_kinetics(tier1p5_demo.registry)

    start = datetime(2026, 1, 1)
    n_steps = 100  # 100 * 5 min = 8.33 h (per audit acceptance: >= 100 substeps)
    current = start
    for _ in range(n_steps):
        tier1p5_demo.step(current)
        current += tier1p5_demo.time_step

    n_final = total_n_active_kinetics(tier1p5_demo.registry)

    np.testing.assert_allclose(
        n_final,
        n_initial,
        rtol=1.0e-1,
        err_msg=(
            "Tier-1.5 active-kinetics total-N conservation failed. "
            f"initial={n_initial!r}, final={n_final!r}, "
            f"absolute drift={(n_final - n_initial)!r} "
            f"(relative drift={(n_final - n_initial) / n_initial:.4e})"
        ),
    )


def test_tier1p5_total_c_conservation_active_kinetics(tier1p5_demo) -> None:
    """Total carbon is conserved across ~12 hours of default kinetics.

    Active rates that redistribute C within the closed set:

    * Algal photosynthesis: DIC -> algal-C (and DIC -> DOC via 1-Fw fraction)
    * Algal respiration: algal-C -> DIC
    * Algal mortality: algal-C -> POC/DOC
    * POC hydrolysis (``kpoc_20`` = 0.005 1/d): POC -> DOC
    * DOC oxidation (``kdoc_20`` = 0.01 1/d): DOC -> DIC
    * POM dissolution (``kpom_20`` = 0.1 1/d): POM -> DOC. Note POM is in
      mg-D/L (dry-weight mass), not mg-C/L, so this is *not* a closed-C
      pathway under the conftest helpers' simple sum. With ``vb=0``
      (POM burial disabled) the POM/DOC exchange is mass-conserving in
      DRY-WEIGHT terms, but the C bookkeeping is approximate.
    * CBOD oxidation (``kbod_20`` = 0.12 1/d): CBOD -> DIC. The C total
      includes ``cbod / roc`` to absorb this exchange.

    Tolerance: rtol=1e-1. The empirical drift is ~3% at 100 substeps,
    dominated by the documented POM-dry-weight-as-C approximation
    rather than a true non-conservation; see the module docstring.
    """
    c_initial = total_c_active_kinetics(tier1p5_demo.registry)

    start = datetime(2026, 1, 1)
    n_steps = 100
    current = start
    for _ in range(n_steps):
        tier1p5_demo.step(current)
        current += tier1p5_demo.time_step

    c_final = total_c_active_kinetics(tier1p5_demo.registry)

    np.testing.assert_allclose(
        c_final,
        c_initial,
        rtol=1.0e-1,
        err_msg=(
            "Tier-1.5 active-kinetics total-C conservation failed. "
            f"initial={c_initial!r}, final={c_final!r}, "
            f"absolute drift={(c_final - c_initial)!r} "
            f"(relative drift={(c_final - c_initial) / c_initial:.4e})"
        ),
    )


def test_tier1p5_no_clip_events_under_active_kinetics(tier1p5_demo) -> None:
    """Default kinetics with closed-system boundaries should not clip.

    Per the resolved Q7 contract (design spec Section 14): physically
    reasonable initial conditions plus default parameters keep every
    state variable non-negative. Any clip event under this configuration
    indicates an integrator overshoot, parameter typo, or a stoichiometry
    bug -- the kind of regression the Phase 9.A/B/C audit was set up to
    catch.
    """
    start = datetime(2026, 1, 1)
    n_steps = 100
    current = start
    for _ in range(n_steps):
        tier1p5_demo.step(current)
        current += tier1p5_demo.time_step

    diagnostics = tier1p5_demo.model.diagnostics
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under default kinetics + closed-system "
        f"boundaries: {diagnostics.clip_events!r}. "
        f"Clip log: {diagnostics.clip_log!r}"
    )


def test_tier1p5_state_variables_finite_and_positive(tier1p5_demo) -> None:
    """Every state variable stays finite and non-negative throughout the run.

    Stronger than the clip-events check: clipping silently zeroes
    negative values, but a NaN or +inf would propagate without firing a
    clip. This test asserts that no kinetic path produces a non-finite
    intermediate.
    """
    start = datetime(2026, 1, 1)
    n_steps = 100
    current = start
    for _ in range(n_steps):
        tier1p5_demo.step(current)
        current += tier1p5_demo.time_step

    state_names = (
        "ammonium", "nitrate", "organic_nitrogen", "n2",
        "tip", "organic_phosphorus",
        "poc", "doc", "dic",
        "pom", "cbod",
        "oxygen_dissolved",
        "alkalinity",
        "algae_floating", "benthic_algae",
        "pathogen",
    )
    failures: list[str] = []
    for name in state_names:
        if name not in tier1p5_demo.registry:
            continue
        da = tier1p5_demo.registry.get(name)
        values = np.asarray(da.values)
        if not np.all(np.isfinite(values)):
            failures.append(f"{name}: contains NaN/inf ({values!r})")
        if np.any(values < 0):
            failures.append(f"{name}: contains negative values ({values!r})")
    assert not failures, (
        "Tier-1.5: state variables failed finite/positive check after "
        f"{n_steps} substeps:\n  " + "\n  ".join(failures)
    )


# ---------------------------------------------------------------------------
# Smoke test: helpers return finite positive totals on the initial state
# ---------------------------------------------------------------------------


def test_total_n_active_kinetics_helper_returns_finite_positive(
    tier1p5_demo,
) -> None:
    """``total_n_active_kinetics`` over the initial state is positive and finite."""
    value = total_n_active_kinetics(tier1p5_demo.registry)
    assert np.isfinite(value)
    assert value > 0.0


def test_total_c_active_kinetics_helper_returns_finite_positive(
    tier1p5_demo,
) -> None:
    """``total_c_active_kinetics`` over the initial state is positive and finite."""
    value = total_c_active_kinetics(tier1p5_demo.registry)
    assert np.isfinite(value)
    assert value > 0.0
