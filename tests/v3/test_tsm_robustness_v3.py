"""v3 TSM robustness tests (M1, M2, M3).

This file exercises the three MAJOR robustness fixes captured in
``design/clearwater_modules_v3_review_findings.md``:

- **M1** Constructor validation of stability parameters
  (``q_net_depth_ramp_ref``, ``dTdt_max_per_hour``). Silent-disable
  through NaN, negative, or zero inputs is rejected at construction time
  rather than producing a frozen or constant-pegged temperature field at
  run time.

- **M2** ``flux_sediment`` and ``sediment_temperature_change``
  degenerate-layer guard. Cells with ``sediment_thickness <= 0`` (zero
  from missing data or hotstart artifact, negative from a transport-
  coupling bug) return 0.0 so inf/NaN cannot poison
  ``water_temperature`` on adjacent wet cells.

- **M3** ``richardson_number`` NaN propagation and
  ``divide-by-zero`` warning suppression. NaN inputs propagate through
  to the stability function (visible defect, by design); ``wind_speed=0``
  no longer emits a ``RuntimeWarning``.

Inputs use Celsius for water/air temperatures and Kelvin via
``conversions.celsius_to_kelvin`` only inside the kernel, matching the
v3 contract.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.temperature import Temperature


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tsm() -> Temperature:
    """Single shared Temperature instance with v3 defaults.

    Module-scoped because none of the M2/M3 method calls below mutate
    instance state.
    """
    return Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)


def _arr(value, n: int = 1) -> xr.DataArray:
    """Wrap a scalar as a 1-D DataArray of length ``n``."""
    return xr.DataArray(np.full(n, value, dtype=float), dims=["nface"])


# ---------------------------------------------------------------------------
# M1: Constructor validation of stability parameters
# ---------------------------------------------------------------------------


class TestM1ConstructorValidation:
    """M1: ``__init__`` rejects silently-disabling values for the
    stability parameters and accepts the documented disable values."""

    # q_net_depth_ramp_ref ---------------------------------------------------

    def test_q_net_depth_ramp_ref_zero_disables(self):
        """``0.0`` is the documented disable value and must be accepted."""
        t = Temperature(
            wind_a=0.3, wind_b=1.5, wind_c=3.0, q_net_depth_ramp_ref=0.0
        )
        assert t.q_net_depth_ramp_ref == 0.0

    def test_q_net_depth_ramp_ref_default_accepted(self):
        """The v1-hardened default of 0.3 m must construct without error."""
        t = Temperature(
            wind_a=0.3, wind_b=1.5, wind_c=3.0, q_net_depth_ramp_ref=0.3
        )
        assert t.q_net_depth_ramp_ref == pytest.approx(0.3)

    def test_q_net_depth_ramp_ref_negative_rejected(self):
        """Negative values would silently disable the ramp (the
        ``if > 0.0`` branch is False); reject them."""
        with pytest.raises(ValueError, match="q_net_depth_ramp_ref"):
            Temperature(
                wind_a=0.3, wind_b=1.5, wind_c=3.0, q_net_depth_ramp_ref=-1.0
            )

    def test_q_net_depth_ramp_ref_nan_rejected(self):
        """NaN comparisons return False, so NaN would silently disable."""
        with pytest.raises(ValueError, match="q_net_depth_ramp_ref"):
            Temperature(
                wind_a=0.3,
                wind_b=1.5,
                wind_c=3.0,
                q_net_depth_ramp_ref=float("nan"),
            )

    def test_q_net_depth_ramp_ref_inf_rejected(self):
        """``+inf`` is not the documented disable value for the ramp; the
        documented disable value is ``0.0``. Reject ``+inf`` so a typo
        does not produce a near-no-op ramp instead of a true disable."""
        with pytest.raises(ValueError, match="q_net_depth_ramp_ref"):
            Temperature(
                wind_a=0.3,
                wind_b=1.5,
                wind_c=3.0,
                q_net_depth_ramp_ref=float("inf"),
            )

    # dTdt_max_per_hour ------------------------------------------------------

    def test_dTdt_max_per_hour_inf_disables(self):
        """``+inf`` is the documented disable value for the rate cap."""
        t = Temperature(
            wind_a=0.3, wind_b=1.5, wind_c=3.0, dTdt_max_per_hour=float("inf")
        )
        assert t.dTdt_max_per_hour == float("inf")

    def test_dTdt_max_per_hour_default_accepted(self):
        """The v1-hardened default of 5.0 K/hr must construct without error."""
        t = Temperature(
            wind_a=0.3, wind_b=1.5, wind_c=3.0, dTdt_max_per_hour=5.0
        )
        assert t.dTdt_max_per_hour == pytest.approx(5.0)

    def test_dTdt_max_per_hour_zero_rejected(self):
        """Zero would freeze the temperature field (every per-substep
        delta T clipped to zero). Reject."""
        with pytest.raises(ValueError, match="dTdt_max_per_hour"):
            Temperature(
                wind_a=0.3, wind_b=1.5, wind_c=3.0, dTdt_max_per_hour=0.0
            )

    def test_dTdt_max_per_hour_negative_rejected(self):
        """Negative cap produces a constant-pegged field via
        ``np.maximum(-cap, np.minimum(cap, ...))`` with cap < 0. Reject."""
        with pytest.raises(ValueError, match="dTdt_max_per_hour"):
            Temperature(
                wind_a=0.3, wind_b=1.5, wind_c=3.0, dTdt_max_per_hour=-1.0
            )

    def test_dTdt_max_per_hour_nan_rejected(self):
        """NaN cap propagates through the per-substep clip via
        ``np.maximum``/``np.minimum``, poisoning every cell. Reject."""
        with pytest.raises(ValueError, match="dTdt_max_per_hour"):
            Temperature(
                wind_a=0.3,
                wind_b=1.5,
                wind_c=3.0,
                dTdt_max_per_hour=float("nan"),
            )


# ---------------------------------------------------------------------------
# M2: flux_sediment and sediment_temperature_change degenerate-layer guard
# ---------------------------------------------------------------------------


class TestM2SedimentThicknessGuard:
    """M2: ``flux_sediment`` and ``sediment_temperature_change`` return
    0.0 on cells where ``sediment_thickness <= 0``."""

    # flux_sediment ----------------------------------------------------------

    def test_flux_sediment_zero_thickness_returns_zero(self, tsm):
        """``sediment_thickness == 0`` -> zero flux (no inf/NaN)."""
        flux = tsm.flux_sediment(
            water_temperature=_arr(10.0),
            sediment_temperature=_arr(20.0),
            sediment_thickness=_arr(0.0),
        )
        assert np.all(np.isfinite(flux.values))
        assert float(flux.values[0]) == 0.0

    def test_flux_sediment_negative_thickness_returns_zero(self, tsm):
        """``sediment_thickness < 0`` -> zero flux."""
        flux = tsm.flux_sediment(
            water_temperature=_arr(10.0),
            sediment_temperature=_arr(20.0),
            sediment_thickness=_arr(-0.1),
        )
        assert np.all(np.isfinite(flux.values))
        assert float(flux.values[0]) == 0.0

    def test_flux_sediment_positive_thickness_returns_formula(self, tsm):
        """``sediment_thickness > 0`` -> the documented formula value.

        Formula: ``rho * cp * alpha / 0.5 / h * (T_sed - T_water) / 86400``.
        Default ``rho=1600``, ``cp=1673``, ``alpha=0.0432`` (m^2/day),
        ``h=0.1``, ``dT=10`` -> positive flux into the water column.
        """
        flux = tsm.flux_sediment(
            water_temperature=_arr(10.0),
            sediment_temperature=_arr(20.0),
            sediment_thickness=_arr(0.1),
        )
        expected = (
            1600.0
            * 1673.0
            * 0.0432
            / 0.5
            / 0.1
            * (20.0 - 10.0)
            / 86400.0
        )
        assert float(flux.values[0]) == pytest.approx(expected, rel=1e-12)
        assert float(flux.values[0]) > 0.0

    def test_flux_sediment_multicell_mixed_thickness(self, tsm):
        """Vector input with mixed [0.0, 0.1, -0.05] thicknesses returns
        [0.0, formula_value, 0.0]; positive cell is unaffected by the
        guard on its neighbors."""
        thickness = xr.DataArray(
            np.array([0.0, 0.1, -0.05], dtype=float), dims=["nface"]
        )
        flux = tsm.flux_sediment(
            water_temperature=_arr(10.0, n=3),
            sediment_temperature=_arr(20.0, n=3),
            sediment_thickness=thickness,
        )
        values = flux.values
        assert np.all(np.isfinite(values))
        assert values[0] == 0.0
        assert values[2] == 0.0
        expected = (
            1600.0
            * 1673.0
            * 0.0432
            / 0.5
            / 0.1
            * (20.0 - 10.0)
            / 86400.0
        )
        assert values[1] == pytest.approx(expected, rel=1e-12)

    # sediment_temperature_change -------------------------------------------

    def test_sediment_temperature_change_zero_thickness_returns_zero(
        self, tsm
    ):
        """``sediment_thickness == 0`` -> zero delta T (no inf/NaN from
        the ``/h^2`` divisor)."""
        delta = tsm.sediment_temperature_change(
            water_temperature=_arr(20.0),
            sediment_temperature=_arr(15.0),
            sediment_thickness=_arr(0.0),
        )
        assert np.all(np.isfinite(delta.values))
        assert float(delta.values[0]) == 0.0

    def test_sediment_temperature_change_negative_thickness_returns_zero(
        self, tsm
    ):
        """``sediment_thickness < 0`` -> zero delta T."""
        delta = tsm.sediment_temperature_change(
            water_temperature=_arr(20.0),
            sediment_temperature=_arr(15.0),
            sediment_thickness=_arr(-0.05),
        )
        assert np.all(np.isfinite(delta.values))
        assert float(delta.values[0]) == 0.0

    def test_sediment_temperature_change_multicell_mixed_thickness(
        self, tsm
    ):
        """Same multi-cell guard pattern verified on the sediment-side
        update."""
        thickness = xr.DataArray(
            np.array([0.0, 0.1, -0.05], dtype=float), dims=["nface"]
        )
        delta = tsm.sediment_temperature_change(
            water_temperature=_arr(20.0, n=3),
            sediment_temperature=_arr(15.0, n=3),
            sediment_thickness=thickness,
        )
        values = delta.values
        assert np.all(np.isfinite(values))
        assert values[0] == 0.0
        assert values[2] == 0.0
        # Middle cell: positive water-sediment temperature gradient, so
        # the sediment warms.
        assert values[1] > 0.0


# ---------------------------------------------------------------------------
# M3: richardson_number NaN propagation and divide-by-zero suppression
# ---------------------------------------------------------------------------


class TestM3RichardsonRobustness:
    """M3: ``richardson_number`` propagates NaN visibly and does not
    emit a ``divide by zero`` ``RuntimeWarning`` on calm-wind cells."""

    def test_zero_wind_no_runtime_warning(self, tsm):
        """``wind_speed = 0`` produces ``-inf`` from the division (which
        is then clamped to ``-1.0``), but the division itself must not
        emit a ``RuntimeWarning``. v1's posture matches this: numpy's
        ``errstate`` suppresses the divide and invalid warnings."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            rn, rf = tsm.richardson_number(
                wind_speed=_arr(0.0),
                density_air_sat=_arr(1.0),
                density_air=_arr(1.0),
            )
        runtime = [w for w in caught if issubclass(w.category, RuntimeWarning)]
        assert runtime == [], (
            f"Unexpected RuntimeWarning(s): "
            f"{[str(w.message) for w in runtime]}"
        )
        # Both densities are equal so the numerator is zero; the
        # division 0/0 yields NaN, which the M3 clamp preserves
        # (visible defect). The NaN survives and propagates to the
        # stability function. Document and assert this behavior.
        assert np.isnan(float(rn.values[0]))
        assert np.isnan(float(rf.values[0]))

    def test_zero_wind_with_buoyancy_clamps_to_minus_one(self, tsm):
        """When ``density_air > density_air_sat`` (positive buoyancy)
        and ``wind_speed = 0``, the bare formula
        ``g * (density_air - density_air_sat) * 2 / (density_air * u^2)``
        diverges. Note that ``constants.GRAVITY = -9.806`` (negative in
        the v3 sign convention), so the numerator is *negative* when
        ``density_air > density_air_sat``, and the limit goes to
        ``-inf`` rather than ``+inf``. The M3 clamp pins the result to
        ``-1.0`` (the lower bound). Verifies the clamp is reached
        without warnings — and pins the sign convention so a future
        change to ``constants.GRAVITY`` would surface here."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            rn, rf = tsm.richardson_number(
                wind_speed=_arr(0.0),
                density_air_sat=_arr(1.0),
                density_air=_arr(1.2),
            )
        runtime = [w for w in caught if issubclass(w.category, RuntimeWarning)]
        assert runtime == []
        # density_air - density_air_sat = +0.2 > 0; gravity negative;
        # numerator -3.92, denominator 0 -> -inf -> clamped to -1.0.
        assert float(rn.values[0]) == pytest.approx(-1.0)
        # rn = -1.0 lies in the unstable regime (rn < 0 AND rn < -0.01):
        # rf = (1 - 22*-1)^0.80 = 23^0.80.
        assert float(rf.values[0]) == pytest.approx(
            (1.0 - 22.0 * -1.0) ** 0.80, rel=1e-12
        )

    def test_nan_wind_propagates_visibly(self, tsm):
        """NaN ``wind_speed`` from missing meteorology forcing must
        produce NaN ``richardson_number`` and NaN ``richardson_function``
        rather than a silently-clamped finite value. This is the
        visible-defect contract from the M3 fix."""
        rn, rf = tsm.richardson_number(
            wind_speed=_arr(float("nan")),
            density_air_sat=_arr(1.0),
            density_air=_arr(1.1),
        )
        assert np.isnan(float(rn.values[0])), (
            "NaN wind must propagate to richardson_number; silent "
            "clamping to a finite value would mask the bad forcing."
        )
        assert np.isnan(float(rf.values[0])), (
            "NaN richardson_number must propagate to richardson_function."
        )

    def test_zero_and_nan_wind_do_not_raise(self, tsm):
        """Neither ``wind_speed = 0`` nor ``wind_speed = NaN`` may raise
        an exception. Both must return arrays of the same shape as the
        inputs."""
        # wind = 0 (non-buoyant, so 0/0 -> NaN survives; non-raising).
        rn0, rf0 = tsm.richardson_number(
            wind_speed=_arr(0.0),
            density_air_sat=_arr(1.0),
            density_air=_arr(1.0),
        )
        assert rn0.shape == (1,)
        assert rf0.shape == (1,)

        # wind = NaN.
        rn_nan, rf_nan = tsm.richardson_number(
            wind_speed=_arr(float("nan")),
            density_air_sat=_arr(1.0),
            density_air=_arr(1.0),
        )
        assert rn_nan.shape == (1,)
        assert rf_nan.shape == (1,)
