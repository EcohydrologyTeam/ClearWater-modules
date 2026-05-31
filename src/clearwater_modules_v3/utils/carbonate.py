"""Carbonate-system equilibrium solver (S4-3).

Pure, vectorized (numpy / xarray) helpers for the NSM2 pH / carbonate
speciation solver. Kept in ``utils`` so the numerics are unit-testable
in isolation (Tier-3 analytical equilibria; Tier-5 NSM2-Fortran parity
at I=0).

Provenance / decisions (``design/nsm2_alkalinity_ph_fortran_extraction.md``):

* **Kw / K1 / K2** are the NSM2 Fortran ``modAlkalinity``
  ``ComputeAlkalinityDerivedVariables`` temperature formulas VERBATIM
  (freshwater, I=0).
* **D-A-3 ionic-strength correction:** a Davies-equation activity
  correction is applied to the apparent constants. The Davies bracket
  ``sqrt(I)/(1+sqrt(I)) - 0.3*I`` is **exactly 0 at I=0**, so the
  correction vanishes and the I=0 path collapses byte-exactly to the
  NSM2 constants (the mandatory carve-out -> Tier-5 parity at I=0 for
  ANY value of the Davies ``A``). ``A=0.5092`` (25 C Debye-Huckel);
  FIXME(s4-3-review): A(T) refinement and the apparent-constant
  gamma-ratio assignment vs Millero / CE-QUAL-W2 are the Tier-3 /
  Tier-6 validation targets (they do not affect the I=0 identity).
* **D-A-4 graceful failure:** Newton-Raphson (the NSM2 analytic step)
  -> Bisection fallback over [3, 13] for non-converged cells -> hold
  the previous pH (or 7.0) for any still-bad cell. Never raises /
  aborts (the NSM2 Fortran ``stop`` is replaced; v3 clip-with-log
  convention). Tier-5 parity is asserted on the convergent path,
  where v3 == NSM2 exactly.
* **f_NH3** uses the Emerson et al. (1975) ammonium pKa(T) =
  0.09018 + 2729.92/Tk (the canonical aquatic-toxicity relation; W2 /
  QUAL2K). Tier-6 W2 cross-check.

Unit conventions: v3 ``dic`` is mg-C/L (Phase 9.E) and ``alkalinity``
is mg-CaCO3/L. The NSM2 charge-balance residual is in mol/L and eq/L;
``dic_mol = dic_mgC / 12000`` and ``alk_eq = alk_mgCaCO3 / 50000``
(the NSM2 ``modAlkalinity`` conventions).
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from clearwater_data.custom_types import ArrayLike

MG_C_PER_MOL = 12000.0          # mg-C per mol-C
MG_CACO3_PER_EQ = 50000.0       # NSM2 modAlkalinity Alk/50000 (eq/L)
_DAVIES_A = 0.5092              # 25 C Debye-Huckel; see module docstring


def _log10(x: ArrayLike) -> ArrayLike:
    return np.log10(x)


def apparent_constants(
    t_water_c: ArrayLike,
    ionic_strength: ArrayLike = 0.0,
) -> tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Return ``(Kw, K1, K2)`` (mol/L scale).

    NSM2 ``modAlkalinity`` T-only formulas, then a Davies
    ionic-strength activity correction that is identically 1 at I=0
    (=> byte-exact NSM2 at freshwater).
    """
    tk = t_water_c + 273.15
    kw0 = 10.0 ** (
        -4787.3 / tk - 7.1321 * _log10(tk) - 0.010365 * tk + 22.80
    )
    k10 = 10.0 ** (
        -356.3094
        - 0.06091964 * tk
        + 21834.37 / tk
        + 126.8339 * _log10(tk)
        - 1684915.0 / (tk ** 2)
    )
    k20 = 10.0 ** (
        -107.8871
        - 0.03252849 * tk
        + 5151.79 / tk
        + 38.92561 * _log10(tk)
        - 563713.9 / (tk ** 2)
    )
    # Davies: log10(gamma_z) = -A z^2 (sqrt(I)/(1+sqrt(I)) - 0.3 I).
    # bracket >= 0, exactly 0 at I=0 -> gamma=1 -> K'=K0 (byte-exact
    # NSM2; the D-A-3 mandatory carve-out, robust to A).
    i = np.maximum(ionic_strength, 0.0)
    s_i = np.sqrt(i)
    bracket = s_i / (1.0 + s_i) - 0.3 * i
    g1 = 10.0 ** (-_DAVIES_A * 1.0 * bracket)   # z=1 (H+, HCO3-, OH-)
    g2 = 10.0 ** (-_DAVIES_A * 4.0 * bracket)   # z=2 (CO3 2-)
    # Apparent (concentration) constants from the thermodynamic K0
    # (Stumm & Morgan; gamma_CO2(neutral) ~ 1):
    #   K1' = K10 / (gH gHCO3) = K10 / g1^2
    #   K2' = K20 gHCO3/(gH gCO3) = K20 / g2
    #   Kw' = Kw0 / (gH gOH) = Kw0 / g1^2
    kw = kw0 / (g1 * g1)
    k1 = k10 / (g1 * g1)
    k2 = k20 / g2
    return kw, k1, k2


def _residual(
    ph: ArrayLike, dic_mol: ArrayLike, alk_eq: ArrayLike,
    kw: ArrayLike, k1: ArrayLike, k2: ArrayLike,
) -> ArrayLike:
    """Charge-balance residual f(pH) (Chapra Eq. 3.58 / NSM2
    ``modAlkalinity`` ``f``)."""
    h = 10.0 ** (-ph)
    return (
        (k1 * h + 2.0 * k1 * k2) / (h * h + k1 * h + k1 * k2) * dic_mol
        + kw / h
        - h
        - alk_eq
    )


def _zeros_like(x: ArrayLike) -> ArrayLike:
    return xr.zeros_like(x) if isinstance(x, xr.DataArray) else np.zeros_like(
        np.asarray(x, dtype=float)
    )


def solve_ph(
    alkalinity_mgcaco3: ArrayLike,
    dic_mgc: ArrayLike,
    kw: ArrayLike,
    k1: ArrayLike,
    k2: ArrayLike,
    prev_ph: ArrayLike | None = None,
    *,
    imax: int = 25,
    tol: float = 1.0e-8,
) -> tuple[ArrayLike, ArrayLike]:
    """Vectorized graceful pH solve (D-A-4).

    Returns ``(pH, fallback_mask)`` where ``fallback_mask`` is True for
    cells that needed the hold-previous fallback (the caller logs a
    clip-style diagnostics event; never raises).
    """
    alk_eq = alkalinity_mgcaco3 / MG_CACO3_PER_EQ
    dic_mol = dic_mgc / MG_C_PER_MOL

    # Newton-Raphson (the NSM2 modAlkalinity analytic step), pH0=7.
    ph = _zeros_like(alk_eq) + 7.0
    for _ in range(imax):
        h = 10.0 ** (-ph)
        denom = h * h + k1 * h + k1 * k2
        # NSM2 NewtonRaphson: pH - f / (ln10 * (K1 h DIC (h^2+4 K2 h+
        # K1 K2)/denom^2 + Kw/h + h)).
        slope = np.log(10.0) * (
            k1 * h * dic_mol * (h * h + 4.0 * k2 * h + k1 * k2)
            / (denom * denom)
            + kw / h
            + h
        )
        ph = ph - _residual(ph, dic_mol, alk_eq, kw, k1, k2) / slope
        ph = np.clip(ph, 0.0, 14.0)

    converged = np.abs(
        _residual(ph, dic_mol, alk_eq, kw, k1, k2)
    ) < tol

    # Bisection fallback over [3, 13] for non-converged cells.
    lo = _zeros_like(alk_eq) + 3.0
    hi = _zeros_like(alk_eq) + 13.0
    f_lo = _residual(lo, dic_mol, alk_eq, kw, k1, k2)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f_mid = _residual(mid, dic_mol, alk_eq, kw, k1, k2)
        same = (f_lo * f_mid) > 0.0
        lo = xr.where(same, mid, lo) if isinstance(mid, xr.DataArray) \
            else np.where(same, mid, lo)
        hi = xr.where(same, hi, mid) if isinstance(mid, xr.DataArray) \
            else np.where(same, hi, mid)
        f_lo = xr.where(same, f_mid, f_lo) if isinstance(mid, xr.DataArray) \
            else np.where(same, f_mid, f_lo)
    ph_bis = 0.5 * (lo + hi)

    use_bis = ~converged
    ph = xr.where(use_bis, ph_bis, ph) if isinstance(ph, xr.DataArray) \
        else np.where(use_bis, ph_bis, ph)

    converged2 = np.abs(
        _residual(ph, dic_mol, alk_eq, kw, k1, k2)
    ) < (tol * 1.0e3)

    # Hold previous pH (or 7.0) for any still-bad cell; never raise.
    hold = prev_ph if prev_ph is not None else (_zeros_like(alk_eq) + 7.0)
    hold = np.clip(hold, 3.0, 13.0)
    fallback_mask = ~converged2
    ph = xr.where(fallback_mask, hold, ph) \
        if isinstance(ph, xr.DataArray) \
        else np.where(fallback_mask, hold, ph)
    # Final physical clamp (D-A-4): natural-water pH is ~4-10; [3, 13]
    # is the bisection bracket. Real inputs converge well inside this
    # range so the clamp is a no-op for them; it only tames a
    # mathematically-valid-but-non-physical NR root from garbage input.
    ph = np.clip(ph, 3.0, 13.0)
    return ph, fallback_mask


def speciation(
    ph: ArrayLike, dic_mgc: ArrayLike, k1: ArrayLike, k2: ArrayLike,
) -> tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Return ``([CO2*], [HCO3-], [CO3 2-])`` in mg-C/L (same units as
    the v3 ``dic`` state). alpha0/1/2 * DIC."""
    h = 10.0 ** (-ph)
    d = h * h + k1 * h + k1 * k2
    a0 = h * h / d
    a1 = k1 * h / d
    a2 = k1 * k2 / d
    return a0 * dic_mgc, a1 * dic_mgc, a2 * dic_mgc


def henrys_k_co2(t_water_c: ArrayLike) -> ArrayLike:
    """CO2 Henry's-law constant (mol/L/atm). Same v1 empirical T-only
    formula as ``carbon.henrys_k_co2`` (duplicated here to keep the
    utils layer free of a processes import; S4-2.5b will add the
    salinity term in one place once a Weiss source is available)."""
    tk = t_water_c + 273.15
    return 10.0 ** (2385.73 / tk + 0.0152642 * tk - 14.0184)


def pco2_uatm(co2star_mgc: ArrayLike, t_water_c: ArrayLike) -> ArrayLike:
    """Partial pressure of CO2 (uatm) from [CO2*] (mg-C/L)."""
    co2_mol = co2star_mgc / MG_C_PER_MOL
    return co2_mol / henrys_k_co2(t_water_c) * 1.0e6


def f_nh3(ph: ArrayLike, t_water_c: ArrayLike) -> ArrayLike:
    """Un-ionized ammonia fraction. pKa(T) = 0.09018 + 2729.92/Tk
    (Emerson et al. 1975; W2 / QUAL2K). f_NH3 = 1/(1+10^(pKa-pH))."""
    tk = t_water_c + 273.15
    pka = 0.09018 + 2729.92 / tk
    return 1.0 / (1.0 + 10.0 ** (pka - ph))
