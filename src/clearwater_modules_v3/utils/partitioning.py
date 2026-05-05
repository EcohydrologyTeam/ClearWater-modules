"""v3 NSM1 sorption partitioning utility.

Stateless partitioning primitive ported from v1 NSM1's shared library
(``clearwater_modules/shared/processes.py``). Provides the dissolved
fraction of total inorganic phosphorus (TIP) given a partition coefficient
and particulate solids concentration.

Phase 9.B audit correction (``design/clearwater_modules_v3_nsm1_audit_c_dox.md``):
the v1 ``shared.processes.fdp`` (and the inline copy in
``nsm1.processes``) writes ``1 / (1 + kdpo4 * Solid / 0.000001)``. The
unit factor is inverted: ``kdpo4 [L/kg] * Solid [mg/L]`` carries units
``mg/kg``; the conversion to dimensionless requires multiplying by
``1e-6 kg/mg`` (equivalently dividing by ``1e6 mg/kg``). Fortran
``modGlobalParam.f90:228`` writes ``kdpo4(i,r) * Solid(i) / 1.0E6`` which
is the correct direction. v3 re-derives the unit factor (``* 1e-6`` ==
``/ 1e6``) so the dissolved fraction reduces toward 0 only at
*physically realistic* particulate-load magnitudes; under the v1 form
even tiny ``kdpo4 * Solid`` products drive ``fdp`` to ~0.
"""

import xarray as xr


def fdp(
    use_TIP: xr.DataArray,
    Solid: xr.DataArray,
    kdpo4: xr.DataArray,
) -> xr.DataArray:
    """Dissolved fraction of total inorganic phosphorus.

    Computes the dissolved fraction via the linear-equilibrium sorption
    isotherm ``fdp = 1 / (1 + kdpo4 * Solid * 1e-6)``. When ``use_TIP`` is
    false, returns 0 to suppress all TIP-related kinetics downstream.

    Dimensional analysis: ``kdpo4`` is in L/kg; ``Solid`` is in mg/L; the
    product ``kdpo4 * Solid`` carries units ``mg/kg``. Multiplying by
    ``1e-6 kg/mg`` yields a dimensionless ratio (the mass of sorbed
    phosphorus per total). v1 has the unit factor inverted (``/ 1e-6``);
    v3 follows Fortran (``modGlobalParam.f90:228``, ``/ 1.0E6``) which
    is the dimensionally correct form.

    Args:
        use_TIP | bool | total inorganic phosphorus module switch.
        Solid | mg/L | inorganic suspended solids concentration.
        kdpo4 | L/kg | solid-water partition coefficient for orthophosphate.

    Returns:
        dimensionless | dissolved fraction of TIP, in [0, 1].
    """
    return xr.where(use_TIP, 1.0 / (1.0 + kdpo4 * Solid * 1.0e-6), 0.0)
