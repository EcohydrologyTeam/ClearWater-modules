"""NSM1-SCI-N1 (MAJOR) regression: denitrification produces exactly
**1 equivalent of alkalinity per mole of NO3-N reduced**, not 4.

Gold-standard spec Workstream A2.

The upstream NSM1 Fortran (``modAlkalinity.f90:54``), v1, and pre-fix
v3 all carried ``r_alkden = 4/14/1000`` eq/mg-N — 4x the stoichiometric
value (CE-QUAL-W2 ``water-quality.f90:3157``; Stumm & Morgan). Because
the error is shared at every stage it is invisible to Fortran<->v1 and
v1<->v3 parity. v3 corrects it to ``1/14/1000`` as a deliberate,
reference-anchored divergence from upstream.

Non-shared-path contract (spec Section 1(4)): the closed-system
expectation here is built from **independently hardcoded literal
constants** (1 eq/mol-N, 14000 mg-N/mol, 50000 mg-CaCO3/eq), NOT by
reading ``alk.r_alkden`` or importing the parameter DEFAULTS.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.alkalinity import Alkalinity
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


# --- Independently hardcoded reference constants (no shared symbol). ---
EQ_PER_MOL_N = 1.0          # denitrification: 1 eq alkalinity / mol NO3-N
N_MOLAR_MASS_MG = 14000.0   # mg-N per mol N
EQ_TO_MG_CACO3 = 50000.0    # mg-CaCO3 per eq
R_ALKDEN_CORRECT = 1.0 / 14.0 / 1000.0   # eq/mg-N (= EQ_PER_MOL_N / N_MOLAR_MASS_MG)
R_ALKDEN_UPSTREAM_DEFECT = 4.0 / 14.0 / 1000.0  # Fortran/v1 (4x too high)

DT = timedelta(minutes=5)
DT_DAYS = DT.total_seconds() / 86400.0


@dataclass
class _MockNitrogen:
    nitrification_flux_rate: xr.DataArray
    denitrification_flux_rate: xr.DataArray


def _registry(alk0: float) -> InMemoryRegistry:
    reg = InMemoryRegistry()
    reg.register("alkalinity", xr.DataArray(np.array([alk0]), dims="cell"))
    reg.register(
        "water_temperature", xr.DataArray(np.array([20.0]), dims="cell")
    )
    reg.register("depth", xr.DataArray(np.array([1.0]), dims="cell"))
    return reg


def test_scin1_param_default_is_one_eq_per_mol_n():
    """r_alkden default == 1/14/1000 (independently hardcoded), NOT the
    upstream-defect 4/14/1000."""
    alk = Alkalinity(time_step=DT)
    assert alk.r_alkden == R_ALKDEN_CORRECT
    assert alk.r_alkden != R_ALKDEN_UPSTREAM_DEFECT
    np.testing.assert_allclose(
        alk.r_alkden, EQ_PER_MOL_N / N_MOLAR_MASS_MG, rtol=1e-15
    )


def test_scin1_closed_system_one_eq_alkalinity_per_mol_n_denitrified():
    """Closed box, only denitrification active: the alkalinity produced
    equals exactly 1 eq per mol NO3-N reduced (not 4)."""
    denit_flux = 0.7  # mg-N/L/d (water-column NO3 denitrification)

    # Independently computed expectation (v1-mirror, hardcoded literals).
    mol_n_denitrified = denit_flux * DT_DAYS / N_MOLAR_MASS_MG  # mol-N/L
    expected_delta = (
        EQ_PER_MOL_N * mol_n_denitrified * EQ_TO_MG_CACO3
    )  # mg-CaCO3/L; denitrification is an alkalinity SOURCE (+)
    upstream_defect_delta = 4.0 * expected_delta

    alk = Alkalinity(time_step=DT)
    alk.diagnostics = Diagnostics()
    alk.use_nitrogen = True
    alk.nitrogen_process = _MockNitrogen(
        nitrification_flux_rate=xr.DataArray(np.array([0.0]), dims="cell"),
        denitrification_flux_rate=xr.DataArray(
            np.array([denit_flux]), dims="cell"
        ),
    )

    reg = _registry(100.0)
    t = datetime(2026, 5, 16)
    alk.run(t, reg)

    alk_final = np.asarray(reg.get_at_time("alkalinity", t))
    delta = alk_final - 100.0

    np.testing.assert_allclose(delta, expected_delta, rtol=1e-12)
    # Hard anti-regression: must NOT be the 4-eq upstream-defect value.
    assert np.all(np.abs(delta - upstream_defect_delta) > 1e-12)

    # The per-day source rate is exactly 1 eq/mol-N too.
    expected_rate = (
        EQ_PER_MOL_N * (denit_flux / N_MOLAR_MASS_MG) * EQ_TO_MG_CACO3
    )
    np.testing.assert_allclose(
        np.asarray(alk.alk_denitrification_rate), expected_rate, rtol=1e-12
    )
