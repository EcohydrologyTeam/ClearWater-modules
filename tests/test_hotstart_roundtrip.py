"""Round-trip hotstart tests for NSM1 (NutrientBudget) and TSM (EnergyBudget).

These tests verify that the new ``hotstart_dataset`` + ``hotstart_timestep``
constructor kwargs let an orchestrator save kernel state mid-run and resume
from it, producing bit-equivalent results to a fresh contiguous run.

Pattern (per kernel):
  1. Build kernel A from dicts; run 10 substeps.
  2. Save A.dataset.to_netcdf(tmp).
  3. Re-open the netCDF, construct kernel B via hotstart_dataset=...,
     hotstart_timestep=0; run 10 more substeps.
  4. Build a fresh kernel C from the same dicts; run 20 contiguous substeps.
  5. Assert that the final-slot state of B matches the final-slot state of C
     for every data variable, within tight float tolerance.

The kernels are deterministic given identical inputs, so np.allclose with
rtol=1e-12 is the appropriate bar.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1.model import NutrientBudget
from clearwater_modules.tsm.model import EnergyBudget


# ---------------------------------------------------------------------------
# Small fixtures: 5 cells laid out as a 1x5 (y, x) grid so dim ordering matches
# what the kernels build via initial_array's (y, x) dims.
# ---------------------------------------------------------------------------

N_CELLS = 5
N_HALF = 10  # substeps in chunk A and chunk B
N_FULL = 2 * N_HALF  # substeps in the reference contiguous run


@pytest.fixture(scope='function')
def small_array() -> xr.DataArray:
    """Return a (1, 5) array of 1.0 with (y, x) dims."""
    return xr.DataArray(
        data=np.ones((1, N_CELLS), dtype=np.float64),
        dims=['y', 'x'],
        coords={'y': [0], 'x': list(range(N_CELLS))},
        attrs={
            'long_name': 'Initial Array',
            'units': 'm',
            'description': 'Test fixture array.',
        },
    )


@pytest.fixture(scope='function')
def tsm_initial_state(small_array: xr.DataArray) -> dict:
    return {
        'water_temp_c': small_array * 20.0,
        'surface_area': small_array * 100.0,
        'volume': small_array * 1000.0,
    }


@pytest.fixture(scope='function')
def nsm1_initial_state(small_array: xr.DataArray) -> dict:
    # State variables required by NutrientBudget.
    return {
        'Ap': small_array * 1.0,
        'Ab': small_array * 1.0,
        'NH4': small_array * 0.5,
        'NO3': small_array * 0.5,
        'OrgN': small_array * 0.3,
        'N2': small_array * 0.1,
        'TIP': small_array * 0.1,
        'OrgP': small_array * 0.1,
        'POC': small_array * 1.0,
        'DOC': small_array * 1.0,
        'DIC': small_array * 1.0,
        'POM': small_array * 1.0,
        'CBOD': small_array * 5.0,
        'DOX': small_array * 8.0,
        'PX': small_array * 0.0,
        'Alk': small_array * 50.0,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _final_slot_data_vars(model, last_slot: int) -> dict[str, np.ndarray]:
    """Return {var_name: ndarray} of the dataset's final-slot temporal vars.

    Picks only data variables that have the time dim, sliced to ``last_slot``.
    Static (non-temporal) vars are skipped because they are constant and not
    a meaningful test of round-trip correctness.
    """
    out: dict[str, np.ndarray] = {}
    time_dim = model.time_dim
    for name, da in model.dataset.data_vars.items():
        if time_dim in da.dims:
            arr = da.isel({time_dim: last_slot}).values
            out[name] = np.asarray(arr)
    return out


def _assert_states_match(
    name_to_b: dict[str, np.ndarray],
    name_to_c: dict[str, np.ndarray],
    rtol: float = 1e-12,
) -> None:
    """Assert that two name->ndarray maps have identical keys and values.

    Uses np.allclose(equal_nan=True) on each variable. Reports the first
    mismatch with the variable name and max abs diff.
    """
    common = set(name_to_b) & set(name_to_c)
    only_b = set(name_to_b) - common
    only_c = set(name_to_c) - common
    assert not only_b, f"vars present only in resumed run: {sorted(only_b)}"
    assert not only_c, f"vars present only in reference run: {sorted(only_c)}"

    mismatches = []
    for var in sorted(common):
        a = name_to_b[var]
        b = name_to_c[var]
        if a.shape != b.shape:
            mismatches.append(f"{var}: shape mismatch {a.shape} vs {b.shape}")
            continue
        if not np.allclose(a, b, rtol=rtol, atol=0.0, equal_nan=True):
            # Compute max abs diff ignoring NaN positions where both are NaN.
            with np.errstate(invalid='ignore'):
                diff = np.abs(a - b)
            mismatches.append(
                f"{var}: max |diff|={np.nanmax(diff)!r} "
                f"(allclose rtol={rtol})"
            )
    assert not mismatches, "Hotstart round-trip mismatch:\n" + "\n".join(mismatches)


def _run_n_substeps(model, n: int) -> None:
    for _ in range(n):
        model.increment_timestep()


# ---------------------------------------------------------------------------
# TSM hotstart round-trip
# ---------------------------------------------------------------------------

def test_tsm_hotstart_roundtrip(tsm_initial_state, tmp_path):
    """Save TSM state mid-run, hotstart, finish — must equal a contiguous run."""

    # --- Reference: contiguous 20-substep run ------------------------------
    ref = EnergyBudget(
        time_steps=N_FULL,
        initial_state_values={k: v.copy() for k, v in tsm_initial_state.items()},
    )
    _run_n_substeps(ref, N_FULL)
    ref_final = _final_slot_data_vars(ref, N_FULL)

    # --- Chunk A: first 10 substeps, then save -----------------------------
    chunk_a = EnergyBudget(
        time_steps=N_HALF,
        initial_state_values={k: v.copy() for k, v in tsm_initial_state.items()},
    )
    _run_n_substeps(chunk_a, N_HALF)

    save_path = tmp_path / "tsm_state.nc"
    chunk_a.dataset.to_netcdf(save_path)

    # Free file handle / reload from disk to mimic process boundary.
    loaded = xr.open_dataset(save_path).load()

    # --- Chunk B: hotstart from saved dataset, run N_HALF more -------------
    # The base class's _init_from_dataset uses ``time_steps`` arg for coord
    # length and the hotstart dataset's actual sizes for data-var shape, so
    # the two must be consistent. The saved dataset has length N_HALF + 1
    # (slot 0 is the IC, slots 1..N_HALF are substep outputs), so we pass
    # ``time_steps = saved length = N_HALF + 1``. After hotstart, slot 0 of
    # the new dataset is set to the LAST slot of the saved one; we can fill
    # slots 1..N_HALF with N_HALF more substeps.
    n_saved_slots = loaded.sizes[chunk_a.time_dim]
    chunk_b = EnergyBudget(
        time_steps=n_saved_slots,
        hotstart_dataset=loaded,
        hotstart_timestep=0,
    )
    _run_n_substeps(chunk_b, N_HALF)
    b_final = _final_slot_data_vars(chunk_b, N_HALF)

    _assert_states_match(b_final, ref_final, rtol=1e-12)


# ---------------------------------------------------------------------------
# NSM1 hotstart round-trip
# ---------------------------------------------------------------------------

def test_nsm1_hotstart_roundtrip(nsm1_initial_state, tmp_path):
    """Save NSM1 state mid-run, hotstart, finish — must equal a contiguous run."""

    ref = NutrientBudget(
        time_steps=N_FULL,
        initial_state_values={k: v.copy() for k, v in nsm1_initial_state.items()},
    )
    _run_n_substeps(ref, N_FULL)
    ref_final = _final_slot_data_vars(ref, N_FULL)

    chunk_a = NutrientBudget(
        time_steps=N_HALF,
        initial_state_values={k: v.copy() for k, v in nsm1_initial_state.items()},
    )
    _run_n_substeps(chunk_a, N_HALF)

    save_path = tmp_path / "nsm1_state.nc"
    chunk_a.dataset.to_netcdf(save_path)
    loaded = xr.open_dataset(save_path).load()

    n_saved_slots = loaded.sizes[chunk_a.time_dim]
    chunk_b = NutrientBudget(
        time_steps=n_saved_slots,
        hotstart_dataset=loaded,
        hotstart_timestep=0,
    )
    _run_n_substeps(chunk_b, N_HALF)
    b_final = _final_slot_data_vars(chunk_b, N_HALF)

    _assert_states_match(b_final, ref_final, rtol=1e-12)


# ---------------------------------------------------------------------------
# Sanity tests on the new kwargs themselves
# ---------------------------------------------------------------------------

def test_tsm_hotstart_kwargs_accepted(tsm_initial_state, tmp_path):
    """The new kwargs must be accepted without error and set timestep correctly."""
    base = EnergyBudget(
        time_steps=2,
        initial_state_values={k: v.copy() for k, v in tsm_initial_state.items()},
    )
    base.increment_timestep()
    save_path = tmp_path / "tsm_kw.nc"
    base.dataset.to_netcdf(save_path)
    loaded = xr.open_dataset(save_path).load()

    n_saved = loaded.sizes[base.time_dim]
    hot = EnergyBudget(
        time_steps=n_saved,
        hotstart_dataset=loaded,
        hotstart_timestep=0,
    )
    assert hot.timestep == 0
    # Time dim sized to match the saved dataset.
    assert hot.dataset.sizes[hot.time_dim] == n_saved


def test_nsm1_hotstart_kwargs_accepted(nsm1_initial_state, tmp_path):
    base = NutrientBudget(
        time_steps=2,
        initial_state_values={k: v.copy() for k, v in nsm1_initial_state.items()},
    )
    base.increment_timestep()
    save_path = tmp_path / "nsm1_kw.nc"
    base.dataset.to_netcdf(save_path)
    loaded = xr.open_dataset(save_path).load()

    n_saved = loaded.sizes[base.time_dim]
    hot = NutrientBudget(
        time_steps=n_saved,
        hotstart_dataset=loaded,
        hotstart_timestep=0,
    )
    assert hot.timestep == 0
    assert hot.dataset.sizes[hot.time_dim] == n_saved
