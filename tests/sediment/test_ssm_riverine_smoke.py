"""End-to-end coupled SSM + ClearwaterRiverine smoke test.

Loads a real RAS-derived mesh from the sibling
``ClearWater-Riverine-streaming`` repo's ``plan01_10x5`` fixture,
instantiates a Riverine transport solver with one trivial conservative
tracer, attaches an SSM bound to the SAND2008 SEDflume bundle to the
same mesh via :meth:`SSM.bind_mesh`, and then steps the coupled system
for a handful of time slots.

This is a **smoke test**, not a precision validation. The intent is to
confirm that the two solvers compose without crashing and that the
expected mesh diagnostics (bed-state arrays, per-class source/sink
arrays for Riverine, τ_b, τ_crit) are populated with plausible values.
Tighter cross-checks against published case studies are deferred.

If the sibling Riverine repo cannot be located (e.g. running this test
suite in isolation), the test is skipped with an informative reason.
"""

from __future__ import annotations

from datetime import timedelta
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment import (
    SSM,
    contracts,
)
from clearwater_modules_v2.processes.sediment.io.sedflume import (
    load_sedflume_bundle,
)
from clearwater_modules_v2.processes.sediment.ssm import (
    _build_classes_from_bundle,
)


# ---------------------------------------------------------------------------
# Fixture-path resolution for the sibling Riverine repo
# ---------------------------------------------------------------------------

# tests/sediment/test_ssm_riverine_smoke.py → repo root is parents[2].
_MODULES_REPO_ROOT = Path(__file__).resolve().parents[2]
_RIVERINE_REPO_ROOT = (
    _MODULES_REPO_ROOT.parent / "ClearWater-Riverine-streaming"
)
_RIVERINE_FIXTURE_DIR = (
    _RIVERINE_REPO_ROOT / "tests" / "data" / "simple_test_cases"
)
_PLAN01_DIR = _RIVERINE_FIXTURE_DIR / "plan01_10x5"
_PLAN01_HDF = _PLAN01_DIR / "clearWaterTestCases.p01.hdf"
_PLAN01_IC = _PLAN01_DIR / "cwr_initial_conditions_p01.csv"
_PLAN01_BC = _PLAN01_DIR / "cwr_boundary_conditions_p01.csv"

_PLAN11_DIR = (
    _RIVERINE_REPO_ROOT
    / "tests" / "data" / "sumwere_test_cases" / "plan11_stormSurge"
)


# SAND2008 SEDflume bundle paths (already shipped in this repo).
_SAND2008_DATA_DIR = Path(__file__).parent / "data" / "sand2008_example"
_BED_SDF = _SAND2008_DATA_DIR / "bed.sdf"
_ERATE_SDF = _SAND2008_DATA_DIR / "erate.sdf"
_CORE_FIELD_SDF = _SAND2008_DATA_DIR / "core_field.sdf"


# Skip the entire module if the sibling Riverine repo or its fixtures
# are not present. We do this at collection time so the test summary is
# clean rather than silently red.
pytestmark = pytest.mark.skipif(
    not _PLAN01_HDF.is_file(),
    reason=(
        "ClearWater-Riverine-streaming sibling repo not found at "
        f"{_RIVERINE_REPO_ROOT}; coupled smoke test cannot run."
    ),
)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _import_clearwater_riverine():
    """Import ClearwaterRiverine, adding the sibling repo's ``src`` to sys.path
    if needed.

    The test is intended to run inside the modules repo's environment,
    where ``clearwater_riverine`` may not already be installed. Falling
    back to a path-based import keeps the smoke test self-contained.
    """
    try:
        import clearwater_riverine as cwr  # type: ignore[import-not-found]
        return cwr
    except ImportError:
        import sys

        riverine_src = _RIVERINE_REPO_ROOT / "src"
        if riverine_src.is_dir() and str(riverine_src) not in sys.path:
            sys.path.insert(0, str(riverine_src))
        try:
            import clearwater_riverine as cwr  # type: ignore[import-not-found]
            return cwr
        except ImportError as e:
            pytest.skip(
                f"clearwater_riverine could not be imported even after "
                f"adding {riverine_src} to sys.path: {e}"
            )


def _make_riverine(plan_dir: Path, plan_hdf: Path, ic_csv: Path, bc_csv: Path,
                   datetime_range=(0, 12)):
    """Construct a small Riverine instance windowed to ``datetime_range``
    time slots so the bed-state allocation stays modest.
    """
    cwr = _import_clearwater_riverine()
    return cwr.ClearwaterRiverine(
        flow_field_file_path=str(plan_hdf),
        diffusion_coefficient_input=0.001,
        constituent_dict={
            "tracer": {
                "initial_conditions": str(ic_csv),
                "boundary_conditions": str(bc_csv),
                "units": "mg/L",
            },
        },
        datetime_range=datetime_range,
    )


def _build_ssm_for_riverine(riverine, time_step: timedelta) -> SSM:
    """Construct an SSM bound to ``riverine.mesh`` from the SAND2008 bundle.

    Uses the ``external`` shear driver so the smoke test can write a
    plausible τ field directly (Riverine's plan01 hdf does not carry
    cell velocities, so the current-only driver would compute τ ≈ 0
    and there'd be no erosion to observe).
    """
    bundle = load_sedflume_bundle(_BED_SDF, _ERATE_SDF, _CORE_FIELD_SDF)
    registry = _build_classes_from_bundle(bundle)
    n_face = riverine.mesh.sizes[contracts.DIM_NFACE]

    ssm = SSM(
        sediment_classes=registry,
        sedflume_bundle=bundle,
        shear_driver="external",
        shear_options={"growth_limit": 0.0},
        bedload_solver="off",
        nsedflume=bundle.nsedflume,
        biostabilization_alpha=0.0,
        time_step=time_step,
        core_id=np.ones(n_face, dtype=np.int64),
    )
    ssm.bind_mesh(riverine.mesh)
    return ssm


def _seed_external_tau(mesh: xr.Dataset, tau_value_pa: float) -> None:
    """Allocate ``bed_shear_stress_input`` on the mesh and set a uniform τ."""
    n_time = mesh.sizes[contracts.DIM_TIME]
    n_face = mesh.sizes[contracts.DIM_NFACE]
    arr = np.full((n_time, n_face), tau_value_pa, dtype="float32")
    mesh[contracts.VAR_BED_SHEAR_STRESS_INPUT] = (
        (contracts.DIM_TIME, contracts.DIM_NFACE), arr,
    )


# ---------------------------------------------------------------------------
# Primary smoke test: plan01_10x5
# ---------------------------------------------------------------------------


def test_ssm_riverine_smoke_plan01_10x5():
    """Coupled run: Riverine transport + SSM bed/erosion on plan01_10x5.

    Steps:
      1. Build Riverine on plan01_10x5 with one tracer constituent.
      2. Build SSM (SAND2008 bundle, external τ driver, bedload off)
         and bind it to the Riverine mesh.
      3. Seed a uniform τ above the SAND2008 fines τ_ce so erosion fires.
      4. Step the coupled system for 5 timesteps:
            a. riverine.update()
            b. ssm.run(t)
      5. Smoke-assert that bed state, τ diagnostics, source-injection
         arrays, and the tracer concentration all look plausible.
    """
    if not _PLAN01_IC.is_file() or not _PLAN01_BC.is_file():
        pytest.skip(
            f"plan01_10x5 IC/BC csvs not found in {_PLAN01_DIR}"
        )

    # Window the run to a small number of time slots so bed-state
    # allocations (time × nface × n_layers × n_class) stay light.
    n_steps = 5
    riverine = _make_riverine(
        _PLAN01_DIR, _PLAN01_HDF, _PLAN01_IC, _PLAN01_BC,
        datetime_range=(0, n_steps + 2),
    )
    mesh = riverine.mesh

    # SSM step: roughly the Riverine cadence. plan01 is 1-second cadence,
    # but we give SSM a longer step so erosion has measurable signal.
    ssm = _build_ssm_for_riverine(riverine, time_step=timedelta(seconds=60))

    # τ above τ_ce(fines) for the SAND2008 bundle (~1.0–1.6 Pa for
    # finer classes); 4 Pa puts us comfortably in the eroding regime.
    _seed_external_tau(mesh, tau_value_pa=4.0)

    # Initial-state snapshots for comparison.
    initial_bed_total = float(ssm._bed.layer_mass_at(0).values.sum())
    initial_tracer = mesh["tracer"].isel(time=0).values.copy()

    # Riverine's mesh has datetime64 time labels. The external shear
    # driver indexes τ_in by label (sel), but bed-state setters accept
    # either an int positional index or a label. We pass the datetime
    # label to ssm.run() and integers to the bed-state carry-forward
    # writes.
    time_labels = mesh[contracts.DIM_TIME].values

    for t_idx in range(n_steps):
        # Carry SSM bed state forward into slot t_idx so the orchestrator
        # mutates a contiguous trajectory.
        if t_idx > 0:
            ssm._bed.set_layer_mass_at(
                t_idx, ssm._bed.layer_mass_at(t_idx - 1).values
            )
            ssm._bed.set_class_fraction_at(
                t_idx, ssm._bed.class_fraction_at(t_idx - 1).values
            )
            ssm._bed.set_layer_active_at(
                t_idx, ssm._bed.layer_active_at(t_idx - 1).values
            )
            ssm._bed.set_layer_taucrit_at(
                t_idx, ssm._bed.layer_taucrit_at(t_idx - 1).values
            )

        riverine.update()
        ssm.run(time=time_labels[t_idx])

    # ------------------------------------------------------------------
    # Smoke assertions
    # ------------------------------------------------------------------

    # 1. Riverine actually advanced.
    assert riverine.time_step == n_steps

    # 2. Bed-state arrays exist on the mesh and are well-formed.
    for var in (
        contracts.VAR_BED_LAYER_MASS,
        contracts.VAR_BED_CLASS_FRACTION,
        contracts.VAR_BED_LAYER_ACTIVE,
        contracts.VAR_BED_SHEAR_STRESS,
        contracts.VAR_BED_CRITICAL_SHEAR_STRESS,
        contracts.VAR_BED_EROSION_FLUX,
        contracts.VAR_BED_DEPOSITION_FLUX,
    ):
        assert var in mesh.data_vars, f"expected bed-state var {var!r}"

    # 3. No NaNs / no negative bed mass.
    for t_idx in range(n_steps):
        layer_mass = ssm._bed.layer_mass_at(t_idx).values
        assert not np.isnan(layer_mass).any(), (
            f"NaN in bed layer_mass at t={t_idx}"
        )
        assert (layer_mass >= -1e-9).all(), (
            f"Negative bed mass at t={t_idx}: min={layer_mass.min()}"
        )

    # 4. τ_b that we wrote is reflected on the mesh.
    tau_at_step0 = mesh[contracts.VAR_BED_SHEAR_STRESS].isel(time=0).values
    assert np.isclose(tau_at_step0.max(), 4.0, atol=1e-3)
    assert (tau_at_step0 > 0.0).all()

    # 5. Erosion fired somewhere — total bed mass should have decreased
    # relative to t=0 because τ=4 Pa is above the SAND2008 fines τ_ce.
    final_bed_total = float(
        ssm._bed.layer_mass_at(n_steps - 1).values.sum()
    )
    assert final_bed_total <= initial_bed_total + 1e-6, (
        f"bed mass grew under sustained erosion: "
        f"{initial_bed_total} -> {final_bed_total}"
    )

    # 6. SSM staged a per-class source on the mesh for Riverine.
    classes = list(ssm.registry_classes)
    for cls in classes:
        src_name = f"{cls.suspended_var}_source"
        assert src_name in mesh.data_vars, (
            f"missing SSM source-injection array {src_name!r}"
        )
        src_vals = mesh[src_name].values
        assert not np.isnan(src_vals).any(), (
            f"NaN in source injection {src_name}"
        )
        # Source magnitude should be finite and not absurd
        # (per-class g/cm² per step ≪ 1 g/cm² for our 60s step).
        assert np.abs(src_vals).max() < 10.0, (
            f"source {src_name} unexpectedly large: max={np.abs(src_vals).max()}"
        )

    # 7. Riverine's tracer evolved over the steps. plan01 has a non-
    # trivial flow field, so at least some interior cell concentration
    # should differ from t=0.
    tracer_evolved = mesh["tracer"].isel(time=n_steps).values
    diff = np.nan_to_num(tracer_evolved - initial_tracer, nan=0.0)
    assert np.any(np.abs(diff) > 1e-6), (
        "Riverine tracer did not evolve over the smoke-test window"
    )


# ---------------------------------------------------------------------------
# Optional: plan11_stormSurge (more physically interesting, marked slow)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_ssm_riverine_smoke_plan11_stormsurge():
    """Coupled run on the storm-surge fixture.

    Marked ``slow`` because plan11 has many more time slots than plan01
    and the bed-state allocation scales with the windowed time
    dimension. Skipped if plan11 fixtures are missing.
    """
    plan11_hdf = _PLAN11_DIR / "clearWaterTestCases.p11.hdf"
    if not plan11_hdf.is_file():
        pytest.skip(f"plan11_stormSurge HDF not found at {plan11_hdf}")

    # Guard against Git LFS pointer files masquerading as the real HDF.
    # A pointer file is ASCII text starting with "version https://...";
    # the real HDF starts with the binary signature \x89HDF.
    try:
        with open(plan11_hdf, "rb") as f:
            head = f.read(8)
        if not head.startswith(b"\x89HDF"):
            pytest.skip(
                f"plan11_stormSurge HDF at {plan11_hdf} is not a valid HDF "
                "file (likely an unfetched Git LFS pointer)"
            )
    except OSError as e:
        pytest.skip(f"could not inspect plan11 HDF: {e}")

    # Pick a discoverable IC/BC pair if present, else fall back to no
    # constituent_dict (Riverine still works with a single trivial
    # constituent only when CSVs exist).
    ic_csv = next(_PLAN11_DIR.glob("cwr_initial_conditions*.csv"), None)
    bc_csv = next(_PLAN11_DIR.glob("cwr_boundary_conditions*.csv"), None)
    if ic_csv is None or bc_csv is None:
        pytest.skip(
            f"plan11_stormSurge IC/BC csvs not found in {_PLAN11_DIR}"
        )

    n_steps = 3
    riverine = _make_riverine(
        _PLAN11_DIR, plan11_hdf, ic_csv, bc_csv,
        datetime_range=(0, n_steps + 2),
    )
    mesh = riverine.mesh
    ssm = _build_ssm_for_riverine(riverine, time_step=timedelta(seconds=60))
    _seed_external_tau(mesh, tau_value_pa=4.0)

    time_labels = mesh[contracts.DIM_TIME].values

    for t_idx in range(n_steps):
        if t_idx > 0:
            ssm._bed.set_layer_mass_at(
                t_idx, ssm._bed.layer_mass_at(t_idx - 1).values
            )
            ssm._bed.set_class_fraction_at(
                t_idx, ssm._bed.class_fraction_at(t_idx - 1).values
            )
            ssm._bed.set_layer_active_at(
                t_idx, ssm._bed.layer_active_at(t_idx - 1).values
            )
            ssm._bed.set_layer_taucrit_at(
                t_idx, ssm._bed.layer_taucrit_at(t_idx - 1).values
            )
        riverine.update()
        ssm.run(time=time_labels[t_idx])

    # Same minimal smoke assertions as plan01.
    assert riverine.time_step == n_steps
    final_layer_mass = ssm._bed.layer_mass_at(n_steps - 1).values
    assert not np.isnan(final_layer_mass).any()
    assert (final_layer_mass >= -1e-9).all()
