"""Tests for SSM SEDflume input parsers and ancillary loaders.

Reference dataset is the worked example in SAND2008-5621 Figures 1, 2,
and 3, reproduced under ``tests/sediment/data/sand2008_example/``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml

from clearwater_modules_v2.processes.sediment import contracts
from clearwater_modules_v2.processes.sediment.io.csv_loader import (
    load_yaml_config,
)
from clearwater_modules_v2.processes.sediment.io.hotstart import (
    read_hotstart,
    write_hotstart,
)
from clearwater_modules_v2.processes.sediment.io.mesh_mapping import (
    load_unstructured_core_map,
)
from clearwater_modules_v2.processes.sediment.io.sedflume import (
    DYNES_CM2_TO_PA,
    SedflumeBundle,
    load_sedflume_bundle,
    parse_bed_sdf,
    parse_core_field_sdf,
    parse_erate_sdf,
)


DATA_DIR = Path(__file__).parent / "data" / "sand2008_example"
BED = DATA_DIR / "bed.sdf"
ERATE = DATA_DIR / "erate.sdf"
CORE_FIELD = DATA_DIR / "core_field.sdf"


# ---------------------------------------------------------------------------
# bed.sdf
# ---------------------------------------------------------------------------


def test_parse_bed_sdf_header():
    """Header card values match SAND2008-5621 Figure 3."""
    bed = parse_bed_sdf(BED)
    assert bed["var_bed"] == 1
    assert bed["icalc_bl"] == 1
    assert bed["n_layers"] == 5            # KB
    assert bed["_itbm"] == 6
    assert bed["_nsicm"] == 8
    assert bed["_n_class"] == 9
    assert bed["max_deposit_limit"] == pytest.approx(1.0)
    assert bed["zb_skin_um"] == pytest.approx(1500.0)
    # TAUCONST is 10 dynes/cm² → 1 Pa.
    assert bed["tau_const_pa"] == pytest.approx(1.0)


def test_parse_bed_sdf_size_interpolants_match_manual():
    """The NSICM size-interpolant SCND values and TAUCRITE in Pa
    should match the Figure 3 example exactly."""
    bed = parse_bed_sdf(BED)
    np.testing.assert_array_equal(
        bed["size_interpolants_um"],
        [125.0, 222.0, 432.0, 1020.0, 2000.0, 2400.0, 3000.0, 6000.0],
    )
    np.testing.assert_allclose(
        bed["taucrit_per_size_pa"],
        [0.12, 0.227, 0.296, 0.417, 0.546, 0.588, 0.642, 0.848],
        rtol=1e-12,
    )


def test_parse_bed_sdf_d50_and_per_class_arrays():
    """Per-class D50 / TCRE / TCRSUS / settling shape and content."""
    bed = parse_bed_sdf(BED)
    np.testing.assert_array_equal(
        bed["d50_um"],
        [237.7, 427.0, 603.5, 853.5, 1070.5, 1570.5, 2415.0, 3415.0, 5450.1],
    )
    # TCRE row: 1.5 2.4 3.3 4.25 7.6 9.5 10.8 16. 24.8 dynes/cm² → /10 Pa
    np.testing.assert_allclose(
        bed["tau_ce_pa"],
        np.array([1.5, 2.4, 3.3, 4.25, 7.6, 9.5, 10.8, 16.0, 24.8]) * DYNES_CM2_TO_PA,
        rtol=1e-12,
    )
    # All settling speeds set to -1 → "compute with Cheng (1997)".
    np.testing.assert_array_equal(bed["settling_cm_s"], np.full(9, -1.0))


def test_parse_bed_sdf_erate_table_shape_and_values():
    """ENRATE table is NSICM × ITBM with the literal manual values."""
    bed = parse_bed_sdf(BED)
    table = bed["erate_active_table"]
    assert table.shape == (8, 6)
    np.testing.assert_allclose(
        table[0],
        [1e-9, 6.60e-5, 4.66e-4, 3.29e-3, 6.17e-3, 4.36e-2],
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        table[-1],
        [1e-9, 2.03e-8, 1.67e-7, 1.44e-6, 2.88e-6, 2.49e-5],
        rtol=1e-6,
    )


# ---------------------------------------------------------------------------
# erate.sdf
# ---------------------------------------------------------------------------


def test_parse_erate_sdf_two_cores_with_expected_shapes():
    erate = parse_erate_sdf(ERATE, n_layers=5, n_class=9, nsedflume=1)
    assert erate["n_cores"] == 2
    assert erate["layer_thickness_cm"].shape == (2, 5)
    assert erate["layer_taucrit_pa"].shape == (2, 5)
    assert erate["bulk_density_g_cm3"].shape == (2, 5)
    assert erate["particle_size_distribution_pct"].shape == (2, 5, 9)
    assert erate["erate_per_core_cm_s"].shape == (2, 5, 6)
    assert erate["water_density_g_cm3"] == pytest.approx(1.0)
    assert erate["solid_density_g_cm3"] == pytest.approx(2.6)


def test_parse_erate_sdf_per_layer_values_match_manual():
    erate = parse_erate_sdf(ERATE, n_layers=5, n_class=9, nsedflume=1)

    # Layer thicknesses (Figure 1): "0.00 0.00 15.00 15.00 15.00"
    np.testing.assert_array_equal(
        erate["layer_thickness_cm"][0],
        [0.0, 0.0, 15.0, 15.0, 15.0],
    )
    # τ_crit per layer = 4.25 dynes/cm² ⇒ 0.425 Pa
    np.testing.assert_allclose(
        erate["layer_taucrit_pa"][0],
        np.full(5, 4.25 * DYNES_CM2_TO_PA),
        rtol=1e-12,
    )
    # Bulk density 1.9 g/cm³
    np.testing.assert_array_equal(
        erate["bulk_density_g_cm3"][0],
        np.full(5, 1.9),
    )
    # PSD row sums (per layer per core) ≈ 100 % (mass percent format).
    # The SAND2008 manual example sums to 100.172 — a documented
    # rounding artifact in Figure 1, not a parser bug; SEDZLJ
    # renormalizes per s_sedic.f90:327-329 before use.
    psd = erate["particle_size_distribution_pct"]
    np.testing.assert_allclose(psd.sum(axis=-1), 100.0, atol=0.5)


def test_parse_erate_sdf_shear_levels_match_manual():
    """Shear levels per Figure 1 caption: 0, 2, 4, 8, 10, 20 Pa."""
    erate = parse_erate_sdf(ERATE, n_layers=5, n_class=9, nsedflume=1)
    np.testing.assert_array_equal(
        erate["tau_levels_pa"],
        [0.0, 2.0, 4.0, 8.0, 10.0, 20.0],
    )


def test_parse_erate_sdf_per_layer_erosion_rates():
    """A few representative entries from the erate.sdf example.

    For shear-level 4 (8 Pa, index 3), all layers report 7.0e-3 cm/s.
    For shear-level 6 (20 Pa, index 5), all layers report 6.49e-2 cm/s.
    """
    erate = parse_erate_sdf(ERATE, n_layers=5, n_class=9, nsedflume=1)
    rates = erate["erate_per_core_cm_s"]
    np.testing.assert_allclose(rates[0, :, 3], np.full(5, 7.0e-3), rtol=1e-12)
    np.testing.assert_allclose(rates[0, :, 5], np.full(5, 6.49e-2), rtol=1e-12)
    # Both cores hold the same data in this example.
    np.testing.assert_array_equal(rates[0], rates[1])


# ---------------------------------------------------------------------------
# core_field.sdf
# ---------------------------------------------------------------------------


def test_parse_core_field_sdf_snl_format():
    """SAND2008-5621 Figure 2: 9 rows × 15 cols, first 3 cols core 2,
    rest core 1."""
    cf = parse_core_field_sdf(CORE_FIELD)
    assert cf.shape == (9, 15)
    assert (cf[:, :3] == 2).all()
    assert (cf[:, 3:] == 1).all()


def test_parse_core_field_sdf_dsi_format(tmp_path: Path):
    """DSI standard: header lines including 'DSI' marker, then INCORE,
    then per-cell I J CORE triples."""
    p = tmp_path / "core_field.sdf"
    p.write_text(
        "C This file is in DSI standard format\n"
        "C ----------------------------\n"
        "C\n"
        "3\n"
        "C cells follow\n"
        "C I  J  CORE\n"
        "1 1 1\n"
        "2 1 2\n"
        "1 2 3\n"
        "2 2 1\n"
    )
    cf = parse_core_field_sdf(p)
    assert cf.shape == (2, 2)
    # File is 1-indexed; (i=1,j=1)→[0,0]=1, (i=2,j=1)→[0,1]=2,
    # (i=1,j=2)→[1,0]=3, (i=2,j=2)→[1,1]=1.
    assert cf[0, 0] == 1
    assert cf[0, 1] == 2
    assert cf[1, 0] == 3
    assert cf[1, 1] == 1


# ---------------------------------------------------------------------------
# Bundle assembly
# ---------------------------------------------------------------------------


def test_load_sedflume_bundle_round_trip():
    bundle = load_sedflume_bundle(BED, ERATE, CORE_FIELD)
    assert isinstance(bundle, SedflumeBundle)
    assert bundle.n_layers == 5
    assert bundle.n_cores == 2
    assert bundle.nsedflume == 1
    assert bundle.var_bed == 1
    assert bundle.icalc_bl == 1
    assert bundle.core_field_ij is not None
    assert bundle.core_field_ij.shape == (9, 15)
    np.testing.assert_array_equal(
        bundle.size_interpolants_um,
        [125.0, 222.0, 432.0, 1020.0, 2000.0, 2400.0, 3000.0, 6000.0],
    )
    np.testing.assert_allclose(
        bundle.taucrit_per_size_pa,
        [0.12, 0.227, 0.296, 0.417, 0.546, 0.588, 0.642, 0.848],
        rtol=1e-12,
    )


# ---------------------------------------------------------------------------
# YAML loader
# ---------------------------------------------------------------------------


def test_load_yaml_config_round_trip(tmp_path: Path):
    """Build a minimal YAML config (2 classes + 1 core) and verify it
    becomes a well-formed SedflumeBundle."""
    cfg = {
        "global": {
            "var_bed": 1,
            "icalc_bl": 1,
            "zb_skin_um": 1500.0,
            "tau_const_pa": 0.5,
            "bedload_cutoff_um": 64.0,
            "max_deposit_limit": 1.0,
        },
        "nsedflume": 1,
        "sediment_classes": [
            {"label": "silt_fine", "d50_um": 32.0, "tau_ce_pa": 0.15,
             "tau_cs_pa": 0.20, "settling_cm_s": -1},
            {"label": "sand_medium", "d50_um": 250.0, "tau_ce_pa": 0.20,
             "tau_cs_pa": 0.30, "settling_cm_s": 2.5},
        ],
        "bed_layers": {
            "n_layers": 4,
            "layer_thickness_cm": [0, 0, 5, 10],
            "bulk_density_g_cm3": [1.6, 1.6, 1.7, 1.8],
        },
        "size_interpolants": {
            "sizes_um": [50, 250],
            "taucrit_pa": [0.1, 0.3],
            "shear_levels_pa": [0, 1, 2],
            "erate_table_cm_s": [
                [0, 1e-5, 1e-4],
                [0, 1e-6, 1e-5],
            ],
        },
        "cores": [
            {
                "id": 1,
                "cell_indices": [0, 1, 2],
                "layer_taucrit_pa": [0.4, 0.4, 0.4, 0.4],
                "layer_thickness_cm": [0, 0, 5, 10],
                "bulk_density_g_cm3": [1.6, 1.6, 1.7, 1.8],
                "water_density_g_cm3": 1.0,
                "solid_density_g_cm3": 2.65,
                "particle_size_distribution_pct": [
                    [60, 40], [55, 45], [40, 60], [30, 70],
                ],
                "erate_cm_s": [
                    [0, 1e-5, 1e-4],
                    [0, 1e-5, 1e-4],
                    [0, 1e-5, 1e-4],
                    [0, 1e-5, 1e-4],
                ],
            }
        ],
    }
    p = tmp_path / "config.yml"
    with open(p, "w") as fh:
        yaml.safe_dump(cfg, fh)

    bundle = load_yaml_config(p)
    assert isinstance(bundle, SedflumeBundle)
    assert bundle.n_layers == 4
    assert bundle.n_cores == 1
    assert bundle.d50_um.tolist() == [32.0, 250.0]
    assert bundle.tau_ce_pa.tolist() == [0.15, 0.20]
    np.testing.assert_array_equal(
        bundle.layer_thickness_cm[0],
        [0, 0, 5, 10],
    )
    np.testing.assert_array_equal(
        bundle.bulk_density_g_cm3[0],
        [1.6, 1.6, 1.7, 1.8],
    )
    assert bundle.tau_levels_pa.tolist() == [0.0, 1.0, 2.0]
    assert bundle.erate_per_core_cm_s.shape == (1, 4, 3)


# ---------------------------------------------------------------------------
# Hotstart NetCDF round trip
# ---------------------------------------------------------------------------


def test_hotstart_netcdf_round_trip(tmp_path: Path):
    """Build a small mesh with a couple of bed-state vars, write/read
    NetCDF, and verify equality."""
    n_face = 4
    n_layer = 3
    n_class = 2
    rng = np.random.default_rng(0)
    mesh = xr.Dataset(
        {
            contracts.VAR_BED_LAYER_MASS: (
                (contracts.DIM_TIME, contracts.DIM_NFACE, contracts.DIM_LAYER),
                rng.random((1, n_face, n_layer), dtype=np.float32),
            ),
            contracts.VAR_BED_LAYER_BULK_DENSITY: (
                (contracts.DIM_NFACE, contracts.DIM_LAYER),
                np.full((n_face, n_layer), 1.7, dtype=np.float32),
            ),
            contracts.VAR_BED_CLASS_FRACTION: (
                (
                    contracts.DIM_TIME,
                    contracts.DIM_NFACE,
                    contracts.DIM_LAYER,
                    contracts.DIM_CLASS,
                ),
                rng.random((1, n_face, n_layer, n_class), dtype=np.float32),
            ),
            # A non-bed variable that should be filtered out.
            "irrelevant_field": (
                (contracts.DIM_NFACE,),
                np.zeros(n_face),
            ),
        },
        coords={contracts.DIM_TIME: [0]},
    )
    out = tmp_path / "ssm_state.nc"
    write_hotstart(mesh, out)

    with read_hotstart(out) as read_back:
        # Filtered: only bed-state vars persist.
        assert "irrelevant_field" not in read_back.data_vars
        for name in (
            contracts.VAR_BED_LAYER_MASS,
            contracts.VAR_BED_LAYER_BULK_DENSITY,
            contracts.VAR_BED_CLASS_FRACTION,
        ):
            assert name in read_back.data_vars
        xr.testing.assert_equal(
            read_back[contracts.VAR_BED_LAYER_MASS],
            mesh[contracts.VAR_BED_LAYER_MASS],
        )
        xr.testing.assert_equal(
            read_back[contracts.VAR_BED_CLASS_FRACTION],
            mesh[contracts.VAR_BED_CLASS_FRACTION],
        )


def test_hotstart_raises_when_no_bed_state(tmp_path: Path):
    mesh = xr.Dataset({"foo": (("nface",), np.zeros(3))})
    with pytest.raises(ValueError, match="bed-state"):
        write_hotstart(mesh, tmp_path / "x.nc")


# ---------------------------------------------------------------------------
# Mesh mapping CSV
# ---------------------------------------------------------------------------


def test_load_unstructured_core_map_basic(tmp_path: Path):
    p = tmp_path / "core_map.csv"
    p.write_text(
        "Cell_Index,Core_ID\n"
        "0,3\n"
        "2,5\n"
    )
    out = load_unstructured_core_map(p, n_face=5)
    np.testing.assert_array_equal(out, [3, 1, 5, 1, 1])


def test_load_unstructured_core_map_empty_defaults_all_to_one(tmp_path: Path):
    p = tmp_path / "core_map.csv"
    p.write_text("Cell_Index,Core_ID\n")
    out = load_unstructured_core_map(p, n_face=4)
    np.testing.assert_array_equal(out, [1, 1, 1, 1])
