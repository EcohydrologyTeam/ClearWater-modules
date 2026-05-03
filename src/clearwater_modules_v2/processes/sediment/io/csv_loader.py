"""YAML/CSV alternative input loader.

For new projects without legacy SEDZLJ inputs, accepts a single YAML
config that defines sediment classes, bed-layer structure, and
SEDflume-equivalent erosion data. Translates to the same
:class:`SedflumeBundle` shape as :mod:`io.sedflume` so downstream
code is loader-agnostic.

Format documented in design spec §9.3.

Expected YAML schema (illustrative, all keys optional unless marked):

.. code-block:: yaml

    sediment_classes:               # required, ordered list
      - label: silt_fine
        d50_um: 32                  # required
        tau_ce_pa: 0.15
        tau_cs_pa: 0.20
        settling_cm_s: -1           # -1 = compute via Cheng (1997)
      - label: sand_medium
        d50_um: 250
        tau_ce_pa: 0.20
        tau_cs_pa: 0.30
        settling_cm_s: 2.5

    bed_layers:
      n_layers: 8                    # required (KB)
      layer_thickness_cm: [0, 0, 5, 5, 5, 10, 20, 50]
      bulk_density_g_cm3: [1.6, 1.6, 1.6, 1.7, 1.7, 1.8, 1.9, 1.9]

    nsedflume: 1                     # 1 (table) or 2 (power-law)

    size_interpolants:
      sizes_um: [125, 222, 432, 1020, 2000, 2400, 3000, 6000]
      taucrit_pa: [0.12, 0.227, 0.296, 0.417, 0.546, 0.588, 0.642, 0.848]
      shear_levels_pa: [0, 2, 4, 8, 10, 20]
      erate_table_cm_s:              # NSICM rows × ITBM cols (nsedflume=1)
        - [1.0e-9, 6.6e-5, 4.66e-4, 3.29e-3, 6.17e-3, 4.36e-2]
        - ...
      # OR for nsedflume=2:
      actdep_a:    [...]
      actdep_n:    [...]
      actdep_max:  [...]

    cores:                           # required, ordered list
      - id: 1
        cell_indices: [0, 1, 2]
        layer_taucrit_pa: [...]      # n_layers
        layer_thickness_cm: [...]
        bulk_density_g_cm3: [...]
        water_density_g_cm3: 1.0
        solid_density_g_cm3: 2.65
        particle_size_distribution_pct:   # n_layers × n_class
          - [...]
          - ...
        # nsedflume == 1:
        erate_cm_s:                  # n_layers × ITBM
          - [...]
        # OR nsedflume == 2:
        ea: [...]
        en: [...]
        max_rate_cm_s: [...]

    global:                          # optional scalars
      var_bed: 1
      icalc_bl: 0
      zb_skin_um: 1500.0
      tau_const_pa: 0.0
      bedload_cutoff_um: 64.0
      max_deposit_limit: 1.0
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from .sedflume import SedflumeBundle


def _arr(values, dtype=float) -> np.ndarray:
    return np.asarray(values, dtype=dtype)


def load_yaml_config(path: Path | str) -> SedflumeBundle:
    """Parse a YAML sediment config and return a :class:`SedflumeBundle`."""
    with open(path, "r") as fh:
        cfg = yaml.safe_load(fh)

    # --- Sediment classes ---------------------------------------------
    classes = cfg["sediment_classes"]
    n_class = len(classes)
    d50_um = _arr([c["d50_um"] for c in classes])
    tau_ce_pa = _arr([c.get("tau_ce_pa", 0.0) for c in classes])
    tau_cs_pa = _arr([c.get("tau_cs_pa", 0.0) for c in classes])
    settling_cm_s = _arr([c.get("settling_cm_s", -1.0) for c in classes])

    # --- Bed layers ---------------------------------------------------
    bl = cfg["bed_layers"]
    n_layers = int(bl["n_layers"])

    # --- Global scalars -----------------------------------------------
    g = cfg.get("global", {})
    var_bed = int(g.get("var_bed", 1))
    icalc_bl = int(g.get("icalc_bl", 0))
    zb_skin_um = float(g.get("zb_skin_um", 0.0))
    tau_const_pa = float(g.get("tau_const_pa", 0.0))
    bedload_cutoff_um = float(g.get("bedload_cutoff_um", 64.0))
    max_deposit_limit = float(g.get("max_deposit_limit", 1.0))

    nsedflume = int(cfg.get("nsedflume", 1))

    # --- Size-interpolant block --------------------------------------
    si = cfg.get("size_interpolants", {})
    size_interpolants_um = _arr(si.get("sizes_um", []))
    taucrit_per_size_pa = _arr(si.get("taucrit_pa", []))
    shear_levels_pa = _arr(si.get("shear_levels_pa", []))

    erate_active_table: np.ndarray | None = None
    actdep_a: np.ndarray | None = None
    actdep_n: np.ndarray | None = None
    actdep_max: np.ndarray | None = None
    if nsedflume == 1:
        if "erate_table_cm_s" in si:
            erate_active_table = _arr(si["erate_table_cm_s"])
    else:
        actdep_a = _arr(si.get("actdep_a", []))
        actdep_n = _arr(si.get("actdep_n", []))
        actdep_max = _arr(si.get("actdep_max", []))

    # --- Cores --------------------------------------------------------
    cores = cfg["cores"]
    n_cores = len(cores)
    itbm = shear_levels_pa.size if shear_levels_pa.size > 0 else 0

    layer_taucrit = np.zeros((n_cores, n_layers), dtype=float)
    layer_thickness = np.zeros((n_cores, n_layers), dtype=float)
    bulk_density = np.zeros((n_cores, n_layers), dtype=float)
    psd = np.zeros((n_cores, n_layers, n_class), dtype=float)

    erate_per_core_cm_s: np.ndarray | None = None
    ea_per_core: np.ndarray | None = None
    en_per_core: np.ndarray | None = None
    max_rate_per_core_cm_s: np.ndarray | None = None

    if nsedflume == 1:
        erate_per_core_cm_s = np.zeros((n_cores, n_layers, itbm), dtype=float)
    else:
        ea_per_core = np.zeros((n_cores, n_layers), dtype=float)
        en_per_core = np.zeros((n_cores, n_layers), dtype=float)
        max_rate_per_core_cm_s = np.zeros((n_cores, n_layers), dtype=float)

    water_density = 1.0
    solid_density = 2.65
    for ic, core in enumerate(cores):
        layer_taucrit[ic, :] = _arr(core.get("layer_taucrit_pa", np.zeros(n_layers)))[:n_layers]
        # Fall back to the bed_layers default if the core omits these.
        thick = core.get("layer_thickness_cm", bl.get("layer_thickness_cm", np.zeros(n_layers)))
        layer_thickness[ic, :] = _arr(thick)[:n_layers]
        bulk = core.get("bulk_density_g_cm3", bl.get("bulk_density_g_cm3", np.zeros(n_layers)))
        bulk_density[ic, :] = _arr(bulk)[:n_layers]
        water_density = float(core.get("water_density_g_cm3", water_density))
        solid_density = float(core.get("solid_density_g_cm3", solid_density))

        psd_block = core.get("particle_size_distribution_pct")
        if psd_block is not None:
            arr = _arr(psd_block)
            psd[ic, : arr.shape[0], : arr.shape[1]] = arr

        if nsedflume == 1:
            tab = core.get("erate_cm_s")
            if tab is not None:
                arr = _arr(tab)
                erate_per_core_cm_s[ic, : arr.shape[0], : arr.shape[1]] = arr
        else:
            ea_per_core[ic, :] = _arr(core.get("ea", np.zeros(n_layers)))[:n_layers]
            en_per_core[ic, :] = _arr(core.get("en", np.zeros(n_layers)))[:n_layers]
            max_rate_per_core_cm_s[ic, :] = _arr(
                core.get("max_rate_cm_s", np.zeros(n_layers))
            )[:n_layers]

    return SedflumeBundle(
        n_layers=n_layers,
        var_bed=var_bed,
        icalc_bl=icalc_bl,
        nsedflume=nsedflume,
        zb_skin_um=zb_skin_um,
        tau_const_pa=tau_const_pa,
        bedload_cutoff_um=bedload_cutoff_um,
        max_deposit_limit=max_deposit_limit,
        d50_um=d50_um,
        tau_ce_pa=tau_ce_pa,
        tau_cs_pa=tau_cs_pa,
        settling_cm_s=settling_cm_s,
        size_interpolants_um=size_interpolants_um,
        taucrit_per_size_pa=taucrit_per_size_pa,
        erate_active_table=erate_active_table,
        actdep_a=actdep_a,
        actdep_n=actdep_n,
        actdep_max=actdep_max,
        n_cores=n_cores,
        layer_thickness_cm=layer_thickness,
        layer_taucrit_pa=layer_taucrit,
        bulk_density_g_cm3=bulk_density,
        water_density_g_cm3=water_density,
        solid_density_g_cm3=solid_density,
        particle_size_distribution_pct=psd,
        tau_levels_pa=shear_levels_pa,
        erate_per_core_cm_s=erate_per_core_cm_s,
        ea_per_core=ea_per_core,
        en_per_core=en_per_core,
        max_rate_per_core_cm_s=max_rate_per_core_cm_s,
        core_field_ij=None,
    )
