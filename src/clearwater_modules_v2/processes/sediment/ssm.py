"""SSM — top-level Process driver.

Orchestrates one sediment time step in the order specified by design
spec §7.4:

1. Read τ_b from the configured shear driver.
2. Compute D50_avg of surface layer; look up τ_crit (active or in-place).
3. Compute T_act; reorganize active layer.
4. Compute deposition rate per class (Gessler / Krone).
5. Compute erosion rate per layer per class (table or power-law).
6. Update layer mass and PERSED with mass conservation.
7. Update bed elevation and write bed_change diagnostics.
8. Step bedload (standalone or Riverine, per config).
9. Inject net (E − D) per class as Riverine source/sink for the next
   transport step.

Reference: design spec; existing v2 process pattern in
``temperature.py:15`` and ``nitrogen.py:19``.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import xarray as xr

from clearwater_data.variables import VariableRegistry

from ..base import Process, ProcessFactory
from . import (
    armoring,
    bed as bed_mod,
    bedload as bedload_mod,
    consolidation as consolidation_mod,
    contracts,
    coupling,
    deposition as deposition_mod,
    erosion as erosion_mod,
    settling as settling_mod,
    shear as shear_mod,
)
from .classes import SedimentClass, SedimentClassRegistry
from .io import csv_loader, mesh_mapping, sedflume as sedflume_io

if TYPE_CHECKING:  # pragma: no cover
    from ...model import Model


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_classes_from_bundle(bundle: sedflume_io.SedflumeBundle) -> SedimentClassRegistry:
    """Construct a :class:`SedimentClassRegistry` from a SEDflume bundle.

    Class labels are auto-generated as ``class_{i}`` since SEDflume
    inputs do not carry symbolic labels. Per-class τ_ce / τ_cs come
    from the bundle (already in Pa). Settling velocities are passed
    through; ``-1`` sentinels are converted to ``None`` so downstream
    code computes them via Cheng (1997).
    """
    classes: list[SedimentClass] = []
    for i, d50 in enumerate(bundle.d50_um):
        ws = float(bundle.settling_cm_s[i])
        settling = ws if ws > 0.0 else None
        classes.append(
            SedimentClass(
                label=f"class_{i}",
                d50_um=float(d50),
                tau_ce_pa=float(bundle.tau_ce_pa[i]),
                tau_cs_pa=float(bundle.tau_cs_pa[i]),
                settling_cm_s=settling,
                solid_density_g_cm3=float(bundle.solid_density_g_cm3),
            )
        )
    return SedimentClassRegistry.from_iterable(classes)


def _surface_layer_index(layer_active_row: np.ndarray) -> int:
    """Return the 0-origin index of the topmost non-absent layer in a single cell.

    Mirrors the SEDZLJ rule (s_sedzlj.f90:229-233): scan top-down from
    layer 1; first ACTIVE or IN_PLACE layer is the surface. Falls back
    to layer 0 if every layer is absent (degenerate; surface D50 will
    then come from a 0-mass blend ⇒ 0).
    """
    for k in range(layer_active_row.size):
        if layer_active_row[k] != bed_mod.LAYER_ABSENT:
            return k
    return 0


def _surface_class_fraction(
    class_fraction: np.ndarray,        # (nface, n_layer, n_class)
    layer_active: np.ndarray,          # (nface, n_layer) int8
) -> np.ndarray:
    """Per-cell surface-layer class fractions, shape (nface, n_class)."""
    n_face, _, n_class = class_fraction.shape
    out = np.zeros((n_face, n_class), dtype="float64")
    for f in range(n_face):
        k = _surface_layer_index(layer_active[f])
        out[f] = class_fraction[f, k]
    return out


def _wentworth_sand_mask(d50_um_array: np.ndarray) -> np.ndarray:
    """Wentworth-classification sand mask: D50 in [62.5, 2000] μm.

    Used to compute the bed-surface sand fraction ``F_s`` consumed by
    the Wilcock-Crowe (2003) closure. The bound 62.5 μm is the classical
    silt/sand boundary; 2000 μm is the sand/granule boundary.
    """
    d = np.asarray(d50_um_array, dtype="float64")
    return (d >= 62.5) & (d <= 2000.0)


def _surface_geometric_mean_um(
    surface_class_fraction: np.ndarray,   # (nface, n_class)
    d50_um_array: np.ndarray,             # (n_class,)
) -> np.ndarray:
    """Per-cell surface geometric-mean grain size :math:`d_{sg}` (μm).

    Computed via the standard mass-weighted log-mean
    :math:`\\log d_{sg} = \\sum_i F_i \\log d_i`, used by Wilcock & Crowe
    (2003) eq. 4 as the reference grain for the surface Shields stress.
    Falls back to the unweighted mean when a cell carries zero surface
    mass (degenerate; the closure will then apply its built-in default).
    """
    log_d = np.log(np.maximum(np.asarray(d50_um_array, dtype="float64"), 1.0e-12))
    fac = np.asarray(surface_class_fraction, dtype="float64")
    fac_sum = fac.sum(axis=-1, keepdims=True)
    safe_sum = np.where(fac_sum > 0.0, fac_sum, 1.0)
    fac_norm = np.where(fac_sum > 0.0, fac / safe_sum, 0.0)
    log_dsg = (fac_norm * log_d[None, :]).sum(axis=-1)
    out = np.exp(log_dsg)
    # Where the surface had zero mass, fall back to the array median.
    fallback = float(np.exp(log_d.mean()))
    return np.where(fac_sum.squeeze(-1) > 0.0, out, fallback)


def _surface_taucrit_pa(
    layer_taucrit: np.ndarray,        # (nface, n_layer) Pa
    layer_active: np.ndarray,         # (nface, n_layer) int8
) -> np.ndarray:
    """Per-cell surface-layer τ_crit (Pa), shape (nface,)."""
    n_face, _ = layer_taucrit.shape
    out = np.zeros(n_face, dtype="float64")
    for f in range(n_face):
        k = _surface_layer_index(layer_active[f])
        out[f] = layer_taucrit[f, k]
    return out


def _ensure_per_class_bed_var(mesh: xr.Dataset, name: str) -> None:
    """Allocate a per-class diagnostic flux var on the mesh if absent."""
    if name in mesh.data_vars:
        return
    n_face = mesh.sizes[contracts.DIM_NFACE]
    n_class = mesh.sizes[contracts.DIM_CLASS]
    if contracts.DIM_TIME in mesh.dims:
        n_time = mesh.sizes[contracts.DIM_TIME]
        mesh[name] = (
            (contracts.DIM_TIME, contracts.DIM_NFACE, contracts.DIM_CLASS),
            np.zeros((n_time, n_face, n_class), dtype="float32"),
        )
    else:
        mesh[name] = (
            (contracts.DIM_NFACE, contracts.DIM_CLASS),
            np.zeros((n_face, n_class), dtype="float32"),
        )


# ---------------------------------------------------------------------------
# Process
# ---------------------------------------------------------------------------


class SSM(Process):
    """Sediment Simulation Module — EFDC SEDZLJ port.

    Parameters
    ----------
    sediment_classes : list[SedimentClass] | SedimentClassRegistry
        Ordered sediment-size classes for this run.
    sedflume_bundle : SedflumeBundle
        Parsed SEDflume input data (from :func:`io.sedflume.load_sedflume_bundle`
        or :func:`io.csv_loader.load_yaml_config`).
    shear_driver : str
        ``"external" | "current_only" | "wave_current"``.
    shear_options : dict | None
        Extra kwargs forwarded to the shear-driver constructor.
    bedload_solver : str
        ``"standalone" | "riverine" | "off"``.
    bedload_transport_function : str
        Name of the bedload transport-rate closure (Stage-1 menu of seven
        formulas; see :data:`bedload.BEDLOAD_TRANSPORT_FUNCTIONS`).
        Default ``"van_rijn"`` preserves backwards compatibility with the
        original SEDZLJ-port behaviour.
    bed_streaming_interval_multiplier : int
        Bed-state Zarr flush every N transport flushes
        (default :data:`contracts.DEFAULT_BED_STREAMING_INTERVAL_MULTIPLIER`).
    nsedflume : int
        SEDflume erosion-rate formulation: 1 (table) or 2 (power-law).
    biostabilization_alpha : float
        ESM-vegetation feedback coefficient on τ_ce.
    time_step : timedelta
        SSM step (typically equal to the transport step).
    core_id : array-like or None
        Optional per-cell core ID vector (length ``nface``). If None,
        defaults to all-ones (single-core) at ``init_process`` time.
    """

    # State and forcing variables this process touches in the registry.
    # The exact set is computed at init_process time once the registry
    # exists and the sediment classes are known (one suspended-conc
    # variable per class), so this list is intentionally minimal here.
    variables = []

    def __init__(
        self,
        sediment_classes,
        sedflume_bundle,
        shear_driver: str = "current_only",
        shear_options: dict[str, Any] | None = None,
        bedload_solver: str = "standalone",
        bedload_transport_function: str = contracts.DEFAULT_BEDLOAD_TRANSPORT_FUNCTION,
        bed_streaming_interval_multiplier: int = contracts.DEFAULT_BED_STREAMING_INTERVAL_MULTIPLIER,
        nsedflume: int | None = None,
        biostabilization_alpha: float = 0.5,
        time_step: timedelta = timedelta(hours=1),
        core_id=None,
        consolidation_model: consolidation_mod.ConsolidationModel | None = None,
    ) -> None:
        Process.__init__(self, time_step)
        self.registry_classes = (
            sediment_classes
            if isinstance(sediment_classes, SedimentClassRegistry)
            else SedimentClassRegistry.from_iterable(sediment_classes)
        )
        self.sedflume_bundle = sedflume_bundle
        self.shear_driver_name = shear_driver
        self.shear_options = dict(shear_options) if shear_options else {}
        self.bedload_solver_name = bedload_solver
        # Validate the transport-function name eagerly so misconfigurations
        # surface at construction time, not deep in the run loop.
        bl_fn_name = str(bedload_transport_function).strip().lower()
        if bl_fn_name not in contracts.BEDLOAD_FUNCTIONS:
            raise ValueError(
                f"Unknown bedload_transport_function {bedload_transport_function!r}; "
                f"expected one of {contracts.BEDLOAD_FUNCTIONS}."
            )
        self.bedload_transport_function_name = bl_fn_name
        self._bedload_transport_function = None  # bound in _instantiate_drivers
        self.bed_streaming_interval_multiplier = bed_streaming_interval_multiplier
        # Default nsedflume from the bundle when caller did not override.
        self.nsedflume = (
            int(nsedflume) if nsedflume is not None
            else int(sedflume_bundle.nsedflume)
        )
        self.biostabilization_alpha = float(biostabilization_alpha)
        # Optional cohesive-bed consolidation model (Sanford-Maa 2001).
        # When None, behaviour matches SEDZLJ (constant per-layer τ_ce);
        # when provided, the per-layer τ_ce gate for cohesive classes is
        # aged via the model's effective_tau_ce(layer_age_s) call.
        self.consolidation_model: consolidation_mod.ConsolidationModel | None = (
            consolidation_model
        )
        self._user_core_id = (
            np.asarray(core_id, dtype=np.int64) if core_id is not None else None
        )

        # Bound at init_process / Stage 2 wire-up.
        self._shear_driver = None
        self._erosion_model = None
        self._bedload_solver = None
        self._bed: bed_mod.BedState | None = None
        self._riverine = None
        self._mesh: xr.Dataset | None = None
        self._core_id_da: xr.DataArray | None = None
        self._settling_velocity_cm_s: np.ndarray | None = None
        self._max_deposit_fraction: float = 1.0

    # -----------------------------------------------------------------
    # from_config
    # -----------------------------------------------------------------

    @ProcessFactory.register("sediment")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "SSM":
        """Build an SSM instance from a YAML-derived config dict.

        Required keys
        -------------
        ``input_format``
            ``"sedflume"`` | ``"yaml"``. Selects how the bed/erosion
            input is loaded.
        ``sedflume_yaml`` (when ``input_format == "yaml"``)
            Path to a single YAML file consumed by
            :func:`io.csv_loader.load_yaml_config`.
        ``bed_sdf``, ``erate_sdf`` (when ``input_format == "sedflume"``)
            Paths to the SEDflume ``bed.sdf`` and ``erate.sdf`` files.
            ``core_field_sdf`` is optional.

        Optional keys
        -------------
        ``shear_driver``, ``shear_options``, ``bedload_solver``,
        ``bedload_transport_function`` (one of
        :data:`contracts.BEDLOAD_FUNCTIONS`; default ``"van_rijn"``),
        ``nsedflume``, ``biostabilization_alpha``, ``time_step``
        (parsed with :func:`pandas.Timedelta`),
        ``bed_streaming_interval_multiplier``, ``core_map_csv``
        (Stage 3: stored for later resolution against the mesh).
        """
        input_format = str(config.get("input_format", "sedflume")).lower()

        if input_format == "yaml":
            yaml_path = Path(config["sedflume_yaml"])
            bundle = csv_loader.load_yaml_config(yaml_path)
        elif input_format == "sedflume":
            bundle = sedflume_io.load_sedflume_bundle(
                config["bed_sdf"],
                config["erate_sdf"],
                config.get("core_field_sdf"),
            )
        else:
            raise ValueError(
                f"Unknown input_format {input_format!r}; expected 'sedflume' or 'yaml'."
            )

        registry = _build_classes_from_bundle(bundle)

        # time_step: accept timedelta or a parseable string.
        ts_raw = config.get("time_step")
        if ts_raw is None:
            time_step = timedelta(hours=1)
        elif isinstance(ts_raw, timedelta):
            time_step = ts_raw
        else:
            time_step = pd.Timedelta(ts_raw).to_pytimedelta()

        # core_map_csv resolved against nface in init_process; stash now.
        core_map_csv = config.get("core_map_csv")

        # ------------------------------------------------------------------
        # Optional consolidation model. Schema:
        #   consolidation:
        #     enabled: true
        #     model: sanford_maa
        #     tau_ce_zero_pa: 0.10
        #     tau_ce_inf_pa:  0.50
        #     consolidation_time_s: 604800
        # If absent or ``enabled: false``, no consolidation is applied
        # (SEDZLJ-equivalent behaviour).
        # ------------------------------------------------------------------
        consolidation_model = None
        consol_cfg = config.get("consolidation")
        if consol_cfg and bool(consol_cfg.get("enabled", False)):
            model_name = str(consol_cfg.get("model", "sanford_maa")).lower()
            if model_name == "sanford_maa":
                consolidation_model = consolidation_mod.SanfordMaaConsolidation(
                    tau_ce_zero_pa=float(
                        consol_cfg.get(
                            "tau_ce_zero_pa",
                            contracts.DEFAULT_CONSOLIDATION_TAU_CE_ZERO_PA,
                        )
                    ),
                    tau_ce_inf_pa=float(
                        consol_cfg.get(
                            "tau_ce_inf_pa",
                            contracts.DEFAULT_CONSOLIDATION_TAU_CE_INF_PA,
                        )
                    ),
                    consolidation_time_s=float(
                        consol_cfg.get(
                            "consolidation_time_s",
                            contracts.DEFAULT_CONSOLIDATION_TIME_S,
                        )
                    ),
                )
            else:
                raise ValueError(
                    f"Unknown consolidation.model {model_name!r}; "
                    "expected 'sanford_maa'."
                )

        ssm = SSM(
            sediment_classes=registry,
            sedflume_bundle=bundle,
            shear_driver=str(config.get("shear_driver", "current_only")),
            shear_options=config.get("shear_options"),
            bedload_solver=str(config.get("bedload_solver", "standalone")),
            bedload_transport_function=str(
                config.get(
                    "bedload_transport_function",
                    contracts.DEFAULT_BEDLOAD_TRANSPORT_FUNCTION,
                )
            ),
            bed_streaming_interval_multiplier=int(
                config.get(
                    "bed_streaming_interval_multiplier",
                    contracts.DEFAULT_BED_STREAMING_INTERVAL_MULTIPLIER,
                )
            ),
            nsedflume=config.get("nsedflume"),
            biostabilization_alpha=float(config.get("biostabilization_alpha", 0.5)),
            time_step=time_step,
            consolidation_model=consolidation_model,
        )
        ssm._core_map_csv_path = core_map_csv  # consumed in init_process
        return ssm

    # -----------------------------------------------------------------
    # init_process
    # -----------------------------------------------------------------

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """Bind to peer processes, allocate bed state, instantiate drivers.

        Tolerant of test stubs: ``model`` may be ``None`` or lack
        ``has_process`` / ``get_process`` (in which case Riverine is
        skipped and the caller must populate ``self._mesh`` and
        ``self._bed`` directly via :meth:`bind_mesh`).
        """
        # Try to find Riverine. Defensive: tests may pass a lightweight
        # stub object that does not implement has_process.
        riverine = None
        has_proc = getattr(model, "has_process", None) if model is not None else None
        if callable(has_proc) and has_proc("Riverine"):
            riverine_proc = model.get_process("Riverine")
            riverine = getattr(riverine_proc, "riverine_instance", None)
            self._riverine = riverine

        # Mesh acquisition: prefer Riverine's mesh, else require the
        # caller to inject one via bind_mesh() before run().
        if riverine is not None and hasattr(riverine, "mesh"):
            self._mesh = riverine.mesh
            self._register_suspended_constituents(riverine)
            self._allocate_bed_state(self._mesh)

        # Always wire up the (Riverine-independent) drivers.
        self._instantiate_drivers()

    def bind_mesh(self, mesh: xr.Dataset) -> None:
        """Test-only / standalone entry point: bind directly to a mesh.

        When SSM runs without a full ClearWater Model and Riverine
        instance (e.g. unit / integration tests), call this after
        ``__init__`` and before ``run`` to allocate bed state on a
        user-provided mesh dataset.
        """
        self._mesh = mesh
        self._allocate_bed_state(mesh)
        if self._shear_driver is None:
            self._instantiate_drivers()

    # -- helpers used by init_process ---------------------------------

    def _register_suspended_constituents(self, riverine) -> None:
        """Register one suspended-class constituent per class with Riverine.

        Stage 3 will set initial / boundary conditions; for now we
        register each class as a zero-decay scalar with no IC/BC
        sources. The constituent_dict format mirrors
        ``transport.py:188-201``.
        """
        if not hasattr(riverine, "constituent_dict") or riverine.constituent_dict is None:
            riverine.constituent_dict = {}

        for cls in self.registry_classes:
            cname = cls.suspended_var
            if cname in riverine.constituent_dict:
                continue
            riverine.constituent_dict[cname] = {
                "initial_conditions": None,
                "boundary_conditions": None,
                "units": "mg/L",
                "decay_rate": 0.0,
            }

    def _allocate_bed_state(self, mesh: xr.Dataset) -> None:
        """Initialize bed-state arrays on the mesh from the SEDflume bundle."""
        bundle = self.sedflume_bundle
        n_face = mesh.sizes[contracts.DIM_NFACE]
        n_layers = bundle.n_layers
        n_class = len(self.registry_classes)

        # Resolve per-cell core_id: explicit kwarg → CSV → all-ones.
        if self._user_core_id is not None:
            core_id_1based = self._user_core_id.copy()
        else:
            csv_path = getattr(self, "_core_map_csv_path", None)
            if csv_path:
                core_id_1based = mesh_mapping.load_unstructured_core_map(
                    csv_path, n_face=n_face
                )
            else:
                core_id_1based = np.ones(n_face, dtype=np.int64)

        if core_id_1based.shape != (n_face,):
            raise ValueError(
                f"core_id length {core_id_1based.shape} != n_face={n_face}"
            )
        # Convert to 0-origin core indices for downstream array gathers.
        core_id_0based = (core_id_1based - 1).astype(np.int64)
        if core_id_0based.min() < 0 or core_id_0based.max() >= bundle.n_cores:
            raise ValueError(
                f"core_id values out of range [1, {bundle.n_cores}]; "
                f"got [{core_id_1based.min()}, {core_id_1based.max()}]"
            )
        self._core_id_da = xr.DataArray(
            core_id_0based, dims=(contracts.DIM_NFACE,), name="ssm_core_id"
        )

        # ---- Initial-condition arrays per cell ------------------------
        # Layer mass (g/cm²) = thickness_cm × bulk_density_g_cm3
        thick = bundle.layer_thickness_cm[core_id_0based]              # (nface, n_layers)
        bulkd = bundle.bulk_density_g_cm3[core_id_0based]              # (nface, n_layers)
        layer_mass = (thick * bulkd).astype("float64")

        # PSD pct → fractions, normalized so layer-sum = 1.
        psd_pct = bundle.particle_size_distribution_pct[core_id_0based]  # (nface, n_layers, n_class)
        psd = psd_pct.astype("float64") / 100.0
        # Per-layer per-cell renormalize. Where a layer has zero mass or
        # all-zero PSD, leave fractions at zero (no division by zero).
        layer_sum = psd.sum(axis=-1, keepdims=True)
        safe_sum = np.where(layer_sum > 0.0, layer_sum, 1.0)
        psd = np.where(layer_sum > 0.0, psd / safe_sum, 0.0)

        # Layer active state: in-place core where mass > 0, else absent.
        # Active and deposition layers (indices 0 and 1) start absent
        # unless the bundle reports nonzero thickness for them
        # (legacy sedflume always sets layers 0/1 thickness to 0).
        layer_active = np.full(
            (n_face, n_layers), bed_mod.LAYER_ABSENT, dtype="int8"
        )
        layer_active[layer_mass > 0.0] = bed_mod.LAYER_IN_PLACE

        # Per-layer τ_crit from bundle (Pa).
        taucor = bundle.layer_taucrit_pa[core_id_0based].astype("float64")

        bed = bed_mod.initialize_bed_state(
            mesh=mesh,
            registry=self.registry_classes,
            n_layers=n_layers,
            initial_layer_mass=layer_mass,
            initial_class_fraction=psd,
            bulk_density=bulkd.astype("float64"),
            initial_layer_active=layer_active,
            taucor_initial=taucor,
        )
        self._bed = bed

        # Pre-allocate per-class diagnostic flux vars in case the bedload
        # solver or the orchestration code references them later.
        _ensure_per_class_bed_var(mesh, contracts.VAR_BED_EROSION_FLUX)
        _ensure_per_class_bed_var(mesh, contracts.VAR_BED_DEPOSITION_FLUX)

        # max-deposit-fraction from bundle (Stage 1 default = 1.0).
        if bundle.max_deposit_limit and bundle.max_deposit_limit > 0.0:
            self._max_deposit_fraction = float(bundle.max_deposit_limit)

    def _instantiate_drivers(self) -> None:
        """Create the configured shear, erosion, and bedload models."""
        # Shear driver --------------------------------------------------
        sd_name = self.shear_driver_name.lower()
        if sd_name == "external":
            self._shear_driver = shear_mod.ExternalShearDriver(**self.shear_options)
        elif sd_name == "current_only":
            self._shear_driver = shear_mod.CurrentOnlyShearDriver(**self.shear_options)
        elif sd_name == "wave_current":
            self._shear_driver = shear_mod.WaveCurrentShearDriver(**self.shear_options)
        else:
            raise ValueError(
                f"Unknown shear_driver {self.shear_driver_name!r}; "
                "expected 'external', 'current_only', or 'wave_current'."
            )

        # Erosion model -------------------------------------------------
        bundle = self.sedflume_bundle
        if self.nsedflume == 1:
            if bundle.erate_per_core_cm_s is None:
                raise ValueError(
                    "nsedflume=1 requires erate_per_core_cm_s in the bundle."
                )
            erate_active = bundle.erate_active_table
            if erate_active is None:
                # Synthesize an active-layer table by repeating the deepest
                # in-place rate per shear bracket. Keeps the constructor
                # happy when the bundle did not include a separate
                # active-layer table.
                erate_active = bundle.erate_per_core_cm_s[0, -1, :].reshape(1, -1)
                erate_active = np.tile(
                    erate_active, (bundle.size_interpolants_um.size, 1)
                )
            self._erosion_model = erosion_mod.SedflumeTableErosionModel(
                tau_levels_pa=bundle.tau_levels_pa,
                erate_per_core=bundle.erate_per_core_cm_s,
                erate_active_per_size=erate_active,
                size_interpolants_um=bundle.size_interpolants_um,
                taucrit_per_size_pa=bundle.taucrit_per_size_pa,
                consolidation_model=self.consolidation_model,
            )
        elif self.nsedflume == 2:
            if (
                bundle.ea_per_core is None
                or bundle.en_per_core is None
                or bundle.max_rate_per_core_cm_s is None
            ):
                raise ValueError(
                    "nsedflume=2 requires (ea_per_core, en_per_core, max_rate_per_core_cm_s)."
                )
            actdep_a = (
                bundle.actdep_a if bundle.actdep_a is not None
                else np.full(bundle.size_interpolants_um.size, 0.0)
            )
            actdep_n = (
                bundle.actdep_n if bundle.actdep_n is not None
                else np.zeros_like(actdep_a)
            )
            actdep_max = (
                bundle.actdep_max if bundle.actdep_max is not None
                else np.full_like(actdep_a, 1.0)
            )
            self._erosion_model = erosion_mod.PowerLawErosionModel(
                ea_per_core=bundle.ea_per_core,
                en_per_core=bundle.en_per_core,
                max_rate_per_core=bundle.max_rate_per_core_cm_s,
                actdep_a=actdep_a,
                actdep_n=actdep_n,
                actdep_max=actdep_max,
                consolidation_model=self.consolidation_model,
            )
        else:
            raise ValueError(
                f"nsedflume must be 1 or 2; got {self.nsedflume}"
            )

        # Bedload transport-rate closure (Stage-1 menu of seven formulas).
        # The standalone / Riverine solvers default to van Rijn for the
        # advection velocity; the closure object is exposed on the SSM
        # instance for downstream code that wants to compute q_b directly
        # (e.g. mass-budget reporting, future implicit bed-evolution).
        self._bedload_transport_function = bedload_mod.get_transport_function(
            self.bedload_transport_function_name
        )

        # Bedload solver ------------------------------------------------
        bl_name = self.bedload_solver_name.lower()
        if bl_name == "off":
            self._bedload_solver = None
        elif bl_name == "standalone":
            self._bedload_solver = bedload_mod.BedloadStandaloneExplicit(
                self.registry_classes,
                bedload_cutoff_um=(
                    self.sedflume_bundle.bedload_cutoff_um
                    if self.sedflume_bundle.bedload_cutoff_um > 0.0
                    else contracts.DEFAULT_BEDLOAD_CUTOFF_UM
                ),
                transport_function=self._bedload_transport_function,
            )
        elif bl_name == "riverine":
            if self._riverine is None:
                raise ValueError(
                    "bedload_solver='riverine' requires a Riverine instance; "
                    "none was bound during init_process."
                )
            self._bedload_solver = bedload_mod.BedloadRiverineConstituent(
                self.registry_classes,
                self._riverine,
                bedload_cutoff_um=(
                    self.sedflume_bundle.bedload_cutoff_um
                    if self.sedflume_bundle.bedload_cutoff_um > 0.0
                    else contracts.DEFAULT_BEDLOAD_CUTOFF_UM
                ),
                transport_function=self._bedload_transport_function,
            )
        else:
            raise ValueError(
                f"Unknown bedload_solver {self.bedload_solver_name!r}; "
                "expected 'standalone', 'riverine', or 'off'."
            )

        # Per-class settling velocities ---------------------------------
        self._settling_velocity_cm_s = settling_mod.resolve_settling_velocities(
            self.registry_classes
        )

    # -----------------------------------------------------------------
    # run
    # -----------------------------------------------------------------

    def run(
        self,
        time,
        registry: VariableRegistry | None = None,
        *,
        tau_pa_override: xr.DataArray | None = None,
        suspended_concentration: xr.DataArray | None = None,
        bottom_water_layer_depth_m: xr.DataArray | None = None,
    ) -> None:
        """Orchestrate one sediment time step.

        Steps follow design spec §7.4 (paraphrased here, see top-of-file
        docstring). The orchestration writes results back to the mesh
        in place; if a Riverine instance is bound, per-class net
        (E − D) is staged on ``mesh[f"{cls.suspended_var}_source"]``
        for the next transport step.

        Test-mode kwargs (``tau_pa_override``, ``suspended_concentration``,
        ``bottom_water_layer_depth_m``) let unit tests bypass the
        Riverine read-back and feed deterministic forcing in.
        """
        if self._bed is None or self._mesh is None:
            raise RuntimeError(
                "SSM.run() called before init_process / bind_mesh established "
                "the bed state. Call bind_mesh(mesh) or include the SSM in a "
                "Model with a Riverine peer."
            )

        bed = self._bed
        mesh = self._mesh
        n_face = mesh.sizes[contracts.DIM_NFACE]
        n_layers = bed.n_layers
        n_classes = bed.n_classes
        dt = float(self.time_step.total_seconds())

        # --- Pull state at the requested time slot -----------------------
        layer_mass = np.asarray(
            bed.layer_mass_at(time).values, dtype="float64"
        ).copy()                                                # (nface, n_layers)
        class_fraction = np.asarray(
            bed.class_fraction_at(time).values, dtype="float64"
        ).copy()                                                # (nface, n_layers, n_class)
        layer_active = np.asarray(
            bed.layer_active_at(time).values, dtype="int8"
        ).copy()                                                # (nface, n_layers)
        layer_taucrit = np.asarray(
            bed.layer_taucrit_at(time).values, dtype="float64"
        )                                                       # (nface, n_layers)
        layer_initial_mass = np.asarray(
            bed.layer_initial_mass.values, dtype="float64"
        )                                                       # (nface, n_layers)
        bulk_density = np.asarray(
            bed.layer_bulk_density.values, dtype="float64"
        )                                                       # (nface, n_layers)

        # --- (2) Surface D50_avg and τ_crit lookup ----------------------
        surface_pers = _surface_class_fraction(class_fraction, layer_active)
        surface_pers_da = xr.DataArray(
            surface_pers,
            dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
        )
        d50_um_array = self.registry_classes.d50_um_array
        d50_surface = armoring.compute_d50_avg(surface_pers_da, d50_um_array)

        tau_crit_surface = armoring.interpolate_taucrit_from_d50(
            d50_surface,
            self.sedflume_bundle.size_interpolants_um,
            self.sedflume_bundle.taucrit_per_size_pa,
        )

        # --- (1) Compute τ via configured driver ------------------------
        if tau_pa_override is not None:
            tau_pa = tau_pa_override
            if not isinstance(tau_pa, xr.DataArray):
                tau_pa = xr.DataArray(
                    np.asarray(tau_pa, dtype="float64"),
                    dims=(contracts.DIM_NFACE,),
                )
        else:
            # Read previous τ for the growth limiter; 0 at t=0.
            prev_tau_da = self._read_previous_tau(mesh, time)
            tau_pa = self._shear_driver.compute(
                mesh=mesh,
                time=time,
                d50_surface_um=d50_surface,
                previous_tau_pa=prev_tau_da,
            )

        # --- (4) Vegetation cohesion: per-class effective τ_ce ----------
        bio, root, _ = coupling.read_vegetation_feedback(mesh, time)
        # Build per-cell, per-class baseline τ_ce array.
        tau_ce_per_class = np.array(
            [c.tau_ce_pa if c.tau_ce_pa is not None else 0.0 for c in self.registry_classes],
            dtype="float64",
        )
        tau_ce_da = xr.DataArray(
            np.broadcast_to(tau_ce_per_class, (n_face, n_classes)).copy(),
            dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
        )
        tau_ce_eff_da = erosion_mod.apply_vegetation_cohesion(
            tau_ce_da,
            biostabilization=bio,
            root_cohesion_pa=root,
            biostabilization_alpha=self.biostabilization_alpha,
        )
        tau_ce_eff = np.asarray(
            tau_ce_eff_da.values, dtype="float64"
        )                                                  # (nface, n_class)

        # --- (3 / 5) Active-layer reorganization ------------------------
        # Use per-cell layer-1 bulk density.
        bulkdens_layer1_da = xr.DataArray(
            bulk_density[:, 0], dims=(contracts.DIM_NFACE,)
        )
        bed_mod.reorganize_active_layer(
            bed=bed, t=time,
            tau_pa=tau_pa,
            tau_crit_pa=tau_crit_surface,
            d50_surface_um=d50_surface,
            bulk_density_layer1=bulkdens_layer1_da,
            tactm=contracts.DEFAULT_TACTM,
        )

        # Refresh local caches after reorganization mutated mesh in place.
        layer_mass = np.asarray(
            bed.layer_mass_at(time).values, dtype="float64"
        ).copy()
        class_fraction = np.asarray(
            bed.class_fraction_at(time).values, dtype="float64"
        ).copy()
        layer_active = np.asarray(
            bed.layer_active_at(time).values, dtype="int8"
        ).copy()

        # --- (6) Deposition flux per class -----------------------------
        if suspended_concentration is None:
            suspended_concentration = self._read_suspended_concentration(mesh, time)
        if bottom_water_layer_depth_m is None:
            bottom_water_layer_depth_m = self._read_bottom_depth(mesh, time)

        deposition_flux = deposition_mod.compute_deposition_flux(
            registry=self.registry_classes,
            suspended_concentration=suspended_concentration,
            tau_pa=tau_pa,
            settling_velocity_cm_s=self._settling_velocity_cm_s,
            bottom_water_layer_depth_m=bottom_water_layer_depth_m,
            dt_seconds=dt,
            max_deposit_fraction=self._max_deposit_fraction,
        )                                                  # (nface, n_class) g/cm²

        # --- (7) Per-layer per-class erosion ---------------------------
        tau_arr = np.asarray(tau_pa.values, dtype="float64")  # (nface,)
        erosion_per_class = np.zeros((n_face, n_classes), dtype="float64")

        # Per-class layer mass remaining (PERSED × layer mass).
        # Updates accumulate as we erode layers top-down.
        layer_class_mass = (
            class_fraction * layer_mass[..., None]
        )                                                  # (nface, n_layers, n_class)

        # Per-(layer, class) effective τ_ce. Default: broadcast the
        # per-class vegetation-aged value across layers (no consolidation).
        # When a consolidation model is configured, age-adjust the cohesive
        # classes per-layer using the current age field.
        tau_ce_layer_class = np.broadcast_to(
            tau_ce_eff[:, None, :], (n_face, n_layers, n_classes)
        ).copy()                                           # (nface, n_layers, n_class)
        if self.consolidation_model is not None:
            layer_age = bed.layer_age_at(time)             # (nface, n_layers) DataArray
            cohesive_mask = np.array(
                [c.is_cohesive for c in self.registry_classes], dtype=bool
            )
            tau_ce_layer_class_da = xr.DataArray(
                tau_ce_layer_class,
                dims=(contracts.DIM_NFACE, contracts.DIM_LAYER, contracts.DIM_CLASS),
            )
            aged_da = erosion_mod.apply_consolidation(
                tau_ce_layer_class_da,
                layer_age,
                cohesive_mask,
                consolidation_model=self.consolidation_model,
            )
            tau_ce_layer_class = np.asarray(aged_da.values, dtype="float64")

        for k in range(n_layers):
            mass_k = layer_mass[:, k]
            if not np.any(mass_k > 0.0):
                continue

            # Bulk erosion rate (g/cm²/s) for this layer.
            mass_k_da = xr.DataArray(mass_k, dims=(contracts.DIM_NFACE,))
            mass_k0_da = xr.DataArray(
                layer_initial_mass[:, k], dims=(contracts.DIM_NFACE,)
            )
            bulkd_k_da = xr.DataArray(
                bulk_density[:, k], dims=(contracts.DIM_NFACE,)
            )
            erate = self._erosion_model.erosion_rate(
                tau_pa=tau_pa,
                layer_index=k + 1,                     # 1-origin inside the model
                layer_mass=mass_k_da,
                layer_initial_mass=mass_k0_da,
                bulk_density=bulkd_k_da,
                core_id=self._core_id_da,
            )                                          # (nface,)
            erate_arr = np.asarray(erate.values, dtype="float64")
            # Bulk mass eroded this step (g/cm²).
            bulk_erode = erate_arr * dt                # (nface,)

            # Fractionate per class via PERSED, gated by
            # τ ≥ τ_ce_eff(layer k, class). The gate is per-layer when
            # consolidation is in use (so freshly-deposited cohesive
            # mass on layer 1 may erode at a lower threshold than the
            # in-place core layers below).
            persed_k = class_fraction[:, k, :]          # (nface, n_class)
            gate_k = (
                tau_arr[:, None] >= tau_ce_layer_class[:, k, :]
            )                                           # (nface, n_class)
            class_eroded = (
                bulk_erode[:, None] * persed_k * gate_k.astype("float64")
            )                                           # (nface, n_class)

            # Cap per-class erosion at the available per-class layer mass
            # (s_sedzlj.f90:584-617 algorithm).
            available = layer_class_mass[:, k, :]       # (nface, n_class)
            class_eroded = np.minimum(class_eroded, available)
            class_eroded = np.maximum(class_eroded, 0.0)

            # Accumulate.
            erosion_per_class += class_eroded
            layer_class_mass[:, k, :] -= class_eroded

        # Re-derive per-layer total mass and PERSED from the per-class
        # remainder, preserving the bed.py mass-conservation invariant.
        new_layer_mass = layer_class_mass.sum(axis=-1)     # (nface, n_layers)
        # PERSED = remaining_class_mass / new_layer_mass (where layer is non-empty).
        safe_lm = np.where(new_layer_mass > 0.0, new_layer_mass, 1.0)
        new_persed = np.where(
            new_layer_mass[..., None] > 0.0,
            layer_class_mass / safe_lm[..., None],
            0.0,
        )
        # Layer-active flag: if a layer is now drained, mark absent.
        new_active = layer_active.copy()
        new_active[new_layer_mass <= 0.0] = bed_mod.LAYER_ABSENT

        # --- Add deposition mass per class to layer 1 ------------------
        deposition_arr = np.asarray(
            deposition_flux.values, dtype="float64"
        )                                                  # (nface, n_class)
        # Sum per-class deposition into layer 1 mass; redistribute fractions.
        deposit_total = deposition_arr.sum(axis=-1)        # (nface,)
        old_layer1_mass = new_layer_mass[:, 0]
        old_layer1_class_mass = layer_class_mass[:, 0, :]
        new_layer1_class_mass = old_layer1_class_mass + deposition_arr
        new_layer1_mass = old_layer1_mass + deposit_total
        # Recompute layer-1 PERSED.
        safe_l1 = np.where(new_layer1_mass > 0.0, new_layer1_mass, 1.0)
        new_layer1_persed = np.where(
            new_layer1_mass[..., None] > 0.0,
            new_layer1_class_mass / safe_l1[..., None],
            0.0,
        )

        # Commit back into the layer-major arrays.
        new_layer_mass[:, 0] = new_layer1_mass
        new_persed[:, 0, :] = new_layer1_persed
        # Layer 1 becomes ACTIVE if it received deposition.
        new_active[:, 0] = np.where(
            new_layer1_mass > 0.0, bed_mod.LAYER_ACTIVE, new_active[:, 0]
        )

        # --- (8) Write back to mesh ------------------------------------
        bed.set_layer_mass_at(time, new_layer_mass)
        bed.set_class_fraction_at(time, new_persed)
        bed.set_layer_active_at(time, new_active)

        # Age dilution on deposition: every g/cm² of fresh mass added to
        # layer 1 enters with age 0, so the mass-weighted layer-1 age is
        # diluted accordingly. We pass (old_layer1_mass, deposit_total)
        # so the helper can compute the new mass-weighted age in place.
        # No-op when deposit_total is zero everywhere.
        if np.any(deposit_total > 0.0):
            bed_mod.dilute_layer1_age_on_deposition(
                bed=bed, t=time,
                layer1_mass_before=old_layer1_mass,
                deposited_mass=deposit_total,
            )

        # Per-class flux diagnostics in g/cm²/s.
        erosion_rate = (erosion_per_class / dt).astype("float32")
        deposition_rate = (deposition_arr / dt).astype("float32")
        self._write_per_class_flux(
            mesh, contracts.VAR_BED_EROSION_FLUX, time, erosion_rate
        )
        self._write_per_class_flux(
            mesh, contracts.VAR_BED_DEPOSITION_FLUX, time, deposition_rate
        )

        # Cache τ for the growth limiter on the next step.
        self._write_tau_diagnostics(mesh, time, tau_pa, tau_crit_surface)

        # --- (9) Bed elevation diagnostics -----------------------------
        bed_mod.update_bed_elevation(bed, t=time, dt_seconds=dt)

        # --- (10) Bedload step -----------------------------------------
        if self._bedload_solver is not None:
            # The standalone solver currently expects integer time-step
            # indexing for time-dimensioned bedload arrays. When given a
            # datetime label we resolve the integer index.
            t_for_bedload = self._resolve_time_index(mesh, time)
            # Build the surface-composition context for the chosen
            # transport function. Use the post-erosion / post-deposition
            # surface fractions (new_persed / new_active) so closures see
            # the bed state consistent with the rest of this step.
            surface_pers_post = _surface_class_fraction(new_persed, new_active)
            registry_context = self._build_bedload_registry_context(
                surface_class_fraction=surface_pers_post,
            )
            self._bedload_solver.step(
                mesh=mesh,
                time=t_for_bedload,
                tau_pa=tau_pa,
                dt_seconds=dt,
                registry_context=registry_context,
            )

        # --- (11) Stage net (E − D) on the mesh for Riverine -----------
        net_per_class = erosion_per_class - deposition_arr  # g/cm² per step
        for s, cls in enumerate(self.registry_classes):
            src_name = f"{cls.suspended_var}_source"
            self._write_per_face_source(mesh, src_name, time, net_per_class[:, s])

    # -----------------------------------------------------------------
    # Internal mesh I/O helpers
    # -----------------------------------------------------------------

    def _build_bedload_registry_context(
        self,
        surface_class_fraction: np.ndarray,
    ) -> dict:
        """Construct the registry-wide ``registry_context`` for the bedload step.

        The dict carries everything the seven shipped transport functions
        might need from the bed surface:

        ``surface_class_fraction``
            xr.DataArray of shape ``(nface, ssm_class)`` — per-class mass
            fractions on the topmost non-absent layer (sums to 1 in cells
            with surface mass).
        ``surface_sand_fraction``
            xr.DataArray of shape ``(nface,)`` — sum of fractions for
            classes whose D50 falls in the Wentworth sand window
            [62.5, 2000] μm. Consumed by Wilcock-Crowe (2003) eq. 4.
        ``surface_geometric_mean_um``
            xr.DataArray of shape ``(nface,)`` — mass-weighted log-mean
            grain size on the surface; consumed by Wilcock-Crowe as
            :math:`d_{sg}` in the hiding/exposure law.
        ``pe_ph_ratio``
            Scalar 1.0 — defaults to the uniform-bed limit for
            Wu-Wang-Jia (2000). A future enhancement will derive
            ``p_e/p_h`` from the full surface size-distribution per
            Wu, Wang & Jia (2000) eqs. 3-4.
        ``registry``
            The :class:`SedimentClassRegistry` itself, in case a future
            closure needs the full per-class metadata.

        Closures that don't read these keys ignore them
        (van Rijn, Engelund-Hansen, Toffaleti, Yang).
        """
        d50_array = np.asarray(self.registry_classes.d50_um_array, dtype="float64")
        sand_mask = _wentworth_sand_mask(d50_array)
        surface_sand_fraction = surface_class_fraction[:, sand_mask].sum(axis=-1)
        d_sg_um = _surface_geometric_mean_um(surface_class_fraction, d50_array)
        return {
            "surface_class_fraction": xr.DataArray(
                surface_class_fraction.astype("float64"),
                dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
            ),
            "surface_sand_fraction": xr.DataArray(
                surface_sand_fraction.astype("float64"),
                dims=(contracts.DIM_NFACE,),
            ),
            "surface_geometric_mean_um": xr.DataArray(
                d_sg_um.astype("float64"),
                dims=(contracts.DIM_NFACE,),
            ),
            "pe_ph_ratio": 1.0,
            "registry": self.registry_classes,
        }

    def _read_previous_tau(self, mesh: xr.Dataset, time) -> xr.DataArray:
        """Return τ at the previous time slot, or zeros at t=0."""
        n_face = mesh.sizes[contracts.DIM_NFACE]
        if contracts.VAR_BED_SHEAR_STRESS not in mesh.data_vars:
            return xr.DataArray(
                np.zeros(n_face, dtype="float64"),
                dims=(contracts.DIM_NFACE,),
            )
        idx = self._resolve_time_index(mesh, time)
        if idx <= 0:
            return xr.DataArray(
                np.zeros(n_face, dtype="float64"),
                dims=(contracts.DIM_NFACE,),
            )
        prev = mesh[contracts.VAR_BED_SHEAR_STRESS].isel(
            {contracts.DIM_TIME: idx - 1}
        )
        return prev.astype("float64")

    def _read_suspended_concentration(self, mesh: xr.Dataset, time) -> xr.DataArray:
        """Build a (nface, ssm_class) DataArray of suspended concentrations.

        Pulls each class's constituent value from the mesh at the
        requested time. Missing constituents (e.g. on a synthetic test
        mesh that hasn't seeded any) default to 0 mg/L.
        """
        n_face = mesh.sizes[contracts.DIM_NFACE]
        n_class = len(self.registry_classes)
        c_arr = np.zeros((n_face, n_class), dtype="float64")
        for s, cls in enumerate(self.registry_classes):
            cname = cls.suspended_var
            if cname not in mesh.data_vars:
                continue
            da = mesh[cname]
            if contracts.DIM_TIME in da.dims:
                # Try label-based selection; fall back to positional.
                try:
                    da = da.sel({contracts.DIM_TIME: time})
                except KeyError:
                    idx = self._resolve_time_index(mesh, time)
                    da = da.isel({contracts.DIM_TIME: idx})
            # Constrain to nface (Riverine mesh has nreal+1 ghost cells).
            vals = np.asarray(da.values, dtype="float64").reshape(-1)
            c_arr[: min(vals.size, n_face), s] = vals[: min(vals.size, n_face)]
        return xr.DataArray(
            c_arr,
            dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
        )

    def _read_bottom_depth(self, mesh: xr.Dataset, time) -> xr.DataArray:
        """Return bottom-water-layer depth in metres on (nface,).

        Uses ``face_hydraulic_depth`` if present; defaults to a 1-m
        column (so deposition's mass-cap is generous) on synthetic
        meshes that omit the field.
        """
        n_face = mesh.sizes[contracts.DIM_NFACE]
        if contracts.VAR_FACE_HYDRAULIC_DEPTH in mesh.data_vars:
            depth = mesh[contracts.VAR_FACE_HYDRAULIC_DEPTH]
            if contracts.DIM_TIME in depth.dims:
                try:
                    depth = depth.sel({contracts.DIM_TIME: time})
                except KeyError:
                    idx = self._resolve_time_index(mesh, time)
                    depth = depth.isel({contracts.DIM_TIME: idx})
            return depth.astype("float64")
        return xr.DataArray(
            np.ones(n_face, dtype="float64"),
            dims=(contracts.DIM_NFACE,),
        )

    def _resolve_time_index(self, mesh: xr.Dataset, time) -> int:
        """Convert a label-based time selector to a positional index."""
        if isinstance(time, (int, np.integer)):
            return int(time)
        if contracts.DIM_TIME in mesh.indexes:
            try:
                return int(mesh.indexes[contracts.DIM_TIME].get_loc(time))
            except (KeyError, TypeError):
                pass
        return 0

    def _write_per_class_flux(
        self,
        mesh: xr.Dataset,
        name: str,
        time,
        values_face_class: np.ndarray,
    ) -> None:
        """Write a (nface, n_class) per-class flux into a mesh diagnostic var."""
        if name not in mesh.data_vars:
            _ensure_per_class_bed_var(mesh, name)
        da = mesh[name]
        if contracts.DIM_TIME in da.dims:
            idx = self._resolve_time_index(mesh, time)
            da.values[idx, :, :] = values_face_class.astype(da.dtype)
        else:
            da.values[:, :] = values_face_class.astype(da.dtype)

    def _write_tau_diagnostics(
        self,
        mesh: xr.Dataset,
        time,
        tau_pa: xr.DataArray,
        tau_crit_pa: xr.DataArray,
    ) -> None:
        """Stamp the applied τ and surface τ_crit into the mesh diagnostics."""
        for var, da_in in (
            (contracts.VAR_BED_SHEAR_STRESS, tau_pa),
            (contracts.VAR_BED_CRITICAL_SHEAR_STRESS, tau_crit_pa),
        ):
            if var not in mesh.data_vars:
                continue
            target = mesh[var]
            arr = np.asarray(da_in.values, dtype=target.dtype)
            if contracts.DIM_TIME in target.dims:
                idx = self._resolve_time_index(mesh, time)
                target.values[idx, :] = arr
            else:
                target.values[:] = arr

    def _write_per_face_source(
        self,
        mesh: xr.Dataset,
        name: str,
        time,
        values_face: np.ndarray,
    ) -> None:
        """Write or overwrite a per-face source/sink array on the mesh."""
        n_face = mesh.sizes[contracts.DIM_NFACE]
        arr = np.asarray(values_face, dtype="float32").reshape(-1)
        if arr.size != n_face:
            raise ValueError(
                f"source values length {arr.size} != n_face {n_face} for {name}"
            )
        if name in mesh.data_vars and contracts.DIM_TIME in mesh[name].dims:
            idx = self._resolve_time_index(mesh, time)
            mesh[name].values[idx, :] = arr
        else:
            mesh[name] = ((contracts.DIM_NFACE,), arr.copy())
