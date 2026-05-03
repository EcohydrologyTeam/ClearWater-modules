"""SSM canonical data contracts (FROZEN — do not modify without coordinating).

This module is the single source of truth for the xarray variable names,
dimensions, units, and dtypes that SSM consumes from and writes to the
shared mesh dataset. All SSM submodules (shear, bed, erosion, deposition,
bedload, etc.) and ESM coupling code MUST reference the constants here
rather than hard-coding strings, so that downstream changes to the schema
require touching only this file.

Reference: SSM design specification, section 6
(`ClearWater-Riverine-streaming/design/ssm_design_spec.md`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final


# ---------------------------------------------------------------------------
# Mesh dimensions (match clearwater_riverine.variables conventions)
# ---------------------------------------------------------------------------

DIM_TIME: Final[str] = "time"
DIM_NFACE: Final[str] = "nface"
DIM_NEDGE: Final[str] = "nedge"
DIM_LAYER: Final[str] = "ssm_layer"     # bed-layer index, 1..K_B (top-down)
DIM_CLASS: Final[str] = "ssm_class"     # sediment-class index, 0..N_class-1


# ---------------------------------------------------------------------------
# Inputs SSM reads from the mesh dataset (populated by RAS reader, ESM, or external)
# ---------------------------------------------------------------------------

# RAS-derived hydraulics (already populated by clearwater_riverine.io.hdf)
VAR_VOLUME: Final[str] = "volume"                                 # (time, nface) m^3
VAR_WATER_SURFACE_ELEVATION: Final[str] = "water_surface_elev"    # (time, nface) m
VAR_FACE_HYDRAULIC_DEPTH: Final[str] = "face_hydraulic_depth"     # (time, nface) m
VAR_EDGE_VELOCITY: Final[str] = "edge_velocity"                   # (time, nedge) m/s
VAR_MANNINGS_N: Final[str] = "mannings_n"                         # (nface,)     s/m^(1/3)

# Optional: direct-import bed shear (Mode A driver)
VAR_BED_SHEAR_STRESS_INPUT: Final[str] = "bed_shear_stress_input"  # (time, nface) Pa

# Optional: wave forcing (Mode C driver)
VAR_WAVE_ORBITAL_VELOCITY: Final[str] = "wave_orbital_velocity"   # (time, nface) m/s
VAR_WAVE_FREQUENCY: Final[str] = "wave_frequency"                 # (time, nface) rad/s
VAR_WAVE_DIRECTION: Final[str] = "wave_direction"                 # (time, nface) rad

# Optional: ESM-supplied vegetation feedback fields
VAR_COMPOSITE_MANNINGS_N: Final[str] = "composite_manning_n"      # (time, nface) s/m^(1/3)
VAR_VEGETATION_BIOSTABILIZATION: Final[str] = "vegetation_biostabilization"  # (time, nface) [0,1]
VAR_VEGETATION_ROOT_COHESION: Final[str] = "vegetation_root_cohesion"        # (time, nface) Pa
VAR_VEGETATION_FRONTAL_AREA: Final[str] = "vegetation_frontal_area"          # (time, nface) m^2/m^2


# ---------------------------------------------------------------------------
# SSM-owned bed state on the mesh dataset
# ---------------------------------------------------------------------------

# Per-layer dynamic state
VAR_BED_LAYER_MASS: Final[str] = "ssm_bed_layer_mass"                       # (time, nface, ssm_layer) g/cm^2  -- TSED
VAR_BED_LAYER_INITIAL_MASS: Final[str] = "ssm_bed_layer_initial_mass"       # (nface, ssm_layer)        g/cm^2  -- TSED0 (for SEDflume depth interp)
VAR_BED_CLASS_FRACTION: Final[str] = "ssm_bed_class_fraction"               # (time, nface, ssm_layer, ssm_class) [0,1] -- PERSED
VAR_BED_LAYER_ACTIVE: Final[str] = "ssm_bed_layer_active"                   # (time, nface, ssm_layer)  int8 {0,1,2} -- LAYERACTIVE
VAR_BED_LAYER_TAUCRIT: Final[str] = "ssm_bed_layer_taucrit"                 # (time, nface, ssm_layer)  Pa  -- TAUCOR
VAR_BED_LAYER_BULK_DENSITY: Final[str] = "ssm_bed_layer_bulk_density"       # (nface, ssm_layer)        g/cm^3 -- BULKDENS (constant per SEDZLJ)
VAR_BED_LAYER_THICKNESS: Final[str] = "ssm_bed_layer_thickness"             # (time, nface, ssm_layer)  m  -- HBED

# Per-cell aggregate bed state
VAR_BED_TOTAL_THICKNESS: Final[str] = "ssm_bed_total_thickness"             # (time, nface) m
VAR_BED_D50_SURFACE: Final[str] = "ssm_bed_d50_surface"                     # (time, nface) μm  -- D50AVG of surface layer
VAR_BED_ELEVATION: Final[str] = "ssm_bed_elevation"                         # (time, nface) m   -- absolute
VAR_BED_CHANGE: Final[str] = "ssm_bed_change"                               # (time, nface) m   -- per-step delta
VAR_BED_CUMULATIVE_CHANGE: Final[str] = "ssm_bed_cumulative_change"         # (time, nface) m   -- running sum

# Bedload (only populated if ICALC_BL > 0)
VAR_BEDLOAD_MASS: Final[str] = "ssm_bedload_mass"                           # (time, nface, ssm_class) g/cm^2  -- CBL

# Computed shear-stress fields (used internally and exposed to ESM)
VAR_BED_SHEAR_STRESS: Final[str] = "ssm_bed_shear_stress"                   # (time, nface) Pa  -- τ_b applied
VAR_BED_CRITICAL_SHEAR_STRESS: Final[str] = "ssm_bed_critical_shear_stress"  # (time, nface) Pa  -- τ_crit at surface

# Per-class diagnostic fluxes (optional but recommended for output / debugging)
VAR_BED_EROSION_FLUX: Final[str] = "ssm_bed_erosion_flux"                   # (time, nface, ssm_class) g/cm^2/s
VAR_BED_DEPOSITION_FLUX: Final[str] = "ssm_bed_deposition_flux"             # (time, nface, ssm_class) g/cm^2/s


# ---------------------------------------------------------------------------
# Suspended-sediment constituent naming
# ---------------------------------------------------------------------------

#: Prefix for suspended-sediment-class constituents registered with Riverine.
#: Full name format: f"{SUSPENDED_PREFIX}{class_label}", e.g. "ssm_suspended_silt_fine".
#: Riverine treats each as an independent transported scalar (mg/L).
SUSPENDED_PREFIX: Final[str] = "ssm_suspended_"


def suspended_var_name(class_label: str) -> str:
    """Canonical suspended-class constituent name for Riverine registration."""
    return f"{SUSPENDED_PREFIX}{class_label}"


# ---------------------------------------------------------------------------
# Per-class advection-coefficient field naming (for Bedload Mode B / Riverine ext)
# ---------------------------------------------------------------------------

#: Prefix for per-constituent advection-coefficient fields. Used when bedload
#: is implemented as a Riverine constituent with its own face-velocity field
#: (van Rijn bedload velocity, not water velocity).
#: Full name format: f"{ADVECTION_COEF_PREFIX}{class_label}", on (time, nedge).
ADVECTION_COEF_PREFIX: Final[str] = "ssm_advection_coef_"


def advection_coef_var_name(class_label: str) -> str:
    """Canonical per-class advection-coefficient edge field for Riverine."""
    return f"{ADVECTION_COEF_PREFIX}{class_label}"


# ---------------------------------------------------------------------------
# Schema descriptors (for validation, documentation, and Zarr encoding)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class VarSpec:
    """Descriptor for a single SSM-managed xarray variable."""
    name: str
    dims: tuple[str, ...]
    dtype: str
    units: str
    description: str
    role: str  # "input" | "bed_state" | "diagnostic" | "constituent" | "advection_coef"


#: Bed-state variables that SSM owns and writes each sediment time step.
#: Used to register with `_release_to_stream` for streaming Zarr output.
BED_STATE_SPECS: Final[tuple[VarSpec, ...]] = (
    VarSpec(VAR_BED_LAYER_MASS, (DIM_TIME, DIM_NFACE, DIM_LAYER),
            "float32", "g cm-2", "Per-layer dry sediment mass (TSED)", "bed_state"),
    VarSpec(VAR_BED_LAYER_INITIAL_MASS, (DIM_NFACE, DIM_LAYER),
            "float32", "g cm-2", "Per-layer initial mass (TSED0); used for SEDflume depth interpolation", "bed_state"),
    VarSpec(VAR_BED_CLASS_FRACTION, (DIM_TIME, DIM_NFACE, DIM_LAYER, DIM_CLASS),
            "float32", "1", "Per-layer per-class mass fraction (PERSED), normalized so layer-sum = 1", "bed_state"),
    VarSpec(VAR_BED_LAYER_ACTIVE, (DIM_TIME, DIM_NFACE, DIM_LAYER),
            "int8", "1", "Layer state: 0=absent, 1=active/deposited, 2=in-place core (LAYERACTIVE)", "bed_state"),
    VarSpec(VAR_BED_LAYER_TAUCRIT, (DIM_TIME, DIM_NFACE, DIM_LAYER),
            "float32", "Pa", "Per-layer critical shear stress for erosion (TAUCOR)", "bed_state"),
    VarSpec(VAR_BED_LAYER_BULK_DENSITY, (DIM_NFACE, DIM_LAYER),
            "float32", "g cm-3", "Per-layer dry bulk density (BULKDENS); SEDZLJ holds this constant", "bed_state"),
    VarSpec(VAR_BED_LAYER_THICKNESS, (DIM_TIME, DIM_NFACE, DIM_LAYER),
            "float32", "m", "Per-layer thickness (HBED) = 0.01 * TSED / BULKDENS", "bed_state"),
    VarSpec(VAR_BED_TOTAL_THICKNESS, (DIM_TIME, DIM_NFACE),
            "float32", "m", "Total bed thickness (sum over layers)", "bed_state"),
    VarSpec(VAR_BED_D50_SURFACE, (DIM_TIME, DIM_NFACE),
            "float32", "um", "Mass-weighted mean D50 of surface layer (D50AVG)", "bed_state"),
    VarSpec(VAR_BED_ELEVATION, (DIM_TIME, DIM_NFACE),
            "float32", "m", "Absolute bed surface elevation", "bed_state"),
    VarSpec(VAR_BED_CHANGE, (DIM_TIME, DIM_NFACE),
            "float32", "m", "Per-step bed-elevation change (negative = erosion)", "diagnostic"),
    VarSpec(VAR_BED_CUMULATIVE_CHANGE, (DIM_TIME, DIM_NFACE),
            "float32", "m", "Running sum of bed-elevation change since simulation start", "diagnostic"),
    VarSpec(VAR_BEDLOAD_MASS, (DIM_TIME, DIM_NFACE, DIM_CLASS),
            "float32", "g cm-2", "Per-class bedload mass per unit bed area (CBL)", "bed_state"),
    VarSpec(VAR_BED_SHEAR_STRESS, (DIM_TIME, DIM_NFACE),
            "float32", "Pa", "Applied bed shear stress (computed or imported)", "diagnostic"),
    VarSpec(VAR_BED_CRITICAL_SHEAR_STRESS, (DIM_TIME, DIM_NFACE),
            "float32", "Pa", "Critical shear stress for erosion at surface (τ_crit)", "diagnostic"),
    VarSpec(VAR_BED_EROSION_FLUX, (DIM_TIME, DIM_NFACE, DIM_CLASS),
            "float32", "g cm-2 s-1", "Per-class erosion flux from bed (positive = up)", "diagnostic"),
    VarSpec(VAR_BED_DEPOSITION_FLUX, (DIM_TIME, DIM_NFACE, DIM_CLASS),
            "float32", "g cm-2 s-1", "Per-class deposition flux to bed (positive = down)", "diagnostic"),
)


#: Lookup table for any bed-state variable name.
BED_STATE_BY_NAME: Final[dict[str, VarSpec]] = {s.name: s for s in BED_STATE_SPECS}


# ---------------------------------------------------------------------------
# Streaming defaults
# ---------------------------------------------------------------------------

#: Bed state evolves more slowly than suspended sediment, so it can flush less
#: often. Default: bed-state Zarr flush every Nth transport-streaming flush.
DEFAULT_BED_STREAMING_INTERVAL_MULTIPLIER: Final[int] = 10

#: Zarr chunk size hints. Time chunk should match the streaming interval; spatial
#: chunk targets ~1M cells per chunk for the Albany-scale (587k) reference case.
DEFAULT_TIME_CHUNK: Final[int] = 100
DEFAULT_FACE_CHUNK: Final[int] = 600_000
DEFAULT_LAYER_CHUNK: Final[int] = 8
DEFAULT_CLASS_CHUNK: Final[int] = 8


# ---------------------------------------------------------------------------
# Physical / unit constants used in equations
# ---------------------------------------------------------------------------

#: Gravitational acceleration (cm/s^2). SEDZLJ uses CGS for sediment math.
G_CGS: Final[float] = 980.0

#: Kinematic viscosity of water at ~20°C (cm^2/s).
NU_CGS: Final[float] = 0.01

#: Default sediment specific gravity (quartz).
DEFAULT_SOLID_SPECIFIC_GRAVITY: Final[float] = 2.65

#: Default water density (g/cm^3).
DEFAULT_WATER_DENSITY_CGS: Final[float] = 1.0

#: Default cutoff between cohesive (Krone deposition) and non-cohesive
#: (Gessler deposition + bedload eligibility) sediment classes (μm).
DEFAULT_BEDLOAD_CUTOFF_UM: Final[float] = 64.0

#: Default skin-friction roughness fallback (m). Used by the current-only
#: shear driver when D50_surface is unavailable.
DEFAULT_ZB_SKIN_M: Final[float] = 1.5e-3

#: Default active-layer thickness multiplier (TACTM in SEDZLJ).
DEFAULT_TACTM: Final[float] = 2.0

#: Default per-step shear-stress growth limiter (fraction). Suppresses shock-
#: induced oscillations near abrupt velocity changes; mirrors SEDZLJ_SHEAR
#: behaviour at s_shear.f90:315.
DEFAULT_SHEAR_GROWTH_LIMIT: Final[float] = 0.10
