"""Sediment-class definitions and registry.

A sediment class is a single discrete grain-size bin with associated
critical shear stresses, settling velocity, and solid density. SSM
transports one suspended-concentration constituent per class through
ClearWater-Riverine; the bed stores per-class mass fractions per layer.

Reference: SAND2008-5621 §"Sediment Bed Layers"; design spec §5.1.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

from . import contracts


@dataclass(frozen=True)
class SedimentClass:
    """One discrete sediment-size class.

    Parameters
    ----------
    label : str
        Short identifier (e.g. ``"silt_fine"``, ``"sand_medium"``). Used as
        the suffix for the suspended-concentration constituent registered
        with Riverine; see :func:`contracts.suspended_var_name`.
    d50_um : float
        Median grain diameter in micrometres.
    tau_ce_pa : float
        Critical shear stress for erosion (τ\\_ce) in pascals. If left as
        the sentinel ``None``, will be computed via Soulsby (1997) at
        registry-build time.
    tau_cs_pa : float
        Critical shear stress for suspension (τ\\_cs) in pascals. If left as
        ``None``, will be computed via van Rijn (1984) eqs 8–9.
    settling_cm_s : float
        Settling velocity in cm/s. If ``-1`` (the SEDZLJ sentinel) or
        ``None``, will be computed via Cheng (1997) from ``d50_um``.
    solid_density_g_cm3 : float
        Particle solid density (default 2.65, quartz).

    Notes
    -----
    Internally SSM stores grain size in μm and critical shear in Pa to
    match modern convention; the SEDflume loader converts the legacy CGS
    units (μm and dynes/cm²) on read.
    """

    label: str
    d50_um: float
    tau_ce_pa: Optional[float] = None
    tau_cs_pa: Optional[float] = None
    settling_cm_s: Optional[float] = None
    solid_density_g_cm3: float = contracts.DEFAULT_SOLID_SPECIFIC_GRAVITY

    @property
    def is_cohesive(self) -> bool:
        """True iff D50 is below the configured bedload cutoff."""
        return self.d50_um < contracts.DEFAULT_BEDLOAD_CUTOFF_UM

    @property
    def is_bedload_eligible(self) -> bool:
        """True iff D50 is at or above the bedload cutoff (sand and coarser)."""
        return not self.is_cohesive

    @property
    def suspended_var(self) -> str:
        """Riverine constituent name for this class's suspended concentration."""
        return contracts.suspended_var_name(self.label)

    @property
    def advection_coef_var(self) -> str:
        """Per-class advection-coefficient field name (used only by the
        Riverine bedload solver, see :class:`bedload.BedloadRiverineConstituent`)."""
        return contracts.advection_coef_var_name(self.label)


@dataclass
class SedimentClassRegistry:
    """Ordered, validated collection of sediment classes for a single SSM run.

    Order matters: the index assigned to each class becomes the
    :data:`contracts.DIM_CLASS` index in all bed-state DataArrays. Once
    constructed and frozen (:meth:`freeze`) the registry is immutable.
    """

    classes: list[SedimentClass] = field(default_factory=list)
    _frozen: bool = field(default=False, init=False, repr=False)

    @classmethod
    def from_iterable(cls, items: Iterable[SedimentClass]) -> "SedimentClassRegistry":
        reg = cls(list(items))
        reg.freeze()
        return reg

    def freeze(self) -> None:
        """Validate uniqueness, sort-stability of labels, and lock the registry."""
        labels = [c.label for c in self.classes]
        if len(labels) != len(set(labels)):
            raise ValueError(f"Duplicate sediment-class labels: {labels}")
        if len(self.classes) == 0:
            raise ValueError("SedimentClassRegistry must contain at least one class")
        self._frozen = True

    def __len__(self) -> int:
        return len(self.classes)

    def __iter__(self):
        return iter(self.classes)

    def __getitem__(self, idx: int) -> SedimentClass:
        return self.classes[idx]

    def by_label(self, label: str) -> SedimentClass:
        for c in self.classes:
            if c.label == label:
                return c
        raise KeyError(label)

    @property
    def labels(self) -> list[str]:
        return [c.label for c in self.classes]

    @property
    def d50_um_array(self):
        """Numpy array of D50 in μm, ordered by class index."""
        import numpy as np
        return np.asarray([c.d50_um for c in self.classes], dtype="float64")
