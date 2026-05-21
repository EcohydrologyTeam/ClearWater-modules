"""v3 ``Process`` and ``ProcessFactory`` abstractions.

Class definitions in-place. Originally inherited verbatim from v2 by
re-export; the v3-self-sufficient refactor moved the class bodies in-tree
so v3 owns its own ``ProcessFactory`` registry. See the architecture
spec section 3 non-goal note: no new framework — these classes mirror
v2's contract exactly, only relocated.

------------------------------------------------------------------
Integration patterns (guideline, not enforced contract)
------------------------------------------------------------------

v3-native processes implement ``run(time, registry)`` to advance state
by one substep. Two integration patterns are common:

(a) Compute a per-second rate of change and apply Forward Euler
    ``state_new = state_old + rate * time_step_seconds``. This is
    suitable for processes with simple linear kinetics (e.g., NSM1
    reaction terms).
(b) Compute the per-substep ``delta_state`` directly and add it to the
    current state. This is suitable for processes whose update depends
    non-linearly on the substep length (e.g., the v3 ``Temperature``
    thin-water guards).

Both patterns are valid; the choice belongs to the process author. M16
(review-findings 2026-05-04) demoted an earlier numbered "5-step
contract" to this guideline because real v3 processes (notably
``Temperature``) follow pattern (b) and the contract wording did not
reflect actual practice.

------------------------------------------------------------------
M5 ordering contract for ``init_process`` and ``from_hotstart``
------------------------------------------------------------------

The v3 ``Model.__init_model`` invokes per-process initialization in a
specific ordering that processes MUST respect when they implement
hotstart support:

1. Data sources (or first chunk) are loaded into the registry.
2. If a hotstart dataset is supplied, the registry is **seeded** from
   the saved dataset (``Model.__seed_from_hotstart``).
3. ``init_process`` is called on every process.
4. If a hotstart dataset is supplied, ``from_hotstart`` is called on
   every process that defines it.

The ordering means ``init_process`` runs **before** ``from_hotstart``,
and so:

- ``init_process`` MUST set the process's substep-internal state
  assuming a *fresh start*. It does not need to know whether a
  hotstart was supplied; it sets defaults appropriate to a fresh run.
- ``from_hotstart`` MUST override those fresh-start defaults with
  values restored from the saved dataset's ``attrs``. If a process
  adds new internal substep state in a future revision but forgets to
  also handle it in ``from_hotstart``, fresh-start runs and
  hotstart-resume runs will silently diverge starting at the first
  substep where the un-restored state matters.

Together, the two methods ensure that fresh-start and
hotstart-resume produce equivalent post-init state when the same
initial conditions are loaded. This invariant is a load-bearing
property of the v3 hotstart design (TSM design spec section 3.2).

The default ``Process.from_hotstart`` is implicit (``getattr`` -based
in the Model). Processes opt in by defining the method; processes
that don't have any internal substep state are free to omit it.
"""

from abc import ABC, abstractmethod
from datetime import datetime, timedelta
import functools
from typing import Callable, TYPE_CHECKING

from clearwater_data.variables import VariableRegistry

if TYPE_CHECKING:
    from clearwater_modules_v3.model import Model


class Process(ABC):
    """
    Base class for processes. Defines a class-level variable registry and a run method to be implemented by subclasses.
    """

    # Class-level definition of variables associated with this process
    variables = []
    time_step_seconds: int

    # Phase H-9 (2026-05-21): per-substep sibling-coupling dependency
    # declaration. Each subclass overrides this tuple with the
    # ``process_name()`` of any sibling process whose step-scoped
    # ``self.<name>`` rate caches it reads inside ``run()``. Model
    # validates the registered process order at init time and raises
    # if a reader is constructed before its writer. Default empty
    # tuple = no sibling dependencies (the legacy contract). Examples:
    #   * DOX reads Nitrogen.nitrification_flux_rate +
    #     FloatingAlgae.algal_growth_rate + BenthicAlgae.balgae_*_rate +
    #     Carbon.doc_dic_oxidation_rate.
    #   * Phosphorus reads FloatingAlgae.algal_growth_rate +
    #     BenthicAlgae.balgae_growth_rate.
    upstream_processes: tuple[str, ...] = ()

    def __init__(self, time_step: timedelta) -> None:
        self.time_step = time_step
        self.time_step_seconds = self.time_step.total_seconds()

    @classmethod
    def from_config(
        cls, config: dict, variable_registry: VariableRegistry
    ) -> "Process":
        return cls(**config)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """
        Initialize of the process.
        """
        # base method assumes no initialization is needed
        pass

    def validate(self, registry: VariableRegistry) -> None:
        """
        Validate the process.
        """
        for variable in self.variables:
            if variable not in registry:
                raise ValueError(
                    f"Variable {variable} not found. Are you sure you provided a valid configuration for {variable}?"
                )

    def finalize_process(self, model: "Model", registry: VariableRegistry) -> None:
        return None

    @abstractmethod
    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """
        Run the process. To be implemented by subclasses.
        """
        raise NotImplementedError

    def process_name(self) -> str:
        return self.__class__.__name__


class ProcessFactory:
    processes: dict[str, Callable] = {}

    @classmethod
    def from_config(
        cls, process_name: str, config: dict, variable_registry: VariableRegistry
    ) -> "Process":
        if process_name not in cls.processes:
            raise ValueError(
                f"Process type {process_name} not registered. Did you register the process at the {__name__}"
            )
        return cls.processes[process_name](config, variable_registry)

    @classmethod
    def register(cls, process_name: str):
        def init_method(from_config_method: Callable) -> Callable:
            cls.processes[process_name] = from_config_method

            @functools.wraps(from_config_method)
            def from_config_method(config: dict) -> "Process":
                return from_config_method(config)

            return from_config_method

        return init_method


__all__ = ["Process", "ProcessFactory"]
