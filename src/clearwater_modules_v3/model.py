"""v3 ``Model`` orchestration.

The v3 ``Model`` is a re-implementation of the upstream v2 ``Model`` that
adds three orchestration-level capabilities ported from v1, and resolves
the four chunking TODOs left in v2's ``__process_loop_chunked``:

1. **Kernel optimization.** The per-substep dispatch precomputes which
   processes fire at each step index, replacing the per-step modulo check
   with an O(1) tuple lookup.
2. **Wet-mask gating.** A registry-level wet-mask is applied at the
   orchestration layer: after each process runs, NaN is written for the
   process's declared variables on dry cells. Processes do not need to
   carry their own ``xr.where(volume > 0, ...)`` masks (the per-process
   mask in ``Temperature.run`` becomes redundant with this Model-level
   mask and is slated for removal in a Phase 4 cleanup; gap-analysis
   row N8).
3. **Hotstart from xarray Dataset.** The constructor accepts a
   ``hotstart_dataset`` and ``hotstart_timestep``; when present, the
   registry is seeded from the saved dataset before processes are
   initialized. Each Process may optionally implement
   ``from_hotstart(state)`` to restore process-internal substep state;
   v3's ``Temperature`` uses this to disable its
   ``__skip_first_time_step`` flag after a hotstart.
4. **Chunking.** Chunk boundaries are precomputed as **integer step
   indices** (``interior_chunk_step_indices = set[int]``), not as
   datetime values. Step-index comparison is exact-integer,
   timezone-independent, and immune to floating-point drift in
   ``current_time += time_step`` arithmetic. The per-chunk write of
   the trailing partial chunk is handled exactly once after the main
   loop exits. (Mirrors the riverine pattern of "precompute once,
   compare in the hot loop"; uses step-index instead of datetime
   identity per the C7 review fix.)

The class is a drop-in for v2 ``Model`` when none of the new constructor
kwargs are passed: behavior is identical (modulo the kernel-optimization
fast path, which is observationally identical and only changes
microseconds-per-step performance).

Subclassing v2's ``Model`` was considered and rejected because v2's
private state uses Python name-mangling (``_Model__processes`` etc.),
which makes inheritance fragile. v3's ``Model`` is a parallel
implementation that mirrors v2's lifecycle and method names so existing
v2 code that knows how to use ``Model`` works against v3 unchanged.
"""

from __future__ import annotations

import os
from datetime import datetime, timedelta
from logging import getLogger
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd
import xarray as xr

from clearwater_data.io.base import ChunkedDataSource, DataSource
from clearwater_data.io.zarr import ChunkedZarrDataStore
from clearwater_data.variables import VariableRegistry

from clearwater_modules_v3.processes.base import Process
from clearwater_modules_v3.utils.numerics import Diagnostics

LOGGER = getLogger(__name__)


class Model:
    """v3 orchestration model. See module docstring for capabilities."""

    def __init__(
        self,
        processes: Iterable[Process],
        variable_registry: VariableRegistry,
        variable_data_sources: dict[str, DataSource | ChunkedDataSource],
        start_time: datetime,
        end_time: datetime,
        time_step: timedelta,
        output_variables: Iterable[str],
        simulation_directory: os.PathLike | None = None,
        chunk_size: timedelta | None = None,
        # ----- v3 additions -----
        wet_mask_variable: str | None = None,
        wet_mask_threshold: float = 0.0,
        wet_mask_provider: Callable[[datetime, VariableRegistry], xr.DataArray] | None = None,
        hotstart_dataset: xr.Dataset | None = None,
        hotstart_timestep: datetime | str | int | None = None,
    ) -> None:
        """Initialize the model.

        Args (v3 additions; everything above is unchanged from v2):
            wet_mask_variable: Registry variable name to threshold for the
                wet/dry mask. When set, the value at each substep is read
                from the registry, compared against ``wet_mask_threshold``,
                and the resulting boolean array is used to write NaN on
                dry cells for every variable each Process declares. ``None``
                disables the mask (legacy v2 behavior).
            wet_mask_threshold: Threshold value paired with
                ``wet_mask_variable``. The comparison is a
                **strict-inequality**: cells with ``value > threshold``
                are wet, cells with ``value <= threshold`` are dry.

                The default of ``0.0`` has the literal semantic "any
                positive value is wet". For variables that represent
                a physical volume, surface area, or depth, this default
                is too permissive: a cell with ``value=1e-300`` (well
                below any meaningful physical threshold and within
                machine round-off) is classified as wet, and any
                downstream divide ``X / value`` produces ``inf`` or
                ``NaN`` that the wet-mask cannot subsequently clean up.

                **Recommendation (caller-side, not enforced by the
                Model):** when the variable has units of m^3, m^2, or m,
                set ``wet_mask_threshold`` to a small positive epsilon
                appropriate to the variable's units, e.g.
                ``1e-6`` m^2 for ``wetted_surface_area``. Values close
                to zero but above the epsilon are still classified as
                dry, avoiding the ``1 / very_small`` numerical hazard.
                The Model itself does not pick an epsilon for you,
                because the appropriate value is variable- and
                application-dependent (M11; review-findings
                2026-05-04).
            wet_mask_provider: Optional callable that takes
                ``(current_time, registry)`` and returns a boolean DataArray
                wet-mask. If provided, takes precedence over
                ``wet_mask_variable``/``wet_mask_threshold``. Useful for
                tests and programmatic mask construction without going
                through the registry.
            hotstart_dataset: xarray Dataset whose state is loaded into
                the registry as the simulation initial condition. ``None``
                means "use ``variable_data_sources`` as initial condition"
                (legacy v2 behavior).
            hotstart_timestep: Slice of ``hotstart_dataset`` to use for the
                initial condition. May be a datetime, a string parseable
                by pandas, or an integer index. ``None`` defaults to the
                last time slice of the saved dataset.
        """
        self.__processes: tuple[Process, ...] = tuple(processes)
        self.__registry: VariableRegistry = variable_registry
        self.__variable_data_sources: dict[
            str, DataSource | ChunkedDataSource
        ] = variable_data_sources
        self.__start_time: datetime = start_time
        self.__end_time: datetime = end_time
        self.__time_step: timedelta = time_step
        self.__output_variables: list[str] = list(output_variables) if output_variables else []
        # C1 fix (review-findings 2026-05-04): wrap with Path so the
        # ``self.__simulation_directory / "model_outputs.zarr"`` operator
        # in __init_output_source works on the default config. Previously
        # this was a bare ``str`` "." which raised TypeError on the
        # default-constructed Model whenever ``output_variables`` was
        # non-empty.
        #
        # m12 fix (review-findings 2026-05-04): the previous truthy check
        # ``if simulation_directory else Path(".")`` rewrote any
        # falsy-but-non-None argument (including ``""`` and ``Path("")``)
        # to ``Path(".")`` silently. The predicate below treats
        # ``None``, ``""``, and ``Path("")`` as "not provided" and falls
        # back to cwd; every other value (including absolute paths,
        # relative paths, and ``Path(".")`` itself) is taken as-is.
        if simulation_directory is None or simulation_directory in ("", Path("")):
            self.__simulation_directory: Path = Path(".")
        else:
            self.__simulation_directory = Path(simulation_directory)

        self.__chunked_mode: bool = chunk_size is not None
        self.__chunk_size: timedelta | None = chunk_size

        self.__wet_mask_variable: str | None = wet_mask_variable
        self.__wet_mask_threshold: float = wet_mask_threshold
        self.__wet_mask_provider = wet_mask_provider
        self.__hotstart_dataset: xr.Dataset | None = hotstart_dataset
        self.__hotstart_timestep = hotstart_timestep

        # Run-level diagnostics (resolved Q7, NSM1 design spec Section 14):
        # the clip-with-log contract for state integrators routes per-step
        # clip events here so Tier 1 conservation tests can assert
        # ``model.diagnostics.clip_events == 0`` directly. Public by design;
        # Process classes capture a reference in ``init_process`` and pass
        # it to ``clip_negative_state`` from each ``Process.run``.
        self.diagnostics: Diagnostics = Diagnostics()

        self.__init_complete: bool = False
        # M10 fix (review-findings 2026-05-04): a second ``run()`` call
        # against the same Model instance silently re-iterates from
        # ``start_time`` against an already-advanced registry. The flag
        # below is set after the substep loop exits and consulted on
        # subsequent ``run()`` calls so the second call raises rather
        # than corrupting state.
        self.__run_complete: bool = False
        self.__output_data_store: ChunkedZarrDataStore | None = None

        # Filled in by __init_model. ``__process_schedule[i]`` is the tuple
        # of processes that fire at step index ``i``. Precomputed so the
        # hot loop avoids modulo arithmetic per step (gap-analysis O3).
        self.__process_schedule: tuple[tuple[Process, ...], ...] = ()

    # ---------- public lifecycle ----------

    def validate(self) -> None:
        if self.__start_time >= self.__end_time:
            raise ValueError("Start time must be before end time.")

    def init_model(self) -> None:
        self.__init_model()
        self.__init_complete = True

    def run(self) -> None:
        # M10 fix (review-findings 2026-05-04): a second ``run()`` against
        # the same Model instance is a footgun — the registry has already
        # been advanced from ``start_time`` to ``end_time``, so a second
        # call would re-iterate against post-final state without the
        # caller noticing. Disallow it explicitly. Callers that genuinely
        # need to re-run the same configuration should construct a fresh
        # Model with the same arguments (which re-loads the data sources
        # into a fresh registry).
        if self.__run_complete:
            raise RuntimeError(
                "Model.run() has already completed for this instance. "
                "Construct a fresh Model to run the configuration again; "
                "the current registry has been advanced from start_time "
                "to end_time and a second run() would re-iterate against "
                "post-final state."
            )
        if not self.__init_complete:
            self.__init_model()
            self.__init_complete = True
        if self.__chunked_mode:
            self.__process_loop_chunked()
        else:
            self.__process_loop_full()
        self.__run_complete = True

    def has_process(self, process_type: type[Process] | str) -> bool:
        if isinstance(process_type, str):
            return any(
                p.process_name().lower() == process_type.lower()
                for p in self.__processes
            )
        return any(isinstance(p, process_type) for p in self.__processes)

    def get_process(self, process_type: type[Process] | str) -> Process:
        if not self.has_process(process_type):
            raise ValueError(f"Process {process_type} not found in model.")
        if isinstance(process_type, str):
            return next(
                p
                for p in self.__processes
                if p.process_name().lower() == process_type.lower()
            )
        return next(p for p in self.__processes if isinstance(p, process_type))

    # ---------- v3 hotstart export hook ----------

    def to_hotstart(self) -> dict[str, dict[str, Any]]:
        """Collect per-process substep state for a hotstart save.

        Returns a dict mapping ``process_name`` to whatever the process's
        ``to_hotstart()`` returns (typically a small dict of internal
        flags). Processes that do not implement ``to_hotstart`` are
        skipped. This is the inverse of the ``hotstart_dataset`` /
        ``from_hotstart`` import path. The Model itself does not save the
        registry's xarray data here — that is the caller's job — only the
        opaque per-process state that doesn't live in the registry.
        """
        snapshot: dict[str, dict[str, Any]] = {}
        for process in self.__processes:
            to_hotstart = getattr(process, "to_hotstart", None)
            if callable(to_hotstart):
                snapshot[process.process_name()] = to_hotstart()
        return snapshot

    # ---------- private lifecycle ----------

    def __init_model(self) -> None:
        """Initialize the registry, processes, schedule, and output store.

        Ordering invariant (M5 contract, review-findings 2026-05-04):

        1. Load every data source (or first chunk) into the registry.
        2. If a hotstart dataset was provided, overwrite registry contents
           with the saved-dataset slice at ``hotstart_timestep``.
        3. Run ``init_process`` on every process. **By contract,
           ``init_process`` sets each process's substep-internal state
           assuming a fresh start.** A process author who adds new
           internal substep state and writes initialization in
           ``init_process`` MUST also restore that state from the
           saved dataset in ``from_hotstart`` (step 4 below).
        4. If a hotstart dataset was provided, call ``from_hotstart`` on
           every process that defines it. **By contract,
           ``from_hotstart`` MUST override the fresh-start defaults set
           by ``init_process`` with the values saved in the hotstart
           dataset's ``attrs``.**

        Steps 3 and 4 together ensure that a fresh-start run and a
        hotstart-resume run produce equivalent post-init state when
        identical initial conditions are loaded. The cost of forgetting
        this contract is silent: the fresh-start path and the
        hotstart-resume path will diverge at the first substep where the
        un-restored state matters.

        Step 5 builds the precomputed process firing schedule.
        Step 6 initializes the output data store.
        Step 7 (M14) validates that ``wet_mask_variable``, if set, is
        registered, so a typo is caught at init rather than at the first
        substep.
        """
        # Step 1: load every data source (or first chunk) into the registry.
        # M7 fix (review-findings 2026-05-04): when chunk_size is None
        # (non-chunked simulation), a ChunkedDataSource was previously
        # asked for a window of just one substep; only the first slice
        # of the dataset would be loaded. The full simulation duration
        # must be read in non-chunked mode so subsequent substeps have
        # data available.
        if self.__chunked_mode:
            chunk_end = self.__start_time + self.__chunk_size
        else:
            chunk_end = self.__end_time
        for variable_name, data_source in self.__variable_data_sources.items():
            if isinstance(data_source, ChunkedDataSource):
                data = data_source.read_chunk(self.__start_time, chunk_end)
            else:
                data = data_source.read(variable_name)
            self.__registry.register(variable_name, data)

        # Step 2: optionally seed the registry from the hotstart dataset.
        # This overwrites whatever was loaded in step 1 for variables that
        # appear in the saved dataset.
        if self.__hotstart_dataset is not None:
            self.__seed_from_hotstart()

        # Step 3 (M5 contract): process-level initialization. Processes
        # see the post-hotstart registry state. ``init_process`` sets
        # default substep-internal state assuming a *fresh start*; if the
        # process has any internal substep state that needs to survive a
        # hotstart, ``from_hotstart`` (step 4) must override the
        # defaults set here.
        for process in self.__processes:
            process.init_process(self, self.__registry)

        # Phase H-9 (2026-05-21): validate per-substep process order
        # against the declared ``upstream_processes`` DAG. Each
        # Process subclass declares the names of sibling processes
        # whose step-scoped ``self.<rate>`` caches it reads inside
        # ``run()``. The Model invokes processes in the order of
        # ``self.__processes`` per substep, so if a reader is
        # constructed before its writer, the reader sees the
        # PREVIOUS substep's cache (a silent one-substep lag).
        # Raise at init time so the user fixes the construction
        # order rather than discovering the lag via a calibration-
        # off-by-N-percent mystery.
        seen: set[str] = set()
        for process in self.__processes:
            name = process.process_name()
            upstream = tuple(
                getattr(process, "upstream_processes", ())
            )
            missing = [u for u in upstream if u not in seen]
            if missing:
                raise ValueError(
                    f"Process order violation: {name!r} declares "
                    f"upstream_processes = {upstream!r} but the "
                    f"processes {missing!r} are NOT registered before "
                    f"it in the Model's process list. {name!r} reads "
                    "step-scoped rate caches from those processes "
                    "inside run(); without them registered earlier, "
                    "the reader sees stale or zero values. Re-order "
                    "the ``processes=`` argument to Model() so the "
                    "dependencies appear before this process."
                )
            seen.add(name)

        # Step 4 (M5 contract): invoke ``from_hotstart`` on any process
        # that defines it, passing whatever process-specific state is
        # available. By contract this MUST override the fresh-start
        # defaults set by ``init_process`` with the saved values.
        if self.__hotstart_dataset is not None:
            self.__restore_process_hotstart()

        # Step 5: build the precomputed process firing schedule.
        self.__process_schedule = self.__build_process_schedule()

        # Step 6: output store.
        self.__init_output_source()

        # Step 7 (M14 fix, review-findings 2026-05-04): validate that
        # ``wet_mask_variable`` exists in the registry now, before the
        # substep loop. Pre-fix, a typo such as ``wetted_surface_aera``
        # only surfaced as a KeyError on the first substep, deep inside
        # ``__compute_wet_mask`` (or worse, was masked by a fallback in
        # the registry). Surfacing it at init time short-circuits long
        # simulations that would otherwise fail mid-run.
        if self.__wet_mask_variable is not None:
            if self.__wet_mask_variable not in self.__registry:
                # Best-effort enumeration of what *is* registered, to help
                # the caller identify the typo. The real
                # ``VariableRegistry`` stores entries in ``_registry``;
                # the test stubs use ``_data``. Fall back to an empty
                # tuple if neither attribute is exposed.
                declared_dict = (
                    getattr(self.__registry, "_registry", None)
                    or getattr(self.__registry, "_data", None)
                    or {}
                )
                declared = sorted(declared_dict.keys()) if hasattr(
                    declared_dict, "keys"
                ) else []
                raise KeyError(
                    f"wet_mask_variable={self.__wet_mask_variable!r} not "
                    f"registered; declared variables: {declared!r}"
                )

        # Step 8 (provider coverage check): every declared process input must
        # have a provider. By now all providers have had their turn to populate
        # the registry — data sources (step 1), hotstart seed (step 2),
        # init_process incl. the riverine bridge (step 3), and from_hotstart
        # (step 4). Any ``process.variables`` entry still absent from the
        # registry has no provider and would otherwise surface as a latent
        # runtime KeyError mid-substep. Turn it into a clear init-time error.
        #
        # Known limitation: constituents read via ``registry.get_at_time(...)``
        # that are NOT declared in any ``process.variables`` (e.g.
        # ``algae_floating``) are not covered here — only declared inputs are.
        # ``depth`` IS declared (e.g. by POM/Temperature), so the riverine
        # depth wiring is covered by this check.
        required_by: dict[str, list[str]] = {}
        for process in self.__processes:
            for variable in process.variables:
                required_by.setdefault(variable, []).append(
                    process.process_name()
                )
        missing = [
            variable
            for variable in sorted(required_by)
            if variable not in self.__registry
        ]
        if missing:
            # Mirror the wet-mask check's enumeration of what *is* registered
            # to help the caller see what providers did supply. The real
            # ``VariableRegistry`` stores entries in ``_registry``; the test
            # stubs use ``_data``.
            declared_dict = (
                getattr(self.__registry, "_registry", None)
                or getattr(self.__registry, "_data", None)
                or {}
            )
            declared = sorted(declared_dict.keys()) if hasattr(
                declared_dict, "keys"
            ) else []
            details = "; ".join(
                f"{variable!r} (declared by {required_by[variable]!r})"
                for variable in missing
            )
            raise KeyError(
                f"Process input(s) have no provider: {details}. Each must be "
                "supplied by a data source, the riverine bridge, a hotstart "
                "dataset, or pre-registration before the substep loop. "
                f"Registered variables: {declared!r}."
            )

    def __finalize_model(self) -> None:
        # C2 fix (review-findings 2026-05-04): the upstream ``Process``
        # base class does not define ``finalize_process``, so calling it
        # unconditionally crashes every chunked run with
        # ``AttributeError`` after the final chunk write. The
        # corresponding ``init_process`` is a no-op default on the base;
        # ``finalize_process`` should follow the same contract but
        # doesn't (yet). Use ``getattr`` so processes that opt in by
        # defining ``finalize_process`` are honored, and processes that
        # don't are silently skipped.
        for process in self.__processes:
            finalize = getattr(process, "finalize_process", None)
            if callable(finalize):
                finalize(self, self.__registry)

    @staticmethod
    def __resolve_space_dimensions(variable) -> list[tuple[str, Any]]:
        """Resolve the (dim_name, coord_values) pairs that size an output
        variable's spatial axes for the model-output zarr store.

        Generalizes what v2 carried as a hardcoded ``water_temperature``
        special-case. Every riverine-bridged constituent self-reports its
        spatial dimension via ``DataArrayVariable(space_dimension="nface")``
        (see ``Riverine._bridge_mesh_to_registry``), so this works for any
        constituent rather than just temperature.

        Two paths, in precedence order:

        1. The variable self-reports ``space_dimension`` (one name or a
           list). Use ``space_dimension_values`` when the variable supplies
           them; otherwise derive the coordinate values from the underlying
           array (``data[dim].values`` -- which returns the default integer
           index when the dim has no explicit coordinate, as the riverine
           mesh constituents do).
        2. The variable does NOT self-report a spatial dimension (e.g. a
           ``water_temperature`` array registered by a data source without a
           ``space_dimension`` kwarg). Fall back to every non-time dimension
           of the underlying array. This preserves the v2 behavior that let
           ``water_temperature`` size the store off its ``nface`` axis,
           without a name-specific check.

        Variables with no array-like payload (e.g. scalar float providers)
        contribute no spatial dimension.
        """
        declared = variable.space_dimension
        names = (
            [declared]
            if isinstance(declared, str)
            else list(declared)
            if declared is not None
            else None
        )

        data = None
        get = getattr(variable, "get", None)
        if callable(get):
            try:
                data = get()
            except Exception:
                data = None

        def coord_values(dim: str):
            if data is not None and hasattr(data, "__getitem__"):
                try:
                    return data[dim].values
                except Exception:
                    return None
            return None

        resolved: list[tuple[str, Any]] = []
        if names is not None:
            declared_values = variable.space_dimension_values
            single = len(names) == 1
            for index, dim in enumerate(names):
                values = None
                if declared_values is not None:
                    values = declared_values if single else declared_values[index]
                if values is None:
                    values = coord_values(dim)
                if values is not None:
                    resolved.append((dim, values))
            return resolved

        # No self-reported spatial dimension: derive from the array's
        # non-time dimensions so externally-registered fields on the mesh
        # (water_temperature, etc.) still size the output store.
        dims = getattr(data, "dims", None)
        if dims is not None:
            for dim in dims:
                if dim == "time":
                    continue
                values = coord_values(dim)
                if values is not None:
                    resolved.append((dim, values))
        return resolved

    def __init_output_source(self) -> None:
        if self.__output_variables is None or len(self.__output_variables) == 0:
            return

        space_dimensions: dict[str, Any] = {}
        for variable_name in self.__output_variables:
            variable = self.__registry.get_variable(variable_name)
            for dim, values in self.__resolve_space_dimensions(variable):
                if dim not in space_dimensions:
                    space_dimensions[dim] = values

        self.__output_data_store = ChunkedZarrDataStore(
            store_path=self.__simulation_directory / "model_outputs.zarr",
            start_date=self.__start_time,
            end_date=self.__end_time,
            time_step=self.__time_step,
            variables=self.__output_variables,
            spatial_field=(
                list(space_dimensions.keys()) if len(space_dimensions) > 0 else None
            ),
            spatial_field_values=(
                list(space_dimensions.values()) if len(space_dimensions) > 0 else None
            ),
            chunk_size=self.__chunk_size,
        )

    def __load_chunk_data(self, chunk_start: datetime, chunk_end: datetime) -> None:
        for variable_name, data_source in self.__variable_data_sources.items():
            if not isinstance(data_source, ChunkedDataSource):
                continue
            data = data_source.read_chunk(chunk_start, chunk_end)
            self.__registry.register(variable_name, data)

    # ---------- precomputed process schedule (kernel-opt) ----------

    def __build_process_schedule(self) -> tuple[tuple[Process, ...], ...]:
        """Precompute, for each step index, the tuple of processes that fire.

        Firing semantic: a process fires every Nth substep starting at
        ``start_time``, where ``N = process.time_step_seconds /
        time_step_seconds``. Equivalently, at substep index ``i`` a
        process fires when ``(i * time_step_seconds) %
        process.time_step_seconds == 0``. This is timezone-independent
        and well-defined for naive ``datetime`` objects.

        C6 fix (review-findings 2026-05-04): the previous implementation
        used ``self.__start_time.timestamp()`` to seed the schedule.
        ``datetime.timestamp()`` on a naive ``datetime`` is interpreted
        in the local timezone (POSIX rule), so the same model
        configuration produced different schedules on hosts in
        different timezones for any process whose ``time_step_seconds``
        did not divide 86400. The new semantic is keyed off
        delta-seconds-from-``start_time`` and is reproducible across
        hosts.

        For the common case where ``process.time_step_seconds`` is an
        integer multiple of ``time_step_seconds`` and the user happens
        to start on a wall-clock-aligned boundary, the new semantic
        agrees with the old one substep-for-substep.

        Validation: ``process.time_step_seconds`` must be an integer
        multiple of the model's ``time_step_seconds`` (mirroring the
        C7 ``chunk_size`` validation). If not, this raises
        ``ValueError``: there is no substep grid on which a non-divisor
        cadence could fire under any deterministic semantic.
        """
        n_steps = self.__count_substeps()
        time_step_seconds = self.__time_step.total_seconds()
        # Validate cadence alignment up front (C6 / mirrors C7 chunk_size check).
        for p in self.__processes:
            interval = p.time_step_seconds
            if interval % time_step_seconds != 0:
                raise ValueError(
                    "process.time_step_seconds must be an integer multiple of "
                    "model time_step; got "
                    f"process={p.process_name()!r}, "
                    f"process.time_step_seconds={interval}, "
                    f"time_step={self.__time_step!r}"
                )
        process_intervals = [(p, p.time_step_seconds) for p in self.__processes]
        schedule: list[tuple[Process, ...]] = []
        for i in range(n_steps + 1):
            delta_seconds = i * time_step_seconds
            firing = tuple(
                p for p, interval in process_intervals if delta_seconds % interval == 0
            )
            schedule.append(firing)
        return tuple(schedule)

    def __count_substeps(self) -> int:
        """Number of model substeps from ``start_time`` (inclusive) to
        ``end_time`` (exclusive)."""
        delta_seconds = (self.__end_time - self.__start_time).total_seconds()
        return int(delta_seconds // self.__time_step.total_seconds())

    # ---------- wet-mask ----------

    def __compute_wet_mask(
        self, current_time: datetime
    ) -> xr.DataArray | None:
        """Resolve the wet-mask for ``current_time``.

        Returns ``None`` when no wet-mask is configured. Otherwise
        returns a boolean DataArray (True == wet).
        """
        if self.__wet_mask_provider is not None:
            mask = self.__wet_mask_provider(current_time, self.__registry)
            return None if mask is None else mask.astype(bool)
        if self.__wet_mask_variable is None:
            return None
        value = self.__registry.get_at_time(self.__wet_mask_variable, current_time)
        return value > self.__wet_mask_threshold

    def __apply_wet_mask(
        self,
        process: Process,
        current_time: datetime,
        wet_mask: xr.DataArray | None,
    ) -> None:
        """Write NaN into masked-dry cells for each output variable the process writes.

        No-op when ``wet_mask is None``. Variables that aren't currently
        in the registry, or whose dtype is non-floating, are skipped.

        C5 fix (review-findings 2026-05-04): the previous implementation
        masked **every** variable declared in ``process.variables``,
        which conflates inputs the process *reads* (forcings like
        ``wind_speed``, ``air_temperature``, ``solar_radiation``) with
        outputs the process *writes* (``water_temperature``,
        ``sediment_temperature``). NaN-masking forcings on dry cells
        silently corrupts the input data for the next substep / process
        / chunk. The fix: only mask variables a process declares as
        *outputs*. Processes opt in by setting a class-level
        ``output_variables: list[str]``; for backward compatibility,
        processes that don't declare it fall back to the prior behavior
        (mask everything in ``.variables``). v3 ``Temperature`` declares
        ``output_variables = ["water_temperature", "sediment_temperature"]``.
        """
        if wet_mask is None:
            return
        # Use output_variables if declared; otherwise fall back to the
        # full variables list for backward compatibility.
        output_names = getattr(process, "output_variables", None)
        if output_names is None:
            output_names = getattr(process, "variables", ()) or ()
        for variable_name in output_names:
            try:
                value = self.__registry.get_at_time(variable_name, current_time)
            except KeyError:
                # Variable not yet in registry; nothing to mask.
                continue
            if not hasattr(value, "dtype") or not np.issubdtype(value.dtype, np.floating):
                continue
            masked = xr.where(wet_mask, value, np.nan)
            self.__registry.set_at_time(variable_name, current_time, masked)

    # ---------- hotstart ----------

    def __seed_from_hotstart(self) -> None:
        """Replace registry contents with the saved-dataset slice at
        ``hotstart_timestep``. Variables already in the registry are
        overwritten; variables only in the saved dataset are added.

        M8 fix (review-findings 2026-05-04): the previous fallback
        ``next(iter(ds.dims))`` silently treated the first available
        dim as the time dim. If the saved dataset only had a spatial
        dim (e.g. ``nface``), ``isel({nface: -1})`` produced single-cell
        scalars rather than per-cell snapshots. The new contract:

        - If the dataset has a recognized time dim (``time``,
          ``time_step``, or ``datetime``), use it.
        - Otherwise, if ``hotstart_timestep`` is None, treat the dataset
          as a single-snapshot dataset and use it as-is (no slicing).
        - Otherwise, raise ``ValueError`` — the caller asked for a
          specific time slice but the dataset has no recognizable time
          axis to slice along.
        """
        ds = self.__hotstart_dataset
        slice_ts = self.__hotstart_timestep
        # Identify the time dimension of the saved dataset.
        time_dim = None
        for candidate in ("time", "time_step", "datetime"):
            if candidate in ds.dims:
                time_dim = candidate
                break

        if time_dim is None:
            # No recognizable time dim. Two cases:
            if slice_ts is not None:
                # Caller asked for a specific slice but there is no
                # time axis — fail loudly rather than silently picking
                # an arbitrary dim (which would, for an ``nface``-only
                # dataset, slice space-as-time and produce scalars).
                raise ValueError(
                    "hotstart_dataset has no recognizable time dimension "
                    "(expected one of 'time', 'time_step', 'datetime') but "
                    f"hotstart_timestep={slice_ts!r} was provided. Either "
                    "provide a dataset with a time dimension or pass "
                    "hotstart_timestep=None to use the single-snapshot "
                    "dataset as-is."
                )
            # Single-snapshot dataset: use as-is.
            sliced = ds
        else:
            if slice_ts is None:
                sliced = ds.isel({time_dim: -1})
            elif isinstance(slice_ts, int):
                sliced = ds.isel({time_dim: slice_ts})
            else:
                sliced = ds.sel({time_dim: slice_ts}, method="nearest")

        for variable_name, da in sliced.data_vars.items():
            self.__registry.register(variable_name, da)

    def __restore_process_hotstart(self) -> None:
        """Invoke ``from_hotstart`` on every process that defines it.

        The saved dataset's ``attrs`` is the substrate for per-process
        state (each process is responsible for its own attr-key
        convention; the v3 ``Temperature`` uses
        ``"temperature.skip_first_time_step"`` for example).
        """
        ds = self.__hotstart_dataset
        if ds is None:
            return
        attrs = dict(ds.attrs)
        for process in self.__processes:
            from_hotstart = getattr(process, "from_hotstart", None)
            if not callable(from_hotstart):
                continue
            from_hotstart(attrs)

    # ---------- substep loops ----------

    def __process_loop_full(self) -> None:
        current_time = self.__start_time
        step_index = 0
        while current_time < self.__end_time:
            LOGGER.info("Running timestep: %s", current_time)
            # Resolves pattern-alignment spec §10 Q1: ``Diagnostics``
            # carries the current substep index so ``clip_negative_state``
            # log records can attribute clips to a step without
            # threading ``step=`` through every Process.run signature.
            self.diagnostics.current_step = step_index
            firing = self.__process_schedule[step_index]
            wet_mask = self.__compute_wet_mask(current_time) if firing else None
            for process in firing:
                process.run(current_time, self.__registry)
                self.__apply_wet_mask(process, current_time, wet_mask)
            current_time += self.__time_step
            step_index += 1
        self.__save_output_model(self.__start_time, self.__end_time)
        # M6 fix (review-findings 2026-05-04): symmetry with the chunked
        # path. Without this, processes with finalize_process bodies
        # silently skip finalization in non-chunked mode.
        self.__finalize_model()

    def __process_loop_chunked(self) -> None:
        """Chunked substep loop. Mirrors the chunking pattern in
        ``clearwater_riverine`` (see TSM design spec §3.2 resolution).

        Resolution of v2's four TODOs:

        - "this need actual chunking logic" — the per-chunk
          input-load / substep-loop / per-chunk-output-write is the
          chunking logic; the prior placeholder comment is removed.
        - "look at riverine's code and mirror where applicable" /
          "align with riverine" — the outer loop precomputes chunk
          boundaries (as step indices, not datetimes — see C7 fix
          below) and tests membership once per substep, the same
          shape as riverine's pattern.
        - "confirm if this is necessary to write out the last chunk or if
          it will be handled in the loop above" — yes, the last partial
          chunk is written exactly once, after the substep loop exits.
          The previous unconditional post-loop write was redundant for
          integer-multiple total durations and double-wrote the last
          chunk; v3 only writes the trailing partial chunk once.

        C7 fix (review-findings 2026-05-04): the chunk-boundary check
        compared ``current_time`` (a ``datetime``) against a
        ``pd.DatetimeIndex`` of ``pd.Timestamp`` (always ns-precision,
        sometimes tz-aware). ``Timestamp.__hash__`` vs
        ``datetime.__hash__`` are not symmetric across all
        configurations, and sub-second ``time_step`` arithmetic
        accumulates floating-point drift that misses every boundary.
        Replaced with **integer step-index** comparison: precompute
        which step indices land on a chunk boundary, then check
        ``step_index in interior_chunk_step_indices``. Exact-integer,
        timezone-independent, drift-immune.
        """
        # Compute chunk boundaries as step indices, not as time values.
        # The first boundary (step_index == 0) is the start; the last
        # boundary (step_index == n_steps) coincides with end_time when
        # the total duration is an integer multiple of chunk_size and
        # is handled by the post-loop write below. Interior boundaries
        # are 1 .. n_chunks - 1.
        chunk_size_seconds = self.__chunk_size.total_seconds()
        time_step_seconds = self.__time_step.total_seconds()
        if chunk_size_seconds % time_step_seconds != 0:
            # Chunk size must be an integer multiple of the substep so
            # boundaries align with substep grid. Otherwise the chunk
            # transition would land between two substeps.
            raise ValueError(
                "chunk_size must be an integer multiple of time_step; got "
                f"chunk_size={self.__chunk_size!r}, time_step={self.__time_step!r}"
            )
        steps_per_chunk = int(round(chunk_size_seconds / time_step_seconds))
        n_steps = self.__count_substeps()
        # Interior boundaries: step indices that mark the start of each
        # new chunk after the first. Excludes 0 (== start_time) and any
        # boundary >= n_steps (== or past end_time).
        interior_chunk_step_indices = {
            i for i in range(steps_per_chunk, n_steps, steps_per_chunk)
        }

        current_chunk_start = self.__start_time
        # Initial chunk's data is already loaded by __init_model; subsequent
        # chunks are loaded on-demand at chunk transitions below.

        current_time = self.__start_time
        step_index = 0
        while current_time < self.__end_time:
            if step_index in interior_chunk_step_indices:
                # Finalize this chunk before crossing the boundary.
                self.__save_output_model(
                    start_time=current_chunk_start,
                    end_time=current_time,
                )
                # Load the next chunk's data; the minus-one-timestep on the
                # start ensures previous-step lookups still resolve at the
                # chunk boundary.
                self.__load_chunk_data(
                    current_time - self.__time_step,
                    current_time + self.__chunk_size,
                )
                current_chunk_start = current_time

            LOGGER.info("Running timestep: %s", current_time)
            # Pattern-alignment spec §10 Q1 — same step-index propagation
            # as in __process_loop_full above.
            self.diagnostics.current_step = step_index
            firing = self.__process_schedule[step_index]
            wet_mask = self.__compute_wet_mask(current_time) if firing else None
            for process in firing:
                process.run(current_time, self.__registry)
                self.__apply_wet_mask(process, current_time, wet_mask)
            current_time += self.__time_step
            step_index += 1

        # Final chunk — written exactly once for the closing partial chunk.
        self.__save_output_model(
            start_time=current_chunk_start,
            end_time=self.__end_time,
        )
        self.__finalize_model()

    def __save_output_model(self, start_time: datetime, end_time: datetime) -> None:
        if self.__output_data_store is None:
            return
        for variable_name in self.__output_variables:
            # Phase I-5 (2026-05-21): migrate to the new-style accessor.
            variable = self.__registry.get_variable(variable_name).get_data()
            self.__output_data_store.write_chunk(
                data=variable,
                parameter_name=variable.name,
                start_time=start_time,
                end_time=end_time,
            )


__all__ = ["Model"]
