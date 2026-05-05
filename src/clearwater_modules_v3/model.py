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

LOGGER = getLogger(__name__)


class Model:
    """v3 orchestration model. See module docstring for capabilities."""

    def __init__(
        self,
        processes: tuple[Process, ...],
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
                ``wet_mask_variable``. Cells with value strictly greater
                than the threshold are considered wet.
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
        self.__simulation_directory: Path = (
            Path(simulation_directory) if simulation_directory else Path(".")
        )

        self.__chunked_mode: bool = chunk_size is not None
        self.__chunk_size: timedelta | None = chunk_size

        self.__wet_mask_variable: str | None = wet_mask_variable
        self.__wet_mask_threshold: float = wet_mask_threshold
        self.__wet_mask_provider = wet_mask_provider
        self.__hotstart_dataset: xr.Dataset | None = hotstart_dataset
        self.__hotstart_timestep = hotstart_timestep

        self.__init_complete: bool = False
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
        if not self.__init_complete:
            self.__init_model()
            self.__init_complete = True
        if self.__chunked_mode:
            self.__process_loop_chunked()
        else:
            self.__process_loop_full()

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
        # Step 1: load every data source (or first chunk) into the registry.
        for variable_name, data_source in self.__variable_data_sources.items():
            if isinstance(data_source, ChunkedDataSource):
                data = data_source.read_chunk(
                    self.__start_time,
                    self.__start_time + (self.__chunk_size or self.__time_step),
                )
            else:
                data = data_source.read(variable_name)
            self.__registry.register(variable_name, data)

        # Step 2: optionally seed the registry from the hotstart dataset.
        # This overwrites whatever was loaded in step 1 for variables that
        # appear in the saved dataset.
        if self.__hotstart_dataset is not None:
            self.__seed_from_hotstart()

        # Step 3: process-level initialization. Processes see the
        # post-hotstart registry state.
        for process in self.__processes:
            process.init_process(self, self.__registry)

        # Step 4: invoke from_hotstart on any process that defines it,
        # passing whatever process-specific state is available.
        if self.__hotstart_dataset is not None:
            self.__restore_process_hotstart()

        # Step 5: build the precomputed process firing schedule.
        self.__process_schedule = self.__build_process_schedule()

        # Step 6: output store.
        self.__init_output_source()

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

    def __init_output_source(self) -> None:
        if self.__output_variables is None or len(self.__output_variables) == 0:
            return

        space_dimensions: dict[str, Any] = {}
        for variable_name in self.__output_variables:
            variable = self.__registry.get_variable(variable_name)
            if (
                variable.space_dimension is not None
                and variable.space_dimension not in space_dimensions
            ):
                space_dimensions[variable.space_dimension] = (
                    variable.space_dimension_values
                )

            # v2 carried this manual override for water_temperature so the
            # output store sees the riverine mesh's ``nface`` dimension.
            # Preserved verbatim until variable.space_dimension is wired
            # up cleanly in riverine.
            if variable_name == "water_temperature":
                data = variable.get()
                space_dimensions["nface"] = data["nface"].values

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

        Replaces the per-step modulo check
        ``current_time_seconds % process.time_step_seconds == 0`` with an
        O(1) tuple lookup. Same firing semantics: a process fires whenever
        the absolute time-since-epoch (seconds) is divisible by its
        ``time_step_seconds``. We enumerate the model's grid of substep
        times once and bind each step to its firing-process tuple.
        """
        n_steps = self.__count_substeps()
        schedule: list[tuple[Process, ...]] = []
        time_step_seconds = self.__time_step.total_seconds()
        process_intervals = [(p, p.time_step_seconds) for p in self.__processes]
        start_seconds = self.__start_time.timestamp()
        for i in range(n_steps + 1):
            current_seconds = start_seconds + i * time_step_seconds
            firing = tuple(
                p for p, interval in process_intervals if current_seconds % interval == 0
            )
            schedule.append(firing)
        return tuple(schedule)

    def __count_substeps(self) -> int:
        """Number of model substeps from ``start_time`` (inclusive) to
        ``end_time`` (exclusive)."""
        delta_seconds = (self.__end_time - self.__start_time).total_seconds()
        return int(delta_seconds // self.__time_step.total_seconds())

    def __step_index(self, current_time: datetime) -> int:
        """Map a substep ``current_time`` to its index in ``__process_schedule``."""
        delta_seconds = (current_time - self.__start_time).total_seconds()
        return int(round(delta_seconds / self.__time_step.total_seconds()))

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
        """Write NaN into masked-dry cells for each variable the process declares.

        No-op when ``wet_mask is None``. Variables that aren't currently
        in the registry, or whose dtype is non-floating, are skipped.
        """
        if wet_mask is None:
            return
        for variable_name in getattr(process, "variables", ()) or ():
            try:
                value = self.__registry.get_at_time(variable_name, current_time)
            except Exception:
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
        """
        ds = self.__hotstart_dataset
        slice_ts = self.__hotstart_timestep
        # Identify the time dimension of the saved dataset. Fall back to
        # "time" if the dataset doesn't expose a clear single time dim.
        time_dim = None
        for candidate in ("time", "time_step", "datetime"):
            if candidate in ds.dims:
                time_dim = candidate
                break
        if time_dim is None and ds.dims:
            time_dim = next(iter(ds.dims))

        if slice_ts is None:
            sliced = ds.isel({time_dim: -1}) if time_dim else ds
        elif isinstance(slice_ts, int):
            sliced = ds.isel({time_dim: slice_ts}) if time_dim else ds
        else:
            sliced = ds.sel({time_dim: slice_ts}, method="nearest") if time_dim else ds

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
            variable = self.__registry.get(variable_name)
            self.__output_data_store.write_chunk(
                data=variable,
                parameter_name=variable.name,
                start_time=start_time,
                end_time=end_time,
            )


__all__ = ["Model"]
