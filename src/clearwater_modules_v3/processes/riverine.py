from __future__ import annotations

from clearwater_modules_v3.processes.base import Process, ProcessFactory
from datetime import datetime, timedelta
from clearwater_data.variables import VariableRegistry, DataArrayVariable
from pathlib import Path
from typing import TYPE_CHECKING

# clearwater_riverine is imported LAZILY (inside the methods that build or use a
# mesh) so importing the universal APL processes package does NOT require CWR to
# be installed. CWR is the riverine (HEC-RAS-2D) transport engine; only
# consumers that actually run the Riverine process need it. The import here is
# under TYPE_CHECKING purely so the ``cwr.*`` annotation below resolves for type
# checkers (with ``from __future__ import annotations`` it is never evaluated at
# runtime).
if TYPE_CHECKING:
    import clearwater_riverine as cwr
    from clearwater_modules_v3.model import Model


class Riverine(Process):
    """
    Riverine process.
    """

    variables = []

    # Mesh-constituent (CW-Riverine fork name) -> v3 canonical registry name.
    _MESH_TO_CANONICAL = {
        "Ap": "algae_floating",
        "Ab": "benthic_algae",
        "NH4": "ammonium",
        "NO3": "nitrate",
        "OrgN": "organic_nitrogen",
        "N2": "n2",
        "TIP": "tip",  # v3 convention; was phosphorus_total_inorganic
        "OrgP": "organic_phosphorus",
        "POC": "poc",
        "DOC": "doc",
        "DIC": "dic",
        "CBOD": "cbod",
        "POM": "pom",
        "DOX": "oxygen_dissolved",
        "Alk": "alkalinity",
        "PX": "pathogen",
    }

    def __init__(
        self,
        riverine_instance: cwr.ClearwaterRiverine,
        time_step: timedelta = timedelta(seconds=30),
    ) -> None:
        self.riverine_instance = riverine_instance
        Process.__init__(self, time_step)

    @staticmethod
    def from_file_path(
        configuration_path: str | Path,
        variable_registry: VariableRegistry,
        start_datetime: str | datetime,
        end_datetime: str | datetime,
        time_step: timedelta = timedelta(seconds=30),
    ) -> "Riverine":
        # TODO: This will be removed once Riverine is updated to use datetime objects
        # or the pandas.to_datetime function is used to convert the start and end datetimes
        if isinstance(start_datetime, datetime):
            start_datetime = start_datetime.strftime("%m-%d-%y %H:%M:%S")
        if isinstance(end_datetime, datetime):
            end_datetime = end_datetime.strftime("%m-%d-%y %H:%M:%S")

        import clearwater_riverine as cwr  # lazy: only riverine consumers need CWR

        return Riverine(
            cwr.ClearwaterRiverine(
                config_filepath=configuration_path,
                start_datetime=start_datetime,
                end_datetime=end_datetime,
                variable_registry=variable_registry,
            ),
            time_step,
        )

    @ProcessFactory.register("riverine")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "Riverine":
        return Riverine.from_file_path(**config, variable_registry=variable_registry)

    def _bridge_mesh_to_registry(self, registry: VariableRegistry) -> None:
        """(Re)point canonical registry names at the current mesh objects.

        Idempotent and cheap: rebinds references (``copy(deep=False)``
        shares buffers). Safe to call every substep; required after
        ``update()`` so a chunk reload's freshly re-registered DataArrays
        are picked up rather than left stranded on the previous chunk's
        buffers.

        The mesh is a ``MeshView`` exposing constituents by item access
        (``mesh["Ap"]``) and membership via ``name in mesh``.
        """
        mesh = self.riverine_instance.mesh
        for mesh_name, canonical in self._MESH_TO_CANONICAL.items():
            if mesh_name in mesh:
                registry.register(
                    canonical,
                    DataArrayVariable(mesh[mesh_name].copy(deep=False)),
                    overwrite=True,  # re-bridge: upsert, not first-insert
                )
        # depth: the cell mean water-column depth, resolved on the riverine
        # side by precedence (RAS Cell Hydraulic Depth -> volume / wetted_
        # surface_area -> WSE - bed) and exposed on demand under
        # 'coupling_depth' once enable_coupling_depth() has been called.
        # Chunk-safe: refreshed per chunk on the riverine side. Bridge it
        # like a constituent.
        if "coupling_depth" not in mesh:
            raise KeyError(
                "Riverine coupling requires 'coupling_depth' (the resolved "
                "cell mean water-column depth) from the transport mesh. Ensure "
                "the ClearWater-Riverine model is coupling-enabled "
                "(enable_coupling_depth() must have been called) or depth could "
                "not be resolved; see design/clearwater_modules_v3_riverine_"
                "process_meshview_compat.md."
            )
        registry.register(
            "depth",
            DataArrayVariable(mesh["coupling_depth"].copy(deep=False)),
            overwrite=True,
        )

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """
        Initialize the riverine process.

        ``model`` is retained for interface compatibility but no longer read.
        """

        def check_variable_in_registry(variable_name: str) -> None:
            if variable_name not in registry:
                raise KeyError(
                    f"Variable '{variable_name}' not found in registry. Did models initialize correctly."
                )

        # register the water temperature, volume, and surface area to the registry
        # verify required varables have initialized correctly
        check_variable_in_registry("water_temperature")
        check_variable_in_registry("volume")
        try:
            check_variable_in_registry("wetted_surface_area")
        except KeyError:
            # Temperature reads wetted_surface_area from the registry. The
            # riverine mesh does not auto-register it, so derive it from the
            # elevation-volume lookups and register it on the shared registry.
            # calculate_wetted_surface_area takes the registry (not the model
            # instance) and RETURNS a DataArrayVariable without registering it.
            import clearwater_riverine.utilities as cwr_utils  # lazy import

            registry.register(
                "wetted_surface_area",
                cwr_utils.calculate_wetted_surface_area(registry),
            )

        # Turn on resolved-depth computation for this coupled run. Idempotent:
        # enables the flag, seed-computes the resolved depth for the current
        # window, and registers it under 'coupling_depth' (refreshed per chunk
        # on the riverine side). Standalone runs never call this. Must precede
        # the first _bridge_mesh_to_registry call below, which reads
        # mesh["coupling_depth"].
        self.riverine_instance.enable_coupling_depth()

        # Seed the canonical names at t0 (for substep 0, where run() skips
        # the first riverine update). The same chunk-safe re-bridge runs in
        # run() after each update() so the names track the current chunk.
        self._bridge_mesh_to_registry(registry)

        # The riverine model use current time_step as the start point and
        # updates the model at the time_step + delta time.
        # The rest of modules uses a definition of time_step as the time
        # to be updated. This boolean allows us to skip the first time_step
        self.__skip_first_time_step = True

    def finalize_process(self, model, registry) -> None:
        self.riverine_instance.finalize()

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """
        Run the riverine process.
        """
        if self.__skip_first_time_step:
            self.__skip_first_time_step = False
            return

        # run the next time step
        self.riverine_instance.update()

        # Chunk-safe re-bridge: re-point the canonical registry names at the
        # current mesh objects. In chunked mode, update() may load a new
        # chunk, re-registering FRESH DataArrays for the constituents and
        # average_depth; re-bridging here picks them up rather than leaving
        # the canonical names stranded on the previous chunk's buffers.
        self._bridge_mesh_to_registry(registry)
