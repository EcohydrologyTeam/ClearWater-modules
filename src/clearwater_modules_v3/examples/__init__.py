"""v3 NSM1 example helpers.

Phase 8.A: a small set of helper modules that construct fully-runnable v3
NSM1 Models without going through the YAML configuration path. The YAML
``init_from_file`` entry point requires a ``data_sources`` block that
maps every process variable to a real data provider (typically a CSV
file or a Riverine mesh). For headless demos, validation runs, and the
Phase 7.C example notebook, that requirement is heavyweight: there is
no Riverine mesh, the time horizon is short, and the inputs are
synthetic.

The helpers in this package construct the same logical Model that
``init_from_file`` would build from ``config/nsm1_default.yml``, but
without the CSV/zarr round-trip. They are the recommended
"NSM1-only, no transport, no I/O" entry point for users who just want
to exercise the kinetics.
"""

from clearwater_modules_v3.examples.nsm1_demo_setup import (
    InMemoryRegistry,
    build_nsm1_demo,
    default_initial_conditions,
    default_process_parameters,
)

__all__ = [
    "InMemoryRegistry",
    "build_nsm1_demo",
    "default_initial_conditions",
    "default_process_parameters",
]
