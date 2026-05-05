"""v3 NSM1 parameter library.

Each submodule defines a single ``DEFAULTS`` dict containing the parameter
name -> default value mapping for one parameter group, ported from v1
``clearwater_modules/nsm1/constants.py`` TypedDicts. Seven critical default-
value corrections are applied at the port (sentinel-999 settling/SOD/reaeration
overrides plus the ``pressure_mb`` magnitude error); see
``parameter_defaults_corrections.md`` for the full list and rationale.

Each ``DEFAULTS`` dict is re-exported here under a ``<GROUP>_DEFAULTS`` alias
for convenient import by Process classes::

    from clearwater_modules_v3.parameters import NITROGEN_DEFAULTS

Process classes typically merge user-provided values over these defaults at
construction (Phase 1.3 work, separate from this parameter library).
"""

from clearwater_modules_v3.parameters.algae import DEFAULTS as ALGAE_DEFAULTS
from clearwater_modules_v3.parameters.alkalinity import DEFAULTS as ALKALINITY_DEFAULTS
from clearwater_modules_v3.parameters.balgae import DEFAULTS as BALGAE_DEFAULTS
from clearwater_modules_v3.parameters.carbon import DEFAULTS as CARBON_DEFAULTS
from clearwater_modules_v3.parameters.cbod import DEFAULTS as CBOD_DEFAULTS
from clearwater_modules_v3.parameters.dox import DEFAULTS as DOX_DEFAULTS
from clearwater_modules_v3.parameters.global_parameters import DEFAULTS as GLOBAL_PARAMETERS_DEFAULTS
from clearwater_modules_v3.parameters.global_vars import DEFAULTS as GLOBAL_VARS_DEFAULTS
from clearwater_modules_v3.parameters.n2 import DEFAULTS as N2_DEFAULTS
from clearwater_modules_v3.parameters.nitrogen import DEFAULTS as NITROGEN_DEFAULTS
from clearwater_modules_v3.parameters.pathogen import DEFAULTS as PATHOGEN_DEFAULTS
from clearwater_modules_v3.parameters.phosphorus import DEFAULTS as PHOSPHORUS_DEFAULTS
from clearwater_modules_v3.parameters.pom import DEFAULTS as POM_DEFAULTS

__all__ = [
    "ALGAE_DEFAULTS",
    "ALKALINITY_DEFAULTS",
    "BALGAE_DEFAULTS",
    "CARBON_DEFAULTS",
    "CBOD_DEFAULTS",
    "DOX_DEFAULTS",
    "GLOBAL_PARAMETERS_DEFAULTS",
    "GLOBAL_VARS_DEFAULTS",
    "N2_DEFAULTS",
    "NITROGEN_DEFAULTS",
    "PATHOGEN_DEFAULTS",
    "PHOSPHORUS_DEFAULTS",
    "POM_DEFAULTS",
]
