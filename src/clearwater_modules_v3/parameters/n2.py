"""v3 NSM1 dissolved nitrogen gas (N2) parameter defaults.

Consumed by: ``N2`` Process.
Source: v1 ``clearwater_modules/nsm1/constants.py`` ``N2StaticVariables`` /
``DEFAULT_N2``.
Corrections applied: none.

v1's ``N2StaticVariables`` is empty: N2 atmospheric exchange is computed via
Henry's law from ``pressure_mb`` and water temperature, not from N2-specific
parameters. This module is provided for symmetry; populate it if v3 adds
N2-specific kinetics.
"""

DEFAULTS: dict[str, float | int | bool] = {}
