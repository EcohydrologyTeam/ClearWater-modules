"""v2-helper-contract test for v3's ``init_from_config``.

v3's ``clearwater_modules_v3.config.init`` reuses two private helpers from
``clearwater_modules_v2.config.init``: ``__init_processes`` and
``__init_model_data``. v3 looks them up by exact attribute name. If v2
renames, removes, or changes the signature of either helper, v3's
production startup path will break.

Purpose of this test: when v2 changes, v3 should fail loudly here, not
silently in a production simulation. CI catches the contract drift before
the next v2 sync lands. See finding C9 in
``design/clearwater_modules_v3_review_findings.md``.
"""

from __future__ import annotations

import inspect

import pytest

from clearwater_modules_v2.config import init as _v2_init


# Exact attribute names v3 depends on. Module-level double-underscore
# names are NOT class-mangled, so these are the real attribute names on
# the v2 init module.
_INIT_PROCESSES = "__init_processes"
_INIT_MODEL_DATA = "__init_model_data"


def test_v2_init_processes_exists_and_is_callable() -> None:
    """``__init_processes`` must exist on v2 init module and be callable."""
    helper = getattr(_v2_init, _INIT_PROCESSES, None)
    assert helper is not None, (
        f"v2 config.init no longer exposes `{_INIT_PROCESSES}`; v3's "
        f"init_from_config will fail at runtime. Update v3's reuse "
        f"contract."
    )
    assert callable(helper), f"`{_INIT_PROCESSES}` must be callable"


def test_v2_init_model_data_exists_and_is_callable() -> None:
    """``__init_model_data`` must exist on v2 init module and be callable."""
    helper = getattr(_v2_init, _INIT_MODEL_DATA, None)
    assert helper is not None, (
        f"v2 config.init no longer exposes `{_INIT_MODEL_DATA}`; v3's "
        f"init_from_config will fail at runtime. Update v3's reuse "
        f"contract."
    )
    assert callable(helper), f"`{_INIT_MODEL_DATA}` must be callable"


def test_v2_init_processes_signature_matches_v3_expectation() -> None:
    """Pin parameter names of ``__init_processes`` so renames break CI.

    v3 calls this helper as
    ``__init_processes(config, variable_registry, default_time_step=time_step)``.
    A rename of any of those parameters upstream would make the v3 call
    fail at runtime; this test surfaces the rename at test time instead.
    """
    helper = getattr(_v2_init, _INIT_PROCESSES)
    sig = inspect.signature(helper)
    param_names = list(sig.parameters.keys())
    assert param_names == ["config", "variable_registry", "default_time_step"], (
        f"v2 `{_INIT_PROCESSES}` parameter names changed: got "
        f"{param_names!r}, expected "
        f"['config', 'variable_registry', 'default_time_step']. "
        f"Update v3's call site and this contract test together."
    )


def test_v2_init_model_data_signature_matches_v3_expectation() -> None:
    """Pin parameter names of ``__init_model_data`` so renames break CI.

    v3 calls this helper as
    ``__init_model_data(config=..., variables=..., start_time=...,
    end_time=..., time_step=...)``. v3 uses keyword arguments, so a
    rename of any of these names upstream would silently raise a
    ``TypeError`` at runtime. This test surfaces the rename at test time.
    """
    helper = getattr(_v2_init, _INIT_MODEL_DATA)
    sig = inspect.signature(helper)
    param_names = list(sig.parameters.keys())
    assert param_names == [
        "config",
        "variables",
        "start_time",
        "end_time",
        "time_step",
    ], (
        f"v2 `{_INIT_MODEL_DATA}` parameter names changed: got "
        f"{param_names!r}, expected "
        f"['config', 'variables', 'start_time', 'end_time', 'time_step']. "
        f"Update v3's call site and this contract test together."
    )


@pytest.mark.parametrize(
    "name",
    [_INIT_PROCESSES, _INIT_MODEL_DATA],
)
def test_no_module_level_name_mangling(name: str) -> None:
    """Confirm the assumption that module-level dunder names are not mangled.

    Python name-mangling (``_ClassName__attr``) applies only inside class
    bodies. v3 relies on this: it does ``getattr(_v2_init,
    "__init_processes")`` directly. If this assumption ever broke (for
    instance, via an import shim that mangles the names), this test
    would fail. The mangled form ``_init__init_processes`` must NOT
    exist as a separate alias.
    """
    mangled = f"_init{name}"
    assert not hasattr(_v2_init, mangled), (
        f"Unexpected mangled alias `{mangled}` on v2 config.init. v3's "
        f"helper resolution assumes the unmangled `{name}` is the "
        f"canonical attribute name."
    )
