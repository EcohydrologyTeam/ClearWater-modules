"""Tests that everything can be imported as we expect."""
import pytest
import clearwater_modules_python

@pytest.fixture
def sub_modules() -> list[str]:
    """Return a list of sub-modules."""
    return [
        'tsm',
    ]

def test_sub_modules(sub_modules) -> None:
    """Test that all sub-modules can be imported."""
    for sub_module in sub_modules:
        assert hasattr(clearwater_modules_python, sub_module)
