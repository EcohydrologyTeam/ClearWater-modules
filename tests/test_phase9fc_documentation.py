"""Phase 9.F.C documentation-fix regression tests.

Pins three documentation-only outcomes from Phase 9.F.C:

1. ``h2`` FIXME cleared from ``parameters/pom.py`` (Section 2.5).
2. ``ksbod_20`` FIXME cleared from ``parameters/cbod.py`` (Section 2.3).
3. SedFlux defensive ``NotImplementedError`` guards in
   ``Nitrogen.__init__`` and ``Phosphorus.__init__`` (Section 2.1).

The Phase 9.F.C resolution markers in the corrections doc itself are
also pinned here so that a future edit cannot quietly downgrade a
RESOLVED section back to "under review".

See ``parameter_defaults_corrections.md`` Sections 2.1, 2.3, 2.5.
"""
from datetime import timedelta
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src" / "clearwater_modules_v3"
CORRECTIONS_DOC = SRC / "parameter_defaults_corrections.md"


# ---------------------------------------------------------------------------
# Section 2.5: h2 FIXME cleared, Di Toro / QUAL2K provenance recorded
# ---------------------------------------------------------------------------

def test_h2_fixme_cleared():
    """Section 2.5: h2 line in parameters/pom.py no longer carries a
    FIXME; inline comment now records the Di Toro / QUAL2K H_2
    provenance."""
    pom_text = (SRC / "parameters" / "pom.py").read_text()

    # Locate the h2 entry line.
    h2_line = next(
        line for line in pom_text.splitlines()
        if line.lstrip().startswith("'h2'")
    )
    assert "FIXME" not in h2_line, (
        f"h2 line should no longer carry a FIXME tag (Phase 9.F.C): "
        f"{h2_line!r}"
    )
    assert "active sediment layer thickness" in h2_line, (
        f"h2 inline comment should record the canonical role: {h2_line!r}"
    )
    assert "Di Toro" in h2_line or "QUAL2K" in h2_line, (
        f"h2 inline comment should cite Di Toro/QUAL2K: {h2_line!r}"
    )

    # Module docstring should also describe the role and cite the
    # Phase 9.F.C correction.
    assert "active sediment layer thickness" in pom_text
    assert "Phase 9.F.C" in pom_text


def test_pom_process_docstring_explains_h2_role():
    """Section 2.5: processes/pom.py module docstring gained a
    Phase 9.F.C "Conceptual note" explaining POM = bed-sediment POM
    (Fortran POM2) and h2 = Di Toro H_2 anaerobic layer thickness."""
    process_text = (SRC / "processes" / "pom.py").read_text()
    assert "Phase 9.F.C" in process_text
    assert "POM2" in process_text, (
        "POM Process docstring should reference Fortran POM2 / Di Toro "
        "layer 2 to clarify the bed-sediment identity"
    )
    assert "Di Toro" in process_text and "QUAL2K" in process_text


# ---------------------------------------------------------------------------
# Section 2.3: ksbod_20 FIXME cleared, citation block recorded
# ---------------------------------------------------------------------------

def test_ksbod_fixme_cleared():
    """Section 2.3: ksbod_20 line in parameters/cbod.py no longer
    carries a FIXME; inline comment now records the QUAL2K /
    Brown & Barnwell / EPA-TMDL citation block."""
    cbod_text = (SRC / "parameters" / "cbod.py").read_text()

    ksbod_line = next(
        line for line in cbod_text.splitlines()
        if line.lstrip().startswith("'ksbod_20'")
    )
    assert "FIXME" not in ksbod_line, (
        f"ksbod_20 line should no longer carry a FIXME tag (Phase 9.F.C): "
        f"{ksbod_line!r}"
    )

    # The module docstring should now carry the QUAL2K / Brown & Barnwell
    # / EPA TMDL citation block.
    assert "Phase 9.F.C" in cbod_text
    assert "QUAL2K" in cbod_text
    assert "QUAL2E" in cbod_text or "Brown" in cbod_text
    assert "EPA TMDL" in cbod_text or "Book II" in cbod_text


# ---------------------------------------------------------------------------
# Section 2.1: SedFlux defensive guards
# ---------------------------------------------------------------------------

def test_nitrogen_use_sedflux_true_raises_notimplementederror():
    """Section 2.1 / Option B: Nitrogen(parameters={'use_SedFlux': True})
    raises NotImplementedError to prevent silent partial behavior."""
    from clearwater_modules_v2.processes.nitrogen import Nitrogen

    with pytest.raises(NotImplementedError) as excinfo:
        Nitrogen(parameters={"use_SedFlux": True})
    msg = str(excinfo.value)
    assert "use_SedFlux" in msg
    assert "NSM2" in msg or "Section 2.1" in msg


def test_nitrogen_use_sedflux_false_constructs_normally():
    """Section 2.1: Nitrogen(parameters={'use_SedFlux': False}) and
    Nitrogen() with no params should both construct without raising
    (the guard fires only on the True opt-in)."""
    from clearwater_modules_v2.processes.nitrogen import Nitrogen

    Nitrogen(parameters={"use_SedFlux": False})
    Nitrogen()


def test_phosphorus_use_sedflux_true_raises_notimplementederror():
    """Section 2.1 / Option B: Phosphorus(parameters={'use_SedFlux': True})
    raises NotImplementedError to prevent silent partial behavior."""
    from clearwater_modules_v3.processes.phosphorus import Phosphorus

    with pytest.raises(NotImplementedError) as excinfo:
        Phosphorus(parameters={"use_SedFlux": True})
    msg = str(excinfo.value)
    assert "use_SedFlux" in msg
    assert "NSM2" in msg or "Section 2.1" in msg


def test_phosphorus_use_sedflux_false_constructs_normally():
    """Section 2.1: Phosphorus(parameters={'use_SedFlux': False}) and
    Phosphorus() with no params should both construct without raising."""
    from clearwater_modules_v3.processes.phosphorus import Phosphorus

    Phosphorus(parameters={"use_SedFlux": False})
    Phosphorus()


# ---------------------------------------------------------------------------
# Corrections-doc resolution markers
# ---------------------------------------------------------------------------

def test_section_2_1_marked_resolved_in_phase9fc():
    """Pin Section 2.1 RESOLVED status so a future edit cannot quietly
    downgrade it."""
    corrections_text = CORRECTIONS_DOC.read_text()
    section_start = corrections_text.find("### 2.1 ")
    assert section_start != -1
    section_end = corrections_text.find("\n### ", section_start + 1)
    section_text = corrections_text[section_start:section_end]
    assert "RESOLVED in Phase 9.F.C" in section_text, (
        "Section 2.1 should be marked RESOLVED in Phase 9.F.C"
    )
    assert "NotImplementedError" in section_text, (
        "Section 2.1 should document the Phase 9.F.C defensive guard"
    )


def test_section_2_3_marked_resolved_in_phase9fc():
    corrections_text = CORRECTIONS_DOC.read_text()
    section_start = corrections_text.find("### 2.3 ")
    assert section_start != -1
    section_end = corrections_text.find("\n### ", section_start + 1)
    section_text = corrections_text[section_start:section_end]
    assert "RESOLVED in Phase 9.F.C" in section_text


def test_section_2_5_marked_resolved_in_phase9fc():
    corrections_text = CORRECTIONS_DOC.read_text()
    section_start = corrections_text.find("### 2.5 ")
    assert section_start != -1
    section_end = corrections_text.find("\n### ", section_start + 1)
    section_text = corrections_text[section_start:section_end]
    assert "RESOLVED in Phase 9.F.C" in section_text


def test_section_2_header_records_resolution_count():
    """Section 2 header should record the 7 RESOLVED / 1 deferred state
    after Phase 9.F.C."""
    corrections_text = CORRECTIONS_DOC.read_text()
    assert "7 RESOLVED" in corrections_text
