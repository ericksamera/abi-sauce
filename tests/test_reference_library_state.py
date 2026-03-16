from __future__ import annotations

from abi_sauce.reference_library_state import (
    ReferenceLibraryState,
    clear_reference_library_state,
    get_stored_reference,
    list_stored_references,
    read_reference_library_state,
    store_reference,
)


def test_read_reference_library_state_defaults_to_empty() -> None:
    session_state: dict[str, object] = {}

    reference_library_state = read_reference_library_state(session_state)

    assert reference_library_state == ReferenceLibraryState()


def test_store_reference_round_trip_and_lookup() -> None:
    session_state: dict[str, object] = {}

    stored_reference = store_reference(
        session_state,
        name="amplicon_001",
        reference_text=">amplicon_001\nAACCGGTT\n",
    )

    reference_library_state = read_reference_library_state(session_state)

    assert tuple(reference_library_state.references_by_id) == (
        stored_reference.reference_id,
    )
    assert (
        get_stored_reference(
            session_state,
            reference_id=stored_reference.reference_id,
        )
        == stored_reference
    )
    assert list_stored_references(session_state) == (stored_reference,)


def test_store_reference_dedupes_exact_name_and_text_matches() -> None:
    session_state: dict[str, object] = {}

    first_reference = store_reference(
        session_state,
        name="amplicon_001",
        reference_text=">amplicon_001\nAACCGGTT\n",
    )
    second_reference = store_reference(
        session_state,
        name="amplicon_001",
        reference_text=">amplicon_001\nAACCGGTT\n",
    )

    assert second_reference == first_reference
    assert list_stored_references(session_state) == (first_reference,)


def test_clear_reference_library_state_removes_session_values() -> None:
    session_state: dict[str, object] = {}

    store_reference(
        session_state,
        name="amplicon_001",
        reference_text=">amplicon_001\nAACCGGTT\n",
    )
    clear_reference_library_state(session_state)

    assert read_reference_library_state(session_state) == ReferenceLibraryState()
