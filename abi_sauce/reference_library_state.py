from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Final, TypeAlias
from uuid import uuid4

SessionStateReader: TypeAlias = Any
SessionStateWriter: TypeAlias = Any

_REFERENCES_BY_ID_SESSION_KEY: Final[str] = (
    "abi_sauce.reference_library.references_by_id"
)


@dataclass(frozen=True, slots=True)
class StoredReference:
    """One reusable reference sequence stored in session state."""

    reference_id: str
    name: str
    reference_text: str


@dataclass(frozen=True, slots=True)
class ReferenceLibraryState:
    """Session-scoped reusable reference library."""

    references_by_id: dict[str, StoredReference] = field(default_factory=dict)


def read_reference_library_state(
    session_state: SessionStateReader,
) -> ReferenceLibraryState:
    """Read the reusable reference library from session storage."""
    return ReferenceLibraryState(
        references_by_id=_coerce_references_by_id(
            session_state.get(_REFERENCES_BY_ID_SESSION_KEY)
        )
    )


def write_reference_library_state(
    session_state: SessionStateWriter,
    reference_library_state: ReferenceLibraryState,
) -> None:
    """Persist the complete reusable reference library into session storage."""
    session_state[_REFERENCES_BY_ID_SESSION_KEY] = dict(
        reference_library_state.references_by_id
    )


def list_stored_references(
    session_state: SessionStateReader,
) -> tuple[StoredReference, ...]:
    """Return stored references in a stable display order."""
    return tuple(
        sorted(
            read_reference_library_state(session_state).references_by_id.values(),
            key=lambda reference: (reference.name.lower(), reference.reference_id),
        )
    )


def get_stored_reference(
    session_state: SessionStateReader,
    *,
    reference_id: str,
) -> StoredReference | None:
    """Return one stored reference by id, if it exists."""
    return read_reference_library_state(session_state).references_by_id.get(
        reference_id
    )


def store_reference(
    session_state: SessionStateWriter,
    *,
    name: str,
    reference_text: str,
) -> StoredReference:
    """Store one reusable reference in session state, deduping exact duplicates."""
    normalized_name = _normalized_required_text(name, field_name="name")
    normalized_reference_text = _normalized_required_text(
        reference_text,
        field_name="reference_text",
    )
    current_state = read_reference_library_state(session_state)

    for stored_reference in current_state.references_by_id.values():
        if (
            stored_reference.name == normalized_name
            and stored_reference.reference_text == normalized_reference_text
        ):
            return stored_reference

    stored_reference = StoredReference(
        reference_id=f"reference-{uuid4().hex[:8]}",
        name=normalized_name,
        reference_text=normalized_reference_text,
    )
    references_by_id = dict(current_state.references_by_id)
    references_by_id[stored_reference.reference_id] = stored_reference
    write_reference_library_state(
        session_state,
        ReferenceLibraryState(references_by_id=references_by_id),
    )
    return stored_reference


def clear_reference_library_state(session_state: SessionStateWriter) -> None:
    """Remove the reusable reference library from session storage."""
    session_state.pop(_REFERENCES_BY_ID_SESSION_KEY, None)


def _normalized_required_text(value: str, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be a string")
    stripped_value = value.strip()
    if not stripped_value:
        raise ValueError(f"{field_name} must not be empty")
    return stripped_value


def _coerce_references_by_id(value: object) -> dict[str, StoredReference]:
    if not isinstance(value, dict):
        return {}
    return {
        reference_id: stored_reference
        for reference_id, stored_reference in value.items()
        if isinstance(reference_id, str)
        and isinstance(stored_reference, StoredReference)
    }
