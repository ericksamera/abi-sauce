from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final, Literal, TypeAlias
from uuid import uuid4

from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.reference_alignment_types import StrandPolicy
from abi_sauce.services.batch_parse import BatchSignature

AlignmentEngineKind = Literal["pairwise", "reference_single", "reference_multi"]

SessionStateReader: TypeAlias = Any
SessionStateWriter: TypeAlias = Any

_ALIGNMENT_BATCH_SIGNATURE_SESSION_KEY: Final[str] = (
    "abi_sauce.alignments.batch_signature"
)
_ALIGNMENTS_BY_ID_SESSION_KEY: Final[str] = "abi_sauce.alignments.alignments_by_id"
_SELECTED_ALIGNMENT_ID_SESSION_KEY: Final[str] = "abi_sauce.alignments.selected_id"


@dataclass(frozen=True, slots=True)
class AlignmentDefinition:
    """One saved alignment workspace definition."""

    alignment_id: str
    name: str
    source_filenames: tuple[str, ...]
    engine_kind: AlignmentEngineKind = "pairwise"
    assembly_config: AssemblyConfig = AssemblyConfig()
    reference_name: str | None = None
    reference_text: str | None = None
    strand_policy: StrandPolicy = "auto"


@dataclass(frozen=True, slots=True)
class AlignmentSessionState:
    """Session-scoped saved alignment definitions for the active batch."""

    batch_signature: BatchSignature | None = None
    alignments_by_id: dict[str, AlignmentDefinition] = field(default_factory=dict)
    selected_alignment_id: str | None = None


def read_alignment_session_state(
    session_state: SessionStateReader,
) -> AlignmentSessionState:
    """Read saved alignment definitions from session state."""
    return AlignmentSessionState(
        batch_signature=_coerce_batch_signature(
            session_state.get(_ALIGNMENT_BATCH_SIGNATURE_SESSION_KEY)
        ),
        alignments_by_id=_coerce_alignments_by_id(
            session_state.get(_ALIGNMENTS_BY_ID_SESSION_KEY)
        ),
        selected_alignment_id=_coerce_selected_alignment_id(
            session_state.get(_SELECTED_ALIGNMENT_ID_SESSION_KEY)
        ),
    )


def sync_alignment_session_state(
    session_state: SessionStateWriter,
    *,
    batch_signature: BatchSignature,
    parsed_record_names: Iterable[str],
) -> AlignmentSessionState:
    """Sanitize saved alignments for the currently active parsed batch."""
    record_names = frozenset(parsed_record_names)
    current_state = read_alignment_session_state(session_state)

    sanitized_alignments_by_id = {
        alignment_id: definition
        for alignment_id, definition in current_state.alignments_by_id.items()
        if definition.source_filenames
        and all(
            source_filename in record_names
            for source_filename in definition.source_filenames
        )
    }
    selected_alignment_id = (
        current_state.selected_alignment_id
        if current_state.selected_alignment_id in sanitized_alignments_by_id
        else next(iter(sanitized_alignments_by_id), None)
    )

    next_state = AlignmentSessionState(
        batch_signature=batch_signature,
        alignments_by_id=sanitized_alignments_by_id,
        selected_alignment_id=selected_alignment_id,
    )
    write_alignment_session_state(session_state, next_state)
    return next_state


def write_alignment_session_state(
    session_state: SessionStateWriter,
    alignment_state: AlignmentSessionState,
) -> None:
    """Persist the complete saved-alignment state to session storage."""
    session_state[_ALIGNMENT_BATCH_SIGNATURE_SESSION_KEY] = (
        alignment_state.batch_signature
    )
    session_state[_ALIGNMENTS_BY_ID_SESSION_KEY] = dict(
        alignment_state.alignments_by_id
    )
    session_state[_SELECTED_ALIGNMENT_ID_SESSION_KEY] = (
        alignment_state.selected_alignment_id
    )


def create_alignment_definition(
    session_state: SessionStateWriter,
    *,
    name: str,
    source_filenames: Iterable[str],
    engine_kind: AlignmentEngineKind,
    assembly_config: AssemblyConfig = AssemblyConfig(),
    reference_name: str | None = None,
    reference_text: str | None = None,
    strand_policy: StrandPolicy = "auto",
) -> AlignmentDefinition:
    """Create and persist one new saved alignment definition."""
    current_state = read_alignment_session_state(session_state)
    source_filenames_tuple = tuple(source_filenames)
    definition = AlignmentDefinition(
        alignment_id=f"alignment-{uuid4().hex[:8]}",
        name=_normalized_alignment_name(
            name,
            source_filenames_tuple,
            engine_kind=engine_kind,
            reference_name=reference_name,
            existing_names=(
                definition.name
                for definition in current_state.alignments_by_id.values()
            ),
        ),
        source_filenames=source_filenames_tuple,
        engine_kind=engine_kind,
        assembly_config=assembly_config,
        reference_name=_normalized_optional_text(reference_name),
        reference_text=_normalized_optional_text(reference_text),
        strand_policy=strand_policy,
    )
    alignments_by_id = dict(current_state.alignments_by_id)
    alignments_by_id[definition.alignment_id] = definition
    write_alignment_session_state(
        session_state,
        AlignmentSessionState(
            batch_signature=current_state.batch_signature,
            alignments_by_id=alignments_by_id,
            selected_alignment_id=definition.alignment_id,
        ),
    )
    return definition


def update_alignment_definition(
    session_state: SessionStateWriter,
    *,
    alignment_id: str,
    name: str,
    source_filenames: Iterable[str],
    engine_kind: AlignmentEngineKind,
    assembly_config: AssemblyConfig = AssemblyConfig(),
    reference_name: str | None = None,
    reference_text: str | None = None,
    strand_policy: StrandPolicy = "auto",
) -> AlignmentDefinition:
    """Update and persist one existing saved alignment definition."""
    current_state = read_alignment_session_state(session_state)
    if alignment_id not in current_state.alignments_by_id:
        raise KeyError(alignment_id)

    source_filenames_tuple = tuple(source_filenames)
    existing_names = [
        definition.name
        for current_id, definition in current_state.alignments_by_id.items()
        if current_id != alignment_id
    ]
    updated_definition = AlignmentDefinition(
        alignment_id=alignment_id,
        name=_normalized_alignment_name(
            name,
            source_filenames_tuple,
            engine_kind=engine_kind,
            reference_name=reference_name,
            existing_names=existing_names,
        ),
        source_filenames=source_filenames_tuple,
        engine_kind=engine_kind,
        assembly_config=assembly_config,
        reference_name=_normalized_optional_text(reference_name),
        reference_text=_normalized_optional_text(reference_text),
        strand_policy=strand_policy,
    )

    alignments_by_id = dict(current_state.alignments_by_id)
    alignments_by_id[alignment_id] = updated_definition
    write_alignment_session_state(
        session_state,
        AlignmentSessionState(
            batch_signature=current_state.batch_signature,
            alignments_by_id=alignments_by_id,
            selected_alignment_id=alignment_id,
        ),
    )
    return updated_definition


def delete_alignment_definition(
    session_state: SessionStateWriter,
    *,
    alignment_id: str,
) -> None:
    """Delete one saved alignment definition."""
    current_state = read_alignment_session_state(session_state)
    alignments_by_id = dict(current_state.alignments_by_id)
    alignments_by_id.pop(alignment_id, None)
    next_selected_id = (
        current_state.selected_alignment_id
        if current_state.selected_alignment_id in alignments_by_id
        else next(iter(alignments_by_id), None)
    )
    write_alignment_session_state(
        session_state,
        AlignmentSessionState(
            batch_signature=current_state.batch_signature,
            alignments_by_id=alignments_by_id,
            selected_alignment_id=next_selected_id,
        ),
    )


def set_selected_alignment_id(
    session_state: SessionStateWriter,
    selected_alignment_id: str | None,
) -> None:
    """Persist the currently selected saved alignment."""
    session_state[_SELECTED_ALIGNMENT_ID_SESSION_KEY] = selected_alignment_id


def clear_alignment_session_state(session_state: SessionStateWriter) -> None:
    """Remove all saved alignment state from session storage."""
    session_state.pop(_ALIGNMENT_BATCH_SIGNATURE_SESSION_KEY, None)
    session_state.pop(_ALIGNMENTS_BY_ID_SESSION_KEY, None)
    session_state.pop(_SELECTED_ALIGNMENT_ID_SESSION_KEY, None)


def suggest_alignment_name(
    source_filenames: Iterable[str],
    *,
    engine_kind: AlignmentEngineKind,
    reference_name: str | None = None,
    existing_names: Iterable[str] = (),
) -> str:
    """Return a stable human-readable default name for one alignment."""
    normalized_source_filenames = tuple(source_filenames)
    return _normalized_alignment_name(
        "",
        normalized_source_filenames,
        engine_kind=engine_kind,
        reference_name=reference_name,
        existing_names=existing_names,
    )


def _normalized_alignment_name(
    value: str,
    source_filenames: tuple[str, ...],
    *,
    engine_kind: AlignmentEngineKind,
    reference_name: str | None,
    existing_names: Iterable[str],
) -> str:
    stripped_value = value.strip()
    base_name = (
        stripped_value
        if stripped_value
        else _derived_name_from_sources(
            source_filenames,
            engine_kind=engine_kind,
            reference_name=reference_name,
        )
    )
    return _dedupe_name(base_name, existing_names=existing_names)


def _derived_name_from_sources(
    source_filenames: tuple[str, ...],
    *,
    engine_kind: AlignmentEngineKind,
    reference_name: str | None,
) -> str:
    if not source_filenames:
        return "Alignment"

    stems = [Path(source_filename).stem for source_filename in source_filenames]
    if engine_kind == "pairwise":
        if len(stems) == 1:
            return stems[0]
        if len(stems) == 2:
            return f"{stems[0]} + {stems[1]}"
        return f"{stems[0]} + {len(stems) - 1} more"

    reference_label = _normalized_optional_text(reference_name) or "reference"
    if len(stems) == 1:
        return f"{stems[0]} vs {reference_label}"
    return f"{stems[0]} + {len(stems) - 1} more vs {reference_label}"


def _dedupe_name(base_name: str, *, existing_names: Iterable[str]) -> str:
    existing_name_set = {name for name in existing_names if isinstance(name, str)}
    if base_name not in existing_name_set:
        return base_name

    suffix = 2
    while True:
        candidate_name = f"{base_name} ({suffix})"
        if candidate_name not in existing_name_set:
            return candidate_name
        suffix += 1


def _normalized_optional_text(value: str | None) -> str | None:
    if not isinstance(value, str):
        return None
    stripped_value = value.strip()
    return stripped_value or None


def _coerce_batch_signature(value: object) -> BatchSignature | None:
    if not isinstance(value, tuple):
        return None
    signature: list[tuple[str, int, str]] = []
    for item in value:
        if (
            not isinstance(item, tuple)
            or len(item) != 3
            or not isinstance(item[0], str)
            or not isinstance(item[1], int)
            or not isinstance(item[2], str)
        ):
            return None
        signature.append((item[0], item[1], item[2]))
    return tuple(signature)


def _coerce_alignments_by_id(value: object) -> dict[str, AlignmentDefinition]:
    if not isinstance(value, dict):
        return {}
    return {
        alignment_id: definition
        for alignment_id, definition in value.items()
        if isinstance(alignment_id, str) and isinstance(definition, AlignmentDefinition)
    }


def _coerce_selected_alignment_id(value: object) -> str | None:
    return value if isinstance(value, str) else None
