from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final, Literal, TypeAlias
from uuid import uuid4

from abi_sauce.assembly import AssemblyConfig
from abi_sauce.services.batch import BatchSignature

AssemblyEngineKind = Literal["pairwise"]

SessionStateReader: TypeAlias = Any
SessionStateWriter: TypeAlias = Any

_ASSEMBLY_BATCH_SIGNATURE_SESSION_KEY: Final[str] = "abi_sauce.assembly.batch_signature"
_ASSEMBLIES_BY_ID_SESSION_KEY: Final[str] = "abi_sauce.assembly.assemblies_by_id"
_SELECTED_ASSEMBLY_ID_SESSION_KEY: Final[str] = "abi_sauce.assembly.selected_id"
_EXPORT_SELECTED_ASSEMBLY_IDS_SESSION_KEY: Final[str] = (
    "abi_sauce.assembly.export_selected_ids"
)
_AUTO_PROMPTED_SIGNATURE_SESSION_KEY: Final[str] = (
    "abi_sauce.assembly.auto_prompted_signature"
)


@dataclass(frozen=True, slots=True)
class AssemblyDefinition:
    """One saved user-defined assembly specification."""

    assembly_id: str
    name: str
    source_filenames: tuple[str, ...]
    config: AssemblyConfig = AssemblyConfig()
    engine_kind: AssemblyEngineKind = "pairwise"


@dataclass(frozen=True, slots=True)
class AssemblySessionState:
    """Session-scoped saved assembly definitions for the active batch."""

    batch_signature: BatchSignature | None = None
    assemblies_by_id: dict[str, AssemblyDefinition] = field(default_factory=dict)
    selected_assembly_id: str | None = None
    export_selected_assembly_ids: frozenset[str] = frozenset()


def read_assembly_session_state(
    session_state: SessionStateReader,
) -> AssemblySessionState:
    """Read saved assembly definitions from session state."""
    return AssemblySessionState(
        batch_signature=_coerce_batch_signature(
            session_state.get(_ASSEMBLY_BATCH_SIGNATURE_SESSION_KEY)
        ),
        assemblies_by_id=_coerce_assemblies_by_id(
            session_state.get(_ASSEMBLIES_BY_ID_SESSION_KEY)
        ),
        selected_assembly_id=_coerce_selected_assembly_id(
            session_state.get(_SELECTED_ASSEMBLY_ID_SESSION_KEY)
        ),
        export_selected_assembly_ids=_coerce_export_selected_ids(
            session_state.get(_EXPORT_SELECTED_ASSEMBLY_IDS_SESSION_KEY)
        ),
    )


def sync_assembly_session_state(
    session_state: SessionStateWriter,
    *,
    batch_signature: BatchSignature,
    parsed_record_names: Iterable[str],
) -> AssemblySessionState:
    """Sanitize saved assemblies for the currently active parsed batch."""
    record_names = frozenset(parsed_record_names)
    current_state = read_assembly_session_state(session_state)

    sanitized_assemblies_by_id = {
        assembly_id: definition
        for assembly_id, definition in current_state.assemblies_by_id.items()
        if definition.source_filenames
        and all(
            source_filename in record_names
            for source_filename in definition.source_filenames
        )
    }
    sanitized_export_selected_ids = frozenset(
        assembly_id
        for assembly_id in current_state.export_selected_assembly_ids
        if assembly_id in sanitized_assemblies_by_id
    )
    selected_assembly_id = (
        current_state.selected_assembly_id
        if current_state.selected_assembly_id in sanitized_assemblies_by_id
        else next(iter(sanitized_assemblies_by_id), None)
    )

    next_state = AssemblySessionState(
        batch_signature=batch_signature,
        assemblies_by_id=sanitized_assemblies_by_id,
        selected_assembly_id=selected_assembly_id,
        export_selected_assembly_ids=sanitized_export_selected_ids,
    )
    write_assembly_session_state(session_state, next_state)
    return next_state


def write_assembly_session_state(
    session_state: SessionStateWriter,
    assembly_state: AssemblySessionState,
) -> None:
    """Persist the complete saved-assembly state to session storage."""
    session_state[_ASSEMBLY_BATCH_SIGNATURE_SESSION_KEY] = (
        assembly_state.batch_signature
    )
    session_state[_ASSEMBLIES_BY_ID_SESSION_KEY] = dict(assembly_state.assemblies_by_id)
    session_state[_SELECTED_ASSEMBLY_ID_SESSION_KEY] = (
        assembly_state.selected_assembly_id
    )
    session_state[_EXPORT_SELECTED_ASSEMBLY_IDS_SESSION_KEY] = tuple(
        sorted(assembly_state.export_selected_assembly_ids)
    )


def create_assembly_definition(
    session_state: SessionStateWriter,
    *,
    name: str,
    source_filenames: Iterable[str],
    config: AssemblyConfig,
    engine_kind: AssemblyEngineKind = "pairwise",
) -> AssemblyDefinition:
    """Create and persist one new saved assembly definition."""
    current_state = read_assembly_session_state(session_state)
    source_filenames_tuple = tuple(source_filenames)
    definition = AssemblyDefinition(
        assembly_id=f"assembly-{uuid4().hex[:8]}",
        name=_normalized_assembly_name(
            name,
            source_filenames_tuple,
            existing_names=(
                definition.name
                for definition in current_state.assemblies_by_id.values()
            ),
        ),
        source_filenames=source_filenames_tuple,
        config=config,
        engine_kind=engine_kind,
    )
    assemblies_by_id = dict(current_state.assemblies_by_id)
    assemblies_by_id[definition.assembly_id] = definition
    write_assembly_session_state(
        session_state,
        AssemblySessionState(
            batch_signature=current_state.batch_signature,
            assemblies_by_id=assemblies_by_id,
            selected_assembly_id=definition.assembly_id,
            export_selected_assembly_ids=current_state.export_selected_assembly_ids,
        ),
    )
    return definition


def update_assembly_definition(
    session_state: SessionStateWriter,
    *,
    assembly_id: str,
    name: str,
    source_filenames: Iterable[str],
    config: AssemblyConfig,
    engine_kind: AssemblyEngineKind = "pairwise",
) -> AssemblyDefinition:
    """Update and persist one existing saved assembly definition."""
    current_state = read_assembly_session_state(session_state)
    if assembly_id not in current_state.assemblies_by_id:
        raise KeyError(assembly_id)

    source_filenames_tuple = tuple(source_filenames)
    existing_names = [
        definition.name
        for current_id, definition in current_state.assemblies_by_id.items()
        if current_id != assembly_id
    ]
    updated_definition = AssemblyDefinition(
        assembly_id=assembly_id,
        name=_normalized_assembly_name(
            name,
            source_filenames_tuple,
            existing_names=existing_names,
        ),
        source_filenames=source_filenames_tuple,
        config=config,
        engine_kind=engine_kind,
    )

    assemblies_by_id = dict(current_state.assemblies_by_id)
    assemblies_by_id[assembly_id] = updated_definition
    write_assembly_session_state(
        session_state,
        AssemblySessionState(
            batch_signature=current_state.batch_signature,
            assemblies_by_id=assemblies_by_id,
            selected_assembly_id=assembly_id,
            export_selected_assembly_ids=current_state.export_selected_assembly_ids,
        ),
    )
    return updated_definition


def delete_assembly_definition(
    session_state: SessionStateWriter,
    *,
    assembly_id: str,
) -> None:
    """Delete one saved assembly definition."""
    current_state = read_assembly_session_state(session_state)
    assemblies_by_id = dict(current_state.assemblies_by_id)
    assemblies_by_id.pop(assembly_id, None)
    next_selected_id = (
        current_state.selected_assembly_id
        if current_state.selected_assembly_id in assemblies_by_id
        else next(iter(assemblies_by_id), None)
    )
    write_assembly_session_state(
        session_state,
        AssemblySessionState(
            batch_signature=current_state.batch_signature,
            assemblies_by_id=assemblies_by_id,
            selected_assembly_id=next_selected_id,
            export_selected_assembly_ids=frozenset(
                selected_id
                for selected_id in current_state.export_selected_assembly_ids
                if selected_id in assemblies_by_id
            ),
        ),
    )


def set_selected_assembly_id(
    session_state: SessionStateWriter,
    selected_assembly_id: str | None,
) -> None:
    """Persist the currently selected saved assembly."""
    session_state[_SELECTED_ASSEMBLY_ID_SESSION_KEY] = selected_assembly_id


def set_export_selected_assembly_ids(
    session_state: SessionStateWriter,
    assembly_ids: Iterable[str],
) -> None:
    """Persist the set of saved assemblies selected for export."""
    session_state[_EXPORT_SELECTED_ASSEMBLY_IDS_SESSION_KEY] = tuple(
        sorted(
            {
                assembly_id
                for assembly_id in assembly_ids
                if isinstance(assembly_id, str)
            }
        )
    )


def clear_assembly_session_state(session_state: SessionStateWriter) -> None:
    """Remove all saved assembly state from session storage."""
    session_state.pop(_ASSEMBLY_BATCH_SIGNATURE_SESSION_KEY, None)
    session_state.pop(_ASSEMBLIES_BY_ID_SESSION_KEY, None)
    session_state.pop(_SELECTED_ASSEMBLY_ID_SESSION_KEY, None)
    session_state.pop(_EXPORT_SELECTED_ASSEMBLY_IDS_SESSION_KEY, None)
    session_state.pop(_AUTO_PROMPTED_SIGNATURE_SESSION_KEY, None)


def suggest_assembly_name(
    source_filenames: Iterable[str],
    *,
    existing_names: Iterable[str] = (),
) -> str:
    """Return a stable human-readable default name for one assembly."""
    normalized_source_filenames = tuple(source_filenames)
    return _normalized_assembly_name(
        "",
        normalized_source_filenames,
        existing_names=existing_names,
    )


def _normalized_assembly_name(
    value: str,
    source_filenames: tuple[str, ...],
    *,
    existing_names: Iterable[str],
) -> str:
    stripped_value = value.strip()
    base_name = (
        stripped_value
        if stripped_value
        else _derived_name_from_sources(source_filenames)
    )
    return _dedupe_name(base_name, existing_names=existing_names)


def _derived_name_from_sources(source_filenames: tuple[str, ...]) -> str:
    if not source_filenames:
        return "Assembly"
    stems = [Path(source_filename).stem for source_filename in source_filenames]
    if len(stems) == 1:
        return stems[0]
    if len(stems) == 2:
        return f"{stems[0]} + {stems[1]}"
    return f"{stems[0]} + {len(stems) - 1} more"


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


def _coerce_assemblies_by_id(value: object) -> dict[str, AssemblyDefinition]:
    if not isinstance(value, dict):
        return {}
    return {
        assembly_id: definition
        for assembly_id, definition in value.items()
        if isinstance(assembly_id, str) and isinstance(definition, AssemblyDefinition)
    }


def _coerce_selected_assembly_id(value: object) -> str | None:
    return value if isinstance(value, str) else None


def _coerce_export_selected_ids(value: object) -> frozenset[str]:
    if not isinstance(value, tuple | list | set | frozenset):
        return frozenset()
    return frozenset(
        selected_id for selected_id in value if isinstance(selected_id, str)
    )
