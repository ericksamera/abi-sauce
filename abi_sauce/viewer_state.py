from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Final, TypeAlias

from abi_sauce.services.batch import BatchSignature
from abi_sauce.trim_state import BatchTrimState, TrimScope
from abi_sauce.trimming import TrimConfig

SessionStateReader: TypeAlias = Any
SessionStateWriter: TypeAlias = Any

_VIEWER_BATCH_SIGNATURE_SESSION_KEY: Final[str] = "abi_sauce.viewer.batch_signature"
_SELECTED_RECORD_NAME_SESSION_KEY: Final[str] = "abi_sauce.viewer.selected_record_name"
_TRIM_SCOPE_SESSION_KEY: Final[str] = "abi_sauce.viewer.trim_scope"
_GLOBAL_TRIM_CONFIG_SESSION_KEY: Final[str] = "abi_sauce.viewer.global_trim_config"
_TRIM_CONFIGS_BY_RECORD_SESSION_KEY: Final[str] = (
    "abi_sauce.viewer.trim_configs_by_record"
)


@dataclass(frozen=True, slots=True)
class ViewerSessionState:
    """Shared sample-selection and trim state for viewer pages."""

    batch_signature: BatchSignature | None = None
    selected_record_name: str | None = None
    trim_state: BatchTrimState = BatchTrimState()


def read_viewer_session_state(session_state: SessionStateReader) -> ViewerSessionState:
    """Read the shared viewer state from session storage."""
    missing = object()
    default_trim_state = BatchTrimState()
    raw_global_trim_config = session_state.get(_GLOBAL_TRIM_CONFIG_SESSION_KEY, missing)
    global_trim_config = (
        default_trim_state.global_trim_config
        if raw_global_trim_config is missing
        else _coerce_trim_config(raw_global_trim_config)
    )
    return ViewerSessionState(
        batch_signature=_coerce_batch_signature(
            session_state.get(_VIEWER_BATCH_SIGNATURE_SESSION_KEY)
        ),
        selected_record_name=_coerce_selected_record_name(
            session_state.get(_SELECTED_RECORD_NAME_SESSION_KEY)
        ),
        trim_state=BatchTrimState(
            trim_scope=_coerce_trim_scope(session_state.get(_TRIM_SCOPE_SESSION_KEY)),
            global_trim_config=global_trim_config,
            trim_configs_by_record=_coerce_trim_configs_by_record(
                session_state.get(_TRIM_CONFIGS_BY_RECORD_SESSION_KEY)
            ),
        ),
    )


def sync_viewer_session_state(
    session_state: SessionStateWriter,
    *,
    batch_signature: BatchSignature,
    parsed_record_names: Iterable[str],
) -> ViewerSessionState:
    """Sanitize the shared viewer state for the current parsed batch."""
    record_names = tuple(parsed_record_names)
    first_record_name = record_names[0] if record_names else None
    current_state = read_viewer_session_state(session_state)

    selected_record_name = (
        current_state.selected_record_name
        if current_state.selected_record_name in record_names
        else first_record_name
    )
    sanitized_trim_configs = {
        record_name: config
        for record_name, config in current_state.trim_state.trim_configs_by_record.items()
        if record_name in record_names
    }
    next_state = ViewerSessionState(
        batch_signature=batch_signature,
        selected_record_name=selected_record_name,
        trim_state=BatchTrimState(
            trim_scope=current_state.trim_state.trim_scope,
            global_trim_config=current_state.trim_state.global_trim_config,
            trim_configs_by_record=sanitized_trim_configs,
        ),
    )
    write_viewer_session_state(session_state, next_state)
    return next_state


def write_viewer_session_state(
    session_state: SessionStateWriter,
    viewer_state: ViewerSessionState,
) -> None:
    """Persist the complete shared viewer state into session storage."""
    session_state[_VIEWER_BATCH_SIGNATURE_SESSION_KEY] = viewer_state.batch_signature
    session_state[_SELECTED_RECORD_NAME_SESSION_KEY] = viewer_state.selected_record_name
    set_batch_trim_state(session_state, viewer_state.trim_state)


def get_selected_record_name(session_state: SessionStateReader) -> str | None:
    """Return the shared currently selected record name, if available."""
    return read_viewer_session_state(session_state).selected_record_name


def set_selected_record_name(
    session_state: SessionStateWriter,
    selected_record_name: str | None,
) -> None:
    """Persist the shared selected-record choice."""
    session_state[_SELECTED_RECORD_NAME_SESSION_KEY] = selected_record_name


def get_batch_trim_state(session_state: SessionStateReader) -> BatchTrimState:
    """Return the shared batch trim state from session storage."""
    return read_viewer_session_state(session_state).trim_state


def set_batch_trim_state(
    session_state: SessionStateWriter,
    trim_state: BatchTrimState,
) -> None:
    """Persist the shared batch trim state into session storage."""
    session_state[_TRIM_SCOPE_SESSION_KEY] = trim_state.trim_scope
    session_state[_GLOBAL_TRIM_CONFIG_SESSION_KEY] = trim_state.global_trim_config
    session_state[_TRIM_CONFIGS_BY_RECORD_SESSION_KEY] = dict(
        trim_state.trim_configs_by_record
    )


def clear_viewer_session_state(session_state: SessionStateWriter) -> None:
    """Remove the shared viewer/session trim state from session storage."""
    session_state.pop(_VIEWER_BATCH_SIGNATURE_SESSION_KEY, None)
    session_state.pop(_SELECTED_RECORD_NAME_SESSION_KEY, None)
    session_state.pop(_TRIM_SCOPE_SESSION_KEY, None)
    session_state.pop(_GLOBAL_TRIM_CONFIG_SESSION_KEY, None)
    session_state.pop(_TRIM_CONFIGS_BY_RECORD_SESSION_KEY, None)


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


def _coerce_selected_record_name(value: object) -> str | None:
    return value if isinstance(value, str) else None


def _coerce_trim_scope(value: object) -> TrimScope:
    if value == "all":
        return "all"
    if value == "selected":
        return "selected"
    return "all"


def _coerce_trim_config(value: object) -> TrimConfig | None:
    return value if isinstance(value, TrimConfig) else None


def _coerce_trim_configs_by_record(value: object) -> dict[str, TrimConfig]:
    if not isinstance(value, dict):
        return {}
    trim_configs_by_record: dict[str, TrimConfig] = {}
    for record_name, config in value.items():
        if isinstance(record_name, str) and isinstance(config, TrimConfig):
            trim_configs_by_record[record_name] = config
    return trim_configs_by_record
