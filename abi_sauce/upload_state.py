from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Final, TypeAlias

from abi_sauce.models import SequenceUpload
from abi_sauce.services.batch import BatchSignature, ParsedBatch, build_batch_signature

_ACTIVE_UPLOADS_SESSION_KEY: Final[str] = "abi_sauce.active_uploads"
_ACTIVE_PARSED_BATCH_SESSION_KEY: Final[str] = "abi_sauce.active_parsed_batch"


SessionStateReader: TypeAlias = Any
SessionStateWriter: TypeAlias = Any


@dataclass(frozen=True, slots=True)
class ActiveBatchState:
    """Session-scoped active uploads plus their parsed batch, if available."""

    uploads: tuple[SequenceUpload, ...] = ()
    parsed_batch: ParsedBatch | None = None

    @property
    def signature(self) -> BatchSignature | None:
        """Return the current batch signature, if one can be resolved."""
        if self.parsed_batch is not None:
            return self.parsed_batch.signature
        if not self.uploads:
            return None
        return build_batch_signature(self.uploads)


def read_active_batch_state(
    session_state: SessionStateReader,
) -> ActiveBatchState:
    """Read the current active-batch session state into one pure dataclass."""
    parsed_batch = get_active_parsed_batch(session_state)
    uploads = (
        parsed_batch.uploads
        if parsed_batch is not None
        else _coerce_uploads_tuple(session_state.get(_ACTIVE_UPLOADS_SESSION_KEY))
    )
    return ActiveBatchState(
        uploads=uploads,
        parsed_batch=parsed_batch,
    )


def get_active_uploads(
    session_state: SessionStateReader,
) -> tuple[SequenceUpload, ...]:
    """Return the current active normalized uploads from session state."""
    return read_active_batch_state(session_state).uploads


def get_active_parsed_batch(
    session_state: SessionStateReader,
) -> ParsedBatch | None:
    """Return the current active parsed batch from session state, if any."""
    parsed_batch = session_state.get(_ACTIVE_PARSED_BATCH_SESSION_KEY)
    return parsed_batch if isinstance(parsed_batch, ParsedBatch) else None


def get_active_batch_signature(
    session_state: SessionStateReader,
) -> BatchSignature | None:
    """Return the signature for the active batch stored in session state."""
    return read_active_batch_state(session_state).signature


def set_active_uploads(
    session_state: SessionStateWriter,
    uploads: Iterable[SequenceUpload],
    *,
    parsed_batch: ParsedBatch | None = None,
) -> None:
    """Persist one active normalized upload batch into session state."""
    uploads_tuple = tuple(uploads)
    if parsed_batch is not None and parsed_batch.uploads != uploads_tuple:
        raise ValueError("parsed_batch.uploads must match uploads")

    session_state[_ACTIVE_UPLOADS_SESSION_KEY] = uploads_tuple
    session_state[_ACTIVE_PARSED_BATCH_SESSION_KEY] = parsed_batch


def set_active_parsed_batch(
    session_state: SessionStateWriter,
    parsed_batch: ParsedBatch,
) -> None:
    """Persist one active parsed batch and its normalized uploads."""
    set_active_uploads(
        session_state,
        parsed_batch.uploads,
        parsed_batch=parsed_batch,
    )


def merge_uploads(
    existing_uploads: Iterable[SequenceUpload],
    incoming_uploads: Iterable[SequenceUpload],
) -> tuple[SequenceUpload, ...]:
    """Return a deterministic filename-keyed merge of existing and incoming uploads."""
    uploads_by_filename = {upload.filename: upload for upload in existing_uploads}
    for upload in incoming_uploads:
        uploads_by_filename[upload.filename] = upload

    return tuple(
        uploads_by_filename[filename] for filename in sorted(uploads_by_filename)
    )


def clear_active_batch(session_state: SessionStateWriter) -> None:
    """Remove all active upload/parse state from the session mapping."""
    session_state.pop(_ACTIVE_UPLOADS_SESSION_KEY, None)
    session_state.pop(_ACTIVE_PARSED_BATCH_SESSION_KEY, None)


def _coerce_uploads_tuple(value: object) -> tuple[SequenceUpload, ...]:
    if not isinstance(value, tuple | list):
        return ()

    uploads: list[SequenceUpload] = []
    for item in value:
        if not isinstance(item, SequenceUpload):
            return ()
        uploads.append(item)

    return tuple(uploads)
