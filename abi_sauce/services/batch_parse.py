from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass, replace
import hashlib
from typing import Protocol

from abi_sauce.exceptions import AbiParseError
from abi_sauce.models import SequenceRecord, SequenceUpload
from abi_sauce.parsers.abi import parse_ab1_upload

BatchSignature = tuple[tuple[str, int, str], ...]


class UploadedFileLike(Protocol):
    """Minimal file-uploader interface needed by the batch service."""

    name: str

    def getvalue(self) -> bytes:
        """Return the raw file contents."""


@dataclass(frozen=True, slots=True)
class ParsedBatch:
    """Pure parsed batch state derived from uploaded files."""

    uploads: tuple[SequenceUpload, ...]
    parsed_records: dict[str, SequenceRecord]
    parse_errors: dict[str, str]
    signature: BatchSignature


def normalize_uploaded_files(
    uploaded_files: Iterable[UploadedFileLike],
) -> tuple[SequenceUpload, ...]:
    """Convert uploaded Streamlit files into sorted framework-agnostic uploads."""
    return tuple(
        SequenceUpload(
            filename=uploaded_file.name,
            content=uploaded_file.getvalue(),
        )
        for uploaded_file in sorted(uploaded_files, key=lambda file: file.name)
    )


def build_batch_signature(uploads: Iterable[SequenceUpload]) -> BatchSignature:
    """Return a stable signature for a batch of uploads."""
    return tuple(
        (
            upload.filename,
            upload.size_bytes,
            _content_digest(upload.content),
        )
        for upload in uploads
    )


def _content_digest(content: bytes) -> str:
    """Return a stable digest for uploaded file contents."""
    return hashlib.blake2b(content, digest_size=16).hexdigest()


def parse_uploaded_batch(
    uploaded_files: Iterable[UploadedFileLike],
    *,
    parse_upload: Callable[[SequenceUpload], SequenceRecord] = parse_ab1_upload,
) -> ParsedBatch:
    """Normalize and parse one batch of uploaded ABI files."""
    return parse_uploads(
        normalize_uploaded_files(uploaded_files),
        parse_upload=parse_upload,
    )


def parse_uploads(
    uploads: Iterable[SequenceUpload],
    *,
    parse_upload: Callable[[SequenceUpload], SequenceRecord] = parse_ab1_upload,
) -> ParsedBatch:
    """Parse normalized uploads into records plus parse failures."""
    uploads_tuple = tuple(uploads)
    parsed_records: dict[str, SequenceRecord] = {}
    parse_errors: dict[str, str] = {}

    for upload in uploads_tuple:
        try:
            parsed_records[upload.filename] = parse_upload(upload)
        except AbiParseError as exc:
            parse_errors[upload.filename] = str(exc)

    return ParsedBatch(
        uploads=uploads_tuple,
        parsed_records=parsed_records,
        parse_errors=parse_errors,
        signature=build_batch_signature(uploads_tuple),
    )


def replace_parsed_batch_record(
    parsed_batch: ParsedBatch,
    *,
    source_filename: str,
    record: SequenceRecord,
) -> ParsedBatch:
    """Return a ParsedBatch with one parsed record replaced by source filename."""
    if source_filename not in parsed_batch.parsed_records:
        raise KeyError(source_filename)

    updated_parsed_records = dict(parsed_batch.parsed_records)
    updated_parsed_records[source_filename] = record
    return replace(
        parsed_batch,
        parsed_records=updated_parsed_records,
        parse_errors=dict(parsed_batch.parse_errors),
    )


__all__ = [
    "BatchSignature",
    "ParsedBatch",
    "UploadedFileLike",
    "build_batch_signature",
    "normalize_uploaded_files",
    "parse_uploaded_batch",
    "parse_uploads",
    "replace_parsed_batch_record",
]
