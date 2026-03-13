from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
import hashlib
from typing import Protocol

from abi_sauce.batch import (
    BatchExportPolicy,
    BatchSummary,
    ExportFormat,
    build_batch_export_policy,
    build_batch_summary,
)
from abi_sauce.export import to_fasta_batch, to_fastq_batch, to_zip_batch
from abi_sauce.exceptions import AbiParseError
from abi_sauce.models import SequenceRecord, SequenceUpload
from abi_sauce.parsers.abi import parse_ab1_upload
from abi_sauce.trimming import TrimConfig, TrimResult, trim_sequence_record

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


@dataclass(frozen=True, slots=True)
class PreparedBatch:
    """Parsed batch plus derived trim/export state."""

    uploads: tuple[SequenceUpload, ...]
    parsed_records: dict[str, SequenceRecord]
    parse_errors: dict[str, str]
    signature: BatchSignature
    trim_results: dict[str, TrimResult]
    batch_summary: BatchSummary
    batch_export_policy: BatchExportPolicy


@dataclass(frozen=True, slots=True)
class BatchExportSelection:
    """Resolved export subset for one export choice."""

    export_format: ExportFormat
    require_min_length: bool
    eligible_records: tuple[SequenceRecord, ...]
    ineligible_reasons: tuple[tuple[str, tuple[str, ...]], ...]


@dataclass(frozen=True, slots=True)
class BatchDownloadArtifact:
    """Resolved export subset plus serialized batch download payload."""

    export_format: ExportFormat
    require_min_length: bool
    concatenate_batch: bool
    eligible_records: tuple[SequenceRecord, ...]
    ineligible_reasons: tuple[tuple[str, tuple[str, ...]], ...]
    data: str | bytes = ""
    filename: str = ""
    mime: str = ""

    @property
    def is_downloadable(self) -> bool:
        """Return whether a concrete download payload is available."""
        return bool(self.filename)


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


def parse_uploaded_batch(uploaded_files: Iterable[UploadedFileLike]) -> ParsedBatch:
    """Normalize and parse one batch of uploaded ABI files."""
    return parse_uploads(normalize_uploaded_files(uploaded_files))


def parse_uploads(uploads: Iterable[SequenceUpload]) -> ParsedBatch:
    """Parse normalized uploads into records plus parse failures."""
    uploads_tuple = tuple(uploads)
    parsed_records: dict[str, SequenceRecord] = {}
    parse_errors: dict[str, str] = {}

    for upload in uploads_tuple:
        try:
            parsed_records[upload.filename] = parse_ab1_upload(upload)
        except AbiParseError as exc:
            parse_errors[upload.filename] = str(exc)

    return ParsedBatch(
        uploads=uploads_tuple,
        parsed_records=parsed_records,
        parse_errors=parse_errors,
        signature=build_batch_signature(uploads_tuple),
    )


def apply_trim_config(
    parsed_batch: ParsedBatch,
    trim_config: TrimConfig,
) -> PreparedBatch:
    """Apply one trim config across an entire parsed batch."""
    return apply_trim_configs(
        parsed_batch,
        default_trim_config=trim_config,
    )


def apply_trim_configs(
    parsed_batch: ParsedBatch,
    *,
    default_trim_config: TrimConfig | None = None,
    trim_configs_by_name: Mapping[str, TrimConfig] | None = None,
) -> PreparedBatch:
    """Apply trimming across a parsed batch with optional per-record configs."""
    resolved_default_trim_config = (
        TrimConfig() if default_trim_config is None else default_trim_config
    )
    resolved_trim_configs_by_name = dict(trim_configs_by_name or {})

    trim_results = {
        name: trim_sequence_record(
            parsed_record,
            resolved_trim_configs_by_name.get(name, resolved_default_trim_config),
        )
        for name, parsed_record in parsed_batch.parsed_records.items()
    }

    return PreparedBatch(
        uploads=parsed_batch.uploads,
        parsed_records=parsed_batch.parsed_records,
        parse_errors=parsed_batch.parse_errors,
        signature=parsed_batch.signature,
        trim_results=trim_results,
        batch_summary=build_batch_summary(
            uploads=parsed_batch.uploads,
            trim_results=trim_results,
            parse_errors=parsed_batch.parse_errors,
        ),
        batch_export_policy=build_batch_export_policy(trim_results=trim_results),
    )


def select_batch_export(
    prepared_batch: PreparedBatch,
    *,
    export_format: ExportFormat,
    require_min_length: bool = False,
) -> BatchExportSelection:
    """Resolve the currently eligible batch export subset."""
    return BatchExportSelection(
        export_format=export_format,
        require_min_length=require_min_length,
        eligible_records=prepared_batch.batch_export_policy.eligible_records(
            export_format=export_format,
            require_min_length=require_min_length,
        ),
        ineligible_reasons=prepared_batch.batch_export_policy.ineligible_reasons_by_filename(
            export_format=export_format,
            require_min_length=require_min_length,
        ),
    )


def prepare_batch_download(
    prepared_batch: PreparedBatch,
    *,
    export_format: ExportFormat,
    concatenate_batch: bool,
    filename_stem: str,
    require_min_length: bool = False,
) -> BatchDownloadArtifact:
    """Build one batch download artifact for the current export selection."""
    export_selection = select_batch_export(
        prepared_batch,
        export_format=export_format,
        require_min_length=require_min_length,
    )

    if not export_selection.eligible_records:
        return BatchDownloadArtifact(
            export_format=export_format,
            require_min_length=require_min_length,
            concatenate_batch=concatenate_batch,
            eligible_records=export_selection.eligible_records,
            ineligible_reasons=export_selection.ineligible_reasons,
        )

    export_data, export_filename, export_mime = _serialize_batch_download(
        records=export_selection.eligible_records,
        export_format=export_format,
        concatenate_batch=concatenate_batch,
        filename_stem=filename_stem,
    )

    return BatchDownloadArtifact(
        export_format=export_format,
        require_min_length=require_min_length,
        concatenate_batch=concatenate_batch,
        eligible_records=export_selection.eligible_records,
        ineligible_reasons=export_selection.ineligible_reasons,
        data=export_data,
        filename=export_filename,
        mime=export_mime,
    )


def _serialize_batch_download(
    *,
    records: tuple[SequenceRecord, ...],
    export_format: ExportFormat,
    concatenate_batch: bool,
    filename_stem: str,
) -> tuple[str | bytes, str, str]:
    """Serialize the current eligible batch export selection."""
    normalized_filename_stem = filename_stem.strip() or "abi-sauce-batch"

    if concatenate_batch and export_format == "fasta":
        return (
            to_fasta_batch(records),
            f"{normalized_filename_stem}.fasta",
            "text/plain",
        )

    if concatenate_batch and export_format == "fastq":
        return (
            to_fastq_batch(records),
            f"{normalized_filename_stem}.fastq",
            "text/plain",
        )

    return (
        to_zip_batch(records, export_format=export_format),
        f"{normalized_filename_stem}.zip",
        "application/zip",
    )
