from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from abi_sauce.models import SequenceRecord
from abi_sauce.orientation import materialize_oriented_record
from abi_sauce.trimming import TrimResult


@dataclass(frozen=True, slots=True)
class PreparedTrimmedRead:
    """One trimmed read prepared in its current display orientation."""

    source_filename: str
    display_name: str
    raw_record: SequenceRecord
    trim_result: TrimResult
    display_sequence: str
    display_qualities: list[int] | None
    trimmed_length: int
    has_trace_data: bool
    has_qualities: bool


def prepare_trimmed_read(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    source_filename: str | None = None,
) -> PreparedTrimmedRead:
    """Return one trimmed read with display-oriented sequence and qualities."""
    resolved_source_filename = _resolved_source_filename(
        raw_record=raw_record,
        trim_result=trim_result,
        source_filename=source_filename,
    )
    display_record = materialize_oriented_record(trim_result.record)
    return PreparedTrimmedRead(
        source_filename=resolved_source_filename,
        display_name=_display_name(trim_result.record, resolved_source_filename),
        raw_record=raw_record,
        trim_result=trim_result,
        display_sequence=display_record.sequence.upper(),
        display_qualities=(
            None
            if display_record.qualities is None
            else [int(value) for value in display_record.qualities]
        ),
        trimmed_length=trim_result.trimmed_length,
        has_trace_data=raw_record.trace_data is not None,
        has_qualities=trim_result.record.qualities is not None,
    )


def prepare_trimmed_reads(
    *,
    source_filenames: Iterable[str],
    raw_records_by_source_filename: Mapping[str, SequenceRecord],
    trim_results_by_source_filename: Mapping[str, TrimResult],
) -> tuple[PreparedTrimmedRead, ...]:
    """Return prepared trimmed reads for one ordered source-filename selection."""
    return tuple(
        prepare_trimmed_read(
            source_filename=source_filename,
            raw_record=raw_records_by_source_filename[source_filename],
            trim_result=trim_results_by_source_filename[source_filename],
        )
        for source_filename in source_filenames
    )


def nonempty_prepared_reads(
    prepared_reads: Iterable[PreparedTrimmedRead],
) -> tuple[PreparedTrimmedRead, ...]:
    """Return only prepared reads with non-empty display sequences."""
    return tuple(
        prepared_read
        for prepared_read in prepared_reads
        if prepared_read.display_sequence
    )


def _resolved_source_filename(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    source_filename: str | None,
) -> str:
    if isinstance(source_filename, str) and source_filename:
        return source_filename
    annotation_source_filename = trim_result.record.annotations.get("source_filename")
    if isinstance(annotation_source_filename, str) and annotation_source_filename:
        return annotation_source_filename
    if raw_record.name:
        return raw_record.name
    if raw_record.record_id:
        return raw_record.record_id
    return "sequence"


def _display_name(record: SequenceRecord, source_filename: str) -> str:
    return record.name or record.record_id or source_filename
