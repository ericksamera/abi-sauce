from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Literal

from abi_sauce.models import SequenceOrientation, SequenceRecord, SequenceUpload
from abi_sauce.trimming import TrimResult

ExportFormat = Literal["fasta", "fastq"]


@dataclass(frozen=True, slots=True)
class BatchRecordSummary:
    """Per-record batch summary derived from a trimmed sequence record."""

    source_filename: str
    record_id: str
    display_name: str
    original_length: int
    trimmed_length: int
    orientation: SequenceOrientation
    total_bases_removed: int
    fixed_bases_removed_left: int
    fixed_bases_removed_right: int
    quality_bases_removed_left: int
    quality_bases_removed_right: int
    passed_min_length: bool
    has_qualities: bool
    has_trace_data: bool
    fastq_exportable: bool

    def to_row(self) -> dict[str, object]:
        """Return the default table row shown in the Streamlit batch summary."""
        return {
            "filename": self.source_filename,
            "record": self.display_name,
            "trimmed_length": self.trimmed_length,
            "orientation": self.orientation,
            "bases_removed": self.total_bases_removed,
            "fixed_left": self.fixed_bases_removed_left,
            "fixed_right": self.fixed_bases_removed_right,
            "quality_left": self.quality_bases_removed_left,
            "quality_right": self.quality_bases_removed_right,
            "passed_min_length": self.passed_min_length,
            "fastq_exportable": self.fastq_exportable,
        }


@dataclass(frozen=True, slots=True)
class BatchExportRecord:
    """Per-record export eligibility derived from a trimmed sequence record."""

    source_filename: str
    display_name: str
    record: SequenceRecord
    passed_min_length: bool
    fasta_exportable: bool
    fastq_exportable: bool
    fastq_ineligible_reasons: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class BatchExportPolicy:
    """Batch-level export eligibility and filter helpers."""

    records: tuple[BatchExportRecord, ...] = ()

    @property
    def fasta_records(self) -> tuple[BatchExportRecord, ...]:
        return self.eligible_policy_records(export_format="fasta")

    @property
    def fastq_records(self) -> tuple[BatchExportRecord, ...]:
        return self.eligible_policy_records(export_format="fastq")

    @property
    def passing_min_length_records(self) -> tuple[BatchExportRecord, ...]:
        return tuple(record for record in self.records if record.passed_min_length)

    @property
    def failing_min_length_records(self) -> tuple[BatchExportRecord, ...]:
        return tuple(record for record in self.records if not record.passed_min_length)

    def eligible_policy_records(
        self,
        *,
        export_format: ExportFormat,
        require_min_length: bool = False,
    ) -> tuple[BatchExportRecord, ...]:
        """Return per-record export entries eligible for the current export rules."""
        return tuple(
            record
            for record in self.records
            if not _record_ineligible_reasons(
                record,
                export_format=export_format,
                require_min_length=require_min_length,
            )
        )

    def eligible_records(
        self,
        *,
        export_format: ExportFormat,
        require_min_length: bool = False,
    ) -> tuple[SequenceRecord, ...]:
        """Return trimmed SequenceRecord objects eligible for export."""
        return tuple(
            record.record
            for record in self.eligible_policy_records(
                export_format=export_format,
                require_min_length=require_min_length,
            )
        )

    def ineligible_reasons_by_filename(
        self,
        *,
        export_format: ExportFormat,
        require_min_length: bool = False,
    ) -> tuple[tuple[str, tuple[str, ...]], ...]:
        """Return excluded filenames with the reasons they were filtered out."""
        return tuple(
            (record.source_filename, reasons)
            for record in self.records
            if (
                reasons := _record_ineligible_reasons(
                    record,
                    export_format=export_format,
                    require_min_length=require_min_length,
                )
            )
        )


@dataclass(frozen=True, slots=True)
class BatchSummary:
    """Aggregate summary for one uploaded ABI batch."""

    total_uploaded_files: int
    parsed_files: int
    failed_files: int
    trimmed_records: int
    records_passing_min_length: int
    records_failing_min_length: int
    fastq_exportable_records: int
    records: tuple[BatchRecordSummary, ...] = ()
    parse_errors: tuple[tuple[str, str], ...] = ()

    def table_rows(self) -> list[dict[str, object]]:
        """Return default rows for rendering a batch summary table."""
        return [record.to_row() for record in self.records]


def build_batch_summary(
    *,
    uploads: Iterable[SequenceUpload],
    trim_results: Mapping[str, TrimResult],
    parse_errors: Mapping[str, str] | None = None,
) -> BatchSummary:
    """Build a pure batch summary from uploads, trim results, and parse errors."""
    uploads_list = list(uploads)
    parse_errors = dict(parse_errors or {})

    records = tuple(
        _build_record_summary(source_filename=source_filename, trim_result=trim_result)
        for source_filename, trim_result in sorted(trim_results.items())
    )

    records_passing_min_length = sum(record.passed_min_length for record in records)
    fastq_exportable_records = sum(record.fastq_exportable for record in records)

    return BatchSummary(
        total_uploaded_files=len(uploads_list),
        parsed_files=len(records),
        failed_files=len(parse_errors),
        trimmed_records=len(records),
        records_passing_min_length=records_passing_min_length,
        records_failing_min_length=len(records) - records_passing_min_length,
        fastq_exportable_records=fastq_exportable_records,
        records=records,
        parse_errors=tuple(sorted(parse_errors.items())),
    )


def build_batch_export_policy(
    *,
    trim_results: Mapping[str, TrimResult],
) -> BatchExportPolicy:
    """Build pure batch export eligibility from trimmed results."""
    records = tuple(
        _build_export_record(source_filename=source_filename, trim_result=trim_result)
        for source_filename, trim_result in sorted(trim_results.items())
    )
    return BatchExportPolicy(records=records)


def _build_record_summary(
    *,
    source_filename: str,
    trim_result: TrimResult,
) -> BatchRecordSummary:
    """Build a per-record batch summary row from one trim result."""
    record = trim_result.record
    has_qualities = record.qualities is not None

    return BatchRecordSummary(
        source_filename=source_filename,
        record_id=record.record_id,
        display_name=record.name or record.record_id or source_filename,
        original_length=trim_result.original_length,
        trimmed_length=trim_result.trimmed_length,
        orientation=record.orientation,
        total_bases_removed=trim_result.original_length - trim_result.trimmed_length,
        fixed_bases_removed_left=trim_result.fixed_bases_removed_left,
        fixed_bases_removed_right=trim_result.fixed_bases_removed_right,
        quality_bases_removed_left=trim_result.quality_bases_removed_left,
        quality_bases_removed_right=trim_result.quality_bases_removed_right,
        passed_min_length=trim_result.passed_min_length,
        has_qualities=has_qualities,
        has_trace_data=record.trace_data is not None,
        fastq_exportable=has_qualities
        and record.qualities is not None
        and len(record.qualities) == len(record.sequence),
    )


def _build_export_record(
    *,
    source_filename: str,
    trim_result: TrimResult,
) -> BatchExportRecord:
    """Build a per-record export eligibility entry from one trim result."""
    record = trim_result.record
    fastq_ineligible_reasons = _fastq_ineligible_reasons(record)

    return BatchExportRecord(
        source_filename=source_filename,
        display_name=record.name or record.record_id or source_filename,
        record=record,
        passed_min_length=trim_result.passed_min_length,
        fasta_exportable=True,
        fastq_exportable=not fastq_ineligible_reasons,
        fastq_ineligible_reasons=fastq_ineligible_reasons,
    )


def _fastq_ineligible_reasons(record: SequenceRecord) -> tuple[str, ...]:
    """Return FASTQ-specific export blockers for a trimmed record."""
    if record.qualities is None:
        return ("missing per-base qualities",)
    if len(record.qualities) != len(record.sequence):
        return ("sequence/quality length mismatch",)
    return ()


def _record_ineligible_reasons(
    record: BatchExportRecord,
    *,
    export_format: ExportFormat,
    require_min_length: bool,
) -> tuple[str, ...]:
    """Return the reasons a record should be excluded from the current export."""
    reasons: list[str] = []

    if export_format == "fastq":
        reasons.extend(record.fastq_ineligible_reasons)

    if require_min_length and not record.passed_min_length:
        reasons.append("below minimum length")

    return tuple(reasons)
