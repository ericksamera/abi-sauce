from __future__ import annotations

from dataclasses import dataclass

from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.chromatogram import ChromatogramView, build_chromatogram_view
from abi_sauce.reference_alignment import align_trimmed_read_to_reference
from abi_sauce.reference_alignment_multi import align_trimmed_reads_to_reference
from abi_sauce.reference_alignment_multi_presenters import (
    reference_multi_columns_to_rows,
    reference_multi_members_to_rows,
)
from abi_sauce.reference_alignment_multi_trace import (
    ReferenceMultiAlignmentTraceView,
    build_reference_multi_alignment_trace_view,
)
from abi_sauce.reference_alignment_presenters import alignment_events_to_rows
from abi_sauce.reference_alignment_trace import (
    ReferenceAlignmentTraceView,
    build_reference_alignment_trace_view,
)
from abi_sauce.reference_alignment_types import (
    AlignmentResult,
    ReferenceMultiAlignmentResult,
    StrandPolicy,
)
from abi_sauce.services.batch_trim import PreparedBatch


@dataclass(frozen=True, slots=True)
class ComputedReferenceAlignment:
    """Resolved reference-alignment result plus page-facing projections."""

    source_filename: str
    alignment_result: AlignmentResult
    event_rows: tuple[dict[str, object], ...]
    chromatogram_view: ChromatogramView
    trace_view: ReferenceAlignmentTraceView | None = None


@dataclass(frozen=True, slots=True)
class ComputedReferenceMultiAlignment:
    """Resolved shared-reference multi-read alignment plus page projections."""

    source_filenames: tuple[str, ...]
    result: ReferenceMultiAlignmentResult
    member_rows: tuple[dict[str, object], ...]
    column_rows: tuple[dict[str, object], ...]
    trace_view: ReferenceMultiAlignmentTraceView | None = None


def compute_reference_alignment(
    prepared_batch: PreparedBatch,
    *,
    source_filename: str,
    reference_text: str,
    reference_name: str | None = None,
    strand_policy: StrandPolicy = "auto",
    include_matches: bool = False,
) -> ComputedReferenceAlignment:
    """Compute one reference alignment from the current prepared batch."""
    try:
        raw_record = prepared_batch.parsed_records[source_filename]
        trim_result = prepared_batch.trim_results[source_filename]
    except KeyError as exc:
        raise KeyError(source_filename) from exc

    alignment_result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text=reference_text,
        reference_name=reference_name,
        strand_policy=strand_policy,
    )
    return ComputedReferenceAlignment(
        source_filename=source_filename,
        alignment_result=alignment_result,
        event_rows=tuple(
            alignment_events_to_rows(
                alignment_result,
                include_matches=include_matches,
            )
        ),
        chromatogram_view=build_chromatogram_view(raw_record, trim_result),
        trace_view=build_reference_alignment_trace_view(
            result=alignment_result,
            source_filename=source_filename,
            raw_record=raw_record,
            trim_result=trim_result,
        ),
    )


def compute_reference_multi_alignment(
    prepared_batch: PreparedBatch,
    *,
    source_filenames: tuple[str, ...] | list[str],
    reference_text: str,
    reference_name: str | None = None,
    strand_policy: StrandPolicy = "auto",
    config: AssemblyConfig | None = None,
    include_reference_matches: bool = False,
) -> ComputedReferenceMultiAlignment:
    """Compute one shared-reference multi-read alignment from the current batch."""
    resolved_source_filenames = tuple(source_filenames)
    missing_filenames = [
        source_filename
        for source_filename in resolved_source_filenames
        if source_filename not in prepared_batch.parsed_records
        or source_filename not in prepared_batch.trim_results
    ]
    if missing_filenames:
        raise KeyError(", ".join(sorted(missing_filenames)))

    result = align_trimmed_reads_to_reference(
        source_filenames=resolved_source_filenames,
        raw_records_by_source_filename={
            source_filename: prepared_batch.parsed_records[source_filename]
            for source_filename in resolved_source_filenames
        },
        trim_results_by_source_filename={
            source_filename: prepared_batch.trim_results[source_filename]
            for source_filename in resolved_source_filenames
        },
        reference_text=reference_text,
        reference_name=reference_name,
        strand_policy=strand_policy,
        config=config,
    )
    return ComputedReferenceMultiAlignment(
        source_filenames=resolved_source_filenames,
        result=result,
        member_rows=tuple(reference_multi_members_to_rows(result)),
        column_rows=tuple(
            reference_multi_columns_to_rows(
                result,
                include_reference_matches=include_reference_matches,
            )
        ),
        trace_view=build_reference_multi_alignment_trace_view(
            result=result,
            raw_records_by_source_filename={
                source_filename: prepared_batch.parsed_records[source_filename]
                for source_filename in resolved_source_filenames
            },
            trim_results_by_source_filename={
                source_filename: prepared_batch.trim_results[source_filename]
                for source_filename in resolved_source_filenames
            },
        ),
    )


__all__ = [
    "ComputedReferenceAlignment",
    "ComputedReferenceMultiAlignment",
    "compute_reference_alignment",
    "compute_reference_multi_alignment",
]
