from __future__ import annotations

from dataclasses import dataclass

from abi_sauce.chromatogram import ChromatogramView, build_chromatogram_view
from abi_sauce.reference_alignment import align_trimmed_read_to_reference
from abi_sauce.reference_alignment_presenters import alignment_events_to_rows
from abi_sauce.reference_alignment_types import AlignmentResult, StrandPolicy
from abi_sauce.services.batch_trim import PreparedBatch


@dataclass(frozen=True, slots=True)
class ComputedReferenceAlignment:
    """Resolved reference-alignment result plus page-facing projections."""

    source_filename: str
    alignment_result: AlignmentResult
    event_rows: tuple[dict[str, object], ...]
    chromatogram_view: ChromatogramView


def compute_reference_alignment(
    prepared_batch: PreparedBatch,
    *,
    source_filename: str,
    reference_text: str,
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
    )


__all__ = [
    "ComputedReferenceAlignment",
    "compute_reference_alignment",
]
