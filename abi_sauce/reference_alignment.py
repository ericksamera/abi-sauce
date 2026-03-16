from __future__ import annotations

from Bio import Align

from abi_sauce.alignment_policy import (
    alignment_strands_for_policy,
    build_semiglobal_aligner,
    select_best_oriented_alignment,
)
from abi_sauce.models import SequenceRecord
from abi_sauce.oriented_reads import prepare_trimmed_read
from abi_sauce.reference_alignment_presenters import (
    alignment_events_to_rows,
    format_alignment_block,
)
from abi_sauce.reference_alignment_types import (
    AlignmentEvent,
    AlignmentResult,
    ChosenStrand,
    StrandPolicy,
)
from abi_sauce.trace_coordinates import trace_position_for_oriented_query_index
from abi_sauce.trimming import TrimResult

_ALLOWED_REFERENCE_CHARACTERS = frozenset("ACGTURYKMSWBDHVN")


def normalize_reference(reference_text: str) -> tuple[str, str]:
    """Return a `(reference_name, sequence)` tuple from FASTA or plain text."""
    lines = [line.strip() for line in reference_text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("Reference sequence is empty.")

    reference_name = "reference"
    if lines[0].startswith(">"):
        reference_name = lines[0][1:].strip() or "reference"
        lines = lines[1:]

    sequence = "".join(lines).replace(" ", "").upper()
    invalid = sorted(set(sequence) - _ALLOWED_REFERENCE_CHARACTERS)
    if not sequence:
        raise ValueError("Reference sequence is empty.")
    if invalid:
        raise ValueError(
            "Reference contains unsupported characters: " + ", ".join(invalid)
        )

    return reference_name, sequence


def build_aligner(
    *,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    open_internal_gap_score: float = -3.0,
    extend_internal_gap_score: float = -1.0,
) -> Align.PairwiseAligner:
    """Build a pairwise aligner with internal penalties and free terminal overhangs."""
    return build_semiglobal_aligner(
        match_score=match_score,
        mismatch_score=mismatch_score,
        open_internal_gap_score=open_internal_gap_score,
        extend_internal_gap_score=extend_internal_gap_score,
    )


def align_trimmed_read_to_reference(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    reference_text: str,
    reference_name: str | None = None,
    strand_policy: StrandPolicy = "auto",
    aligner: Align.PairwiseAligner | None = None,
    context_window: int = 5,
) -> AlignmentResult:
    """Align one trimmed read against a reference and return summary + event details."""
    normalized_reference_name, reference_sequence = normalize_reference(reference_text)
    resolved_reference_name = (
        normalized_reference_name if reference_name is None else reference_name
    )

    prepared_query_read = prepare_trimmed_read(
        raw_record=raw_record,
        trim_result=trim_result,
    )
    if not prepared_query_read.display_sequence:
        raise ValueError("Trimmed query sequence is empty.")

    resolved_aligner = build_aligner() if aligner is None else aligner
    allowed_strands = alignment_strands_for_policy(strand_policy)
    best_query_alignment = select_best_oriented_alignment(
        target_sequence=reference_sequence,
        display_query_sequence=prepared_query_read.display_sequence,
        display_query_qualities=prepared_query_read.display_qualities,
        aligner=resolved_aligner,
        allowed_strands=allowed_strands,
    )

    if best_query_alignment is None:
        raise ValueError("No alignment could be generated.")

    (
        aligned_reference,
        match_line,
        aligned_query,
        events,
        matches,
        mismatches,
        insertions,
        deletions,
        reference_start,
        reference_end,
        query_start,
        query_end,
    ) = _extract_alignment_details(
        alignment=best_query_alignment.alignment,
        raw_record=raw_record,
        trim_result=trim_result,
        reference_sequence=reference_sequence,
        oriented_query_sequence=best_query_alignment.sequence,
        oriented_query_qualities=best_query_alignment.qualities,
        strand=best_query_alignment.strand,
        context_window=context_window,
    )

    aligned_columns = matches + mismatches + insertions + deletions
    percent_identity = (matches / aligned_columns * 100.0) if aligned_columns else 0.0

    return AlignmentResult(
        sample_name=prepared_query_read.display_name,
        reference_name=resolved_reference_name,
        strand=best_query_alignment.strand,
        score=best_query_alignment.score,
        reference_start=reference_start,
        reference_end=reference_end,
        query_start=query_start,
        query_end=query_end,
        percent_identity=percent_identity,
        mismatch_count=mismatches,
        insertion_count=insertions,
        deletion_count=deletions,
        aligned_reference=aligned_reference,
        match_line=match_line,
        aligned_query=aligned_query,
        events=tuple(events),
    )


def _extract_alignment_details(
    *,
    alignment,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    reference_sequence: str,
    oriented_query_sequence: str,
    oriented_query_qualities: list[int] | None,
    strand: ChosenStrand,
    context_window: int,
) -> tuple[
    str,
    str,
    str,
    list[AlignmentEvent],
    int,
    int,
    int,
    int,
    int | None,
    int | None,
    int | None,
    int | None,
]:
    target_indices = alignment.indices[0]
    query_indices = alignment.indices[1]

    aligned_reference_chars: list[str] = []
    aligned_query_chars: list[str] = []
    match_line_chars: list[str] = []

    events: list[AlignmentEvent] = []
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0

    reference_positions: list[int] = []
    query_positions: list[int] = []

    for column_index, (target_index, query_index) in enumerate(
        zip(target_indices, query_indices, strict=True)
    ):
        resolved_target_index = int(target_index)
        resolved_query_index = int(query_index)

        ref_base = (
            "-"
            if resolved_target_index < 0
            else reference_sequence[resolved_target_index]
        )
        query_base = (
            "-"
            if resolved_query_index < 0
            else oriented_query_sequence[resolved_query_index]
        )

        aligned_reference_chars.append(ref_base)
        aligned_query_chars.append(query_base)

        if resolved_target_index >= 0:
            reference_positions.append(resolved_target_index)
        if resolved_query_index >= 0:
            query_positions.append(resolved_query_index)

        if resolved_target_index >= 0 and resolved_query_index >= 0:
            if ref_base == query_base:
                event_type = "match"
                match_line_chars.append("|")
                matches += 1
            else:
                event_type = "mismatch"
                match_line_chars.append(".")
                mismatches += 1
        elif resolved_target_index < 0 and resolved_query_index >= 0:
            event_type = "insertion"
            match_line_chars.append(" ")
            insertions += 1
        elif resolved_target_index >= 0 and resolved_query_index < 0:
            event_type = "deletion"
            match_line_chars.append(" ")
            deletions += 1
        else:
            continue

        qscore = (
            None
            if resolved_query_index < 0 or oriented_query_qualities is None
            else int(oriented_query_qualities[resolved_query_index])
        )
        flank_q_left, flank_q_right = _flanking_qualities_for_column(
            oriented_query_qualities,
            query_indices,
            column_index=column_index,
        )
        trace_x = trace_position_for_oriented_query_index(
            raw_record=raw_record,
            trim_result=trim_result,
            oriented_query_index=(
                None if resolved_query_index < 0 else resolved_query_index
            ),
            strand=strand,
        )
        context_ref = _context_window(
            reference_sequence,
            _context_index_for_column(target_indices, column_index),
            window=context_window,
        )
        context_query = _context_window(
            oriented_query_sequence,
            _context_index_for_column(query_indices, column_index),
            window=context_window,
        )

        events.append(
            AlignmentEvent(
                ref_pos=(
                    None if resolved_target_index < 0 else resolved_target_index + 1
                ),
                query_pos=(
                    None if resolved_query_index < 0 else resolved_query_index + 1
                ),
                event_type=event_type,
                ref_base=ref_base,
                query_base=query_base,
                qscore=qscore,
                flank_q_left=flank_q_left,
                flank_q_right=flank_q_right,
                trace_x=trace_x,
                context_ref=context_ref,
                context_query=context_query,
            )
        )

    reference_start = None if not reference_positions else reference_positions[0] + 1
    reference_end = None if not reference_positions else reference_positions[-1] + 1
    query_start = None if not query_positions else query_positions[0] + 1
    query_end = None if not query_positions else query_positions[-1] + 1

    return (
        "".join(aligned_reference_chars),
        "".join(match_line_chars),
        "".join(aligned_query_chars),
        events,
        matches,
        mismatches,
        insertions,
        deletions,
        reference_start,
        reference_end,
        query_start,
        query_end,
    )


def _flanking_qualities_for_column(
    qualities: list[int] | None,
    indices,
    *,
    column_index: int,
) -> tuple[int | None, int | None]:
    if qualities is None:
        return (None, None)

    current_index = int(indices[column_index])
    if current_index >= 0:
        left = qualities[current_index - 1] if current_index - 1 >= 0 else None
        right = (
            qualities[current_index + 1] if current_index + 1 < len(qualities) else None
        )
        return (
            None if left is None else int(left),
            None if right is None else int(right),
        )

    previous_index = None
    for scan_index in range(column_index - 1, -1, -1):
        value = int(indices[scan_index])
        if value >= 0:
            previous_index = value
            break

    next_index = None
    for scan_index in range(column_index + 1, len(indices)):
        value = int(indices[scan_index])
        if value >= 0:
            next_index = value
            break

    left = None if previous_index is None else int(qualities[previous_index])
    right = None if next_index is None else int(qualities[next_index])
    return (left, right)


def _context_index_for_column(indices, column_index: int) -> int | None:
    current_index = int(indices[column_index])
    if current_index >= 0:
        return current_index

    for scan_index in range(column_index + 1, len(indices)):
        value = int(indices[scan_index])
        if value >= 0:
            return value

    for scan_index in range(column_index - 1, -1, -1):
        value = int(indices[scan_index])
        if value >= 0:
            return value

    return None


def _context_window(sequence: str, index: int | None, *, window: int) -> str:
    if index is None:
        return ""

    left = max(index - window, 0)
    right = min(index + window + 1, len(sequence))
    return sequence[left:right]


__all__ = [
    "StrandPolicy",
    "ChosenStrand",
    "AlignmentEvent",
    "AlignmentResult",
    "normalize_reference",
    "build_aligner",
    "align_trimmed_read_to_reference",
    "alignment_events_to_rows",
    "format_alignment_block",
]
