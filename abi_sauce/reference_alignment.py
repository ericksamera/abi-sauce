from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from Bio import Align

from abi_sauce.models import SequenceOrientation, SequenceRecord
from abi_sauce.orientation import (
    materialize_oriented_record,
    reverse_complement_sequence,
)
from abi_sauce.trimming import TrimResult

StrandPolicy = Literal["auto", "forward", "reverse-complement"]
ChosenStrand = Literal["forward", "reverse-complement"]

_ALLOWED_REFERENCE_CHARACTERS = frozenset("ACGTURYKMSWBDHVN")


@dataclass(frozen=True, slots=True)
class AlignmentEvent:
    ref_pos: int | None
    query_pos: int | None
    event_type: str
    ref_base: str
    query_base: str
    qscore: int | None
    flank_q_left: int | None
    flank_q_right: int | None
    trace_x: int | None
    context_ref: str
    context_query: str


@dataclass(frozen=True, slots=True)
class AlignmentResult:
    sample_name: str
    reference_name: str
    strand: ChosenStrand
    score: float
    reference_start: int | None
    reference_end: int | None
    query_start: int | None
    query_end: int | None
    percent_identity: float
    mismatch_count: int
    insertion_count: int
    deletion_count: int
    aligned_reference: str
    match_line: str
    aligned_query: str
    events: tuple[AlignmentEvent, ...]


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
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score

    aligner.open_internal_insertion_score = open_internal_gap_score
    aligner.extend_internal_insertion_score = extend_internal_gap_score
    aligner.open_internal_deletion_score = open_internal_gap_score
    aligner.extend_internal_deletion_score = extend_internal_gap_score

    aligner.open_left_insertion_score = 0.0
    aligner.extend_left_insertion_score = 0.0
    aligner.open_right_insertion_score = 0.0
    aligner.extend_right_insertion_score = 0.0
    aligner.open_left_deletion_score = 0.0
    aligner.extend_left_deletion_score = 0.0
    aligner.open_right_deletion_score = 0.0
    aligner.extend_right_deletion_score = 0.0

    return aligner


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
    if strand_policy not in {"auto", "forward", "reverse-complement"}:
        raise ValueError(f"Unsupported strand policy: {strand_policy}")

    normalized_reference_name, reference_sequence = normalize_reference(reference_text)
    resolved_reference_name = (
        normalized_reference_name if reference_name is None else reference_name
    )

    display_trimmed_record = materialize_oriented_record(trim_result.record)
    trimmed_sequence = display_trimmed_record.sequence.upper()
    if not trimmed_sequence:
        raise ValueError("Trimmed query sequence is empty.")

    resolved_aligner = build_aligner() if aligner is None else aligner

    candidate_queries: list[tuple[ChosenStrand, str, list[int] | None]] = []
    if strand_policy in {"auto", "forward"}:
        candidate_queries.append(
            (
                "forward",
                trimmed_sequence,
                _candidate_query_qualities(
                    display_trimmed_record.qualities,
                    strand="forward",
                ),
            )
        )
    if strand_policy in {"auto", "reverse-complement"}:
        candidate_queries.append(
            (
                "reverse-complement",
                reverse_complement_sequence(trimmed_sequence),
                _candidate_query_qualities(
                    display_trimmed_record.qualities,
                    strand="reverse-complement",
                ),
            )
        )

    best_alignment = None
    best_query_sequence = ""
    best_query_qualities = None
    chosen_strand: ChosenStrand = "forward"
    best_score = float("-inf")

    for strand, oriented_query_sequence, oriented_query_qualities in candidate_queries:
        alignment = resolved_aligner.align(reference_sequence, oriented_query_sequence)[
            0
        ]
        if alignment.score > best_score:
            best_alignment = alignment
            best_query_sequence = oriented_query_sequence
            best_query_qualities = oriented_query_qualities
            chosen_strand = strand
            best_score = float(alignment.score)

    if best_alignment is None:
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
        alignment=best_alignment,
        raw_record=raw_record,
        trim_result=trim_result,
        reference_sequence=reference_sequence,
        oriented_query_sequence=best_query_sequence,
        oriented_query_qualities=best_query_qualities,
        strand=chosen_strand,
        context_window=context_window,
    )

    aligned_columns = matches + mismatches + insertions + deletions
    percent_identity = (matches / aligned_columns * 100.0) if aligned_columns else 0.0

    return AlignmentResult(
        sample_name=trim_result.record.name or trim_result.record.record_id,
        reference_name=resolved_reference_name,
        strand=chosen_strand,
        score=best_score,
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


def alignment_events_to_rows(
    result: AlignmentResult,
    *,
    include_matches: bool = False,
) -> list[dict[str, object]]:
    """Convert alignment events into Streamlit-table-friendly rows."""
    rows: list[dict[str, object]] = []
    for event in result.events:
        if not include_matches and event.event_type == "match":
            continue
        rows.append(
            {
                "ref_pos": event.ref_pos,
                "query_pos": event.query_pos,
                "type": event.event_type,
                "ref_base": event.ref_base,
                "query_base": event.query_base,
                "qscore": event.qscore,
                "flank_q_left": event.flank_q_left,
                "flank_q_right": event.flank_q_right,
                "trace_x": event.trace_x,
                "context_ref": event.context_ref,
                "context_query": event.context_query,
            }
        )
    return rows


def format_alignment_block(
    result: AlignmentResult,
    *,
    line_width: int = 80,
) -> str:
    """Render a three-line gapped alignment block for display."""
    lines: list[str] = []
    for start in range(0, len(result.aligned_reference), line_width):
        end = start + line_width
        lines.append(result.aligned_reference[start:end])
        lines.append(result.match_line[start:end])
        lines.append(result.aligned_query[start:end])
        lines.append("")
    return "\n".join(lines).rstrip()


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
        trace_x = _trace_position_for_query_index(
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


def _candidate_query_qualities(
    display_qualities: list[int] | None,
    *,
    strand: ChosenStrand,
) -> list[int] | None:
    if display_qualities is None:
        return None
    if strand == "forward":
        return [int(value) for value in display_qualities]
    return [int(value) for value in reversed(display_qualities)]


def _trace_position_for_query_index(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    oriented_query_index: int | None,
    strand: ChosenStrand,
) -> int | None:
    if oriented_query_index is None:
        return None

    trace_data = raw_record.trace_data
    if trace_data is None:
        return None

    raw_start = trim_result.bases_removed_left
    trimmed_length = trim_result.trimmed_length
    if trimmed_length <= 0:
        return None

    display_query_index = _display_query_index_for_aligned_query_index(
        oriented_query_index,
        trimmed_length=trimmed_length,
        strand=strand,
    )
    raw_trimmed_index = _raw_trimmed_index_for_display_query_index(
        display_query_index,
        trimmed_length=trimmed_length,
        display_orientation=raw_record.orientation,
    )
    raw_base_index = raw_start + raw_trimmed_index

    if raw_base_index < 0 or raw_base_index >= len(trace_data.base_positions):
        return None

    raw_trace_position = int(trace_data.base_positions[raw_base_index])
    trace_length = _sanitized_trace_length(trace_data)
    if trace_length <= 0:
        return raw_trace_position
    if raw_trace_position < 0 or raw_trace_position >= trace_length:
        return None
    return _display_trace_position(
        raw_trace_position,
        trace_length=trace_length,
        display_orientation=raw_record.orientation,
    )


def _display_query_index_for_aligned_query_index(
    aligned_query_index: int,
    *,
    trimmed_length: int,
    strand: ChosenStrand,
) -> int:
    if strand == "forward":
        return aligned_query_index
    return trimmed_length - 1 - aligned_query_index


def _raw_trimmed_index_for_display_query_index(
    display_query_index: int,
    *,
    trimmed_length: int,
    display_orientation: SequenceOrientation,
) -> int:
    if display_orientation == "forward":
        return display_query_index
    return trimmed_length - 1 - display_query_index


def _sanitized_trace_length(trace_data) -> int:
    channel_lengths = [
        len(signal) for signal in trace_data.channels.values() if len(signal) > 0
    ]
    if not channel_lengths:
        return 0
    return min(channel_lengths)


def _display_trace_position(
    raw_trace_position: int,
    *,
    trace_length: int,
    display_orientation: SequenceOrientation,
) -> int:
    if display_orientation == "forward":
        return raw_trace_position
    return int(float(trace_length - 1) - float(raw_trace_position))


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
