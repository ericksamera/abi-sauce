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

AssemblyStrand = Literal["forward", "reverse-complement"]
ConflictResolution = Literal[
    "concordant",
    "single_read",
    "quality_resolved",
    "ambiguous",
]


@dataclass(frozen=True, slots=True)
class AssemblyConfig:
    """Scoring and acceptance thresholds for pairwise read assembly."""

    match_score: float = 1.0
    mismatch_score: float = -1.0
    open_internal_gap_score: float = -3.0
    extend_internal_gap_score: float = -1.0
    min_overlap_length: int = 25
    min_percent_identity: float = 70.0
    quality_margin: int = 3


@dataclass(frozen=True, slots=True)
class AssemblyConflict:
    """One consensus-support column worth surfacing in the UI."""

    column_index: int
    left_base: str
    right_base: str
    consensus_base: str
    resolution: ConflictResolution
    left_query_pos: int | None
    right_query_pos: int | None
    left_quality: int | None
    right_quality: int | None
    left_trace_x: int | None
    right_trace_x: int | None


@dataclass(frozen=True, slots=True)
class AssemblyResult:
    """Derived pairwise assembly result for two trimmed reads."""

    left_source_filename: str
    right_source_filename: str
    left_display_name: str
    right_display_name: str
    chosen_right_orientation: AssemblyStrand
    score: float
    accepted: bool
    rejection_reason: str | None
    overlap_length: int
    percent_identity: float
    mismatch_count: int
    aligned_left: str
    match_line: str
    aligned_right: str
    gapped_consensus: str
    consensus_sequence: str
    conflicts: tuple[AssemblyConflict, ...] = ()

    @property
    def conflict_count(self) -> int:
        return len(self.conflicts)


def build_assembly_aligner(
    *,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    open_internal_gap_score: float = -3.0,
    extend_internal_gap_score: float = -1.0,
) -> Align.PairwiseAligner:
    """Build a semi-global style aligner for pairwise assembly."""
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


def assemble_trimmed_pair(
    *,
    left_source_filename: str,
    left_raw_record: SequenceRecord,
    left_trim_result: TrimResult,
    right_source_filename: str,
    right_raw_record: SequenceRecord,
    right_trim_result: TrimResult,
    config: AssemblyConfig | None = None,
    aligner: Align.PairwiseAligner | None = None,
) -> AssemblyResult:
    """Assemble two currently trimmed reads into one consensus result."""
    resolved_config = AssemblyConfig() if config is None else config
    resolved_aligner = (
        build_assembly_aligner(
            match_score=resolved_config.match_score,
            mismatch_score=resolved_config.mismatch_score,
            open_internal_gap_score=resolved_config.open_internal_gap_score,
            extend_internal_gap_score=resolved_config.extend_internal_gap_score,
        )
        if aligner is None
        else aligner
    )

    left_display_record = materialize_oriented_record(left_trim_result.record)
    right_display_record = materialize_oriented_record(right_trim_result.record)
    left_sequence = left_display_record.sequence.upper()
    right_display_sequence = right_display_record.sequence.upper()

    if not left_sequence:
        return _empty_assembly_result(
            left_source_filename=left_source_filename,
            right_source_filename=right_source_filename,
            left_display_name=_display_name(
                left_trim_result.record, left_source_filename
            ),
            right_display_name=_display_name(
                right_trim_result.record, right_source_filename
            ),
            chosen_right_orientation="forward",
            rejection_reason="left read trims down to an empty sequence",
        )
    if not right_display_sequence:
        return _empty_assembly_result(
            left_source_filename=left_source_filename,
            right_source_filename=right_source_filename,
            left_display_name=_display_name(
                left_trim_result.record, left_source_filename
            ),
            right_display_name=_display_name(
                right_trim_result.record, right_source_filename
            ),
            chosen_right_orientation="forward",
            rejection_reason="right read trims down to an empty sequence",
        )

    candidates: list[tuple[AssemblyStrand, str, list[int] | None]] = [
        (
            "forward",
            right_display_sequence,
            _candidate_member_qualities(
                right_display_record.qualities,
                strand="forward",
            ),
        ),
        (
            "reverse-complement",
            reverse_complement_sequence(right_display_sequence),
            _candidate_member_qualities(
                right_display_record.qualities,
                strand="reverse-complement",
            ),
        ),
    ]

    best_alignment = None
    best_right_sequence = ""
    best_right_qualities = None
    best_right_orientation: AssemblyStrand = "forward"
    best_score = float("-inf")

    for (
        right_orientation,
        oriented_right_sequence,
        oriented_right_qualities,
    ) in candidates:
        alignment = resolved_aligner.align(left_sequence, oriented_right_sequence)[0]
        if alignment.score > best_score:
            best_alignment = alignment
            best_right_sequence = oriented_right_sequence
            best_right_qualities = oriented_right_qualities
            best_right_orientation = right_orientation
            best_score = float(alignment.score)

    if best_alignment is None:
        return _empty_assembly_result(
            left_source_filename=left_source_filename,
            right_source_filename=right_source_filename,
            left_display_name=_display_name(
                left_trim_result.record, left_source_filename
            ),
            right_display_name=_display_name(
                right_trim_result.record, right_source_filename
            ),
            chosen_right_orientation="forward",
            rejection_reason="no alignment could be generated",
        )

    result = _extract_assembly_result(
        alignment=best_alignment,
        left_source_filename=left_source_filename,
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        left_oriented_sequence=left_sequence,
        left_oriented_qualities=_candidate_member_qualities(
            left_display_record.qualities,
            strand="forward",
        ),
        right_source_filename=right_source_filename,
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        right_oriented_sequence=best_right_sequence,
        right_oriented_qualities=best_right_qualities,
        chosen_right_orientation=best_right_orientation,
        quality_margin=resolved_config.quality_margin,
        score=best_score,
    )

    if result.overlap_length <= 0:
        return _replace_result_acceptance(
            result,
            accepted=False,
            rejection_reason="no overlapping aligned columns were found",
        )
    if result.overlap_length < resolved_config.min_overlap_length:
        return _replace_result_acceptance(
            result,
            accepted=False,
            rejection_reason=(
                "overlap length below threshold "
                f"({result.overlap_length} < {resolved_config.min_overlap_length})"
            ),
        )
    if result.percent_identity < resolved_config.min_percent_identity:
        return _replace_result_acceptance(
            result,
            accepted=False,
            rejection_reason=(
                "percent identity below threshold "
                f"({result.percent_identity:.1f}% < "
                f"{resolved_config.min_percent_identity:.1f}%)"
            ),
        )

    return result


def assembly_conflicts_to_rows(
    result: AssemblyResult,
) -> list[dict[str, object]]:
    """Convert assembly conflicts into Streamlit-table-friendly rows."""
    return [
        {
            "column": conflict.column_index,
            "left_base": conflict.left_base,
            "right_base": conflict.right_base,
            "consensus_base": conflict.consensus_base,
            "resolution": conflict.resolution,
            "left_query_pos": conflict.left_query_pos,
            "right_query_pos": conflict.right_query_pos,
            "left_quality": conflict.left_quality,
            "right_quality": conflict.right_quality,
            "left_trace_x": conflict.left_trace_x,
            "right_trace_x": conflict.right_trace_x,
        }
        for conflict in result.conflicts
    ]


def format_assembly_block(
    result: AssemblyResult,
    *,
    line_width: int = 80,
) -> str:
    """Render a four-line gapped assembly block for display."""
    lines: list[str] = []
    for start in range(0, len(result.aligned_left), line_width):
        end = start + line_width
        lines.append(result.aligned_left[start:end])
        lines.append(result.match_line[start:end])
        lines.append(result.aligned_right[start:end])
        lines.append(result.gapped_consensus[start:end])
        lines.append("")
    return "\n".join(lines).rstrip()


def consensus_record_from_result(
    result: AssemblyResult,
    *,
    name: str | None = None,
) -> SequenceRecord:
    """Project one assembly result into an exportable consensus SequenceRecord."""
    consensus_name = (
        name.strip()
        if isinstance(name, str) and name.strip()
        else f"{result.left_display_name}__{result.right_display_name}_consensus"
    )
    return SequenceRecord(
        record_id=consensus_name,
        name=consensus_name,
        description=(
            f"assembly consensus from {result.left_display_name} and "
            f"{result.right_display_name}"
        ),
        sequence=result.consensus_sequence,
        source_format="assembly",
        orientation="forward",
        qualities=None,
        trace_data=None,
        annotations={
            "assembly_left_source_filename": result.left_source_filename,
            "assembly_right_source_filename": result.right_source_filename,
            "assembly_right_orientation": result.chosen_right_orientation,
            "assembly_overlap_length": result.overlap_length,
            "assembly_percent_identity": result.percent_identity,
            "assembly_conflict_count": result.conflict_count,
            "assembly_accepted": result.accepted,
            "assembly_rejection_reason": result.rejection_reason,
        },
    )


def _extract_assembly_result(
    *,
    alignment,
    left_source_filename: str,
    left_raw_record: SequenceRecord,
    left_trim_result: TrimResult,
    left_oriented_sequence: str,
    left_oriented_qualities: list[int] | None,
    right_source_filename: str,
    right_raw_record: SequenceRecord,
    right_trim_result: TrimResult,
    right_oriented_sequence: str,
    right_oriented_qualities: list[int] | None,
    chosen_right_orientation: AssemblyStrand,
    quality_margin: int,
    score: float,
) -> AssemblyResult:
    target_indices = alignment.indices[0]
    query_indices = alignment.indices[1]

    aligned_left_chars: list[str] = []
    aligned_right_chars: list[str] = []
    match_line_chars: list[str] = []
    gapped_consensus_chars: list[str] = []
    conflicts: list[AssemblyConflict] = []

    overlap_length = 0
    match_count = 0
    mismatch_count = 0

    for column_index, (target_index, query_index) in enumerate(
        zip(target_indices, query_indices, strict=True),
        start=1,
    ):
        resolved_target_index = int(target_index)
        resolved_query_index = int(query_index)

        left_base = (
            "-"
            if resolved_target_index < 0
            else left_oriented_sequence[resolved_target_index]
        )
        right_base = (
            "-"
            if resolved_query_index < 0
            else right_oriented_sequence[resolved_query_index]
        )
        aligned_left_chars.append(left_base)
        aligned_right_chars.append(right_base)

        left_quality = (
            None
            if resolved_target_index < 0 or left_oriented_qualities is None
            else int(left_oriented_qualities[resolved_target_index])
        )
        right_quality = (
            None
            if resolved_query_index < 0 or right_oriented_qualities is None
            else int(right_oriented_qualities[resolved_query_index])
        )
        left_query_pos = (
            None if resolved_target_index < 0 else resolved_target_index + 1
        )
        right_query_pos = None if resolved_query_index < 0 else resolved_query_index + 1

        consensus_base, resolution = _resolve_consensus_column(
            left_base=left_base,
            right_base=right_base,
            left_quality=left_quality,
            right_quality=right_quality,
            quality_margin=quality_margin,
        )
        gapped_consensus_chars.append(
            "-" if left_base == "-" and right_base == "-" else consensus_base
        )

        if left_base != "-" and right_base != "-":
            overlap_length += 1
            if left_base == right_base:
                match_count += 1
                match_line_chars.append("|")
            else:
                mismatch_count += 1
                match_line_chars.append(".")
        else:
            match_line_chars.append(" ")

        if resolution != "concordant":
            conflicts.append(
                AssemblyConflict(
                    column_index=column_index,
                    left_base=left_base,
                    right_base=right_base,
                    consensus_base=consensus_base,
                    resolution=resolution,
                    left_query_pos=left_query_pos,
                    right_query_pos=right_query_pos,
                    left_quality=left_quality,
                    right_quality=right_quality,
                    left_trace_x=_trace_position_for_member_query_index(
                        raw_record=left_raw_record,
                        trim_result=left_trim_result,
                        oriented_query_index=(
                            None if resolved_target_index < 0 else resolved_target_index
                        ),
                        strand="forward",
                    ),
                    right_trace_x=_trace_position_for_member_query_index(
                        raw_record=right_raw_record,
                        trim_result=right_trim_result,
                        oriented_query_index=(
                            None if resolved_query_index < 0 else resolved_query_index
                        ),
                        strand=chosen_right_orientation,
                    ),
                )
            )

    percent_identity = (match_count / overlap_length) * 100.0 if overlap_length else 0.0

    return AssemblyResult(
        left_source_filename=left_source_filename,
        right_source_filename=right_source_filename,
        left_display_name=_display_name(left_trim_result.record, left_source_filename),
        right_display_name=_display_name(
            right_trim_result.record, right_source_filename
        ),
        chosen_right_orientation=chosen_right_orientation,
        score=score,
        accepted=True,
        rejection_reason=None,
        overlap_length=overlap_length,
        percent_identity=percent_identity,
        mismatch_count=mismatch_count,
        aligned_left="".join(aligned_left_chars),
        match_line="".join(match_line_chars),
        aligned_right="".join(aligned_right_chars),
        gapped_consensus="".join(gapped_consensus_chars),
        consensus_sequence="".join(
            base for base in gapped_consensus_chars if base != "-"
        ),
        conflicts=tuple(conflicts),
    )


def _resolve_consensus_column(
    *,
    left_base: str,
    right_base: str,
    left_quality: int | None,
    right_quality: int | None,
    quality_margin: int,
) -> tuple[str, ConflictResolution]:
    if left_base != "-" and right_base != "-":
        if left_base == right_base:
            return left_base, "concordant"
        if (
            left_quality is not None
            and right_quality is not None
            and abs(left_quality - right_quality) >= quality_margin
        ):
            return (
                (left_base, "quality_resolved")
                if left_quality > right_quality
                else (right_base, "quality_resolved")
            )
        return "N", "ambiguous"

    if left_base != "-":
        return left_base, "single_read"
    if right_base != "-":
        return right_base, "single_read"
    return "N", "ambiguous"


def _replace_result_acceptance(
    result: AssemblyResult,
    *,
    accepted: bool,
    rejection_reason: str | None,
) -> AssemblyResult:
    return AssemblyResult(
        left_source_filename=result.left_source_filename,
        right_source_filename=result.right_source_filename,
        left_display_name=result.left_display_name,
        right_display_name=result.right_display_name,
        chosen_right_orientation=result.chosen_right_orientation,
        score=result.score,
        accepted=accepted,
        rejection_reason=rejection_reason,
        overlap_length=result.overlap_length,
        percent_identity=result.percent_identity,
        mismatch_count=result.mismatch_count,
        aligned_left=result.aligned_left,
        match_line=result.match_line,
        aligned_right=result.aligned_right,
        gapped_consensus=result.gapped_consensus,
        consensus_sequence=result.consensus_sequence,
        conflicts=result.conflicts,
    )


def _empty_assembly_result(
    *,
    left_source_filename: str,
    right_source_filename: str,
    left_display_name: str,
    right_display_name: str,
    chosen_right_orientation: AssemblyStrand,
    rejection_reason: str,
) -> AssemblyResult:
    return AssemblyResult(
        left_source_filename=left_source_filename,
        right_source_filename=right_source_filename,
        left_display_name=left_display_name,
        right_display_name=right_display_name,
        chosen_right_orientation=chosen_right_orientation,
        score=float("-inf"),
        accepted=False,
        rejection_reason=rejection_reason,
        overlap_length=0,
        percent_identity=0.0,
        mismatch_count=0,
        aligned_left="",
        match_line="",
        aligned_right="",
        gapped_consensus="",
        consensus_sequence="",
        conflicts=(),
    )


def _display_name(record: SequenceRecord, source_filename: str) -> str:
    return record.name or record.record_id or source_filename


def _candidate_member_qualities(
    display_qualities: list[int] | None,
    *,
    strand: AssemblyStrand,
) -> list[int] | None:
    if display_qualities is None:
        return None
    if strand == "forward":
        return [int(value) for value in display_qualities]
    return [int(value) for value in reversed(display_qualities)]


def _trace_position_for_member_query_index(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    oriented_query_index: int | None,
    strand: AssemblyStrand,
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
    strand: AssemblyStrand,
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
