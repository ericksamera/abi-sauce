from __future__ import annotations

from collections import Counter
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
    "majority_resolved",
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
class AssemblyColumn:
    """One gapped alignment column carrying enough state for downstream views."""

    column_index: int
    left_base: str
    right_base: str
    consensus_base: str
    resolution: ConflictResolution
    left_query_index: int | None
    right_query_index: int | None
    left_query_pos: int | None
    right_query_pos: int | None
    left_quality: int | None
    right_quality: int | None
    left_trace_x: int | None
    right_trace_x: int | None
    is_overlap: bool
    is_match: bool


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
    columns: tuple[AssemblyColumn, ...] = ()
    conflicts: tuple[AssemblyConflict, ...] = ()

    @property
    def conflict_count(self) -> int:
        return len(self.conflicts)


@dataclass(frozen=True, slots=True)
class MultiAssemblyMember:
    """One member selected for a multi-read assembly attempt."""

    member_index: int
    source_filename: str
    display_name: str
    chosen_orientation: AssemblyStrand
    included: bool
    inclusion_reason: str | None
    trimmed_length: int
    has_trace_data: bool
    has_qualities: bool
    is_seed: bool = False
    alignment_score: float | None = None
    overlap_length: int | None = None
    percent_identity: float | None = None


@dataclass(frozen=True, slots=True)
class MultiAssemblyMemberCell:
    """One included-member cell projected onto the multi-read alignment grid."""

    member_index: int
    base: str
    query_index: int | None
    query_pos: int | None
    quality: int | None
    trace_x: int | None
    is_gap: bool


@dataclass(frozen=True, slots=True)
class MultiAssemblyColumn:
    """One global alignment column for a multi-read assembly."""

    column_index: int
    consensus_base: str
    resolution: ConflictResolution
    member_cells: tuple[MultiAssemblyMemberCell, ...]
    support_counts: tuple[tuple[str, int], ...] = ()
    quality_sums: tuple[tuple[str, int], ...] = ()
    non_gap_member_count: int = 0
    gap_member_count: int = 0
    ambiguous: bool = False


@dataclass(frozen=True, slots=True)
class MultiAssemblyResult:
    """Derived anchored multi-read assembly result."""

    members: tuple[MultiAssemblyMember, ...]
    seed_member_index: int
    accepted: bool
    rejection_reason: str | None
    included_member_indices: tuple[int, ...]
    aligned_member_sequences: tuple[str, ...]
    gapped_consensus: str
    consensus_sequence: str
    columns: tuple[MultiAssemblyColumn, ...] = ()
    included_member_count: int = 0
    excluded_member_count: int = 0
    ambiguous_column_count: int = 0

    @property
    def conflict_count(self) -> int:
        return self.ambiguous_column_count


AssemblyComputationResult = AssemblyResult | MultiAssemblyResult


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


def format_assembly_alignment_fasta(
    result: AssemblyComputationResult,
    *,
    consensus_name: str | None = None,
) -> str:
    """Render one assembly alignment plus gapped consensus as FASTA."""
    lines: list[str] = []

    if isinstance(result, MultiAssemblyResult):
        members_by_index = {
            member.member_index: member for member in result.members if member.included
        }
        for member_index, aligned_sequence in zip(
            result.included_member_indices,
            result.aligned_member_sequences,
            strict=True,
        ):
            member = members_by_index[member_index]
            lines.append(f">{member.display_name}")
            lines.append(aligned_sequence)

        seed_member = next(
            (
                member
                for member in result.members
                if member.member_index == result.seed_member_index
            ),
            None,
        )
        resolved_consensus_name = (
            consensus_name.strip()
            if isinstance(consensus_name, str) and consensus_name.strip()
            else (
                f"{seed_member.display_name}_multi_consensus"
                if seed_member is not None
                else "multi_consensus"
            )
        )
    else:
        lines.append(f">{result.left_display_name}")
        lines.append(result.aligned_left)
        lines.append(f">{result.right_display_name}")
        lines.append(result.aligned_right)
        resolved_consensus_name = (
            consensus_name.strip()
            if isinstance(consensus_name, str) and consensus_name.strip()
            else f"{result.left_display_name}__{result.right_display_name}_consensus"
        )

    lines.append(f">{resolved_consensus_name}__gapped_consensus")
    lines.append(result.gapped_consensus)
    return "\n".join(lines) + "\n"


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
    columns: list[AssemblyColumn] = []
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

        left_trace_x = _trace_position_for_member_query_index(
            raw_record=left_raw_record,
            trim_result=left_trim_result,
            oriented_query_index=(
                None if resolved_target_index < 0 else resolved_target_index
            ),
            strand="forward",
        )
        right_trace_x = _trace_position_for_member_query_index(
            raw_record=right_raw_record,
            trim_result=right_trim_result,
            oriented_query_index=(
                None if resolved_query_index < 0 else resolved_query_index
            ),
            strand=chosen_right_orientation,
        )
        column = AssemblyColumn(
            column_index=column_index,
            left_base=left_base,
            right_base=right_base,
            consensus_base=consensus_base,
            resolution=resolution,
            left_query_index=(
                None if resolved_target_index < 0 else resolved_target_index
            ),
            right_query_index=(
                None if resolved_query_index < 0 else resolved_query_index
            ),
            left_query_pos=left_query_pos,
            right_query_pos=right_query_pos,
            left_quality=left_quality,
            right_quality=right_quality,
            left_trace_x=left_trace_x,
            right_trace_x=right_trace_x,
            is_overlap=left_base != "-" and right_base != "-",
            is_match=left_base != "-" and right_base != "-" and left_base == right_base,
        )
        columns.append(column)
        if resolution != "concordant":
            conflicts.append(_conflict_from_column(column))

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
        columns=tuple(columns),
        conflicts=tuple(conflicts),
    )


def _conflict_from_column(column: AssemblyColumn) -> AssemblyConflict:
    return AssemblyConflict(
        column_index=column.column_index,
        left_base=column.left_base,
        right_base=column.right_base,
        consensus_base=column.consensus_base,
        resolution=column.resolution,
        left_query_pos=column.left_query_pos,
        right_query_pos=column.right_query_pos,
        left_quality=column.left_quality,
        right_quality=column.right_quality,
        left_trace_x=column.left_trace_x,
        right_trace_x=column.right_trace_x,
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
        columns=result.columns,
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
        columns=(),
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


@dataclass(frozen=True, slots=True)
class _ResolvedMultiMemberInput:
    member_index: int
    source_filename: str
    raw_record: SequenceRecord
    trim_result: TrimResult
    display_name: str
    display_sequence: str
    display_qualities: list[int] | None
    trimmed_length: int
    has_trace_data: bool
    has_qualities: bool


@dataclass(frozen=True, slots=True)
class _ProjectedMultiCell:
    base: str
    query_index: int | None
    query_pos: int | None
    quality: int | None
    trace_x: int | None
    is_gap: bool


@dataclass(frozen=True, slots=True)
class _PlacedMultiMember:
    member: MultiAssemblyMember
    seed_column_cells: tuple[_ProjectedMultiCell, ...]
    insertion_cells_by_bucket: dict[int, tuple[_ProjectedMultiCell, ...]]


def assemble_trimmed_multi(
    *,
    source_filenames: tuple[str, ...] | list[str],
    raw_records_by_source_filename: dict[str, SequenceRecord],
    trim_results_by_source_filename: dict[str, TrimResult],
    config: AssemblyConfig | None = None,
    aligner: Align.PairwiseAligner | None = None,
) -> MultiAssemblyResult:
    """Assemble multiple currently trimmed reads onto one seed-aligned grid."""
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

    member_inputs = tuple(
        _resolved_multi_member_input(
            member_index=member_index,
            source_filename=source_filename,
            raw_record=raw_records_by_source_filename[source_filename],
            trim_result=trim_results_by_source_filename[source_filename],
        )
        for member_index, source_filename in enumerate(source_filenames)
    )
    non_empty_member_inputs = tuple(
        member_input for member_input in member_inputs if member_input.display_sequence
    )
    if not non_empty_member_inputs:
        return MultiAssemblyResult(
            members=tuple(
                MultiAssemblyMember(
                    member_index=member_input.member_index,
                    source_filename=member_input.source_filename,
                    display_name=member_input.display_name,
                    chosen_orientation="forward",
                    included=False,
                    inclusion_reason="trimmed read is empty",
                    trimmed_length=member_input.trimmed_length,
                    has_trace_data=member_input.has_trace_data,
                    has_qualities=member_input.has_qualities,
                )
                for member_input in member_inputs
            ),
            seed_member_index=0,
            accepted=False,
            rejection_reason="all selected reads trim down to empty sequences",
            included_member_indices=(),
            aligned_member_sequences=(),
            gapped_consensus="",
            consensus_sequence="",
            columns=(),
            included_member_count=0,
            excluded_member_count=len(member_inputs),
            ambiguous_column_count=0,
        )

    seed_member_input = _select_multi_seed_member(non_empty_member_inputs)
    seed_sequence = seed_member_input.display_sequence

    placements: list[_PlacedMultiMember] = [_build_seed_multi_member(seed_member_input)]
    for member_input in member_inputs:
        if member_input.member_index == seed_member_input.member_index:
            continue
        placements.append(
            _place_member_against_multi_seed(
                seed_member_input=seed_member_input,
                member_input=member_input,
                config=resolved_config,
                aligner=resolved_aligner,
            )
        )

    placements_by_index = {
        placement.member.member_index: placement for placement in placements
    }
    ordered_placements = tuple(
        placements_by_index.get(member_input.member_index)
        or _excluded_multi_member_for_empty_sequence(member_input)
        for member_input in member_inputs
    )
    included_placements = tuple(
        placement for placement in ordered_placements if placement.member.included
    )
    included_member_indices = tuple(
        placement.member.member_index for placement in included_placements
    )

    bucket_columns_by_index = {
        bucket_index: _resolve_multi_insertion_bucket(
            bucket_cells_by_member={
                placement.member.member_index: placement.insertion_cells_by_bucket.get(
                    bucket_index, ()
                )
                for placement in included_placements
            },
            included_member_indices=included_member_indices,
            aligner=resolved_aligner,
        )
        for bucket_index in range(len(seed_sequence) + 1)
    }

    aligned_sequences_by_member_index: dict[int, list[str]] = {
        member_index: [] for member_index in included_member_indices
    }
    columns: list[MultiAssemblyColumn] = []
    column_index = 1

    for bucket_index in range(len(seed_sequence) + 1):
        for bucket_column in bucket_columns_by_index[bucket_index]:
            column = _build_multi_assembly_column(
                column_index=column_index,
                column_cells_by_member=bucket_column,
                included_member_indices=included_member_indices,
                quality_margin=resolved_config.quality_margin,
            )
            if column is None:
                continue
            columns.append(column)
            for member_cell in column.member_cells:
                aligned_sequences_by_member_index[member_cell.member_index].append(
                    member_cell.base
                )
            column_index += 1
        if bucket_index >= len(seed_sequence):
            continue
        column = _build_multi_assembly_column(
            column_index=column_index,
            column_cells_by_member={
                placement.member.member_index: placement.seed_column_cells[bucket_index]
                for placement in included_placements
            },
            included_member_indices=included_member_indices,
            quality_margin=resolved_config.quality_margin,
        )
        if column is None:
            continue
        columns.append(column)
        for member_cell in column.member_cells:
            aligned_sequences_by_member_index[member_cell.member_index].append(
                member_cell.base
            )
        column_index += 1

    gapped_consensus = "".join(column.consensus_base for column in columns)
    consensus_sequence = "".join(base for base in gapped_consensus if base != "-")
    ambiguous_column_count = sum(column.ambiguous for column in columns)
    included_member_count = len(included_member_indices)
    rejected_reason = None
    accepted = True
    if included_member_count < 2:
        accepted = False
        rejected_reason = "fewer than 2 members passed placement thresholds"
    elif not columns or not consensus_sequence:
        accepted = False
        rejected_reason = "no consensus columns could be derived"

    return MultiAssemblyResult(
        members=tuple(placement.member for placement in ordered_placements),
        seed_member_index=seed_member_input.member_index,
        accepted=accepted,
        rejection_reason=rejected_reason,
        included_member_indices=included_member_indices,
        aligned_member_sequences=tuple(
            "".join(aligned_sequences_by_member_index[member_index])
            for member_index in included_member_indices
        ),
        gapped_consensus=gapped_consensus,
        consensus_sequence=consensus_sequence,
        columns=tuple(columns),
        included_member_count=included_member_count,
        excluded_member_count=len(member_inputs) - included_member_count,
        ambiguous_column_count=ambiguous_column_count,
    )


def consensus_record_from_multi_result(
    result: MultiAssemblyResult,
    *,
    name: str | None = None,
) -> SequenceRecord:
    """Project one multi-read assembly result into an exportable consensus record."""
    seed_member = next(
        (
            member
            for member in result.members
            if member.member_index == result.seed_member_index
        ),
        None,
    )
    default_name = (
        f"{seed_member.display_name}_multi_consensus"
        if seed_member is not None
        else "multi_consensus"
    )
    consensus_name = (
        name.strip() if isinstance(name, str) and name.strip() else default_name
    )
    return SequenceRecord(
        record_id=consensus_name,
        name=consensus_name,
        description=(
            f"multi-read assembly consensus from {result.included_member_count} "
            "included reads"
        ),
        sequence=result.consensus_sequence,
        source_format="assembly",
        orientation="forward",
        qualities=None,
        trace_data=None,
        annotations={
            "assembly_engine_kind": "multi",
            "assembly_member_filenames": [
                member.source_filename for member in result.members
            ],
            "assembly_seed_source_filename": (
                None if seed_member is None else seed_member.source_filename
            ),
            "assembly_included_member_count": result.included_member_count,
            "assembly_excluded_member_count": result.excluded_member_count,
            "assembly_ambiguous_column_count": result.ambiguous_column_count,
            "assembly_accepted": result.accepted,
            "assembly_rejection_reason": result.rejection_reason,
        },
    )


def _resolved_multi_member_input(
    *,
    member_index: int,
    source_filename: str,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
) -> _ResolvedMultiMemberInput:
    display_record = materialize_oriented_record(trim_result.record)
    return _ResolvedMultiMemberInput(
        member_index=member_index,
        source_filename=source_filename,
        raw_record=raw_record,
        trim_result=trim_result,
        display_name=_display_name(trim_result.record, source_filename),
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


def _select_multi_seed_member(
    member_inputs: tuple[_ResolvedMultiMemberInput, ...],
) -> _ResolvedMultiMemberInput:
    return sorted(
        member_inputs,
        key=lambda member_input: (
            -len(member_input.display_sequence),
            -_mean_quality(member_input.display_qualities),
            member_input.source_filename,
        ),
    )[0]


def _build_seed_multi_member(
    seed_member_input: _ResolvedMultiMemberInput,
) -> _PlacedMultiMember:
    oriented_qualities = _candidate_member_qualities(
        seed_member_input.display_qualities,
        strand="forward",
    )
    seed_column_cells = tuple(
        _projected_multi_cell(
            raw_record=seed_member_input.raw_record,
            trim_result=seed_member_input.trim_result,
            strand="forward",
            base=base,
            query_index=query_index,
            oriented_qualities=oriented_qualities,
        )
        for query_index, base in enumerate(seed_member_input.display_sequence)
    )
    return _PlacedMultiMember(
        member=MultiAssemblyMember(
            member_index=seed_member_input.member_index,
            source_filename=seed_member_input.source_filename,
            display_name=seed_member_input.display_name,
            chosen_orientation="forward",
            included=True,
            inclusion_reason=None,
            trimmed_length=seed_member_input.trimmed_length,
            has_trace_data=seed_member_input.has_trace_data,
            has_qualities=seed_member_input.has_qualities,
            is_seed=True,
            alignment_score=None,
            overlap_length=len(seed_member_input.display_sequence),
            percent_identity=100.0,
        ),
        seed_column_cells=seed_column_cells,
        insertion_cells_by_bucket={},
    )


def _excluded_multi_member_for_empty_sequence(
    member_input: _ResolvedMultiMemberInput,
) -> _PlacedMultiMember:
    return _PlacedMultiMember(
        member=MultiAssemblyMember(
            member_index=member_input.member_index,
            source_filename=member_input.source_filename,
            display_name=member_input.display_name,
            chosen_orientation="forward",
            included=False,
            inclusion_reason="trimmed read is empty",
            trimmed_length=member_input.trimmed_length,
            has_trace_data=member_input.has_trace_data,
            has_qualities=member_input.has_qualities,
        ),
        seed_column_cells=(),
        insertion_cells_by_bucket={},
    )


def _place_member_against_multi_seed(
    *,
    seed_member_input: _ResolvedMultiMemberInput,
    member_input: _ResolvedMultiMemberInput,
    config: AssemblyConfig,
    aligner: Align.PairwiseAligner,
) -> _PlacedMultiMember:
    if not member_input.display_sequence:
        return _excluded_multi_member_for_empty_sequence(member_input)

    candidates: list[tuple[AssemblyStrand, str, list[int] | None]] = [
        (
            "forward",
            member_input.display_sequence,
            _candidate_member_qualities(
                member_input.display_qualities,
                strand="forward",
            ),
        ),
        (
            "reverse-complement",
            reverse_complement_sequence(member_input.display_sequence),
            _candidate_member_qualities(
                member_input.display_qualities,
                strand="reverse-complement",
            ),
        ),
    ]

    best_alignment = None
    best_orientation: AssemblyStrand = "forward"
    best_sequence = ""
    best_qualities: list[int] | None = None
    best_score = float("-inf")

    for orientation, oriented_sequence, oriented_qualities in candidates:
        alignment = aligner.align(
            seed_member_input.display_sequence, oriented_sequence
        )[0]
        if alignment.score > best_score:
            best_alignment = alignment
            best_orientation = orientation
            best_sequence = oriented_sequence
            best_qualities = oriented_qualities
            best_score = float(alignment.score)

    if best_alignment is None:
        return _PlacedMultiMember(
            member=MultiAssemblyMember(
                member_index=member_input.member_index,
                source_filename=member_input.source_filename,
                display_name=member_input.display_name,
                chosen_orientation="forward",
                included=False,
                inclusion_reason="no alignment could be generated",
                trimmed_length=member_input.trimmed_length,
                has_trace_data=member_input.has_trace_data,
                has_qualities=member_input.has_qualities,
            ),
            seed_column_cells=(),
            insertion_cells_by_bucket={},
        )

    overlap_length, percent_identity = _alignment_overlap_metrics(
        best_alignment,
        target_sequence=seed_member_input.display_sequence,
        query_sequence=best_sequence,
    )
    rejection_reason = _multi_member_rejection_reason(
        overlap_length=overlap_length,
        percent_identity=percent_identity,
        config=config,
    )
    if rejection_reason is not None:
        return _PlacedMultiMember(
            member=MultiAssemblyMember(
                member_index=member_input.member_index,
                source_filename=member_input.source_filename,
                display_name=member_input.display_name,
                chosen_orientation=best_orientation,
                included=False,
                inclusion_reason=rejection_reason,
                trimmed_length=member_input.trimmed_length,
                has_trace_data=member_input.has_trace_data,
                has_qualities=member_input.has_qualities,
                alignment_score=best_score,
                overlap_length=overlap_length,
                percent_identity=percent_identity,
            ),
            seed_column_cells=(),
            insertion_cells_by_bucket={},
        )

    seed_column_cells, insertion_cells_by_bucket = _project_member_against_seed(
        alignment=best_alignment,
        member_input=member_input,
        oriented_sequence=best_sequence,
        oriented_qualities=best_qualities,
        strand=best_orientation,
        seed_length=len(seed_member_input.display_sequence),
    )
    return _PlacedMultiMember(
        member=MultiAssemblyMember(
            member_index=member_input.member_index,
            source_filename=member_input.source_filename,
            display_name=member_input.display_name,
            chosen_orientation=best_orientation,
            included=True,
            inclusion_reason=None,
            trimmed_length=member_input.trimmed_length,
            has_trace_data=member_input.has_trace_data,
            has_qualities=member_input.has_qualities,
            alignment_score=best_score,
            overlap_length=overlap_length,
            percent_identity=percent_identity,
        ),
        seed_column_cells=seed_column_cells,
        insertion_cells_by_bucket=insertion_cells_by_bucket,
    )


def _multi_member_rejection_reason(
    *,
    overlap_length: int,
    percent_identity: float,
    config: AssemblyConfig,
) -> str | None:
    if overlap_length <= 0:
        return "no overlapping aligned columns were found"
    if overlap_length < config.min_overlap_length:
        return (
            "overlap length below threshold "
            f"({overlap_length} < {config.min_overlap_length})"
        )
    if percent_identity < config.min_percent_identity:
        return (
            "percent identity below threshold "
            f"({percent_identity:.1f}% < {config.min_percent_identity:.1f}%)"
        )
    return None


def _alignment_overlap_metrics(
    alignment,
    *,
    target_sequence: str,
    query_sequence: str,
) -> tuple[int, float]:
    overlap_length = 0
    match_count = 0
    for target_index, query_index in zip(
        alignment.indices[0],
        alignment.indices[1],
        strict=True,
    ):
        resolved_target_index = int(target_index)
        resolved_query_index = int(query_index)
        if resolved_target_index < 0 or resolved_query_index < 0:
            continue
        overlap_length += 1
        if (
            target_sequence[resolved_target_index]
            == query_sequence[resolved_query_index]
        ):
            match_count += 1
    percent_identity = (match_count / overlap_length) * 100.0 if overlap_length else 0.0
    return overlap_length, percent_identity


def _project_member_against_seed(
    *,
    alignment,
    member_input: _ResolvedMultiMemberInput,
    oriented_sequence: str,
    oriented_qualities: list[int] | None,
    strand: AssemblyStrand,
    seed_length: int,
) -> tuple[tuple[_ProjectedMultiCell, ...], dict[int, tuple[_ProjectedMultiCell, ...]]]:
    seed_column_cells: list[_ProjectedMultiCell] = [
        _gap_multi_projected_cell() for _ in range(seed_length)
    ]
    insertion_cells_by_bucket: dict[int, list[_ProjectedMultiCell]] = {}
    last_seed_index = -1

    for target_index, query_index in zip(
        alignment.indices[0],
        alignment.indices[1],
        strict=True,
    ):
        resolved_target_index = int(target_index)
        resolved_query_index = int(query_index)

        if resolved_target_index >= 0:
            last_seed_index = resolved_target_index
            if resolved_query_index < 0:
                continue
            seed_column_cells[resolved_target_index] = _projected_multi_cell(
                raw_record=member_input.raw_record,
                trim_result=member_input.trim_result,
                strand=strand,
                base=oriented_sequence[resolved_query_index],
                query_index=resolved_query_index,
                oriented_qualities=oriented_qualities,
            )
            continue

        if resolved_query_index < 0:
            continue
        bucket_index = last_seed_index + 1
        insertion_cells_by_bucket.setdefault(bucket_index, []).append(
            _projected_multi_cell(
                raw_record=member_input.raw_record,
                trim_result=member_input.trim_result,
                strand=strand,
                base=oriented_sequence[resolved_query_index],
                query_index=resolved_query_index,
                oriented_qualities=oriented_qualities,
            )
        )

    return (
        tuple(seed_column_cells),
        {
            bucket_index: tuple(cells)
            for bucket_index, cells in insertion_cells_by_bucket.items()
        },
    )


def _projected_multi_cell(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    strand: AssemblyStrand,
    base: str,
    query_index: int,
    oriented_qualities: list[int] | None,
) -> _ProjectedMultiCell:
    return _ProjectedMultiCell(
        base=base,
        query_index=query_index,
        query_pos=query_index + 1,
        quality=(
            None if oriented_qualities is None else int(oriented_qualities[query_index])
        ),
        trace_x=_trace_position_for_member_query_index(
            raw_record=raw_record,
            trim_result=trim_result,
            oriented_query_index=query_index,
            strand=strand,
        ),
        is_gap=False,
    )


def _gap_multi_projected_cell() -> _ProjectedMultiCell:
    return _ProjectedMultiCell(
        base="-",
        query_index=None,
        query_pos=None,
        quality=None,
        trace_x=None,
        is_gap=True,
    )


def _resolve_multi_insertion_bucket(
    *,
    bucket_cells_by_member: dict[int, tuple[_ProjectedMultiCell, ...]],
    included_member_indices: tuple[int, ...],
    aligner: Align.PairwiseAligner,
) -> tuple[dict[int, _ProjectedMultiCell], ...]:
    non_empty_bucket_cells = [
        (member_index, cells)
        for member_index, cells in bucket_cells_by_member.items()
        if cells
    ]
    if not non_empty_bucket_cells:
        return ()

    non_empty_bucket_cells.sort(key=lambda item: (-len(item[1]), item[0]))
    anchor_member_index, anchor_cells = non_empty_bucket_cells[0]
    resolved_columns = [{anchor_member_index: cell} for cell in anchor_cells]
    for member_index, member_cells in non_empty_bucket_cells[1:]:
        resolved_columns = _merge_multi_bucket_member_cells(
            resolved_columns,
            member_index=member_index,
            member_cells=member_cells,
            aligner=aligner,
        )

    for resolved_column in resolved_columns:
        for member_index in included_member_indices:
            resolved_column.setdefault(member_index, _gap_multi_projected_cell())

    return tuple(resolved_columns)


def _merge_multi_bucket_member_cells(
    resolved_columns: list[dict[int, _ProjectedMultiCell]],
    *,
    member_index: int,
    member_cells: tuple[_ProjectedMultiCell, ...],
    aligner: Align.PairwiseAligner,
) -> list[dict[int, _ProjectedMultiCell]]:
    representative_sequence = "".join(
        _multi_bucket_column_representative(resolved_column)
        for resolved_column in resolved_columns
    )
    member_sequence = "".join(cell.base for cell in member_cells)
    if not representative_sequence:
        return [{member_index: cell} for cell in member_cells]

    alignment = aligner.align(representative_sequence, member_sequence)[0]
    merged_columns: list[dict[int, _ProjectedMultiCell]] = []
    for target_index, query_index in zip(
        alignment.indices[0],
        alignment.indices[1],
        strict=True,
    ):
        resolved_target_index = int(target_index)
        resolved_query_index = int(query_index)
        if resolved_target_index >= 0:
            merged_column = dict(resolved_columns[resolved_target_index])
            merged_column[member_index] = (
                _gap_multi_projected_cell()
                if resolved_query_index < 0
                else member_cells[resolved_query_index]
            )
            merged_columns.append(merged_column)
            continue
        if resolved_query_index < 0:
            continue
        merged_columns.append({member_index: member_cells[resolved_query_index]})
    return merged_columns


def _multi_bucket_column_representative(
    resolved_column: dict[int, _ProjectedMultiCell],
) -> str:
    bases = [
        cell.base
        for cell in resolved_column.values()
        if not cell.is_gap and cell.base != "-"
    ]
    if not bases:
        return "N"
    base_counts = Counter(bases)
    return sorted(
        base_counts.items(),
        key=lambda item: (-item[1], item[0]),
    )[
        0
    ][0]


def _build_multi_assembly_column(
    *,
    column_index: int,
    column_cells_by_member: dict[int, _ProjectedMultiCell],
    included_member_indices: tuple[int, ...],
    quality_margin: int,
) -> MultiAssemblyColumn | None:
    ordered_projected_cells = tuple(
        column_cells_by_member.get(member_index, _gap_multi_projected_cell())
        for member_index in included_member_indices
    )
    non_gap_projected_cells = tuple(
        cell for cell in ordered_projected_cells if not cell.is_gap and cell.base != "-"
    )
    if not non_gap_projected_cells:
        return None

    support_counts: Counter[str] = Counter(
        cell.base for cell in non_gap_projected_cells
    )
    quality_sums: dict[str, int] = {}
    for cell in non_gap_projected_cells:
        quality_sums[cell.base] = quality_sums.get(cell.base, 0) + int(
            0 if cell.quality is None else cell.quality
        )

    consensus_base, resolution = _resolve_multi_consensus_column(
        support_counts=support_counts,
        quality_sums=quality_sums,
        quality_margin=quality_margin,
    )
    return MultiAssemblyColumn(
        column_index=column_index,
        consensus_base=consensus_base,
        resolution=resolution,
        member_cells=tuple(
            MultiAssemblyMemberCell(
                member_index=member_index,
                base=projected_cell.base,
                query_index=projected_cell.query_index,
                query_pos=projected_cell.query_pos,
                quality=projected_cell.quality,
                trace_x=projected_cell.trace_x,
                is_gap=projected_cell.is_gap,
            )
            for member_index, projected_cell in zip(
                included_member_indices,
                ordered_projected_cells,
                strict=True,
            )
        ),
        support_counts=tuple(
            sorted(
                support_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
        ),
        quality_sums=tuple(sorted(quality_sums.items())),
        non_gap_member_count=len(non_gap_projected_cells),
        gap_member_count=len(ordered_projected_cells) - len(non_gap_projected_cells),
        ambiguous=consensus_base == "N",
    )


def _resolve_multi_consensus_column(
    *,
    support_counts: Counter[str],
    quality_sums: dict[str, int],
    quality_margin: int,
) -> tuple[str, ConflictResolution]:
    ranked_bases = sorted(
        support_counts,
        key=lambda base: (-support_counts[base], -quality_sums.get(base, 0), base),
    )
    if len(ranked_bases) == 1:
        base = ranked_bases[0]
        return (
            (base, "single_read") if support_counts[base] == 1 else (base, "concordant")
        )

    top_base = ranked_bases[0]
    next_base = ranked_bases[1]
    top_support = support_counts[top_base]
    next_support = support_counts[next_base]
    top_quality = quality_sums.get(top_base, 0)
    next_quality = quality_sums.get(next_base, 0)

    if top_support > next_support:
        return top_base, "majority_resolved"
    if (top_quality - next_quality) >= quality_margin:
        return top_base, "quality_resolved"
    return "N", "ambiguous"


def _mean_quality(qualities: list[int] | None) -> float:
    if not qualities:
        return 0.0
    return float(sum(qualities)) / float(len(qualities))
