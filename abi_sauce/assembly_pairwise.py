from __future__ import annotations

from Bio import Align

from abi_sauce.alignment_policy import (
    build_semiglobal_aligner,
    oriented_qualities,
    select_best_oriented_alignment,
)
from abi_sauce.assembly_types import (
    AssemblyColumn,
    AssemblyConfig,
    AssemblyConflict,
    AssemblyResult,
    AssemblyStrand,
    ConflictResolution,
)
from abi_sauce.models import SequenceRecord
from abi_sauce.oriented_reads import prepare_trimmed_read
from abi_sauce.trace_coordinates import trace_position_for_oriented_query_index
from abi_sauce.trimming import TrimResult


def build_assembly_aligner(
    *,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    open_internal_gap_score: float = -3.0,
    extend_internal_gap_score: float = -1.0,
) -> Align.PairwiseAligner:
    """Build a semi-global style aligner for pairwise assembly."""
    return build_semiglobal_aligner(
        match_score=match_score,
        mismatch_score=mismatch_score,
        open_internal_gap_score=open_internal_gap_score,
        extend_internal_gap_score=extend_internal_gap_score,
    )


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

    left_prepared_read = prepare_trimmed_read(
        source_filename=left_source_filename,
        raw_record=left_raw_record,
        trim_result=left_trim_result,
    )
    right_prepared_read = prepare_trimmed_read(
        source_filename=right_source_filename,
        raw_record=right_raw_record,
        trim_result=right_trim_result,
    )
    left_sequence = left_prepared_read.display_sequence
    right_display_sequence = right_prepared_read.display_sequence

    if not left_sequence:
        return _empty_assembly_result(
            left_source_filename=left_source_filename,
            right_source_filename=right_source_filename,
            left_display_name=left_prepared_read.display_name,
            right_display_name=right_prepared_read.display_name,
            chosen_right_orientation="forward",
            rejection_reason="left read trims down to an empty sequence",
        )
    if not right_display_sequence:
        return _empty_assembly_result(
            left_source_filename=left_source_filename,
            right_source_filename=right_source_filename,
            left_display_name=left_prepared_read.display_name,
            right_display_name=right_prepared_read.display_name,
            chosen_right_orientation="forward",
            rejection_reason="right read trims down to an empty sequence",
        )

    best_right_alignment = select_best_oriented_alignment(
        target_sequence=left_sequence,
        display_query_sequence=right_display_sequence,
        display_query_qualities=right_prepared_read.display_qualities,
        aligner=resolved_aligner,
    )

    if best_right_alignment is None:
        return _empty_assembly_result(
            left_source_filename=left_source_filename,
            right_source_filename=right_source_filename,
            left_display_name=left_prepared_read.display_name,
            right_display_name=right_prepared_read.display_name,
            chosen_right_orientation="forward",
            rejection_reason="no alignment could be generated",
        )

    result = _extract_assembly_result(
        alignment=best_right_alignment.alignment,
        left_source_filename=left_source_filename,
        left_display_name=left_prepared_read.display_name,
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        left_oriented_sequence=left_sequence,
        left_oriented_qualities=oriented_qualities(
            left_prepared_read.display_qualities,
            strand="forward",
        ),
        right_source_filename=right_source_filename,
        right_display_name=right_prepared_read.display_name,
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        right_oriented_sequence=best_right_alignment.sequence,
        right_oriented_qualities=best_right_alignment.qualities,
        chosen_right_orientation=best_right_alignment.strand,
        quality_margin=resolved_config.quality_margin,
        score=best_right_alignment.score,
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


def _extract_assembly_result(
    *,
    alignment,
    left_source_filename: str,
    left_display_name: str,
    left_raw_record: SequenceRecord,
    left_trim_result: TrimResult,
    left_oriented_sequence: str,
    left_oriented_qualities: list[int] | None,
    right_source_filename: str,
    right_display_name: str,
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

        left_trace_x = trace_position_for_oriented_query_index(
            raw_record=left_raw_record,
            trim_result=left_trim_result,
            oriented_query_index=(
                None if resolved_target_index < 0 else resolved_target_index
            ),
            strand="forward",
        )
        right_trace_x = trace_position_for_oriented_query_index(
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
        left_display_name=left_display_name,
        right_display_name=right_display_name,
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
