from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from Bio import Align

from abi_sauce.alignment_policy import (
    alignment_overlap_metrics,
    alignment_strands_for_policy,
    select_best_oriented_alignment,
)
from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.models import SequenceRecord
from abi_sauce.oriented_reads import PreparedTrimmedRead, prepare_trimmed_reads
from abi_sauce.reference_alignment import build_aligner, normalize_reference
from abi_sauce.reference_alignment_types import (
    ChosenStrand,
    ReferenceConsensusResolution,
    ReferenceMultiAnchorKind,
    ReferenceMultiAlignmentColumn,
    ReferenceMultiAlignmentMember,
    ReferenceMultiAlignmentMemberCell,
    ReferenceMultiAlignmentResult,
    StrandPolicy,
)
from abi_sauce.trace_coordinates import trace_position_for_oriented_query_index
from abi_sauce.trimming import TrimResult


@dataclass(frozen=True, slots=True)
class _ResolvedReferenceMultiMemberInput:
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
class _ProjectedReferenceMultiCell:
    base: str
    query_index: int | None
    query_pos: int | None
    qscore: int | None
    trace_x: int | None
    is_gap: bool


@dataclass(frozen=True, slots=True)
class _PlacedReferenceMultiMember:
    member: ReferenceMultiAlignmentMember
    reference_cells: tuple[_ProjectedReferenceMultiCell, ...]
    insertion_cells_by_bucket: dict[int, tuple[_ProjectedReferenceMultiCell, ...]]


def align_trimmed_reads_to_reference(
    *,
    source_filenames: tuple[str, ...] | list[str],
    raw_records_by_source_filename: dict[str, SequenceRecord],
    trim_results_by_source_filename: dict[str, TrimResult],
    reference_text: str,
    reference_name: str | None = None,
    strand_policy: StrandPolicy = "auto",
    config: AssemblyConfig | None = None,
    aligner: Align.PairwiseAligner | None = None,
) -> ReferenceMultiAlignmentResult:
    """Align multiple trimmed reads independently to one shared reference grid."""
    normalized_reference_name, reference_sequence = normalize_reference(reference_text)
    resolved_reference_name = (
        normalized_reference_name if reference_name is None else reference_name
    )
    resolved_config = AssemblyConfig() if config is None else config
    resolved_aligner = (
        build_aligner(
            match_score=resolved_config.match_score,
            mismatch_score=resolved_config.mismatch_score,
            open_internal_gap_score=resolved_config.open_internal_gap_score,
            extend_internal_gap_score=resolved_config.extend_internal_gap_score,
        )
        if aligner is None
        else aligner
    )

    prepared_reads = prepare_trimmed_reads(
        source_filenames=source_filenames,
        raw_records_by_source_filename=raw_records_by_source_filename,
        trim_results_by_source_filename=trim_results_by_source_filename,
    )
    member_inputs = tuple(
        _resolved_reference_multi_member_input(
            member_index=member_index,
            prepared_read=prepared_read,
        )
        for member_index, prepared_read in enumerate(prepared_reads)
    )

    allowed_strands = alignment_strands_for_policy(strand_policy)
    placements = tuple(
        _place_member_against_reference(
            reference_sequence=reference_sequence,
            member_input=member_input,
            config=resolved_config,
            aligner=resolved_aligner,
            allowed_strands=allowed_strands,
        )
        for member_input in member_inputs
    )
    included_placements = tuple(
        placement for placement in placements if placement.member.included
    )
    included_member_indices = tuple(
        placement.member.member_index for placement in included_placements
    )

    bucket_columns_by_index = {
        bucket_index: _resolve_reference_multi_insertion_bucket(
            bucket_cells_by_member={
                placement.member.member_index: placement.insertion_cells_by_bucket.get(
                    bucket_index,
                    (),
                )
                for placement in included_placements
            },
            included_member_indices=included_member_indices,
            aligner=resolved_aligner,
        )
        for bucket_index in range(len(reference_sequence) + 1)
    }

    aligned_sequences_by_member_index: dict[int, list[str]] = {
        member_index: [] for member_index in included_member_indices
    }
    columns: list[ReferenceMultiAlignmentColumn] = []
    column_index = 1

    for bucket_index in range(len(reference_sequence) + 1):
        for bucket_column in bucket_columns_by_index[bucket_index]:
            column = _build_reference_multi_column(
                column_index=column_index,
                anchor_kind="insertion",
                anchor_index=bucket_index,
                ref_index=None,
                ref_base="-",
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

        if bucket_index >= len(reference_sequence):
            continue

        column = _build_reference_multi_column(
            column_index=column_index,
            anchor_kind="reference",
            anchor_index=bucket_index,
            ref_index=bucket_index,
            ref_base=reference_sequence[bucket_index],
            column_cells_by_member={
                placement.member.member_index: placement.reference_cells[bucket_index]
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

    gapped_reference = "".join(column.ref_base for column in columns)
    gapped_consensus = "".join(column.consensus_base for column in columns)
    consensus_sequence = "".join(base for base in gapped_consensus if base != "-")
    ambiguous_column_count = sum(column.ambiguous for column in columns)
    included_member_count = len(included_member_indices)

    accepted = True
    rejection_reason = None
    if included_member_count < 2:
        accepted = False
        rejection_reason = (
            "fewer than 2 reads passed shared-reference placement thresholds"
        )
    elif not columns:
        accepted = False
        rejection_reason = "no shared-reference columns could be derived"

    return ReferenceMultiAlignmentResult(
        reference_name=resolved_reference_name,
        reference_sequence=reference_sequence,
        accepted=accepted,
        rejection_reason=rejection_reason,
        members=tuple(placement.member for placement in placements),
        included_member_indices=included_member_indices,
        aligned_member_sequences=tuple(
            "".join(aligned_sequences_by_member_index[member_index])
            for member_index in included_member_indices
        ),
        gapped_reference=gapped_reference,
        gapped_consensus=gapped_consensus,
        consensus_sequence=consensus_sequence,
        columns=tuple(columns),
        included_member_count=included_member_count,
        excluded_member_count=len(member_inputs) - included_member_count,
        ambiguous_column_count=ambiguous_column_count,
    )


def _resolved_reference_multi_member_input(
    *,
    member_index: int,
    prepared_read: PreparedTrimmedRead,
) -> _ResolvedReferenceMultiMemberInput:
    return _ResolvedReferenceMultiMemberInput(
        member_index=member_index,
        source_filename=prepared_read.source_filename,
        raw_record=prepared_read.raw_record,
        trim_result=prepared_read.trim_result,
        display_name=prepared_read.display_name,
        display_sequence=prepared_read.display_sequence,
        display_qualities=prepared_read.display_qualities,
        trimmed_length=prepared_read.trimmed_length,
        has_trace_data=prepared_read.has_trace_data,
        has_qualities=prepared_read.has_qualities,
    )


def _place_member_against_reference(
    *,
    reference_sequence: str,
    member_input: _ResolvedReferenceMultiMemberInput,
    config: AssemblyConfig,
    aligner: Align.PairwiseAligner,
    allowed_strands: tuple[ChosenStrand, ...],
) -> _PlacedReferenceMultiMember:
    if not member_input.display_sequence:
        return _excluded_reference_multi_member(
            member_input,
            inclusion_reason="trimmed read is empty",
        )

    best_alignment = select_best_oriented_alignment(
        target_sequence=reference_sequence,
        display_query_sequence=member_input.display_sequence,
        display_query_qualities=member_input.display_qualities,
        aligner=aligner,
        allowed_strands=allowed_strands,
    )
    if best_alignment is None:
        return _excluded_reference_multi_member(
            member_input,
            inclusion_reason="no alignment could be generated",
        )

    overlap_length, percent_identity = alignment_overlap_metrics(
        best_alignment.alignment,
        target_sequence=reference_sequence,
        query_sequence=best_alignment.sequence,
    )
    rejection_reason = _reference_multi_rejection_reason(
        overlap_length=overlap_length,
        percent_identity=percent_identity,
        config=config,
    )
    mismatch_count, insertion_count, deletion_count = _alignment_change_counts(
        best_alignment.alignment,
        reference_sequence=reference_sequence,
        oriented_sequence=best_alignment.sequence,
    )
    if rejection_reason is not None:
        return _excluded_reference_multi_member(
            member_input,
            inclusion_reason=rejection_reason,
            chosen_strand=best_alignment.strand,
            alignment_score=best_alignment.score,
            overlap_length=overlap_length,
            percent_identity=percent_identity,
            mismatch_count=mismatch_count,
            insertion_count=insertion_count,
            deletion_count=deletion_count,
        )

    reference_cells, insertion_cells_by_bucket = _project_member_against_reference(
        alignment=best_alignment.alignment,
        member_input=member_input,
        oriented_sequence=best_alignment.sequence,
        oriented_qualities=best_alignment.qualities,
        strand=best_alignment.strand,
        reference_length=len(reference_sequence),
    )
    return _PlacedReferenceMultiMember(
        member=ReferenceMultiAlignmentMember(
            member_index=member_input.member_index,
            source_filename=member_input.source_filename,
            display_name=member_input.display_name,
            chosen_strand=best_alignment.strand,
            included=True,
            inclusion_reason=None,
            trimmed_length=member_input.trimmed_length,
            has_trace_data=member_input.has_trace_data,
            has_qualities=member_input.has_qualities,
            alignment_score=best_alignment.score,
            overlap_length=overlap_length,
            percent_identity=percent_identity,
            mismatch_count=mismatch_count,
            insertion_count=insertion_count,
            deletion_count=deletion_count,
        ),
        reference_cells=reference_cells,
        insertion_cells_by_bucket=insertion_cells_by_bucket,
    )


def _excluded_reference_multi_member(
    member_input: _ResolvedReferenceMultiMemberInput,
    *,
    inclusion_reason: str,
    chosen_strand: ChosenStrand = "forward",
    alignment_score: float | None = None,
    overlap_length: int | None = None,
    percent_identity: float | None = None,
    mismatch_count: int | None = None,
    insertion_count: int | None = None,
    deletion_count: int | None = None,
) -> _PlacedReferenceMultiMember:
    return _PlacedReferenceMultiMember(
        member=ReferenceMultiAlignmentMember(
            member_index=member_input.member_index,
            source_filename=member_input.source_filename,
            display_name=member_input.display_name,
            chosen_strand=chosen_strand,
            included=False,
            inclusion_reason=inclusion_reason,
            trimmed_length=member_input.trimmed_length,
            has_trace_data=member_input.has_trace_data,
            has_qualities=member_input.has_qualities,
            alignment_score=alignment_score,
            overlap_length=overlap_length,
            percent_identity=percent_identity,
            mismatch_count=mismatch_count,
            insertion_count=insertion_count,
            deletion_count=deletion_count,
        ),
        reference_cells=(),
        insertion_cells_by_bucket={},
    )


def _reference_multi_rejection_reason(
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


def _alignment_change_counts(
    alignment,
    *,
    reference_sequence: str,
    oriented_sequence: str,
) -> tuple[int, int, int]:
    mismatch_count = 0
    insertion_count = 0
    deletion_count = 0
    for target_index, query_index in zip(
        alignment.indices[0],
        alignment.indices[1],
        strict=True,
    ):
        resolved_target_index = int(target_index)
        resolved_query_index = int(query_index)
        if resolved_target_index >= 0 and resolved_query_index >= 0:
            if (
                reference_sequence[resolved_target_index]
                != oriented_sequence[resolved_query_index]
            ):
                mismatch_count += 1
        elif resolved_target_index < 0 and resolved_query_index >= 0:
            insertion_count += 1
        elif resolved_target_index >= 0 and resolved_query_index < 0:
            deletion_count += 1
    return mismatch_count, insertion_count, deletion_count


def _project_member_against_reference(
    *,
    alignment,
    member_input: _ResolvedReferenceMultiMemberInput,
    oriented_sequence: str,
    oriented_qualities: list[int] | None,
    strand: ChosenStrand,
    reference_length: int,
) -> tuple[
    tuple[_ProjectedReferenceMultiCell, ...],
    dict[int, tuple[_ProjectedReferenceMultiCell, ...]],
]:
    reference_cells: list[_ProjectedReferenceMultiCell] = [
        _gap_reference_multi_projected_cell() for _ in range(reference_length)
    ]
    insertion_cells_by_bucket: dict[int, list[_ProjectedReferenceMultiCell]] = {}
    last_reference_index = -1

    for target_index, query_index in zip(
        alignment.indices[0],
        alignment.indices[1],
        strict=True,
    ):
        resolved_target_index = int(target_index)
        resolved_query_index = int(query_index)

        if resolved_target_index >= 0:
            last_reference_index = resolved_target_index
            if resolved_query_index < 0:
                continue
            reference_cells[resolved_target_index] = _projected_reference_multi_cell(
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
        bucket_index = last_reference_index + 1
        insertion_cells_by_bucket.setdefault(bucket_index, []).append(
            _projected_reference_multi_cell(
                raw_record=member_input.raw_record,
                trim_result=member_input.trim_result,
                strand=strand,
                base=oriented_sequence[resolved_query_index],
                query_index=resolved_query_index,
                oriented_qualities=oriented_qualities,
            )
        )

    return (
        tuple(reference_cells),
        {
            bucket_index: tuple(cells)
            for bucket_index, cells in insertion_cells_by_bucket.items()
        },
    )


def _projected_reference_multi_cell(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    strand: ChosenStrand,
    base: str,
    query_index: int,
    oriented_qualities: list[int] | None,
) -> _ProjectedReferenceMultiCell:
    return _ProjectedReferenceMultiCell(
        base=base,
        query_index=query_index,
        query_pos=query_index + 1,
        qscore=(
            None if oriented_qualities is None else int(oriented_qualities[query_index])
        ),
        trace_x=trace_position_for_oriented_query_index(
            raw_record=raw_record,
            trim_result=trim_result,
            oriented_query_index=query_index,
            strand=strand,
        ),
        is_gap=False,
    )


def _gap_reference_multi_projected_cell() -> _ProjectedReferenceMultiCell:
    return _ProjectedReferenceMultiCell(
        base="-",
        query_index=None,
        query_pos=None,
        qscore=None,
        trace_x=None,
        is_gap=True,
    )


def _resolve_reference_multi_insertion_bucket(
    *,
    bucket_cells_by_member: dict[int, tuple[_ProjectedReferenceMultiCell, ...]],
    included_member_indices: tuple[int, ...],
    aligner: Align.PairwiseAligner,
) -> tuple[dict[int, _ProjectedReferenceMultiCell], ...]:
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
        resolved_columns = _merge_reference_multi_bucket_member_cells(
            resolved_columns,
            member_index=member_index,
            member_cells=member_cells,
            aligner=aligner,
        )

    for resolved_column in resolved_columns:
        for member_index in included_member_indices:
            resolved_column.setdefault(
                member_index, _gap_reference_multi_projected_cell()
            )

    return tuple(resolved_columns)


def _merge_reference_multi_bucket_member_cells(
    resolved_columns: list[dict[int, _ProjectedReferenceMultiCell]],
    *,
    member_index: int,
    member_cells: tuple[_ProjectedReferenceMultiCell, ...],
    aligner: Align.PairwiseAligner,
) -> list[dict[int, _ProjectedReferenceMultiCell]]:
    representative_sequence = "".join(
        _reference_multi_bucket_column_representative(resolved_column)
        for resolved_column in resolved_columns
    )
    member_sequence = "".join(cell.base for cell in member_cells)
    if not representative_sequence:
        return [{member_index: cell} for cell in member_cells]

    alignment = aligner.align(representative_sequence, member_sequence)[0]
    merged_columns: list[dict[int, _ProjectedReferenceMultiCell]] = []
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
                _gap_reference_multi_projected_cell()
                if resolved_query_index < 0
                else member_cells[resolved_query_index]
            )
            merged_columns.append(merged_column)
            continue
        if resolved_query_index < 0:
            continue
        merged_columns.append({member_index: member_cells[resolved_query_index]})
    return merged_columns


def _reference_multi_bucket_column_representative(
    resolved_column: dict[int, _ProjectedReferenceMultiCell],
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


def _build_reference_multi_column(
    *,
    column_index: int,
    anchor_kind: ReferenceMultiAnchorKind,
    anchor_index: int,
    ref_index: int | None,
    ref_base: str,
    column_cells_by_member: dict[int, _ProjectedReferenceMultiCell],
    included_member_indices: tuple[int, ...],
    quality_margin: int,
) -> ReferenceMultiAlignmentColumn | None:
    ordered_projected_cells = tuple(
        column_cells_by_member.get(member_index, _gap_reference_multi_projected_cell())
        for member_index in included_member_indices
    )
    non_gap_projected_cells = tuple(
        cell for cell in ordered_projected_cells if not cell.is_gap and cell.base != "-"
    )
    if anchor_kind == "insertion" and not non_gap_projected_cells:
        return None

    support_counts: Counter[str] = Counter(
        cell.base for cell in non_gap_projected_cells
    )
    quality_sums: dict[str, int] = {}
    for cell in non_gap_projected_cells:
        quality_sums[cell.base] = quality_sums.get(cell.base, 0) + int(
            0 if cell.qscore is None else cell.qscore
        )

    consensus_base, resolution = _resolve_reference_multi_consensus_column(
        support_counts=support_counts,
        quality_sums=quality_sums,
        quality_margin=quality_margin,
        ref_base=ref_base,
        allow_deleted=anchor_kind == "reference",
    )
    matches_reference = (
        None if anchor_kind == "insertion" else consensus_base == ref_base
    )
    return ReferenceMultiAlignmentColumn(
        column_index=column_index,
        anchor_kind=anchor_kind,
        anchor_index=anchor_index,
        ref_index=ref_index,
        ref_pos=None if ref_index is None else ref_index + 1,
        ref_base=ref_base,
        consensus_base=consensus_base,
        resolution=resolution,
        member_cells=tuple(
            ReferenceMultiAlignmentMemberCell(
                member_index=member_index,
                base=projected_cell.base,
                query_index=projected_cell.query_index,
                query_pos=projected_cell.query_pos,
                qscore=projected_cell.qscore,
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
        matches_reference=matches_reference,
    )


def _resolve_reference_multi_consensus_column(
    *,
    support_counts: Counter[str],
    quality_sums: dict[str, int],
    quality_margin: int,
    ref_base: str,
    allow_deleted: bool,
) -> tuple[str, ReferenceConsensusResolution]:
    if not support_counts:
        if allow_deleted:
            return "-", "deleted"
        return "N", "ambiguous"

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
    if ref_base != "-" and ref_base in support_counts:
        return "N", "ambiguous"
    return "N", "ambiguous"


__all__ = ["align_trimmed_reads_to_reference"]
