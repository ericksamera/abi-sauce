from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from Bio import Align

from abi_sauce.alignment_policy import (
    alignment_overlap_metrics,
    oriented_qualities,
    select_best_oriented_alignment,
)
from abi_sauce.assembly_pairwise import build_assembly_aligner
from abi_sauce.assembly_types import (
    AssemblyConfig,
    AssemblyStrand,
    ConflictResolution,
    MultiAssemblyColumn,
    MultiAssemblyMember,
    MultiAssemblyMemberCell,
    MultiAssemblyResult,
)
from abi_sauce.models import SequenceRecord
from abi_sauce.oriented_reads import PreparedTrimmedRead, prepare_trimmed_reads
from abi_sauce.trace_coordinates import trace_position_for_oriented_query_index
from abi_sauce.trimming import TrimResult


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

    prepared_reads = prepare_trimmed_reads(
        source_filenames=source_filenames,
        raw_records_by_source_filename=raw_records_by_source_filename,
        trim_results_by_source_filename=trim_results_by_source_filename,
    )
    member_inputs = tuple(
        _resolved_multi_member_input(
            member_index=member_index,
            prepared_read=prepared_read,
        )
        for member_index, prepared_read in enumerate(prepared_reads)
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


def _resolved_multi_member_input(
    *,
    member_index: int,
    prepared_read: PreparedTrimmedRead,
) -> _ResolvedMultiMemberInput:
    return _ResolvedMultiMemberInput(
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
    seed_oriented_qualities = oriented_qualities(
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
            oriented_qualities=seed_oriented_qualities,
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

    best_member_alignment = select_best_oriented_alignment(
        target_sequence=seed_member_input.display_sequence,
        display_query_sequence=member_input.display_sequence,
        display_query_qualities=member_input.display_qualities,
        aligner=aligner,
    )

    if best_member_alignment is None:
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

    overlap_length, percent_identity = alignment_overlap_metrics(
        best_member_alignment.alignment,
        target_sequence=seed_member_input.display_sequence,
        query_sequence=best_member_alignment.sequence,
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
                chosen_orientation=best_member_alignment.strand,
                included=False,
                inclusion_reason=rejection_reason,
                trimmed_length=member_input.trimmed_length,
                has_trace_data=member_input.has_trace_data,
                has_qualities=member_input.has_qualities,
                alignment_score=best_member_alignment.score,
                overlap_length=overlap_length,
                percent_identity=percent_identity,
            ),
            seed_column_cells=(),
            insertion_cells_by_bucket={},
        )

    seed_column_cells, insertion_cells_by_bucket = _project_member_against_seed(
        alignment=best_member_alignment.alignment,
        member_input=member_input,
        oriented_sequence=best_member_alignment.sequence,
        oriented_qualities=best_member_alignment.qualities,
        strand=best_member_alignment.strand,
        seed_length=len(seed_member_input.display_sequence),
    )
    return _PlacedMultiMember(
        member=MultiAssemblyMember(
            member_index=member_input.member_index,
            source_filename=member_input.source_filename,
            display_name=member_input.display_name,
            chosen_orientation=best_member_alignment.strand,
            included=True,
            inclusion_reason=None,
            trimmed_length=member_input.trimmed_length,
            has_trace_data=member_input.has_trace_data,
            has_qualities=member_input.has_qualities,
            alignment_score=best_member_alignment.score,
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
        trace_x=trace_position_for_oriented_query_index(
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
