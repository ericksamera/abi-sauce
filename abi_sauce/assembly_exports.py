from __future__ import annotations

from abi_sauce.assembly_types import (
    AssemblyComputationResult,
    AssemblyResult,
    MultiAssemblyResult,
)
from abi_sauce.models import SequenceRecord


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
