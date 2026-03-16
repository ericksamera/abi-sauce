from __future__ import annotations

from abi_sauce.models import SequenceRecord
from abi_sauce.reference_alignment_types import ReferenceMultiAlignmentResult


def format_reference_multi_alignment_fasta(
    result: ReferenceMultiAlignmentResult,
    *,
    consensus_name: str | None = None,
) -> str:
    """Render a shared-reference multi-read alignment as gapped FASTA."""
    lines = [f">{result.reference_name}", result.gapped_reference]

    members_by_index = {
        member.member_index: member for member in result.members if member.included
    }
    for member_index, aligned_sequence in zip(
        result.included_member_indices,
        result.aligned_member_sequences,
        strict=True,
    ):
        member = members_by_index[member_index]
        lines.append(f">{member.display_name}__{member.chosen_strand}")
        lines.append(aligned_sequence)

    resolved_consensus_name = (
        consensus_name.strip()
        if isinstance(consensus_name, str) and consensus_name.strip()
        else f"{result.reference_name}_consensus"
    )
    lines.append(f">{resolved_consensus_name}__gapped_consensus")
    lines.append(result.gapped_consensus)
    return "\n".join(lines) + "\n"


def consensus_record_from_reference_multi_result(
    result: ReferenceMultiAlignmentResult,
    *,
    name: str | None = None,
) -> SequenceRecord:
    """Project a shared-reference multi-read alignment consensus into a SequenceRecord."""
    consensus_name = (
        name.strip()
        if isinstance(name, str) and name.strip()
        else f"{result.reference_name}_consensus"
    )
    return SequenceRecord(
        record_id=consensus_name,
        name=consensus_name,
        description=(
            f"reference-guided consensus from {result.included_member_count} included reads"
        ),
        sequence=result.consensus_sequence,
        source_format="alignment",
        orientation="forward",
        qualities=None,
        trace_data=None,
        annotations={
            "alignment_engine_kind": "reference_multi",
            "reference_name": result.reference_name,
            "member_filenames": [member.source_filename for member in result.members],
            "included_member_count": result.included_member_count,
            "excluded_member_count": result.excluded_member_count,
            "ambiguous_column_count": result.ambiguous_column_count,
            "alignment_accepted": result.accepted,
            "alignment_rejection_reason": result.rejection_reason,
        },
    )


__all__ = [
    "format_reference_multi_alignment_fasta",
    "consensus_record_from_reference_multi_result",
]
