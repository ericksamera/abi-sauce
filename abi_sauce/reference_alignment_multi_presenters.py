from __future__ import annotations

from abi_sauce.reference_alignment_types import ReferenceMultiAlignmentResult


def reference_multi_columns_to_rows(
    result: ReferenceMultiAlignmentResult,
    *,
    include_reference_matches: bool = False,
) -> list[dict[str, object]]:
    """Convert shared-reference columns into table-friendly support rows."""
    rows: list[dict[str, object]] = []
    for column in result.columns:
        if (
            not include_reference_matches
            and column.anchor_kind == "reference"
            and column.matches_reference is True
            and column.consensus_base == column.ref_base
            and not column.ambiguous
        ):
            continue
        rows.append(
            {
                "column": column.column_index,
                "anchor": column.anchor_kind,
                "anchor_index": column.anchor_index,
                "ref_pos": column.ref_pos,
                "ref_base": column.ref_base,
                "consensus_base": column.consensus_base,
                "resolution": column.resolution,
                "support": ", ".join(
                    f"{base}:{count}" for base, count in column.support_counts
                ),
                "non_gap_members": column.non_gap_member_count,
                "gap_members": column.gap_member_count,
                "matches_reference": column.matches_reference,
            }
        )
    return rows


def reference_multi_members_to_rows(
    result: ReferenceMultiAlignmentResult,
) -> list[dict[str, object]]:
    """Convert shared-reference member summaries into table-friendly rows."""
    return [
        {
            "member": member.display_name,
            "source_filename": member.source_filename,
            "included": member.included,
            "reason": member.inclusion_reason,
            "strand": member.chosen_strand,
            "score": member.alignment_score,
            "overlap_length": member.overlap_length,
            "percent_identity": member.percent_identity,
            "mismatches": member.mismatch_count,
            "insertions": member.insertion_count,
            "deletions": member.deletion_count,
            "trimmed_length": member.trimmed_length,
        }
        for member in result.members
    ]


def format_reference_multi_alignment_block(
    result: ReferenceMultiAlignmentResult,
    *,
    line_width: int = 80,
) -> str:
    """Render a gapped block with reference, members, and consensus."""
    lines: list[str] = []
    member_lines = list(result.aligned_member_sequences)
    for start in range(0, len(result.gapped_reference), line_width):
        end = start + line_width
        lines.append(result.gapped_reference[start:end])
        for aligned_member_sequence in member_lines:
            lines.append(aligned_member_sequence[start:end])
        lines.append(result.gapped_consensus[start:end])
        lines.append("")
    return "\n".join(lines).rstrip()


__all__ = [
    "reference_multi_columns_to_rows",
    "reference_multi_members_to_rows",
    "format_reference_multi_alignment_block",
]
