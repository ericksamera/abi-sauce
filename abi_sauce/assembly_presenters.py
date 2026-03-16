from __future__ import annotations

from abi_sauce.assembly_types import AssemblyResult


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
