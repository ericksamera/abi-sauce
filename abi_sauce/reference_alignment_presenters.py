from __future__ import annotations

from abi_sauce.reference_alignment_types import AlignmentResult


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
                "column": event.column_index,
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


__all__ = [
    "alignment_events_to_rows",
    "format_alignment_block",
]
