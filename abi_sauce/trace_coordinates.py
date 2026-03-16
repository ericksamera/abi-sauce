from __future__ import annotations

from typing import Literal

from abi_sauce.models import SequenceOrientation, SequenceRecord, TraceData
from abi_sauce.trimming import TrimResult

TraceStrand = Literal["forward", "reverse_complement"]


def display_trim_start(
    trim_result: TrimResult,
    *,
    display_orientation: SequenceOrientation,
) -> int:
    """Return the display-space trim start for one trimmed read."""
    if display_orientation == "forward":
        return trim_result.bases_removed_left
    return trim_result.bases_removed_right


def oriented_trim_interval(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    strand: TraceStrand,
) -> tuple[int, int]:
    """Return the oriented full-sequence trim interval for one assembly row."""
    display_start = display_trim_start(
        trim_result,
        display_orientation=raw_record.orientation,
    )
    display_end = display_start + trim_result.trimmed_length
    if strand == "forward":
        return (display_start, display_end)

    full_length = len(raw_record.sequence)
    return (
        full_length - display_end,
        full_length - display_start,
    )


def display_query_index_for_aligned_query_index(
    aligned_query_index: int,
    *,
    trimmed_length: int,
    strand: TraceStrand,
) -> int:
    """Convert one oriented query index into display-space trimmed coordinates."""
    if strand == "forward":
        return aligned_query_index
    return trimmed_length - 1 - aligned_query_index


def raw_trimmed_index_for_display_query_index(
    display_query_index: int,
    *,
    trimmed_length: int,
    display_orientation: SequenceOrientation,
) -> int:
    """Convert one display-space trimmed index back into raw trimmed coordinates."""
    if display_orientation == "forward":
        return display_query_index
    return trimmed_length - 1 - display_query_index


def sanitized_trace_length(trace_data: TraceData | None) -> int:
    """Return the shared usable trace length across non-empty channels."""
    if trace_data is None:
        return 0

    channel_lengths = [
        len(signal) for signal in trace_data.channels.values() if len(signal) > 0
    ]
    if not channel_lengths:
        return 0
    return min(channel_lengths)


def display_trace_position(
    raw_trace_position: int,
    *,
    trace_length: int,
    display_orientation: SequenceOrientation,
) -> int:
    """Project one raw trace position into the record's display orientation."""
    if display_orientation == "forward":
        return raw_trace_position
    return int(float(trace_length - 1) - float(raw_trace_position))


def trace_position_for_oriented_query_index(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    oriented_query_index: int | None,
    strand: TraceStrand,
) -> int | None:
    """Map one oriented trimmed query index onto a display-space trace x value."""
    if oriented_query_index is None:
        return None

    trace_data = raw_record.trace_data
    if trace_data is None:
        return None

    trimmed_length = trim_result.trimmed_length
    if trimmed_length <= 0:
        return None

    raw_start = trim_result.bases_removed_left
    display_query_index = display_query_index_for_aligned_query_index(
        oriented_query_index,
        trimmed_length=trimmed_length,
        strand=strand,
    )
    raw_trimmed_index = raw_trimmed_index_for_display_query_index(
        display_query_index,
        trimmed_length=trimmed_length,
        display_orientation=raw_record.orientation,
    )
    raw_base_index = raw_start + raw_trimmed_index

    if raw_base_index < 0 or raw_base_index >= len(trace_data.base_positions):
        return None

    raw_trace_position = int(trace_data.base_positions[raw_base_index])
    trace_length = sanitized_trace_length(trace_data)
    if trace_length <= 0:
        return raw_trace_position
    if raw_trace_position < 0 or raw_trace_position >= trace_length:
        return None
    return display_trace_position(
        raw_trace_position,
        trace_length=trace_length,
        display_orientation=raw_record.orientation,
    )
