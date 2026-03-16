from __future__ import annotations

from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.trace_coordinates import (
    display_query_index_for_aligned_query_index,
    display_trace_position,
    display_trim_start,
    oriented_trim_interval,
    raw_trimmed_index_for_display_query_index,
    sanitized_trace_length,
    trace_position_for_oriented_query_index,
)
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    name: str,
    sequence: str,
    base_positions: list[int],
    trace_length: int = 100,
) -> SequenceRecord:
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic coordinate record",
        sequence=sequence,
        source_format="abi",
        qualities=[30] * len(sequence),
        trace_data=TraceData(
            channels={
                "DATA9": [1] * trace_length,
                "DATA10": [2] * trace_length,
                "DATA11": [3] * trace_length,
                "DATA12": [4] * trace_length,
            },
            base_positions=base_positions,
            channel_order="GATC",
        ),
    )


def test_display_trim_start_respects_display_orientation() -> None:
    forward_record = make_record(
        name="forward",
        sequence="AACCGGTT",
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75],
    )
    reverse_record = make_record(
        name="reverse",
        sequence="AACCGGTT",
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75],
    )
    reverse_record.orientation = "reverse_complement"

    trim_result = trim_sequence_record(
        forward_record,
        TrimConfig(left_trim=2, right_trim=1),
    )
    reverse_trim_result = trim_sequence_record(
        reverse_record,
        TrimConfig(left_trim=2, right_trim=1),
    )

    assert (
        display_trim_start(
            trim_result,
            display_orientation="forward",
        )
        == 2
    )
    assert (
        display_trim_start(
            reverse_trim_result,
            display_orientation="reverse_complement",
        )
        == 1
    )


def test_display_query_index_and_raw_trimmed_index_handle_reverse_cases() -> None:
    assert (
        display_query_index_for_aligned_query_index(
            1,
            trimmed_length=5,
            strand="forward",
        )
        == 1
    )
    assert (
        display_query_index_for_aligned_query_index(
            1,
            trimmed_length=5,
            strand="reverse_complement",
        )
        == 3
    )

    assert (
        raw_trimmed_index_for_display_query_index(
            3,
            trimmed_length=5,
            display_orientation="forward",
        )
        == 3
    )
    assert (
        raw_trimmed_index_for_display_query_index(
            3,
            trimmed_length=5,
            display_orientation="reverse_complement",
        )
        == 1
    )


def test_display_trace_position_and_sanitized_trace_length() -> None:
    trace_data = TraceData(
        channels={
            "DATA9": [1, 2, 3],
            "DATA10": [1, 2, 3, 4],
            "DATA11": [],
            "DATA12": [1, 2],
        },
        base_positions=[5, 15],
        channel_order="GATC",
    )

    assert sanitized_trace_length(trace_data) == 2
    assert (
        display_trace_position(
            5,
            trace_length=100,
            display_orientation="forward",
        )
        == 5
    )
    assert (
        display_trace_position(
            5,
            trace_length=100,
            display_orientation="reverse_complement",
        )
        == 94
    )


def test_oriented_trim_interval_respects_row_strand_and_display_orientation() -> None:
    raw_record = make_record(
        name="reverse_display",
        sequence="AACCGGTT",
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75],
    )
    raw_record.orientation = "reverse_complement"
    trim_result = trim_sequence_record(
        raw_record,
        TrimConfig(left_trim=2, right_trim=1),
    )

    assert oriented_trim_interval(
        raw_record=raw_record,
        trim_result=trim_result,
        strand="forward",
    ) == (1, 6)
    assert oriented_trim_interval(
        raw_record=raw_record,
        trim_result=trim_result,
        strand="reverse_complement",
    ) == (2, 7)


def test_trace_position_for_oriented_query_index_maps_through_display_coordinates() -> (
    None
):
    raw_record = make_record(
        name="reverse_display_trace",
        sequence="GTTTT",
        base_positions=[5, 15, 25, 35, 45],
    )
    raw_record.orientation = "reverse_complement"
    trim_result = trim_sequence_record(raw_record, TrimConfig())

    assert (
        trace_position_for_oriented_query_index(
            raw_record=raw_record,
            trim_result=trim_result,
            oriented_query_index=4,
            strand="forward",
        )
        == 94
    )
    assert (
        trace_position_for_oriented_query_index(
            raw_record=raw_record,
            trim_result=trim_result,
            oriented_query_index=0,
            strand="reverse_complement",
        )
        == 94
    )
    assert (
        trace_position_for_oriented_query_index(
            raw_record=raw_record,
            trim_result=trim_result,
            oriented_query_index=None,
            strand="forward",
        )
        is None
    )
