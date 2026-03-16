from __future__ import annotations

from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.reference_alignment import align_trimmed_read_to_reference
from abi_sauce.reference_alignment_trace import build_reference_alignment_trace_view
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    base_positions: list[int] | None = None,
    trace_data: TraceData | None = None,
) -> SequenceRecord:
    resolved_trace_data = trace_data
    if resolved_trace_data is None and base_positions is not None:
        trace_length = max(base_positions[-1] + 10, 200)
        resolved_trace_data = TraceData(
            channels={
                "DATA9": list(range(trace_length)),
                "DATA10": [value * 2 for value in range(trace_length)],
                "DATA11": [value * 3 for value in range(trace_length)],
                "DATA12": [value * 4 for value in range(trace_length)],
            },
            base_positions=base_positions,
            channel_order="GATC",
        )
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic reference-alignment trace record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        trace_data=resolved_trace_data,
    )


def test_build_reference_alignment_trace_view_projects_columns_and_trace_spans() -> (
    None
):
    raw_record = make_record(
        name="trace",
        sequence="AACCGGTT",
        qualities=[10, 20, 30, 40, 50, 40, 30, 20],
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75],
    )
    trim_result = trim_sequence_record(
        raw_record,
        TrimConfig(left_trim=2, right_trim=2),
    )
    result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text=">ref\nCCGA\n",
        strand_policy="forward",
    )

    trace_view = build_reference_alignment_trace_view(
        result=result,
        source_filename="trace.ab1",
        raw_record=raw_record,
        trim_result=trim_result,
    )

    assert trace_view.reference_name == "ref"
    assert trace_view.sample_name == "trace"
    assert trace_view.alignment_length == 4
    assert trace_view.x_range == (0.0, 4.0)
    assert len(trace_view.rows) == 1

    row = trace_view.rows[0]
    assert row.display_name == "trace"
    assert row.aligned_sequence == "CCGG"
    assert row.has_trace_signal is True
    assert row.cells[0].raw_left == 20.0
    assert row.cells[0].raw_right == 30.0
    assert row.cells[0].trace_x == 25
    assert row.cells[-1].event_type == "mismatch"
    assert row.cells[-1].quality == 40
    assert row.cells[-1].trace_x == 55


def test_build_reference_alignment_trace_view_handles_reverse_complement_query_rows() -> (
    None
):
    raw_record = make_record(
        name="trace_reverse",
        sequence="GTTTT",
        qualities=[10, 20, 30, 40, 50],
        base_positions=[5, 15, 25, 35, 45],
    )
    raw_record.orientation = "reverse_complement"
    trim_result = trim_sequence_record(raw_record, TrimConfig())
    result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text=">ref\nAAAAT\n",
        strand_policy="forward",
    )

    trace_view = build_reference_alignment_trace_view(
        result=result,
        source_filename="trace_reverse.ab1",
        raw_record=raw_record,
        trim_result=trim_result,
    )

    row = trace_view.rows[0]
    assert result.aligned_query == "AAAAC"
    assert row.aligned_sequence == "AAAAC"
    assert row.cells[-1].event_type == "mismatch"
    assert row.cells[-1].trace_x == 194
    assert row.cells[-1].raw_center == 194.0


def test_build_reference_alignment_trace_view_degrades_without_trace_data() -> None:
    raw_record = make_record(
        name="trace_no_signal",
        sequence="AACCGGTT",
        qualities=[10, 20, 30, 40, 50, 40, 30, 20],
        trace_data=None,
    )
    trim_result = trim_sequence_record(
        raw_record,
        TrimConfig(left_trim=2, right_trim=2),
    )
    result = align_trimmed_read_to_reference(
        raw_record=raw_record,
        trim_result=trim_result,
        reference_text=">ref\nCCGG\n",
        strand_policy="forward",
    )

    trace_view = build_reference_alignment_trace_view(
        result=result,
        source_filename="trace_no_signal.ab1",
        raw_record=raw_record,
        trim_result=trim_result,
    )

    row = trace_view.rows[0]
    assert row.has_trace_signal is False
    assert all(cell.channels == () for cell in row.cells)
