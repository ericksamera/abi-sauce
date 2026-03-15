from __future__ import annotations

from abi_sauce.assembly import (
    AssemblyConfig,
    assemble_trimmed_multi,
    assemble_trimmed_pair,
)
from abi_sauce.assembly_trace import (
    build_multi_assembly_trace_view,
    build_pairwise_assembly_trace_view,
    resolve_assembly_trace_samples_per_cell,
)
from abi_sauce.models import SequenceRecord, TraceData
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
        description="synthetic assembly trace record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        trace_data=resolved_trace_data,
    )


def test_build_pairwise_assembly_trace_view_builds_row_bands_and_alignment_cells() -> (
    None
):
    left_raw_record = make_record(
        name="left",
        sequence="CCCCAAAAC",
        qualities=[40] * 9,
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75, 85],
    )
    right_raw_record = make_record(
        name="right",
        sequence="GTTTT",
        qualities=[40] * 5,
        base_positions=[10, 20, 30, 40, 50],
    )
    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig())
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
    )

    trace_view = build_pairwise_assembly_trace_view(
        result=result,
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
    )

    assert trace_view.alignment_length == len(result.columns)
    assert trace_view.x_range == (0.0, float(len(result.columns)))
    assert len(trace_view.rows) == 2

    left_row, right_row = trace_view.rows
    assert left_row.y_bottom == 3.0
    assert left_row.y_top == 6.0
    assert right_row.y_bottom == 0.0
    assert right_row.y_top == 3.0
    assert left_row.strand == "forward"
    assert right_row.strand == "reverse-complement"
    assert left_row.aligned_sequence == result.aligned_left
    assert right_row.aligned_sequence == result.aligned_right
    assert len(left_row.cells) == len(result.columns)
    assert len(right_row.cells) == len(result.columns)

    right_gap_cells = [cell for cell in right_row.cells if cell.is_gap]
    assert right_gap_cells
    assert all(cell.channels == () for cell in right_gap_cells)

    populated_left_cells = [cell for cell in left_row.cells if cell.has_trace_signal]
    assert populated_left_cells
    first_cell = populated_left_cells[0]
    assert len(first_cell.channels) == 4
    assert all(len(channel.x_values) == 16 for channel in first_cell.channels)
    assert all(len(channel.normalized_signal) == 16 for channel in first_cell.channels)


def test_build_pairwise_assembly_trace_view_carries_raw_base_spans_into_cells() -> None:
    left_raw_record = make_record(
        name="left",
        sequence="AAAACCCCTTTT",
        qualities=[40] * 9 + [35] + [40] * 2,
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115],
    )
    right_raw_record = make_record(
        name="right",
        sequence="CCCCTGTT",
        qualities=[40, 40, 40, 40, 40, 10, 40, 40],
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75],
    )
    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig(left_trim=4))
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=8, min_percent_identity=70.0),
    )

    trace_view = build_pairwise_assembly_trace_view(
        result=result,
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
    )

    left_row, right_row = trace_view.rows
    mismatch_column_index = 5
    left_cell = left_row.cells[mismatch_column_index]
    right_cell = right_row.cells[mismatch_column_index]

    assert left_cell.column_index == 6
    assert left_cell.base == "T"
    assert left_cell.raw_left == 90.0
    assert left_cell.raw_right == 100.0
    assert left_cell.raw_center == 95.0
    assert left_cell.trace_x == 95

    assert right_cell.column_index == 6
    assert right_cell.base == "G"
    assert right_cell.raw_left == 50.0
    assert right_cell.raw_right == 60.0
    assert right_cell.raw_center == 55.0
    assert right_cell.trace_x == 55


def test_build_pairwise_assembly_trace_view_handles_reverse_complement_display_rows() -> (
    None
):
    left_raw_record = make_record(
        name="left",
        sequence="GTTTT",
        qualities=[10, 20, 30, 40, 50],
        base_positions=[5, 15, 25, 35, 45],
    )
    left_raw_record.orientation = "reverse_complement"
    right_raw_record = make_record(
        name="right",
        sequence="AAAAC",
        qualities=[50, 40, 30, 20, 10],
        base_positions=[5, 15, 25, 35, 45],
    )
    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig())
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=5, min_percent_identity=90.0),
    )

    trace_view = build_pairwise_assembly_trace_view(
        result=result,
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
    )

    left_row, right_row = trace_view.rows
    assert result.aligned_left == "AAAAC"
    assert result.aligned_right == "AAAAC"
    assert left_row.aligned_sequence == "AAAAC"
    assert right_row.aligned_sequence == "AAAAC"
    assert left_row.cells[-1].trace_x == 194
    assert left_row.cells[-1].raw_center == 194.0


def test_build_pairwise_assembly_trace_view_degrades_without_renderable_trace_data() -> (
    None
):
    left_raw_record = make_record(
        name="left",
        sequence="CCCCAAAAC",
        qualities=[40] * 9,
        trace_data=None,
    )
    right_raw_record = make_record(
        name="right",
        sequence="GTTTT",
        qualities=[40] * 5,
        trace_data=None,
    )
    left_trim_result = trim_sequence_record(left_raw_record, TrimConfig())
    right_trim_result = trim_sequence_record(right_raw_record, TrimConfig())

    result = assemble_trimmed_pair(
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=90.0),
    )

    trace_view = build_pairwise_assembly_trace_view(
        result=result,
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
    )

    left_row, right_row = trace_view.rows
    assert left_row.has_trace_signal is False
    assert right_row.has_trace_signal is False
    assert all(cell.channels == () for cell in left_row.cells)
    assert all(cell.channels == () for cell in right_row.cells)
    assert left_row.aligned_sequence == result.aligned_left
    assert right_row.aligned_sequence == result.aligned_right


def test_build_multi_assembly_trace_view_builds_one_row_per_included_member() -> None:
    seed_raw_record = make_record(
        name="seed",
        sequence="AACCGGTTA",
        qualities=[40] * 9,
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75, 85],
    )
    shifted_raw_record = make_record(
        name="shifted",
        sequence="ACCGGATTA",
        qualities=[35] * 9,
        base_positions=[15, 25, 35, 45, 55, 65, 75, 85, 95],
    )
    reverse_raw_record = make_record(
        name="reverse",
        sequence="TAACCGGTT",
        qualities=[30] * 9,
        base_positions=[5, 15, 25, 35, 45, 55, 65, 75, 85],
    )

    trim_results = {
        "seed.ab1": trim_sequence_record(seed_raw_record, TrimConfig()),
        "shifted.ab1": trim_sequence_record(shifted_raw_record, TrimConfig()),
        "reverse.ab1": trim_sequence_record(reverse_raw_record, TrimConfig()),
    }
    raw_records = {
        "seed.ab1": seed_raw_record,
        "shifted.ab1": shifted_raw_record,
        "reverse.ab1": reverse_raw_record,
    }
    result = assemble_trimmed_multi(
        source_filenames=("seed.ab1", "shifted.ab1", "reverse.ab1"),
        raw_records_by_source_filename=raw_records,
        trim_results_by_source_filename=trim_results,
        config=AssemblyConfig(min_overlap_length=4, min_percent_identity=70.0),
    )

    trace_view = build_multi_assembly_trace_view(
        result=result,
        raw_records_by_source_filename=raw_records,
        trim_results_by_source_filename=trim_results,
    )

    assert trace_view.alignment_length == len(result.columns)
    assert trace_view.x_range == (0.0, float(len(result.columns)))
    assert len(trace_view.rows) == 3
    assert tuple(row.y_bottom for row in trace_view.rows) == (6.0, 3.0, 0.0)
    assert tuple(row.y_top for row in trace_view.rows) == (9.0, 6.0, 3.0)
    assert tuple(row.strand for row in trace_view.rows) == (
        "forward",
        "forward",
        "reverse-complement",
    )
    assert tuple(row.aligned_sequence for row in trace_view.rows) == (
        "AACCGG-TTA",
        "-ACCGGATTA",
        "AACCGG-TTA",
    )

    insertion_column_index = 6
    seed_cell = trace_view.rows[0].cells[insertion_column_index]
    shifted_cell = trace_view.rows[1].cells[insertion_column_index]
    reverse_cell = trace_view.rows[2].cells[insertion_column_index]

    assert seed_cell.base == "-"
    assert shifted_cell.base == "A"
    assert reverse_cell.base == "-"
    assert seed_cell.channels == ()
    assert reverse_cell.channels == ()
    assert shifted_cell.has_trace_signal is True


def test_resolve_assembly_trace_samples_per_cell_keeps_default_for_small_views() -> (
    None
):
    assert (
        resolve_assembly_trace_samples_per_cell(
            alignment_length=24,
            row_count=2,
            requested_samples_per_cell=None,
        )
        == 16
    )
    assert (
        resolve_assembly_trace_samples_per_cell(
            alignment_length=24,
            row_count=2,
            requested_samples_per_cell=12,
        )
        == 12
    )


def test_resolve_assembly_trace_samples_per_cell_reduces_for_large_views() -> None:
    assert (
        resolve_assembly_trace_samples_per_cell(
            alignment_length=400,
            row_count=2,
            requested_samples_per_cell=None,
        )
        == 8
    )
    assert (
        resolve_assembly_trace_samples_per_cell(
            alignment_length=200,
            row_count=6,
            requested_samples_per_cell=None,
        )
        == 8
    )
