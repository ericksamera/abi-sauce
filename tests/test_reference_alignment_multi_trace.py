from __future__ import annotations

from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.reference_alignment_multi import align_trimmed_reads_to_reference
from abi_sauce.reference_alignment_multi_trace import (
    build_reference_multi_alignment_trace_view,
)
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    base_positions: list[int] | None = None,
) -> SequenceRecord:
    resolved_base_positions = base_positions or list(
        range(5, 5 + (10 * len(sequence)), 10)
    )
    trace_length = max(resolved_base_positions[-1] + 10, 200)
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic shared-reference alignment trace record",
        sequence=sequence,
        source_format="abi",
        qualities=qualities,
        trace_data=TraceData(
            channels={
                "DATA9": list(range(trace_length)),
                "DATA10": [value * 2 for value in range(trace_length)],
                "DATA11": [value * 3 for value in range(trace_length)],
                "DATA12": [value * 4 for value in range(trace_length)],
            },
            base_positions=resolved_base_positions,
            channel_order="GATC",
        ),
    )


def make_trace_view():
    raw_records_by_source_filename = {
        "read_1.ab1": make_record(
            name="read_1",
            sequence="AACCGGTT",
            qualities=[40] * 8,
        ),
        "read_2.ab1": make_record(
            name="read_2",
            sequence="AACCGGTT",
            qualities=[38] * 8,
        ),
        "read_ins.ab1": make_record(
            name="read_ins",
            sequence="AACCGGGTT",
            qualities=[30] * 9,
        ),
    }
    trim_results_by_source_filename = {
        source_filename: trim_sequence_record(record, TrimConfig())
        for source_filename, record in raw_records_by_source_filename.items()
    }
    result = align_trimmed_reads_to_reference(
        source_filenames=("read_1.ab1", "read_2.ab1", "read_ins.ab1"),
        raw_records_by_source_filename=raw_records_by_source_filename,
        trim_results_by_source_filename=trim_results_by_source_filename,
        reference_text=">ref\nAACCGGTT\n",
        strand_policy="forward",
        config=AssemblyConfig(
            min_overlap_length=6,
            min_percent_identity=80.0,
            quality_margin=3,
        ),
    )
    trace_view = build_reference_multi_alignment_trace_view(
        result=result,
        raw_records_by_source_filename=raw_records_by_source_filename,
        trim_results_by_source_filename=trim_results_by_source_filename,
    )
    return result, trace_view


def test_build_reference_multi_alignment_trace_view_builds_stacked_rows() -> None:
    result, trace_view = make_trace_view()

    assert trace_view.reference_name == "ref"
    assert trace_view.alignment_length == len(result.columns)
    assert trace_view.x_range == (0.0, float(len(result.columns)))
    assert len(trace_view.rows) == 3

    first_row = trace_view.rows[0]
    assert first_row.display_name == "read_1"
    assert first_row.aligned_sequence == result.aligned_member_sequences[0]
    assert first_row.has_trace_signal is True
    assert len(first_row.cells) == trace_view.alignment_length

    insertion_column_index = next(
        column.column_index
        for column in result.columns
        if column.anchor_kind == "insertion"
    )
    insertion_cell = first_row.cells[insertion_column_index - 1]
    assert insertion_cell.anchor_kind == "insertion"
    assert insertion_cell.ref_base == "-"
    assert insertion_cell.query_base == "-"

    insertion_row = trace_view.rows[-1]
    inserted_cell = insertion_row.cells[insertion_column_index - 1]
    assert inserted_cell.query_base != "-"
    assert inserted_cell.has_trace_signal is True
    assert inserted_cell.raw_center is not None
