from __future__ import annotations

from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.reference_alignment_multi import align_trimmed_reads_to_reference
from abi_sauce.reference_alignment_multi_trace import (
    build_reference_multi_alignment_trace_view,
)
from abi_sauce.reference_alignment_multi_trace_figure import (
    build_reference_multi_alignment_trace_figure,
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
        description="synthetic shared-reference trace figure record",
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
    return build_reference_multi_alignment_trace_view(
        result=result,
        raw_records_by_source_filename=raw_records_by_source_filename,
        trim_results_by_source_filename=trim_results_by_source_filename,
    )


def test_build_reference_multi_alignment_trace_figure_renders_reference_consensus_and_rows() -> (
    None
):
    trace_view = make_trace_view()

    figure = build_reference_multi_alignment_trace_figure(trace_view)

    trace_names = [trace.name for trace in figure.data]
    assert trace_names[:6] == [
        "Consensus background",
        "Consensus",
        "Consensus hover",
        "Reference background",
        "Reference",
        "Reference hover",
    ]
    assert figure.layout.xaxis.title.text == "Alignment column"
    assert tuple(figure.layout.yaxis.ticktext) == ("read_1", "read_2", "read_ins")
    assert figure.layout.xaxis.rangeslider.visible is True
    assert figure.layout.showlegend is False


def test_build_reference_multi_alignment_trace_figure_adds_selected_column_highlight() -> (
    None
):
    trace_view = make_trace_view()

    figure = build_reference_multi_alignment_trace_figure(
        trace_view,
        selected_column_index=4,
    )

    assert len(figure.layout.shapes) > 0
    selected_shape = figure.layout.shapes[0]
    assert selected_shape.type == "rect"
    assert selected_shape.x0 == 3.0
    assert selected_shape.x1 == 4.0


def test_build_reference_multi_alignment_trace_figure_starts_at_first_trace_column_for_long_reference() -> (
    None
):
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
        reference_text=f">ref\n{'T' * 60}AACCGGTT\n",
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

    figure = build_reference_multi_alignment_trace_figure(trace_view)

    left, right = tuple(figure.layout.xaxis.range)
    first_trace_column_index = min(
        cell.column_index
        for row in trace_view.rows
        for cell in row.cells
        if cell.has_trace_signal
    )

    assert left > 0.0
    assert left <= float(first_trace_column_index)
    assert right <= float(trace_view.alignment_length)


def test_build_reference_multi_alignment_trace_figure_uses_dark_theme_colors() -> None:
    trace_view = make_trace_view()

    figure = build_reference_multi_alignment_trace_figure(
        trace_view,
        theme_type="dark",
    )

    assert figure.layout.plot_bgcolor == "#0E1117"
    assert figure.layout.font.color == "#FAFAFA"
    assert figure.data[6].name == "read_1 background"
