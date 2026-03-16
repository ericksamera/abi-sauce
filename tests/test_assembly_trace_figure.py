from __future__ import annotations

from abi_sauce.assembly_multi import assemble_trimmed_multi
from abi_sauce.assembly_pairwise import assemble_trimmed_pair
from abi_sauce.assembly_types import AssemblyConfig
from abi_sauce.assembly_trace import (
    build_multi_assembly_trace_view,
    build_pairwise_assembly_trace_view,
)
from abi_sauce.assembly_trace_figure import build_assembly_trace_figure
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


def make_trace_view():
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

    return build_pairwise_assembly_trace_view(
        result=result,
        left_source_filename="left.ab1",
        left_raw_record=left_raw_record,
        left_trim_result=left_trim_result,
        right_source_filename="right.ab1",
        right_raw_record=right_raw_record,
        right_trim_result=right_trim_result,
    )


def test_build_assembly_trace_figure_renders_consensus_rows_and_hover_layers() -> None:
    trace_view = make_trace_view()

    figure = build_assembly_trace_figure(trace_view)

    trace_names = [trace.name for trace in figure.data]
    assert trace_names == [
        "Consensus background",
        "Consensus",
        "Consensus hover",
        "left background",
        "left:G",
        "left:A",
        "left:T",
        "left:C",
        "left bases",
        "left hover",
        "right background",
        "right:G",
        "right:A",
        "right:T",
        "right:C",
        "right bases",
        "right hover",
    ]
    assert figure.layout.xaxis.title.text == "Alignment column"
    assert tuple(figure.layout.yaxis.ticktext) == ("left", "right")
    assert tuple(figure.layout.xaxis.range) == (0.0, float(trace_view.alignment_length))
    assert figure.layout.xaxis.rangeslider.visible is True
    assert figure.layout.showlegend is False
    assert figure.layout.dragmode == "pan"


def test_build_assembly_trace_figure_adds_selected_column_highlight() -> None:
    trace_view = make_trace_view()

    figure = build_assembly_trace_figure(
        trace_view,
        selected_column_index=6,
    )

    assert len(figure.layout.shapes) > 0
    selected_shape = figure.layout.shapes[0]
    assert selected_shape.type == "rect"
    assert selected_shape.x0 == 5.0
    assert selected_shape.x1 == 6.0


def test_build_assembly_trace_figure_uses_dark_theme_colors() -> None:
    trace_view = make_trace_view()

    figure = build_assembly_trace_figure(
        trace_view,
        theme_type="dark",
    )

    assert figure.layout.plot_bgcolor == "#0E1117"
    assert figure.layout.font.color == "#FAFAFA"
    assert figure.data[4].line.color == "#E5E7EB"
    assert figure.data[5].line.color == "#22C55E"
    assert figure.data[6].line.color == "#EF4444"
    assert figure.data[7].line.color == "#3B82F6"


def make_multi_trace_view():
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
    return build_multi_assembly_trace_view(
        result=result,
        raw_records_by_source_filename=raw_records,
        trim_results_by_source_filename=trim_results,
    )


def test_build_assembly_trace_figure_supports_multi_trace_views() -> None:
    trace_view = make_multi_trace_view()

    figure = build_assembly_trace_figure(trace_view)

    trace_names = [trace.name for trace in figure.data]
    assert trace_names[:3] == [
        "Consensus background",
        "Consensus",
        "Consensus hover",
    ]
    assert "seed background" in trace_names
    assert "shifted background" in trace_names
    assert "reverse background" in trace_names
    assert tuple(figure.layout.yaxis.ticktext) == ("seed", "shifted", "reverse")

    consensus_hover = next(
        trace for trace in figure.data if trace.name == "Consensus hover"
    )
    assert all("support=" in value for value in consensus_hover.customdata)
    assert all("left=" not in value for value in consensus_hover.customdata)
