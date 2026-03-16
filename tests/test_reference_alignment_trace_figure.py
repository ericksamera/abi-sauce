from __future__ import annotations

from abi_sauce.models import SequenceRecord, TraceData
from abi_sauce.reference_alignment import align_trimmed_read_to_reference
from abi_sauce.reference_alignment_trace import build_reference_alignment_trace_view
from abi_sauce.reference_alignment_trace_figure import (
    build_reference_alignment_trace_figure,
)
from abi_sauce.trimming import TrimConfig, trim_sequence_record


def make_record(
    *,
    name: str,
    sequence: str,
    qualities: list[int] | None = None,
    base_positions: list[int] | None = None,
) -> SequenceRecord:
    trace_length = max(base_positions[-1] + 10, 200) if base_positions else 200
    return SequenceRecord(
        record_id=f"{name}_id",
        name=name,
        description="synthetic reference-alignment trace record",
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
            base_positions=base_positions or [5, 15, 25, 35],
            channel_order="GATC",
        ),
    )


def make_trace_view():
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
    return build_reference_alignment_trace_view(
        result=result,
        source_filename="trace.ab1",
        raw_record=raw_record,
        trim_result=trim_result,
    )


def test_build_reference_alignment_trace_figure_renders_reference_band_query_row_and_hover_layers() -> (
    None
):
    trace_view = make_trace_view()

    figure = build_reference_alignment_trace_figure(trace_view)

    trace_names = [trace.name for trace in figure.data]
    assert trace_names == [
        "Reference background",
        "Reference",
        "Reference hover",
        "trace background",
        "trace:G",
        "trace:A",
        "trace:T",
        "trace:C",
        "trace bases",
        "trace hover",
    ]
    assert figure.layout.xaxis.title.text == "Alignment column"
    assert tuple(figure.layout.yaxis.ticktext) == ("trace",)
    assert tuple(figure.layout.xaxis.range) == (0.0, float(trace_view.alignment_length))
    assert figure.layout.xaxis.rangeslider.visible is True
    assert figure.layout.showlegend is False
    assert figure.layout.dragmode == "pan"


def test_build_reference_alignment_trace_figure_adds_selected_column_highlight() -> (
    None
):
    trace_view = make_trace_view()

    figure = build_reference_alignment_trace_figure(
        trace_view,
        selected_column_index=4,
    )

    assert len(figure.layout.shapes) > 0
    selected_shape = figure.layout.shapes[0]
    assert selected_shape.type == "rect"
    assert selected_shape.x0 == 3.0
    assert selected_shape.x1 == 4.0


def test_build_reference_alignment_trace_figure_uses_dark_theme_colors() -> None:
    trace_view = make_trace_view()

    figure = build_reference_alignment_trace_figure(
        trace_view,
        theme_type="dark",
    )

    assert figure.layout.plot_bgcolor == "#0E1117"
    assert figure.layout.font.color == "#FAFAFA"
    assert figure.data[4].line.color == "#E5E7EB"
    assert figure.data[5].line.color == "#22C55E"
    assert figure.data[6].line.color == "#EF4444"
    assert figure.data[7].line.color == "#3B82F6"
