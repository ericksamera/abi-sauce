from __future__ import annotations

import pytest

from abi_sauce.chromatogram import (
    ChromatogramBaseCall,
    ChromatogramChannel,
    ChromatogramQualitySegment,
    ChromatogramTrimBoundaries,
    ChromatogramView,
)
from abi_sauce.chromatogram_figure import build_chromatogram_figure


def make_view(
    *,
    include_quality: bool = True,
    trim_boundaries: ChromatogramTrimBoundaries | None = None,
    retained_sample_range: tuple[float, float] | None = None,
    has_any_retained_samples: bool = True,
) -> ChromatogramView:
    channels = (
        ChromatogramChannel(
            data_key="DATA9",
            base="G",
            color="black",
            signal=(0, 5, 12, 4, 1),
        ),
        ChromatogramChannel(
            data_key="DATA10",
            base="A",
            color="green",
            signal=(1, 2, 6, 2, 1),
        ),
    )
    base_calls = (
        ChromatogramBaseCall(base_index=0, base="A", position=1, color="green"),
        ChromatogramBaseCall(base_index=1, base="C", position=3, color="blue"),
    )
    quality_segments = (
        ChromatogramQualitySegment(
            base_index=0,
            left=0.0,
            right=2.0,
            center=1.0,
            width=2.0,
            quality=40,
        ),
        ChromatogramQualitySegment(
            base_index=1,
            left=2.0,
            right=4.0,
            center=3.0,
            width=2.0,
            quality=28,
        ),
    )
    return ChromatogramView(
        is_renderable=True,
        x_values=(0, 1, 2, 3, 4),
        channels=channels,
        base_calls=base_calls,
        quality_segments=quality_segments if include_quality else (),
        trim_boundaries=trim_boundaries or ChromatogramTrimBoundaries(),
        retained_sample_range=retained_sample_range,
        has_any_retained_samples=has_any_retained_samples,
    )


def test_build_chromatogram_figure_adds_channel_base_and_quality_traces() -> None:
    figure = build_chromatogram_figure(make_view())

    assert [trace.name for trace in figure.data] == [
        "G trace",
        "A trace",
        "Called peaks",
        "Base calls",
        "Quality scores",
    ]
    assert figure.data[0].type == "scattergl"
    assert figure.data[1].type == "scattergl"
    assert figure.data[2].type == "scattergl"
    assert figure.data[2].mode == "markers"
    assert figure.data[3].type == "scattergl"
    assert figure.data[3].mode == "text"
    assert figure.data[4].type == "bar"
    assert tuple(figure.data[4].x) == (1.0, 3.0)
    assert tuple(figure.data[4].width) == (2.0, 2.0)


def test_build_chromatogram_figure_omits_quality_trace_when_overlay_is_unavailable() -> (
    None
):
    figure = build_chromatogram_figure(make_view(include_quality=False))

    assert [trace.name for trace in figure.data] == [
        "G trace",
        "A trace",
        "Called peaks",
        "Base calls",
    ]
    assert "yaxis2" not in figure.layout


def test_build_chromatogram_figure_adds_trim_marker_shapes() -> None:
    figure = build_chromatogram_figure(
        make_view(
            trim_boundaries=ChromatogramTrimBoundaries(
                left=1.5,
                right=3.5,
            )
        )
    )

    assert len(figure.layout.shapes) == 2
    assert [shape.x0 for shape in figure.layout.shapes] == [1.5, 3.5]
    assert all(shape.type == "line" for shape in figure.layout.shapes)
    assert all(shape.yref == "paper" for shape in figure.layout.shapes)


def test_build_chromatogram_figure_sets_initial_x_range_to_a_50_base_window() -> None:
    base_calls = tuple(
        ChromatogramBaseCall(
            base_index=index,
            base="A",
            position=(index + 1) * 10,
            color="green",
        )
        for index in range(80)
    )
    quality_segments = tuple(
        ChromatogramQualitySegment(
            base_index=index,
            left=(index * 10) + 5.0,
            right=(index * 10) + 15.0,
            center=(index * 10) + 10.0,
            width=10.0,
            quality=40,
        )
        for index in range(80)
    )
    figure = build_chromatogram_figure(
        ChromatogramView(
            is_renderable=True,
            x_values=tuple(range(900)),
            channels=(
                ChromatogramChannel(
                    data_key="DATA9",
                    base="G",
                    color="black",
                    signal=tuple(range(900)),
                ),
            ),
            base_calls=base_calls,
            quality_segments=quality_segments,
            trim_boundaries=ChromatogramTrimBoundaries(left=195.0, right=795.0),
            retained_sample_range=(195.0, 795.0),
        )
    )

    assert tuple(figure.layout.xaxis.range) == (185.0, 705.0)
    assert figure.layout.xaxis.rangeslider.visible is True


def test_build_chromatogram_figure_hides_legend_and_uses_pan_dragmode() -> None:
    figure = build_chromatogram_figure(make_view())

    assert figure.layout.showlegend is False
    assert figure.layout.dragmode == "pan"


def test_build_chromatogram_figure_rejects_non_renderable_views() -> None:
    with pytest.raises(ValueError, match="not renderable"):
        build_chromatogram_figure(
            ChromatogramView(
                is_renderable=False,
                render_failure_reason="missing_trace_data",
            )
        )


def test_build_chromatogram_figure_colors_quality_segments_by_retained_region() -> None:
    figure = build_chromatogram_figure(make_view(retained_sample_range=(0.0, 2.0)))

    quality_trace = next(trace for trace in figure.data if trace.type == "bar")

    assert tuple(quality_trace.marker.color) == (
        "rgba(223, 240, 250, 0.95)",
        "rgba(190, 190, 190, 0.45)",
    )


def test_build_chromatogram_figure_adds_retained_trace_overlays_and_mutes_trimmed_regions() -> (
    None
):
    figure = build_chromatogram_figure(make_view(retained_sample_range=(0.0, 2.0)))

    assert [trace.name for trace in figure.data] == [
        "G trace",
        "G trace (retained)",
        "A trace",
        "A trace (retained)",
        "Called peaks",
        "Base calls",
        "Quality scores",
    ]
    assert figure.data[0].line.color == "rgba(0, 0, 0, 0.40)"
    assert figure.data[1].line.color == "black"
    assert figure.data[2].line.color == "rgba(0, 128, 0, 0.45)"
    assert figure.data[3].line.color == "green"
    assert tuple(figure.data[5].textfont.color) == (
        "green",
        "rgba(0, 0, 255, 0.45)",
    )


def test_build_chromatogram_figure_mutes_everything_when_no_samples_are_retained() -> (
    None
):
    figure = build_chromatogram_figure(make_view(has_any_retained_samples=False))

    assert [trace.name for trace in figure.data] == [
        "G trace",
        "A trace",
        "Called peaks",
        "Base calls",
        "Quality scores",
    ]
    assert figure.data[0].line.color == "rgba(0, 0, 0, 0.40)"
    assert figure.data[1].line.color == "rgba(0, 128, 0, 0.45)"
    quality_trace = next(trace for trace in figure.data if trace.type == "bar")
    assert tuple(quality_trace.marker.color) == (
        "rgba(190, 190, 190, 0.45)",
        "rgba(190, 190, 190, 0.45)",
    )


def test_build_chromatogram_figure_uses_dark_theme_colors_when_requested() -> None:
    figure = build_chromatogram_figure(
        make_view(retained_sample_range=(0.0, 2.0)),
        theme_type="dark",
    )

    assert figure.layout.plot_bgcolor == "#0E1117"
    assert figure.layout.font.color == "#FAFAFA"
    assert figure.layout.xaxis.rangeslider.bgcolor == "#0B0F14"
    assert figure.data[0].line.color == "rgba(229, 231, 235, 0.40)"
    assert figure.data[1].line.color == "#E5E7EB"
    assert tuple(figure.data[5].textfont.color) == (
        "#22C55E",
        "rgba(59, 130, 246, 0.45)",
    )
    quality_trace = next(trace for trace in figure.data if trace.type == "bar")
    assert tuple(quality_trace.marker.color) == (
        "rgba(223, 240, 250, 0.55)",
        "rgba(148, 163, 184, 0.28)",
    )
