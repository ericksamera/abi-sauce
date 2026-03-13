from __future__ import annotations

import pytest

from abi_sauce.chromatogram import (
    ChromatogramBaseCall,
    ChromatogramChannel,
    ChromatogramQualityPoint,
    ChromatogramTrimBoundaries,
    ChromatogramView,
)
from abi_sauce.chromatogram_figure import build_chromatogram_figure


def make_view(
    *,
    include_quality: bool = True,
    trim_boundaries: ChromatogramTrimBoundaries | None = None,
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
    quality_points = (
        ChromatogramQualityPoint(base_index=0, position=1, quality=40),
        ChromatogramQualityPoint(base_index=1, position=3, quality=28),
    )
    return ChromatogramView(
        is_renderable=True,
        x_values=(0, 1, 2, 3, 4),
        channels=channels,
        base_calls=base_calls,
        quality_points=quality_points if include_quality else (),
        trim_boundaries=trim_boundaries or ChromatogramTrimBoundaries(),
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


def test_build_chromatogram_figure_sets_trim_aware_initial_x_range() -> None:
    figure = build_chromatogram_figure(
        ChromatogramView(
            is_renderable=True,
            x_values=tuple(range(100)),
            channels=(
                ChromatogramChannel(
                    data_key="DATA9",
                    base="G",
                    color="black",
                    signal=tuple(range(100)),
                ),
            ),
            base_calls=(
                ChromatogramBaseCall(
                    base_index=0,
                    base="A",
                    position=25,
                    color="green",
                ),
                ChromatogramBaseCall(
                    base_index=1,
                    base="C",
                    position=75,
                    color="blue",
                ),
            ),
            quality_points=(),
            trim_boundaries=ChromatogramTrimBoundaries(left=25.0, right=75.0),
        )
    )

    assert tuple(figure.layout.xaxis.range) == (15.0, 85.0)
    assert figure.layout.xaxis.rangeslider.visible is True


def test_build_chromatogram_figure_rejects_non_renderable_views() -> None:
    with pytest.raises(ValueError, match="not renderable"):
        build_chromatogram_figure(
            ChromatogramView(
                is_renderable=False,
                render_failure_reason="missing_trace_data",
            )
        )
