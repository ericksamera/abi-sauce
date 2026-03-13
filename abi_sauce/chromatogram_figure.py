from __future__ import annotations

from typing import Final

import plotly.graph_objects as go

from abi_sauce.chromatogram import ChromatogramView

_DEFAULT_TRIM_MARKER_COLOR: Final[str] = "#666666"
_DEFAULT_QUALITY_COLOR: Final[str] = "rgba(120, 120, 120, 0.35)"


def build_chromatogram_figure(view: ChromatogramView) -> go.Figure:
    """Build a deterministic Plotly figure from a normalized chromatogram view."""
    if not view.is_renderable:
        raise ValueError("Chromatogram view is not renderable")

    max_signal = _max_signal(view)
    base_label_y = _base_label_y(max_signal)

    figure = go.Figure()

    for channel in view.channels:
        figure.add_trace(
            go.Scatter(
                x=view.x_values,
                y=channel.signal,
                mode="lines",
                name=f"{channel.base} trace",
                line={"color": channel.color},
                hovertemplate=(
                    "sample=%{x}<br>signal=%{y}<extra>" + channel.base + "</extra>"
                ),
            )
        )

    figure.add_trace(
        go.Scatter(
            x=[base_call.position for base_call in view.base_calls],
            y=[base_label_y] * len(view.base_calls),
            text=[base_call.base for base_call in view.base_calls],
            mode="text",
            name="Base calls",
            hovertemplate="base=%{text}<br>sample=%{x}<extra></extra>",
            cliponaxis=False,
        )
    )

    if view.has_quality_overlay:
        figure.add_trace(
            go.Bar(
                x=[point.position for point in view.quality_points],
                y=[point.quality for point in view.quality_points],
                name="Quality scores",
                yaxis="y2",
                marker={"color": _DEFAULT_QUALITY_COLOR},
                hovertemplate="sample=%{x}<br>quality=%{y}<extra></extra>",
            )
        )

    for boundary in (view.trim_boundaries.left, view.trim_boundaries.right):
        if boundary is None:
            continue
        figure.add_shape(
            type="line",
            xref="x",
            yref="paper",
            x0=boundary,
            x1=boundary,
            y0=0,
            y1=1,
            line={
                "color": _DEFAULT_TRIM_MARKER_COLOR,
                "width": 2,
                "dash": "dash",
            },
        )

    figure.update_layout(
        xaxis={
            "title": "Sample index",
            "range": list(_initial_x_range(view)),
            "rangeslider": {"visible": True},
        },
        yaxis={
            "title": "Signal intensity",
            "range": [0.0, base_label_y * 1.08],
        },
        hovermode="x unified",
        bargap=0.0,
        showlegend=True,
    )

    if view.has_quality_overlay:
        max_quality = max(point.quality for point in view.quality_points)
        figure.update_layout(
            yaxis2={
                "title": "PHRED quality",
                "overlaying": "y",
                "side": "right",
                "range": [0.0, max(max_quality * 1.1, 1.0)],
                "showgrid": False,
                "zeroline": False,
            }
        )

    return figure


def _max_signal(view: ChromatogramView) -> float:
    if not view.channels:
        return 1.0
    return float(max(max(channel.signal) for channel in view.channels))


def _base_label_y(max_signal: float) -> float:
    if max_signal <= 0:
        return 1.0
    return max_signal * 1.05


def _initial_x_range(view: ChromatogramView) -> tuple[float, float]:
    if view.trace_length <= 1:
        return (0.0, float(max(view.trace_length - 1, 0)))

    left_boundary = view.trim_boundaries.left
    right_boundary = view.trim_boundaries.right

    if (
        left_boundary is not None
        and right_boundary is not None
        and right_boundary > left_boundary
    ):
        trimmed_width = right_boundary - left_boundary
        padding = max(10.0, trimmed_width * 0.2)
        return (
            _clamp_x(left_boundary - padding, view=view),
            _clamp_x(right_boundary + padding, view=view),
        )

    return (0.0, float(view.trace_length - 1))


def _clamp_x(value: float, *, view: ChromatogramView) -> float:
    return max(0.0, min(value, float(view.trace_length - 1)))
