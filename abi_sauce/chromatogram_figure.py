from __future__ import annotations

from typing import Final

import plotly.graph_objects as go

from abi_sauce.chromatogram import ChromatogramBaseCall, ChromatogramView

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
            go.Scattergl(
                x=view.x_values,
                y=channel.signal,
                mode="lines",
                name=f"{channel.base} trace",
                line={"color": channel.color},
                hoverinfo="skip",
                hovertemplate=None,
            )
        )

    called_peak_positions: list[int] = []
    called_peak_heights: list[int] = []
    called_peak_bases: list[str] = []
    called_peak_colors: list[str] = []
    for base_call in view.base_calls:
        peak_height = _resolve_called_peak_height(view, base_call)
        if peak_height is None:
            continue
        called_peak_positions.append(base_call.position)
        called_peak_heights.append(peak_height)
        called_peak_bases.append(base_call.base)
        called_peak_colors.append(base_call.color)

    if called_peak_positions:
        figure.add_trace(
            go.Scattergl(
                x=called_peak_positions,
                y=called_peak_heights,
                mode="markers",
                name="Called peaks",
                showlegend=False,
                customdata=called_peak_bases,
                hovertemplate=(
                    "base=%{customdata}<br>sample=%{x}<br>signal=%{y}<extra></extra>"
                ),
                marker={
                    "size": 10,
                    "opacity": 0,
                },
            )
        )

    figure.add_trace(
        go.Scattergl(
            x=[base_call.position for base_call in view.base_calls],
            y=[base_label_y] * len(view.base_calls),
            text=[base_call.base for base_call in view.base_calls],
            mode="text",
            name="Base calls",
            hovertemplate="base=%{text}<br>sample=%{x}<extra></extra>",
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


def _resolve_called_peak_height(
    view: ChromatogramView,
    base_call: ChromatogramBaseCall,
) -> int | None:
    for channel in view.channels:
        if channel.base == base_call.base and base_call.position < len(channel.signal):
            return channel.signal[base_call.position]

    fallback_heights = [
        channel.signal[base_call.position]
        for channel in view.channels
        if base_call.position < len(channel.signal)
    ]
    if not fallback_heights:
        return None
    return max(fallback_heights)
