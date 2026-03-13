from __future__ import annotations

from typing import Final

import plotly.graph_objects as go

from abi_sauce.chromatogram import (
    ChromatogramBaseCall,
    ChromatogramChannel,
    ChromatogramQualitySegment,
    ChromatogramView,
)

_DEFAULT_TRIM_MARKER_COLOR: Final[str] = "#666666"
_DEFAULT_RETAINED_QUALITY_COLOR: Final[str] = "rgba(223, 240, 250, 0.95)"
_DEFAULT_TRIMMED_QUALITY_COLOR: Final[str] = "rgba(190, 190, 190, 0.45)"
_MUTED_TRACE_COLORS: Final[dict[str, str]] = {
    "green": "rgba(0, 128, 0, 0.45)",
    "blue": "rgba(0, 0, 255, 0.45)",
    "black": "rgba(0, 0, 0, 0.40)",
    "red": "rgba(255, 0, 0, 0.45)",
    "magenta": "rgba(255, 58, 255, 0.45)",
}
_DEFAULT_INITIAL_VISIBLE_BASES: Final[int] = 50
_DEFAULT_INITIAL_BASE_PADDING_MULTIPLIER: Final[float] = 1.0


def build_chromatogram_figure(view: ChromatogramView) -> go.Figure:
    """Build a deterministic Plotly figure from a normalized chromatogram view."""
    if not view.is_renderable:
        raise ValueError("Chromatogram view is not renderable")

    max_signal = _max_signal(view)
    base_label_y = _base_label_y(max_signal)

    figure = go.Figure()

    for channel in view.channels:
        _add_channel_traces(figure, view=view, channel=channel)

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
            textfont={
                "color": [
                    _base_call_text_color(view, base_call)
                    for base_call in view.base_calls
                ]
            },
            hovertemplate="base=%{text}<br>sample=%{x}<extra></extra>",
        )
    )

    if view.has_quality_overlay:
        figure.add_trace(
            go.Bar(
                x=[segment.center for segment in view.quality_segments],
                y=[segment.quality for segment in view.quality_segments],
                width=[segment.width for segment in view.quality_segments],
                customdata=[
                    [
                        segment.left,
                        segment.right,
                        _quality_segment_status(view, segment),
                    ]
                    for segment in view.quality_segments
                ],
                name="Quality scores",
                yaxis="y2",
                marker={
                    "color": [
                        _quality_segment_color(view, segment)
                        for segment in view.quality_segments
                    ]
                },
                hovertemplate=(
                    "left=%{customdata[0]:.1f}"
                    "<br>right=%{customdata[1]:.1f}"
                    "<br>status=%{customdata[2]}"
                    "<br>quality=%{y}"
                    "<extra></extra>"
                ),
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
                "width": 1,
                "dash": "longdash",
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
        dragmode="pan",
        showlegend=False,
    )

    if view.has_quality_overlay:
        max_quality = max(segment.quality for segment in view.quality_segments)
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

    if not view.base_calls:
        return (0.0, float(view.trace_length - 1))

    average_base_width = _average_base_width(view)
    first_visible_base_index = _first_visible_base_index(view)
    last_visible_base_index = min(
        first_visible_base_index + _DEFAULT_INITIAL_VISIBLE_BASES - 1,
        len(view.base_calls) - 1,
    )
    padding = average_base_width * _DEFAULT_INITIAL_BASE_PADDING_MULTIPLIER

    left_edge = _base_left_edge(
        view,
        base_call_index=first_visible_base_index,
        average_base_width=average_base_width,
    )
    right_edge = _base_right_edge(
        view,
        base_call_index=last_visible_base_index,
        average_base_width=average_base_width,
    )
    left = _clamp_x(left_edge - padding, view=view)
    right = _clamp_x(right_edge + padding, view=view)

    if right <= left:
        return (left, _clamp_x(left + max(average_base_width, 1.0), view=view))
    return (left, right)


def _average_base_width(view: ChromatogramView) -> float:
    segment_widths = [
        segment.width for segment in view.quality_segments if segment.width > 0
    ]
    if segment_widths:
        return sum(segment_widths) / len(segment_widths)

    base_steps = [
        float(right.position - left.position)
        for left, right in zip(view.base_calls, view.base_calls[1:])
        if right.position > left.position
    ]
    if base_steps:
        return sum(base_steps) / len(base_steps)

    if view.trace_length > 1 and view.base_calls:
        return float(view.trace_length - 1) / len(view.base_calls)

    return 1.0


def _first_visible_base_index(view: ChromatogramView) -> int:
    retained_range = view.retained_sample_range
    if retained_range is None or not view.has_any_retained_samples:
        return 0

    for index, base_call in enumerate(view.base_calls):
        if float(base_call.position) >= retained_range[0]:
            return index
    return max(len(view.base_calls) - 1, 0)


def _base_left_edge(
    view: ChromatogramView,
    *,
    base_call_index: int,
    average_base_width: float,
) -> float:
    base_call = view.base_calls[base_call_index]
    quality_segment = _quality_segment_for_base(
        view,
        base_index=base_call.base_index,
    )
    if quality_segment is not None:
        return quality_segment.left

    if base_call_index <= 0:
        return float(base_call.position) - (average_base_width / 2.0)

    previous_base_call = view.base_calls[base_call_index - 1]
    return (previous_base_call.position + base_call.position) / 2.0


def _base_right_edge(
    view: ChromatogramView,
    *,
    base_call_index: int,
    average_base_width: float,
) -> float:
    base_call = view.base_calls[base_call_index]
    quality_segment = _quality_segment_for_base(
        view,
        base_index=base_call.base_index,
    )
    if quality_segment is not None:
        return quality_segment.right

    if base_call_index >= len(view.base_calls) - 1:
        return float(base_call.position) + (average_base_width / 2.0)

    next_base_call = view.base_calls[base_call_index + 1]
    return (base_call.position + next_base_call.position) / 2.0


def _quality_segment_for_base(
    view: ChromatogramView,
    *,
    base_index: int,
) -> ChromatogramQualitySegment | None:
    for quality_segment in view.quality_segments:
        if quality_segment.base_index == base_index:
            return quality_segment
    return None


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


def _add_channel_traces(
    figure: go.Figure,
    *,
    view: ChromatogramView,
    channel: ChromatogramChannel,
) -> None:
    retained_range = view.retained_sample_range
    if retained_range is None and view.has_any_retained_samples:
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
        return

    figure.add_trace(
        go.Scattergl(
            x=view.x_values,
            y=channel.signal,
            mode="lines",
            name=f"{channel.base} trace",
            line={"color": _muted_trace_color(channel.color)},
            hoverinfo="skip",
            hovertemplate=None,
        )
    )

    if retained_range is None:
        return

    retained_x = list(view.x_values)
    retained_y = [
        signal if retained_range[0] <= sample <= retained_range[1] else None
        for sample, signal in zip(view.x_values, channel.signal, strict=True)
    ]
    if not any(value is not None for value in retained_y):
        return

    figure.add_trace(
        go.Scattergl(
            x=retained_x,
            y=retained_y,
            mode="lines",
            name=f"{channel.base} trace (retained)",
            showlegend=False,
            line={"color": channel.color},
            hoverinfo="skip",
            hovertemplate=None,
        )
    )


def _muted_trace_color(color: str) -> str:
    return _MUTED_TRACE_COLORS.get(color, color)


def _is_sample_retained(view: ChromatogramView, sample: float) -> bool:
    retained_range = view.retained_sample_range
    if retained_range is None:
        return view.has_any_retained_samples
    return retained_range[0] <= sample <= retained_range[1]


def _quality_segment_status(
    view: ChromatogramView,
    segment: ChromatogramQualitySegment,
) -> str:
    return "retained" if _is_sample_retained(view, segment.center) else "trimmed"


def _quality_segment_color(
    view: ChromatogramView,
    segment: ChromatogramQualitySegment,
) -> str:
    if _quality_segment_status(view, segment) == "retained":
        return _DEFAULT_RETAINED_QUALITY_COLOR
    return _DEFAULT_TRIMMED_QUALITY_COLOR


def _base_call_text_color(
    view: ChromatogramView,
    base_call: ChromatogramBaseCall,
) -> str:
    if _is_sample_retained(view, float(base_call.position)):
        return base_call.color
    return _muted_trace_color(base_call.color)
