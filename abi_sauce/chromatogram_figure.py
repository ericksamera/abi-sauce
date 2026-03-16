from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import plotly.graph_objects as go

from abi_sauce.chromatogram import (
    ChromatogramBaseCall,
    ChromatogramChannel,
    ChromatogramColumn,
    ChromatogramColumnChannel,
    ChromatogramColumnView,
    ChromatogramQualitySegment,
    ChromatogramView,
)

_DEFAULT_INITIAL_VISIBLE_BASES: Final[int] = 50
_DEFAULT_INITIAL_BASE_PADDING_MULTIPLIER: Final[float] = 1.0


@dataclass(frozen=True, slots=True)
class ChromatogramFigureTheme:
    trace_colors: dict[str, str]
    muted_trace_colors: dict[str, str]
    trim_marker_color: str
    retained_quality_color: str
    trimmed_quality_color: str
    plot_bgcolor: str
    paper_bgcolor: str
    font_color: str
    axis_color: str
    grid_color: str
    rangeslider_bgcolor: str
    rangeslider_bordercolor: str


_LIGHT_THEME: Final[ChromatogramFigureTheme] = ChromatogramFigureTheme(
    trace_colors={
        "green": "green",
        "blue": "blue",
        "black": "black",
        "red": "red",
        "magenta": "magenta",
    },
    muted_trace_colors={
        "green": "rgba(0, 128, 0, 0.45)",
        "blue": "rgba(0, 0, 255, 0.45)",
        "black": "rgba(0, 0, 0, 0.40)",
        "red": "rgba(255, 0, 0, 0.45)",
        "magenta": "rgba(255, 58, 255, 0.45)",
    },
    trim_marker_color="#666666",
    retained_quality_color="rgba(223, 240, 250, 0.95)",
    trimmed_quality_color="rgba(190, 190, 190, 0.45)",
    plot_bgcolor="#FFFFFF",
    paper_bgcolor="rgba(0, 0, 0, 0)",
    font_color="#111111",
    axis_color="#444444",
    grid_color="rgba(0, 0, 0, 0.10)",
    rangeslider_bgcolor="#F4F6F8",
    rangeslider_bordercolor="rgba(0, 0, 0, 0.10)",
)

_DARK_THEME: Final[ChromatogramFigureTheme] = ChromatogramFigureTheme(
    trace_colors={
        "green": "#22C55E",
        "blue": "#3B82F6",
        "black": "#E5E7EB",
        "red": "#EF4444",
        "magenta": "#F472B6",
    },
    muted_trace_colors={
        "green": "rgba(34, 197, 94, 0.45)",
        "blue": "rgba(59, 130, 246, 0.45)",
        "black": "rgba(229, 231, 235, 0.40)",
        "red": "rgba(239, 68, 68, 0.45)",
        "magenta": "rgba(244, 114, 182, 0.45)",
    },
    trim_marker_color="#94A3B8",
    retained_quality_color="rgba(223, 240, 250, 0.55)",
    trimmed_quality_color="rgba(148, 163, 184, 0.28)",
    plot_bgcolor="#0E1117",
    paper_bgcolor="rgba(0, 0, 0, 0)",
    font_color="#FAFAFA",
    axis_color="#CBD5E1",
    grid_color="rgba(203, 213, 225, 0.16)",
    rangeslider_bgcolor="#0B0F14",
    rangeslider_bordercolor="rgba(203, 213, 225, 0.16)",
)


def build_chromatogram_figure(
    view: ChromatogramView,
    *,
    theme_type: str = "light",
) -> go.Figure:
    """Build a deterministic Plotly figure from a normalized chromatogram view."""
    if not view.is_renderable:
        raise ValueError("Chromatogram view is not renderable")

    max_signal = _max_signal(view)
    base_label_y = _base_label_y(max_signal)
    theme = _resolve_figure_theme(theme_type)

    figure = go.Figure()

    for channel in view.channels:
        _add_channel_traces(figure, view=view, channel=channel, theme=theme)

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
                    _base_call_text_color(view, base_call, theme)
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
                        _quality_segment_color(view, segment, theme)
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
                "color": theme.trim_marker_color,
                "width": 1,
                "dash": "longdash",
            },
        )

    figure.update_xaxes(
        title_text="Sample index",
        title_font={"color": theme.font_color},
        tickfont={"color": theme.font_color},
        range=list(_initial_x_range(view)),
        showgrid=False,
        linecolor=theme.axis_color,
        rangeslider={
            "visible": True,
            "bgcolor": theme.rangeslider_bgcolor,
            "bordercolor": theme.rangeslider_bordercolor,
        },
    )
    figure.update_yaxes(
        title_text="Signal intensity",
        title_font={"color": theme.font_color},
        tickfont={"color": theme.font_color},
        range=[0.0, base_label_y * 1.08],
        showgrid=True,
        gridcolor=theme.grid_color,
        linecolor=theme.axis_color,
        zeroline=False,
    )
    figure.update_layout(
        hovermode="x unified",
        bargap=0.0,
        dragmode="pan",
        showlegend=False,
        plot_bgcolor=theme.plot_bgcolor,
        paper_bgcolor=theme.paper_bgcolor,
        font={"color": theme.font_color},
    )

    if view.has_quality_overlay:
        max_quality = max(segment.quality for segment in view.quality_segments)
        figure.update_layout(
            yaxis2={
                "title": {
                    "text": "PHRED quality",
                    "font": {"color": theme.font_color},
                },
                "overlaying": "y",
                "side": "right",
                "range": [0.0, max(max_quality * 1.1, 1.0)],
                "showgrid": False,
                "zeroline": False,
                "tickfont": {"color": theme.font_color},
                "linecolor": theme.axis_color,
            }
        )

    return figure


def build_chromatogram_column_figure(
    view: ChromatogramColumnView,
    *,
    theme_type: str = "light",
) -> go.Figure:
    """Build a fixed-width base-column chromatogram figure."""
    if not view.is_renderable:
        raise ValueError("Chromatogram column view is not renderable")

    theme = _resolve_figure_theme(theme_type)
    max_signal = _max_column_signal(view)
    base_label_y = _base_label_y(max_signal)
    figure = go.Figure()

    for base in _column_trace_bases(view):
        _add_column_channel_traces(
            figure,
            view=view,
            base=base,
            theme=theme,
        )

    peak_columns = [column for column in view.columns if column.peak_height is not None]
    if peak_columns:
        figure.add_trace(
            go.Scattergl(
                x=[column.cell_center for column in peak_columns],
                y=[column.peak_height for column in peak_columns],
                mode="markers",
                name="Called peaks",
                showlegend=False,
                customdata=[
                    [
                        column.base,
                        column.query_pos,
                        column.trace_x,
                        _column_status(column),
                    ]
                    for column in peak_columns
                ],
                hovertemplate=(
                    "base=%{customdata[0]}"
                    "<br>base_index=%{customdata[1]}"
                    "<br>sample=%{customdata[2]}"
                    "<br>status=%{customdata[3]}"
                    "<br>signal=%{y}"
                    "<extra></extra>"
                ),
                marker={"size": 10, "opacity": 0},
            )
        )

    figure.add_trace(
        go.Scattergl(
            x=[column.cell_center for column in view.columns],
            y=[base_label_y] * len(view.columns),
            text=[column.base for column in view.columns],
            mode="text",
            name="Base calls",
            textfont={
                "color": [_column_text_color(column, theme) for column in view.columns]
            },
            customdata=[
                [column.query_pos, column.trace_x, _column_status(column)]
                for column in view.columns
            ],
            hovertemplate=(
                "base=%{text}"
                "<br>base_index=%{customdata[0]}"
                "<br>sample=%{customdata[1]}"
                "<br>status=%{customdata[2]}"
                "<extra></extra>"
            ),
        )
    )

    quality_columns = [column for column in view.columns if column.quality is not None]
    if quality_columns:
        figure.add_trace(
            go.Bar(
                x=[column.cell_center for column in quality_columns],
                y=[column.quality for column in quality_columns],
                width=[view.cell_width] * len(quality_columns),
                customdata=[
                    [
                        column.query_pos,
                        column.raw_left,
                        column.raw_right,
                        _column_status(column),
                    ]
                    for column in quality_columns
                ],
                name="Quality scores",
                yaxis="y2",
                marker={
                    "color": [
                        _column_quality_color(column, theme)
                        for column in quality_columns
                    ]
                },
                hovertemplate=(
                    "base_index=%{customdata[0]}"
                    "<br>left=%{customdata[1]:.1f}"
                    "<br>right=%{customdata[2]:.1f}"
                    "<br>status=%{customdata[3]}"
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
                "color": theme.trim_marker_color,
                "width": 1,
                "dash": "longdash",
            },
        )

    tickvals, ticktext = _column_ticks(view)
    figure.update_xaxes(
        title_text="Base index",
        title_font={"color": theme.font_color},
        tickfont={"color": theme.font_color},
        range=list(_initial_column_x_range(view)),
        tickmode="array",
        tickvals=tickvals,
        ticktext=ticktext,
        showgrid=False,
        linecolor=theme.axis_color,
        rangeslider={
            "visible": True,
            "bgcolor": theme.rangeslider_bgcolor,
            "bordercolor": theme.rangeslider_bordercolor,
        },
    )
    figure.update_yaxes(
        title_text="Signal intensity",
        title_font={"color": theme.font_color},
        tickfont={"color": theme.font_color},
        range=[0.0, base_label_y * 1.08],
        showgrid=True,
        gridcolor=theme.grid_color,
        linecolor=theme.axis_color,
        zeroline=False,
    )
    figure.update_layout(
        hovermode="x unified",
        bargap=0.0,
        dragmode="pan",
        showlegend=False,
        plot_bgcolor=theme.plot_bgcolor,
        paper_bgcolor=theme.paper_bgcolor,
        font={"color": theme.font_color},
    )

    if quality_columns:
        max_quality = max(
            column.quality for column in quality_columns if column.quality is not None
        )
        figure.update_layout(
            yaxis2={
                "title": {
                    "text": "PHRED quality",
                    "font": {"color": theme.font_color},
                },
                "overlaying": "y",
                "side": "right",
                "range": [0.0, max(max_quality * 1.1, 1.0)],
                "showgrid": False,
                "zeroline": False,
                "tickfont": {"color": theme.font_color},
                "linecolor": theme.axis_color,
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
    theme: ChromatogramFigureTheme,
) -> None:
    retained_range = view.retained_sample_range
    if retained_range is None and view.has_any_retained_samples:
        figure.add_trace(
            go.Scattergl(
                x=view.x_values,
                y=channel.signal,
                mode="lines",
                name=f"{channel.base} trace",
                line={"color": _trace_color(channel.color, theme)},
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
            line={"color": _trace_color(channel.color, theme, muted=True)},
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
            line={"color": _trace_color(channel.color, theme)},
            hoverinfo="skip",
            hovertemplate=None,
        )
    )


def _max_column_signal(view: ChromatogramColumnView) -> float:
    signal_values = [
        max(channel.signal)
        for column in view.columns
        for channel in column.channels
        if channel.signal
    ]
    if not signal_values:
        return 1.0
    return float(max(signal_values))


def _column_trace_bases(view: ChromatogramColumnView) -> tuple[str, ...]:
    for column in view.columns:
        if column.channels:
            return tuple(channel.base for channel in column.channels)
    return ("G", "A", "T", "C")


def _column_channel_segment(
    column: ChromatogramColumn,
    *,
    base: str,
) -> ChromatogramColumnChannel | None:
    for channel in column.channels:
        if channel.base == base:
            return channel
    return None


def _column_channel_points(
    view: ChromatogramColumnView,
    *,
    base: str,
    retained_only: bool | None,
) -> tuple[list[float | None], list[float | None]]:
    x_values: list[float | None] = []
    y_values: list[float | None] = []
    for column in view.columns:
        include_column = retained_only is None or column.is_retained is retained_only
        channel = _column_channel_segment(column, base=base)
        if (
            not include_column
            or channel is None
            or not channel.x_values
            or not channel.signal
        ):
            x_values.append(None)
            y_values.append(None)
            continue
        x_values.extend(channel.x_values)
        y_values.extend(channel.signal)
        x_values.append(None)
        y_values.append(None)
    return x_values, y_values


def _add_column_channel_traces(
    figure: go.Figure,
    *,
    view: ChromatogramColumnView,
    base: str,
    theme: ChromatogramFigureTheme,
) -> None:
    any_retained = view.has_any_retained_bases
    any_trimmed = any(not column.is_retained for column in view.columns)

    x_values, y_values = _column_channel_points(
        view,
        base=base,
        retained_only=None,
    )
    if any(value is not None for value in y_values):
        figure.add_trace(
            go.Scattergl(
                x=x_values,
                y=y_values,
                mode="lines",
                name=f"{base} trace",
                line={
                    "color": _trace_color(
                        _channel_color_for_base(base),
                        theme,
                        muted=any_trimmed or not any_retained,
                    )
                },
                hoverinfo="skip",
                hovertemplate=None,
            )
        )

    if not any_retained or not any_trimmed:
        return

    retained_x_values, retained_y_values = _column_channel_points(
        view,
        base=base,
        retained_only=True,
    )
    if not any(value is not None for value in retained_y_values):
        return

    figure.add_trace(
        go.Scattergl(
            x=retained_x_values,
            y=retained_y_values,
            mode="lines",
            name=f"{base} trace (retained)",
            showlegend=False,
            line={"color": _trace_color(_channel_color_for_base(base), theme)},
            hoverinfo="skip",
            hovertemplate=None,
        )
    )


def _channel_color_for_base(base: str) -> str:
    return {
        "A": "green",
        "C": "blue",
        "G": "black",
        "T": "red",
        "N": "magenta",
    }.get(base.upper(), "magenta")


def _column_status(column: ChromatogramColumn) -> str:
    return "retained" if column.is_retained else "trimmed"


def _column_quality_color(
    column: ChromatogramColumn,
    theme: ChromatogramFigureTheme,
) -> str:
    if column.is_retained:
        return theme.retained_quality_color
    return theme.trimmed_quality_color


def _column_text_color(
    column: ChromatogramColumn,
    theme: ChromatogramFigureTheme,
) -> str:
    if column.is_retained:
        return _trace_color(column.color, theme)
    return _trace_color(column.color, theme, muted=True)


def _initial_column_x_range(view: ChromatogramColumnView) -> tuple[float, float]:
    full_left, full_right = view.x_range
    if view.base_count <= 0:
        return (0.0, 1.0)

    first_visible_column_index = _first_visible_column_index(view)
    last_visible_column_index = min(
        first_visible_column_index + _DEFAULT_INITIAL_VISIBLE_BASES - 1,
        view.base_count - 1,
    )
    padding = view.cell_width * _DEFAULT_INITIAL_BASE_PADDING_MULTIPLIER
    left = max(
        full_left,
        view.columns[first_visible_column_index].cell_left - padding,
    )
    right = min(
        full_right,
        view.columns[last_visible_column_index].cell_right + padding,
    )
    if right <= left:
        return (left, min(full_right, left + max(view.cell_width, 1.0)))
    return (left, right)


def _first_visible_column_index(view: ChromatogramColumnView) -> int:
    for index, column in enumerate(view.columns):
        if column.is_retained:
            return index
    return 0


def _column_ticks(view: ChromatogramColumnView) -> tuple[list[float], list[str]]:
    if not view.columns:
        return ([], [])

    tick_step = 10
    tick_indices = list(range(0, len(view.columns), tick_step))
    if len(view.columns) <= tick_step:
        tick_indices = list(range(len(view.columns)))
    return (
        [view.columns[index].cell_center for index in tick_indices],
        [str(view.columns[index].query_pos) for index in tick_indices],
    )


def _resolve_figure_theme(theme_type: str) -> ChromatogramFigureTheme:
    return _DARK_THEME if theme_type == "dark" else _LIGHT_THEME


def _trace_color(
    color: str,
    theme: ChromatogramFigureTheme,
    *,
    muted: bool = False,
) -> str:
    palette = theme.muted_trace_colors if muted else theme.trace_colors
    return palette.get(color, color)


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
    theme: ChromatogramFigureTheme,
) -> str:
    if _quality_segment_status(view, segment) == "retained":
        return theme.retained_quality_color
    return theme.trimmed_quality_color


def _base_call_text_color(
    view: ChromatogramView,
    base_call: ChromatogramBaseCall,
    theme: ChromatogramFigureTheme,
) -> str:
    if _is_sample_retained(view, float(base_call.position)):
        return _trace_color(base_call.color, theme)
    return _trace_color(base_call.color, theme, muted=True)
