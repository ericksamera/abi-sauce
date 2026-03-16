from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import plotly.graph_objects as go

from abi_sauce.reference_alignment_multi_trace import (
    ReferenceMultiAlignmentTraceCell,
    ReferenceMultiAlignmentTraceChannelSegment,
    ReferenceMultiAlignmentTraceRow,
    ReferenceMultiAlignmentTraceView,
)


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentTraceFigureTheme:
    trace_colors: dict[str, str]
    row_background_color: str
    row_border_color: str
    reference_text_color: str
    consensus_text_color: str
    gap_fill_color: str
    insertion_fill_color: str
    selected_column_color: str
    plot_bgcolor: str
    paper_bgcolor: str
    font_color: str
    axis_color: str
    grid_color: str
    rangeslider_bgcolor: str
    rangeslider_bordercolor: str
    resolution_fill_colors: dict[str, str]


_LIGHT_THEME: Final[ReferenceMultiAlignmentTraceFigureTheme] = (
    ReferenceMultiAlignmentTraceFigureTheme(
        trace_colors={
            "green": "green",
            "blue": "blue",
            "black": "black",
            "red": "red",
            "magenta": "magenta",
        },
        row_background_color="rgba(255, 255, 255, 0.0)",
        row_border_color="rgba(0, 0, 0, 0.12)",
        reference_text_color="#111111",
        consensus_text_color="#111111",
        gap_fill_color="rgba(148, 163, 184, 0.10)",
        insertion_fill_color="rgba(59, 130, 246, 0.10)",
        selected_column_color="rgba(234, 179, 8, 0.22)",
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        font_color="#111111",
        axis_color="#444444",
        grid_color="rgba(0, 0, 0, 0.08)",
        rangeslider_bgcolor="#F4F6F8",
        rangeslider_bordercolor="rgba(0, 0, 0, 0.10)",
        resolution_fill_colors={
            "concordant": "rgba(34, 197, 94, 0.08)",
            "single_read": "rgba(59, 130, 246, 0.10)",
            "quality_resolved": "rgba(245, 158, 11, 0.14)",
            "majority_resolved": "rgba(168, 85, 247, 0.14)",
            "ambiguous": "rgba(239, 68, 68, 0.14)",
            "deleted": "rgba(120, 113, 108, 0.16)",
        },
    )
)

_DARK_THEME: Final[ReferenceMultiAlignmentTraceFigureTheme] = (
    ReferenceMultiAlignmentTraceFigureTheme(
        trace_colors={
            "green": "#22C55E",
            "blue": "#3B82F6",
            "black": "#E5E7EB",
            "red": "#EF4444",
            "magenta": "#F472B6",
        },
        row_background_color="rgba(255, 255, 255, 0.0)",
        row_border_color="rgba(203, 213, 225, 0.20)",
        reference_text_color="#FAFAFA",
        consensus_text_color="#FAFAFA",
        gap_fill_color="rgba(148, 163, 184, 0.18)",
        insertion_fill_color="rgba(59, 130, 246, 0.16)",
        selected_column_color="rgba(234, 179, 8, 0.24)",
        plot_bgcolor="#0E1117",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        font_color="#FAFAFA",
        axis_color="#CBD5E1",
        grid_color="rgba(203, 213, 225, 0.16)",
        rangeslider_bgcolor="#0B0F14",
        rangeslider_bordercolor="rgba(203, 213, 225, 0.16)",
        resolution_fill_colors={
            "concordant": "rgba(34, 197, 94, 0.12)",
            "single_read": "rgba(59, 130, 246, 0.14)",
            "quality_resolved": "rgba(245, 158, 11, 0.18)",
            "majority_resolved": "rgba(168, 85, 247, 0.18)",
            "ambiguous": "rgba(239, 68, 68, 0.18)",
            "deleted": "rgba(120, 113, 108, 0.20)",
        },
    )
)

_REFERENCE_BAND_HEIGHT: Final[float] = 0.9
_CONSENSUS_BAND_HEIGHT: Final[float] = 0.9
_BAND_GAP: Final[float] = 0.3
_DEFAULT_INITIAL_VISIBLE_COLUMNS: Final[int] = 50


def build_reference_multi_alignment_trace_figure(
    view: ReferenceMultiAlignmentTraceView,
    *,
    theme_type: str = "light",
    selected_column_index: int | None = None,
) -> go.Figure:
    """Build a stacked shared-reference electropherogram figure."""
    theme = _resolve_figure_theme(theme_type)
    figure = go.Figure()

    consensus_bottom = 0.0
    consensus_top = consensus_bottom + _CONSENSUS_BAND_HEIGHT
    consensus_mid = (consensus_bottom + consensus_top) / 2.0
    rows_y_offset = consensus_top + _BAND_GAP
    reference_bottom = rows_y_offset + view.total_height + _BAND_GAP
    reference_top = reference_bottom + _REFERENCE_BAND_HEIGHT
    reference_mid = (reference_bottom + reference_top) / 2.0
    total_height = reference_top + 0.2

    if (
        selected_column_index is not None
        and 1 <= selected_column_index <= view.alignment_length
    ):
        selected_left, selected_right = _column_bounds(
            selected_column_index,
            cell_width=view.cell_width,
        )
        figure.add_shape(
            type="rect",
            xref="x",
            yref="y",
            x0=selected_left,
            x1=selected_right,
            y0=0.0,
            y1=total_height,
            line={"width": 0},
            fillcolor=theme.selected_column_color,
            layer="below",
        )

    _add_consensus_band(
        figure,
        view=view,
        theme=theme,
        y_bottom=consensus_bottom,
        y_top=consensus_top,
        y_text=consensus_mid,
    )
    _add_reference_band(
        figure,
        view=view,
        theme=theme,
        y_bottom=reference_bottom,
        y_top=reference_top,
        y_text=reference_mid,
    )

    for row in view.rows:
        _add_row_background(
            figure,
            row=row,
            view=view,
            theme=theme,
            y_offset=rows_y_offset,
        )
        _add_row_channel_segments(
            figure,
            row=row,
            theme=theme,
            y_offset=rows_y_offset,
        )
        _add_row_base_labels(
            figure,
            row=row,
            theme=theme,
            y_offset=rows_y_offset,
        )
        _add_row_hover_markers(
            figure,
            row=row,
            y_offset=rows_y_offset,
        )

    figure.update_xaxes(
        title_text="Alignment column",
        title_font={"color": theme.font_color},
        tickfont={"color": theme.font_color},
        range=list(_initial_x_range(view)),
        showgrid=True,
        gridcolor=theme.grid_color,
        linecolor=theme.axis_color,
        zeroline=False,
        rangeslider={
            "visible": True,
            "bgcolor": theme.rangeslider_bgcolor,
            "bordercolor": theme.rangeslider_bordercolor,
        },
    )
    figure.update_yaxes(
        range=[0.0, total_height],
        tickmode="array",
        tickvals=[_row_midpoint(row, y_offset=rows_y_offset) for row in view.rows],
        ticktext=[row.display_name for row in view.rows],
        tickfont={"color": theme.font_color},
        showgrid=False,
        linecolor=theme.axis_color,
        zeroline=False,
    )
    figure.update_layout(
        hovermode="closest",
        dragmode="pan",
        showlegend=False,
        barmode="overlay",
        bargap=0.0,
        plot_bgcolor=theme.plot_bgcolor,
        paper_bgcolor=theme.paper_bgcolor,
        font={"color": theme.font_color},
        margin={"l": 24, "r": 24, "t": 24, "b": 56},
    )
    return figure


def _add_reference_band(
    figure: go.Figure,
    *,
    view: ReferenceMultiAlignmentTraceView,
    theme: ReferenceMultiAlignmentTraceFigureTheme,
    y_bottom: float,
    y_top: float,
    y_text: float,
) -> None:
    _add_background_bars(
        figure,
        x_values=[
            _column_center(column.column_index, cell_width=view.cell_width)
            for column in view.columns
        ],
        y_base=y_bottom,
        height=y_top - y_bottom,
        width=view.cell_width,
        colors=[_column_fill_color(column, theme) for column in view.columns],
        name="Reference background",
    )
    figure.add_trace(
        go.Scattergl(
            x=[
                _column_center(column.column_index, cell_width=view.cell_width)
                for column in view.columns
            ],
            y=[y_text] * len(view.columns),
            text=[column.ref_base for column in view.columns],
            mode="text",
            name="Reference",
            textfont={"color": theme.reference_text_color},
            hoverinfo="skip",
            hovertemplate=None,
        )
    )
    figure.add_trace(
        go.Scattergl(
            x=[
                _column_center(column.column_index, cell_width=view.cell_width)
                for column in view.columns
            ],
            y=[y_text] * len(view.columns),
            mode="markers",
            name="Reference hover",
            showlegend=False,
            customdata=[
                (
                    f"row=reference"
                    f"<br>column={column.column_index}"
                    f"<br>anchor={column.anchor_kind}"
                    f"<br>ref_pos={column.ref_pos}"
                    f"<br>ref_base={column.ref_base}"
                    f"<br>consensus={column.consensus_base}"
                    f"<br>resolution={column.resolution}"
                    f"<br>support={_support_summary(column.support_counts)}"
                )
                for column in view.columns
            ],
            marker={"size": 14, "opacity": 0},
            hovertemplate="%{customdata}<extra></extra>",
        )
    )


def _add_consensus_band(
    figure: go.Figure,
    *,
    view: ReferenceMultiAlignmentTraceView,
    theme: ReferenceMultiAlignmentTraceFigureTheme,
    y_bottom: float,
    y_top: float,
    y_text: float,
) -> None:
    _add_background_bars(
        figure,
        x_values=[
            _column_center(column.column_index, cell_width=view.cell_width)
            for column in view.columns
        ],
        y_base=y_bottom,
        height=y_top - y_bottom,
        width=view.cell_width,
        colors=[_column_fill_color(column, theme) for column in view.columns],
        name="Consensus background",
    )
    figure.add_trace(
        go.Scattergl(
            x=[
                _column_center(column.column_index, cell_width=view.cell_width)
                for column in view.columns
            ],
            y=[y_text] * len(view.columns),
            text=[column.consensus_base for column in view.columns],
            mode="text",
            name="Consensus",
            textfont={"color": theme.consensus_text_color},
            hoverinfo="skip",
            hovertemplate=None,
        )
    )
    figure.add_trace(
        go.Scattergl(
            x=[
                _column_center(column.column_index, cell_width=view.cell_width)
                for column in view.columns
            ],
            y=[y_text] * len(view.columns),
            mode="markers",
            name="Consensus hover",
            showlegend=False,
            customdata=[
                (
                    f"row=consensus"
                    f"<br>column={column.column_index}"
                    f"<br>consensus={column.consensus_base}"
                    f"<br>resolution={column.resolution}"
                    f"<br>support={_support_summary(column.support_counts)}"
                    f"<br>gap_members={column.gap_member_count}"
                )
                for column in view.columns
            ],
            marker={"size": 14, "opacity": 0},
            hovertemplate="%{customdata}<extra></extra>",
        )
    )


def _add_row_background(
    figure: go.Figure,
    *,
    row: ReferenceMultiAlignmentTraceRow,
    view: ReferenceMultiAlignmentTraceView,
    theme: ReferenceMultiAlignmentTraceFigureTheme,
    y_offset: float,
) -> None:
    full_left, full_right = view.x_range
    figure.add_shape(
        type="rect",
        xref="x",
        yref="y",
        x0=full_left,
        x1=full_right,
        y0=row.y_bottom + y_offset,
        y1=row.y_top + y_offset,
        line={"width": 1, "color": theme.row_border_color},
        fillcolor=theme.row_background_color,
        layer="below",
    )
    _add_background_bars(
        figure,
        x_values=[cell.cell_center for cell in row.cells],
        y_base=row.y_bottom + y_offset,
        height=row.y_top - row.y_bottom,
        width=view.cell_width,
        colors=[_cell_fill_color(cell, theme) for cell in row.cells],
        name=f"{row.display_name} background",
    )


def _add_background_bars(
    figure: go.Figure,
    *,
    x_values: list[float],
    y_base: float,
    height: float,
    width: float,
    colors: list[str],
    name: str,
) -> None:
    if not x_values:
        return
    figure.add_trace(
        go.Bar(
            x=x_values,
            y=[height] * len(x_values),
            base=[y_base] * len(x_values),
            width=[width] * len(x_values),
            name=name,
            showlegend=False,
            marker={"color": colors, "line": {"width": 0}},
            hoverinfo="skip",
            hovertemplate=None,
        )
    )


def _add_row_channel_segments(
    figure: go.Figure,
    *,
    row: ReferenceMultiAlignmentTraceRow,
    theme: ReferenceMultiAlignmentTraceFigureTheme,
    y_offset: float,
) -> None:
    if not row.has_trace_signal:
        return

    baseline = row.y_bottom + y_offset + 0.15
    trace_height = max((row.y_top - row.y_bottom) - 1.0, 0.5)
    for base in _row_channel_order(row):
        x_values: list[float | None] = []
        y_values: list[float | None] = []
        for cell in row.cells:
            segment = _cell_channel_segment(cell, base=base)
            if segment is None or not segment.x_values or not segment.normalized_signal:
                x_values.append(None)
                y_values.append(None)
                continue
            x_values.extend(segment.x_values)
            y_values.extend(
                baseline + (trace_height * signal)
                for signal in segment.normalized_signal
            )
            x_values.append(None)
            y_values.append(None)
        if not x_values:
            continue
        figure.add_trace(
            go.Scattergl(
                x=x_values,
                y=y_values,
                mode="lines",
                name=f"{row.display_name}:{base}",
                line={
                    "color": _trace_color_for_segment_base(
                        base=base,
                        fallback_color="magenta",
                        theme=theme,
                    ),
                    "width": 1.2,
                },
                hoverinfo="skip",
                hovertemplate=None,
            )
        )


def _add_row_base_labels(
    figure: go.Figure,
    *,
    row: ReferenceMultiAlignmentTraceRow,
    theme: ReferenceMultiAlignmentTraceFigureTheme,
    y_offset: float,
) -> None:
    text_y = row.y_top + y_offset - 0.35
    figure.add_trace(
        go.Scattergl(
            x=[cell.cell_center for cell in row.cells],
            y=[text_y] * len(row.cells),
            text=[cell.query_base for cell in row.cells],
            mode="text",
            name=f"{row.display_name} bases",
            textfont={
                "color": [
                    (
                        _trace_color_for_segment_base(
                            base=cell.query_base,
                            fallback_color=theme.font_color,
                            theme=theme,
                        )
                        if cell.query_base != "-"
                        else theme.font_color
                    )
                    for cell in row.cells
                ]
            },
            hoverinfo="skip",
            hovertemplate=None,
        )
    )


def _add_row_hover_markers(
    figure: go.Figure,
    *,
    row: ReferenceMultiAlignmentTraceRow,
    y_offset: float,
) -> None:
    figure.add_trace(
        go.Scattergl(
            x=[cell.cell_center for cell in row.cells],
            y=[_row_midpoint(row, y_offset=y_offset)] * len(row.cells),
            mode="markers",
            name=f"{row.display_name} hover",
            showlegend=False,
            customdata=[
                [
                    cell.column_index,
                    cell.anchor_kind,
                    cell.ref_base,
                    cell.query_base,
                    cell.consensus_base,
                    cell.resolution,
                    cell.ref_pos,
                    cell.query_pos,
                    cell.quality,
                    cell.trace_x,
                    cell.raw_left,
                    cell.raw_right,
                    row.display_name,
                    row.strand,
                ]
                for cell in row.cells
            ],
            marker={"size": 16, "opacity": 0},
            hovertemplate=(
                "row=%{customdata[12]}"
                "<br>strand=%{customdata[13]}"
                "<br>column=%{customdata[0]}"
                "<br>anchor=%{customdata[1]}"
                "<br>reference_base=%{customdata[2]}"
                "<br>query_base=%{customdata[3]}"
                "<br>consensus=%{customdata[4]}"
                "<br>resolution=%{customdata[5]}"
                "<br>ref_pos=%{customdata[6]}"
                "<br>query_pos=%{customdata[7]}"
                "<br>quality=%{customdata[8]}"
                "<br>trace_x=%{customdata[9]}"
                "<br>raw_left=%{customdata[10]}"
                "<br>raw_right=%{customdata[11]}"
                "<extra></extra>"
            ),
        )
    )


def _column_fill_color(column, theme: ReferenceMultiAlignmentTraceFigureTheme) -> str:
    if column.anchor_kind == "insertion":
        return theme.insertion_fill_color
    return theme.resolution_fill_colors.get(column.resolution, theme.gap_fill_color)


def _cell_fill_color(
    cell: ReferenceMultiAlignmentTraceCell,
    theme: ReferenceMultiAlignmentTraceFigureTheme,
) -> str:
    if cell.is_gap:
        return theme.gap_fill_color
    return theme.resolution_fill_colors.get(cell.resolution, theme.gap_fill_color)


def _support_summary(support_counts: tuple[tuple[str, int], ...]) -> str:
    if not support_counts:
        return "NA"
    return ", ".join(f"{base}:{count}" for base, count in support_counts)


def _row_channel_order(row: ReferenceMultiAlignmentTraceRow) -> tuple[str, ...]:
    for cell in row.cells:
        if cell.channels:
            return tuple(segment.base for segment in cell.channels)
    return ("G", "A", "T", "C")


def _cell_channel_segment(
    cell: ReferenceMultiAlignmentTraceCell,
    *,
    base: str,
) -> ReferenceMultiAlignmentTraceChannelSegment | None:
    for segment in cell.channels:
        if segment.base == base:
            return segment
    return None


def _trace_color_for_segment_base(
    *,
    base: str,
    fallback_color: str,
    theme: ReferenceMultiAlignmentTraceFigureTheme,
) -> str:
    base_to_color = {
        "A": "green",
        "C": "blue",
        "G": "black",
        "T": "red",
        "N": "magenta",
    }
    color_key = base_to_color.get(base.upper())
    if color_key is None:
        return fallback_color
    return theme.trace_colors.get(color_key, fallback_color)


def _resolve_figure_theme(theme_type: str) -> ReferenceMultiAlignmentTraceFigureTheme:
    return _DARK_THEME if theme_type == "dark" else _LIGHT_THEME


def _initial_x_range(view: ReferenceMultiAlignmentTraceView) -> tuple[float, float]:
    full_left, full_right = view.x_range
    visible_width = min(
        float(_DEFAULT_INITIAL_VISIBLE_COLUMNS) * view.cell_width,
        full_right - full_left,
    )
    if visible_width <= 0:
        return (full_left, full_right)
    return (full_left, full_left + visible_width)


def _row_midpoint(row: ReferenceMultiAlignmentTraceRow, *, y_offset: float) -> float:
    return ((row.y_bottom + y_offset) + (row.y_top + y_offset)) / 2.0


def _column_bounds(
    column_index: int,
    *,
    cell_width: float,
) -> tuple[float, float]:
    left = float(column_index - 1) * cell_width
    return (left, left + cell_width)


def _column_center(
    column_index: int,
    *,
    cell_width: float,
) -> float:
    left, right = _column_bounds(column_index, cell_width=cell_width)
    return (left + right) / 2.0


__all__ = ["build_reference_multi_alignment_trace_figure"]
