#!/usr/bin/env python3
from __future__ import annotations
from typing import Dict, List, Optional, Sequence, Tuple
import plotly.graph_objects as go
from plotly.subplots import make_subplots

BASE_COLORS: Dict[str, str] = {"A": "green", "C": "blue", "G": "black", "T": "red"}
UNKNOWN_COLOR = "#888888"


def merge_gap_spans(columns: List[Tuple[Optional[int], Optional[int]]]):
    """Compress consecutive gap columns into spans for shading."""
    spans, i, n = [], 0, len(columns)
    while i < n:
        iA, iB = columns[i]
        if iA is None or iB is None:
            which = (
                "both" if (iA is None and iB is None) else ("A" if iA is None else "B")
            )
            j = i + 1
            while j < n:
                jA, jB = columns[j]
                same = (
                    (jA is None and jB is None)
                    if which == "both"
                    else ((jA is None) if which == "A" else (jB is None))
                )
                if not same:
                    break
                j += 1
            spans.append((i, j - 1, which))
            i = j
        else:
            i += 1
    return spans


def letters_trace(gapped: str, ncols: int, row: int):
    """Render per-column letters just below the baseline (SVG text; reliable)."""
    if not gapped:
        return None
    xs = list(range(1, ncols + 1))
    colors = [
        BASE_COLORS.get(ch.upper(), UNKNOWN_COLOR) if ch != "-" else "rgba(0,0,0,0)"
        for ch in gapped
    ]
    text = [ch if ch != "-" else "" for ch in gapped]
    y = [-0.10] * ncols  # slightly under baseline; don't clip
    return go.Scatter(
        x=xs,
        y=y,
        mode="text",
        text=text,
        textposition="middle center",
        textfont=dict(size=ninth_font_size(ncols), color=colors),
        hoverinfo="skip",
        showlegend=False,
        name=f"row{row}-letters",
        cliponaxis=False,
    )


def ninth_font_size(ncols: int) -> int:
    # Tiny helper to keep letters readable on dense alignments
    if ncols <= 600:
        return 9
    if ncols <= 1200:
        return 8
    if ncols <= 1800:
        return 7
    return 6


def _per_base_hover(
    cols: List[Tuple[Optional[int], Optional[int]]],
    wins: List[Tuple[int, int]],
    ch: Dict[str, Sequence[float]],
    gapped: str,
    is_A: bool,
) -> tuple[List[float], List[float], List[str]]:
    xs, ys, tt = [], [], []
    A, C, G, T = ch["A"], ch["C"], ch["G"], ch["T"]

    def _dom(sA, sC, sG, sT):
        vals = {
            "A": max(sA) if sA else 0.0,
            "C": max(sC) if sC else 0.0,
            "G": max(sG) if sG else 0.0,
            "T": max(sT) if sT else 0.0,
        }
        return max(vals, key=vals.get) if max(vals.values()) > 0 else None

    for k, (iA, iB) in enumerate(cols, start=1):
        idx = iA if is_A else iB
        if idx is None:
            continue
        L, R = wins[idx]
        L = max(0, int(L))
        R = max(L, int(R))
        sA, sC, sG, sT = A[L : R + 1], C[L : R + 1], G[L : R + 1], T[L : R + 1]
        local_max = max(
            (max(sA) if sA else 0.0),
            (max(sC) if sC else 0.0),
            (max(sG) if sG else 0.0),
            (max(sT) if sT else 0.0),
            0.0,
        )
        dom = _dom(sA, sC, sG, sT)
        called = gapped[k - 1] if k - 1 < len(gapped) else "-"
        extra = (
            f"<br>call: {called}"
            if (called and called != "-" and dom and called != dom)
            else ""
        )
        xs.append(float(k))
        ys.append(float(local_max))
        tt.append(f"Base {k}: {dom or called or '?'}{extra}")
    return xs, ys, tt


def plot_aligned_traces(
    *,
    columns: List[Tuple[Optional[int], Optional[int]]],
    a_series: Dict[str, Tuple[Sequence[float], Sequence[float]]],
    a_label: str,
    a_gapped: str,
    a_windows: List[Tuple[int, int]],
    a_channels: Dict[str, Sequence[float]],
    b_series: Dict[str, Tuple[Sequence[float], Sequence[float]]],
    b_label: str,
    b_gapped: str,
    b_windows: List[Tuple[int, int]],
    b_channels: Dict[str, Sequence[float]],
    letter_mode: str = "axis",  # "axis" | "overlay"
    force_svg: bool = False,
    show_rangeslider: bool = True,
    row_height: int = 260,
) -> go.Figure:
    # If overlay letters are requested, force SVG for reliable text rendering
    force_svg = bool(force_svg or (letter_mode == "overlay"))

    ncols = len(columns)
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.05,
        row_heights=[0.52, 0.48],
        subplot_titles=(a_label, b_label),
    )

    for s, e, which in merge_gap_spans(columns):
        alpha = 0.20 if which == "both" else 0.10
        fig.add_vrect(
            x0=s + 0.5,
            x1=e + 1.5,
            fillcolor=f"rgba(140,140,140,{alpha})",
            line_width=0,
            layer="below",
            row="all",
            col=1,
        )

    Trace = go.Scatter if force_svg else go.Scattergl

    def _add_row(series, row, legendgroup):
        for base in "ACGT":
            xs, ys = series[base]
            fig.add_trace(
                Trace(
                    x=xs,
                    y=ys,
                    mode="lines",
                    line=dict(width=1.2, color=BASE_COLORS[base]),
                    name=f"{legendgroup}: {base}",
                    legendgroup=legendgroup,
                    hoverinfo="skip",
                    connectgaps=False,
                ),
                row=row,
                col=1,
            )

    _add_row(a_series, 1, a_label)
    _add_row(b_series, 2, b_label)

    # Invisible per-base markers (hover text) placed at local maxima
    invisible_marker = dict(
        size=16, color="rgba(0,0,0,0)", line=dict(width=0), symbol="circle"
    )
    ax, ay, atxt = _per_base_hover(columns, a_windows, a_channels, a_gapped, is_A=True)
    bx, by, btxt = _per_base_hover(columns, b_windows, b_channels, b_gapped, is_A=False)
    fig.add_trace(
        Trace(
            x=ax,
            y=ay,
            mode="markers",
            marker=invisible_marker,
            hoverinfo="text",
            hovertext=atxt,
            showlegend=False,
            name="A-hover",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        Trace(
            x=bx,
            y=by,
            mode="markers",
            marker=invisible_marker,
            hoverinfo="text",
            hovertext=btxt,
            showlegend=False,
            name="B-hover",
        ),
        row=2,
        col=1,
    )

    # Letters over the waves (overlay mode uses SVG text)
    if letter_mode == "overlay":
        lt_a = letters_trace(a_gapped, ncols, row=1)
        lt_b = letters_trace(b_gapped, ncols, row=2)
        if lt_a:
            fig.add_trace(lt_a, row=1, col=1)
        if lt_b:
            fig.add_trace(lt_b, row=2, col=1)

    # Axes + shared spike
    fig.update_xaxes(
        showgrid=False,
        zeroline=False,
        ticks="outside",
        showspikes=True,
        spikemode="across+toaxis",
        spikesnap="cursor",
        spikethickness=2,
        spikecolor="#444",
        spikedash="solid",
        row="all",
        col=1,
    )

    # ✅ Letters on the axis (now for BOTH rows)
    if letter_mode == "axis":
        max_dense = 2400  # guardrail for huge alignments
        if ncols <= max_dense:
            tickvals = list(range(1, ncols + 1))
            top_text = [ch if ch != "-" else "" for ch in a_gapped]
            bot_text = [ch if ch != "-" else "" for ch in b_gapped]
            font = dict(size=ninth_font_size(ncols), color="#666")
            common = dict(
                showticklabels=True,
                tickmode="array",
                tickvals=tickvals,
                tickangle=0,
                ticks="",
                tickfont=font,
                automargin=True,
            )
            fig.update_xaxes(row=1, col=1, ticktext=top_text, **common)
            fig.update_xaxes(row=2, col=1, ticktext=bot_text, **common)
        else:
            # too many columns — degrade both rows to coarse numeric ticks
            fig.update_xaxes(tickmode="linear", dtick=10, row=1, col=1)
            fig.update_xaxes(tickmode="linear", dtick=10, row=2, col=1)
    else:
        # overlay letters are drawn as text, so keep light numeric ticks on BOTH rows
        fig.update_xaxes(tickmode="linear", dtick=5, row=1, col=1)
        fig.update_xaxes(tickmode="linear", dtick=5, row=2, col=1)

    # y-axis ranges per row
    def _row_ymax(series):
        m = 0.0
        for _, ys in series.values():
            if ys:
                m = max(m, max(v for v in ys if v == v))
        return m or 1.0

    a_ymax = _row_ymax(a_series)
    b_ymax = _row_ymax(b_series)
    fig.update_yaxes(
        title_text=f"{a_label} • signal (raw)",
        range=[-0.15, a_ymax * 1.02],
        row=1,
        col=1,
    )
    fig.update_yaxes(
        title_text=f"{b_label} • signal (raw)",
        range=[-0.15, b_ymax * 1.02],
        row=2,
        col=1,
    )

    # Layout, rangeslider on the bottom shared x-axis
    fig.update_layout(
        height=max(220, int(row_height)) * 2 + 80,
        margin=dict(l=60, r=28, t=54, b=42),
        hovermode="x",
        hoverdistance=-1,
        spikedistance=-1,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        xaxis_title="Alignment column (includes gaps)",
        xaxis2_title="Alignment column (includes gaps)",
    )
    fig.update_xaxes(
        rangeslider=dict(visible=bool(show_rangeslider), thickness=0.12), row=2, col=1
    )
    return fig
