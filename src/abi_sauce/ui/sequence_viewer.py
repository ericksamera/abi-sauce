#!/usr/bin/env python3
from __future__ import annotations
from typing import Optional, List, Dict, Any, Tuple
import html
import streamlit as st
import plotly.graph_objects as go  # NEW: for the optional Plotly strip
from abi_sauce.models import SequenceAsset


def _build_fasta(name: str, seq: str) -> str:
    if not seq:
        return ""
    wrapped = "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70))
    return f">{name}\n{wrapped}\n"


def _colored_sequence_html(
    seq: str,
    *,
    keep_range_0b: Optional[Tuple[int, int]] = None,  # 0-based inclusive (L, R)
    wrap: int = 70,
    crosshatch_trim: bool = True,
    color_bases: bool = False,  # optional A/C/G/T coloring for kept region
) -> str:
    """Return HTML that renders `seq` in monospace with trimmed flanks dimmed (inline styles)."""
    if not seq:
        return "<em>(empty)</em>"

    # 0b inclusive keep window
    L = R = None
    if keep_range_0b is not None:
        L, R = keep_range_0b
        L = max(0, int(L))
        R = min(len(seq) - 1, int(R))
        if R < L:
            L = R = None

    # Base colors (optional)
    base_colors = {"A": "#16a34a", "C": "#2563eb", "G": "#111827", "T": "#dc2626"}

    mono_div = (
        "font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "
        "'Liberation Mono','Courier New', monospace; font-size: 13px; "
        "line-height: 1.4; white-space: pre-wrap;"
    )
    trim_bg = (
        "background: repeating-linear-gradient(45deg, rgba(120,120,120,0.16) 0px, "
        "rgba(120,120,120,0.16) 3px, rgba(120,120,120,0.28) 3px, rgba(120,120,120,0.28) 6px);"
        if crosshatch_trim
        else ""
    )
    trim_style = (
        f"opacity:0.55; color:#6b7280; {trim_bg} padding:0 1px; border-radius:2px;"
    )

    def keep_style_for(ch: str) -> str:
        if not color_bases:
            return ""
        return f"color:{base_colors.get(ch.upper(), '#111827')};"

    spans: List[str] = []
    n = len(seq)
    for i, ch in enumerate(seq):
        safe = html.escape(ch)
        if L is not None and R is not None and (i < L or i > R):
            spans.append(f'<span style="{trim_style}">{safe}</span>')
        else:
            spans.append(f'<span style="{keep_style_for(ch)}">{safe}</span>')

    lines = ["".join(spans[i : i + wrap]) for i in range(0, n, wrap)]
    return f"<div style='{mono_div}'>{'<br>'.join(lines)}</div>"


def _sequence_plotly_strip(
    seq: str,
    *,
    keep_range_0b: Optional[Tuple[int, int]] = None,
    color_bases: bool = True,
):
    """Small Plotly figure with colored per-base letters + shaded trimmed flanks."""
    n = len(seq)
    if n == 0:
        return None
    x = list(range(1, n + 1))
    text = list(seq)
    # per-point colors
    base_colors = {"A": "green", "C": "blue", "G": "black", "T": "red"}

    def color_for(i: int, ch: str) -> str:
        if keep_range_0b is not None:
            L, R = keep_range_0b
            if i < L or i > R:
                return "#9aa0a6"
        return base_colors.get(ch.upper(), "#111827") if color_bases else "#111827"

    colors = [color_for(i, ch) for i, ch in enumerate(seq)]

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=[0] * n,
            mode="text",
            text=text,
            textposition="middle center",
            textfont=dict(size=12, color=colors),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    # Shade trimmed flanks
    if keep_range_0b is not None:
        L, R = keep_range_0b
        if L > 0:
            fig.add_vrect(
                x0=0.5,
                x1=L + 0.5,
                fillcolor="rgba(150,150,150,0.18)",
                line_width=0,
                layer="below",
            )
        if R < n - 1:
            fig.add_vrect(
                x0=R + 1.5,
                x1=n + 0.5,
                fillcolor="rgba(150,150,150,0.18)",
                line_width=0,
                layer="below",
            )

    fig.update_xaxes(visible=False, range=[0.5, n + 0.5])
    fig.update_yaxes(visible=False, range=[-1, 1])
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        height=56,
    )
    return fig


def render_sequence_block(
    *,
    name: str,
    sequence: str,
    features: Optional[List[Dict[str, Any]]] = None,
    fasta_text: Optional[str] = None,
    download_filename: Optional[str] = None,
    title: str = "Sequence",
    # highlight outside the kept region (0-based inclusive keep)
    keep_range_0b: Optional[Tuple[int, int]] = None,
    crosshatch_trim: bool = True,
    colorize: bool = False,  # when True, use the HTML renderer
    color_bases: bool = False,  # optionally color A/C/G/T in kept region
    plotly_strip: bool = False,  # when True, also show a tiny Plotly strip
) -> None:
    """Generic sequence viewer used by both FASTA assets and AB1 traces."""
    st.subheader(title)
    seq = sequence or ""

    if colorize and (keep_range_0b is not None):
        html_block = _colored_sequence_html(
            seq,
            keep_range_0b=keep_range_0b,
            crosshatch_trim=crosshatch_trim,
            color_bases=color_bases,
        )
        st.markdown(html_block, unsafe_allow_html=True)
        st.caption(
            "Trimmed bases (outside keep window) are dimmed"
            + (" with crosshatch." if crosshatch_trim else ".")
        )
    else:
        st.code(seq[:5000] + ("…" if len(seq) > 5000 else ""), wrap_lines=True)

    if plotly_strip:
        fig = _sequence_plotly_strip(seq, keep_range_0b=keep_range_0b, color_bases=True)
        if fig:
            st.plotly_chart(fig, use_container_width=True)

    with st.expander("Features", expanded=False):
        if features:
            st.json(features)
        else:
            st.caption("(No features)")

    fa = fasta_text or _build_fasta(name, seq)
    if fa:
        st.download_button(
            "Download FASTA", data=fa, file_name=download_filename or f"{name}.fasta"
        )


def render_sequence_asset(asset: SequenceAsset) -> None:
    render_sequence_block(
        name=asset.name,
        sequence=asset.sequence or "",
        features=getattr(asset, "features", None),
        fasta_text=asset.to_fasta(),
        download_filename=f"{asset.name}.fasta",
        title="Sequence",
        colorize=False,
    )
