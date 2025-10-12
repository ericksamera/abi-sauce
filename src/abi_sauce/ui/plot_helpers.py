# src/abi_sauce/ui/plot_helpers.py
from typing import Dict, List, Optional, Tuple, Literal
import math
import plotly.graph_objects as go


# --- Defensive helpers to avoid Plotly/pandas import-order issues in Streamlit ---

def _to_simple_list(seq):
    """
    Convert many array-like types to a plain Python list.
    Handles: None, list/tuple, numpy arrays, pandas Series/Index (even if pandas
    is partially initialized), objects with .tolist(), and generic iterables.
    """
    if seq is None:
        return None
    if isinstance(seq, (list, tuple)):
        return list(seq)
    try:
        tolist = getattr(seq, "tolist", None)
        if callable(tolist):
            out = tolist()
            if isinstance(out, (list, tuple)):
                return list(out)
            try:
                return list(out)
            except Exception:
                return [out]
    except Exception:
        pass
    try:
        return list(seq)
    except Exception:
        return [seq]


def _safe_add_scatter(fig, xs, ys, **scatter_kwargs):
    """Normalize xs/ys to lists and add a Scatter trace."""
    xs = _to_simple_list(xs)
    ys = _to_simple_list(ys)
    trace = go.Scatter(x=xs, y=ys, **scatter_kwargs)
    fig.add_trace(trace)
    return fig


def _safe_add_bar(fig, xs, ys, **bar_kwargs):
    """Normalize xs/ys to lists and add a Bar trace."""
    xs = _to_simple_list(xs)
    ys = _to_simple_list(ys)
    trace = go.Bar(x=xs, y=ys, **bar_kwargs)
    fig.add_trace(trace)
    return fig


def _clamp(v: int, lo: int, hi: int) -> int:
    return max(lo, min(hi, v))


def _compute_base_windows(ploc: List[int], total_samples: int) -> List[Tuple[int, int]]:
    """
    Turn peak locations into per-base sample windows using midpoints between peaks.
    Returns inclusive (start_idx, end_idx) pairs per base.
    """
    n = len(ploc)
    if n == 0:
        return []
    windows: List[Tuple[int, int]] = []
    mids = [(ploc[i] + ploc[i + 1]) / 2.0 for i in range(n - 1)]
    for i in range(n):
        start = 0 if i == 0 else (math.floor(mids[i - 1]) + 1)
        end   = (total_samples - 1) if i == n - 1 else math.floor(mids[i])
        start = _clamp(start, 0, total_samples - 1)
        end   = _clamp(end,   0, total_samples - 1)
        if end < start:
            end = start
        windows.append((start, end))
    return windows


def _resample_channel_to_bases(channel: List[float],
                               windows: List[Tuple[int, int]],
                               method: Literal['max', 'mean', 'median'] = 'max') -> List[float]:
    """Aggregate `channel` samples into one value per base-window."""
    import statistics
    res: List[float] = []
    n = len(channel)
    for s, e in windows:
        s = _clamp(s, 0, n - 1)
        e = _clamp(e, 0, n - 1)
        if e < s:
            e = s
        block = channel[s:e + 1]
        if not block:
            res.append(0.0)
            continue
        if method == 'max':
            res.append(max(block))
        elif method == 'mean':
            res.append(sum(block) / len(block))
        elif method == 'median':
            res.append(statistics.median(block))
        else:
            res.append(max(block))
    return res


# --- Peak labeling / positioning helpers ---

def _dominant_base_at_peak(p: int, channels: Dict[str, List[float]], window: int = 2) -> Optional[str]:
    """Return the base letter with strongest signal near sample index p (±window)."""
    best_base, best_val = None, float("-inf")
    for b in "ACGT":
        arr = channels.get(b) or []
        if not arr:
            continue
        s = _clamp(p - window, 0, len(arr) - 1)
        e = _clamp(p + window, 0, len(arr) - 1)
        if e < s:
            continue
        val = max(arr[s:e + 1])
        if val > best_val:
            best_val = val
            best_base = b
    return best_base


def _peak_heights_samples(ploc: List[int], channels: Dict[str, List[float]], window: int = 2) -> List[float]:
    """For each peak location, return the max across A/C/G/T within ±window samples."""
    heights: List[float] = []
    for p in ploc:
        best = 0.0
        for arr in channels.values():
            if not arr:
                continue
            s = _clamp(p - window, 0, len(arr) - 1)
            e = _clamp(p + window, 0, len(arr) - 1)
            if e >= s:
                best = max(best, max(arr[s:e + 1]))
        heights.append(best)
    return heights


def _peak_heights_bases(windows: List[Tuple[int, int]], channels: Dict[str, List[float]]) -> List[float]:
    """Per-base peak heights (for bases mode): max across A/C/G/T within each base window."""
    heights: List[float] = []
    for s, e in windows:
        best = 0.0
        for arr in channels.values():
            if not arr:
                continue
            n = len(arr)
            ss = _clamp(s, 0, n - 1)
            ee = _clamp(e, 0, n - 1)
            if ee >= ss:
                best = max(best, max(arr[ss:ee + 1]))
        heights.append(best)
    return heights


def _peak_hover_texts(ploc: List[int],
                      seq_len: int,
                      seq: Optional[str],
                      channels: Dict[str, List[float]],
                      quals: Optional[List[float]]) -> List[str]:
    """Hover strings for each base: 'Base N: {dominant}' (+ called base if different) + sample + Q."""
    texts: List[str] = []
    for i in range(seq_len):
        p = ploc[i] if i < len(ploc) else (ploc[-1] if ploc else 0)
        peak_base = _dominant_base_at_peak(p, channels, window=2)
        called = None
        if seq and i < len(seq):
            c = seq[i].upper()
            if c in "ACGTN":
                called = c
        main = peak_base or (called if called and called != "N" else "?")
        extra = f"<br>call: {called}" if (called and called != "N" and peak_base and called != peak_base) else ""
        qtxt = f"<br>Q{int(quals[i])}" if (quals and i < len(quals)) else ""
        texts.append(f"Base {i+1}: {main}{extra}<br>sample {p}{qtxt}")
    return texts


# --- Main API ---

def build_trace_fig(
    channels: Dict[str, List[float]],
    ploc: List[int],
    seq_len: Optional[int] = None,
    quals: Optional[List[float]] = None,
    seq: Optional[str] = None,
    mode: Literal['samples', 'bases'] = 'samples',
    resample_method: Literal['max', 'mean', 'median'] = 'max',
    show_trim: bool = True,
    trim_range: Optional[Tuple[int, int]] = None,  # (left_base_idx, right_base_idx), 1-based inclusive
    uirevision: str = "chrom_trim_base_ticks_v1",
    height: int = 340,
    show_grid: bool = True,
    peak_only_hover: bool = True,  # if True, only peaks show hover text
    show_rangeslider: bool = True, # NEW: show the x-axis range slider
) -> go.Figure:
    """
    Build a Plotly Figure for chromatogram-like signal plotting.
    - Invisible (but hoverable) peak markers positioned at peak tops.
    - PHRED bars fill the entire inter-peak region in 'samples' mode (midpoint-to-midpoint).
      In 'bases' mode, bars fill each base cell.
    - X-axis range slider enabled (toggle via show_rangeslider).
    """

    # sizes and defaults
    total_samples = next((len(arr) for arr in channels.values()), 0)
    if seq_len is None:
        seq_len = len(ploc)

    # align ploc length
    if len(ploc) < seq_len:
        ploc = (ploc[:] + [total_samples - 1] * (seq_len - len(ploc)))[:seq_len]
    elif len(ploc) > seq_len:
        ploc = ploc[:seq_len]

    windows = _compute_base_windows(ploc, total_samples)

    # y scaling
    max_signal = max((max(arr) for arr in channels.values() if arr), default=0.0)
    y_max = max_signal * 1.08 if max_signal > 0 else 1.0

    fig = go.Figure()

    # hover settings
    channel_hoverinfo = "skip" if peak_only_hover else "x+y+name"
    bar_hoverinfo = "skip" if peak_only_hover else "x+y"

    # precompute hover texts
    peak_htexts = _peak_hover_texts(ploc, seq_len, seq, channels, quals)

    # transparent marker style (hoverable but invisible)
    invisible_marker = dict(
        size=16,                      # large hitbox
        color='rgba(0,0,0,0)',        # fully transparent
        line=dict(width=0),
        symbol='circle',
    )

    if mode == 'samples':
        xs = list(range(total_samples))
        for base_letter, arr in channels.items():
            _safe_add_scatter(
                fig, xs, arr,
                mode='lines',
                name=base_letter,
                hoverinfo=channel_hoverinfo,
                hovertemplate=None if peak_only_hover else "%{y}<extra>%{fullData.name}</extra>",
                line=dict(width=1.5),
                opacity=0.95,
            )

        # invisible markers at actual peak tops
        peak_ys = _peak_heights_samples(ploc, channels, window=2)
        _safe_add_scatter(
            fig, ploc, peak_ys,
            mode='markers',
            marker=invisible_marker,
            name='peaks',
            hoverinfo='text',
            hovertext=peak_htexts,
            showlegend=False
        )

        xaxis_title = "Sample index (base ticks shown below)"
        x_range = [0, max(total_samples - 1, 0)]
        tick_step = 1 if seq_len <= 50 else max(1, seq_len // 50)
        tickvals = [ploc[i] for i in range(0, seq_len, tick_step)]
        ticktext = [str(i + 1) for i in range(0, seq_len, tick_step)]

    else:  # mode == 'bases'
        xs = list(range(1, seq_len + 1))

        if not windows:
            for base_letter, _ in channels.items():
                _safe_add_scatter(
                    fig, xs, [0] * seq_len,
                    mode='lines',
                    name=base_letter,
                    hoverinfo=channel_hoverinfo,
                )
            peak_ys_bases = [0.0] * seq_len
        else:
            resampled_by_base: Dict[str, List[float]] = {}
            for base_letter, arr in channels.items():
                resampled_by_base[base_letter] = _resample_channel_to_bases(arr, windows, method=resample_method)
                _safe_add_scatter(
                    fig, xs, resampled_by_base[base_letter],
                    mode='lines',
                    name=base_letter,
                    hoverinfo=channel_hoverinfo,
                    hovertemplate=None if peak_only_hover else "%{y}<extra>%{fullData.name}</extra>",
                    line=dict(width=1.5)
                )
            peak_ys_bases = [
                max(resampled_by_base.get('A', [0]*seq_len)[i],
                    resampled_by_base.get('C', [0]*seq_len)[i],
                    resampled_by_base.get('G', [0]*seq_len)[i],
                    resampled_by_base.get('T', [0]*seq_len)[i])
                for i in range(seq_len)
            ]

        _safe_add_scatter(
            fig, xs[:len(peak_htexts)], peak_ys_bases[:len(peak_htexts)],
            mode='markers',
            marker=invisible_marker,  # invisible but hoverable
            name='peaks',
            hoverinfo='text',
            hovertext=peak_htexts,
            showlegend=False
        )

        xaxis_title = "Base number"
        x_range = [0.5, seq_len + 0.5]
        if seq_len <= 50:
            tickvals = xs
            ticktext = [str(i) for i in xs]
        else:
            step = max(1, seq_len // 50)
            tickvals = [i for i in xs if (i - 1) % step == 0]
            ticktext = [str(i) for i in tickvals]

    # quality overlay (PHRED) — with full inter-peak width in samples mode
    if quals:
        qlen = len(quals)
        if qlen < seq_len:
            quals = quals + [0] * (seq_len - qlen)
        elif qlen > seq_len:
            quals = quals[:seq_len]

        if mode == 'samples':
            if windows:
                q_centers = [ (s + e) / 2.0 for (s, e) in windows ][:len(quals)]
                q_widths  = [ (e - s + 1)    for (s, e) in windows ][:len(quals)]
                q_x = q_centers
                q_y = quals[:len(q_x)]
                _safe_add_bar(
                    fig, q_x, q_y,
                    name='quality',
                    width=q_widths,       # variable widths (midpoint-to-midpoint)
                    marker=dict(opacity=0.35),
                    yaxis='y2',
                    hoverinfo=bar_hoverinfo,
                    showlegend=True
                )
            else:
                q_x = ploc[:len(quals)]
                q_y = quals[:len(q_x)]
                _safe_add_bar(
                    fig, q_x, q_y,
                    name='quality',
                    marker=dict(opacity=0.35),
                    yaxis='y2',
                    hoverinfo=bar_hoverinfo,
                    showlegend=True
                )
        else:
            q_x = list(range(1, seq_len + 1))
            q_y = quals
            _safe_add_bar(
                fig, q_x, q_y,
                name='quality',
                width=0.98,            # fill the base cell
                marker=dict(opacity=0.35),
                yaxis='y2',
                hoverinfo=bar_hoverinfo,
                showlegend=True
            )

    # layout and axes
    fig.update_layout(
        height=height,
        uirevision=uirevision,
        margin=dict(l=40, r=40, t=40, b=40),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='left', x=0),
        hovermode="closest",  # only point under cursor gets hover
        modebar=dict(remove=['lasso2d', 'select2d'])
    )

    # X-axis + range slider
    fig.update_xaxes(
        title_text=xaxis_title,
        showgrid=False,
        range=x_range,
        tickmode='array',
        tickvals=tickvals if 'tickvals' in locals() else None,
        ticktext=ticktext if 'ticktext' in locals() else None,
        tickangle=0,
        zeroline=False,
        rangeslider=dict(visible=show_rangeslider, thickness=0.12)  # <-- range slider
    )

    # primary yaxis (signals)
    fig.update_yaxes(
        title_text='Signal intensity',
        range=[0, y_max],
        showgrid=show_grid,
        zeroline=False,
    )

    # secondary y-axis for quality scores (if present)
    if quals:
        qmax = max(quals) if quals else 1
        fig.update_layout(
            yaxis2=dict(
                title='PHRED',
                overlaying='y',
                side='right',
                range=[0, max(1, qmax * 1.05)],
                showgrid=False,
            )
        )

    # optional trim highlight
    if show_trim and trim_range:
        left_base, right_base = map(int, trim_range)
        if mode == 'samples':
            if windows and 1 <= left_base <= len(windows) and 1 <= right_base <= len(windows):
                x0 = windows[left_base - 1][0]
                x1 = windows[right_base - 1][1]
                fig.add_shape(
                    type="rect",
                    x0=x0, x1=x1, y0=0, y1=y_max,
                    xref='x', yref='y',
                    fillcolor="LightSalmon",
                    opacity=0.18, layer="below", line_width=0,
                )
        else:
            fig.add_shape(
                type="rect",
                x0=left_base - 0.5, x1=right_base + 0.5, y0=0, y1=1,
                xref='x', yref='paper',
                fillcolor="LightSalmon",
                opacity=0.18, layer="below", line_width=0,
            )

    return fig
