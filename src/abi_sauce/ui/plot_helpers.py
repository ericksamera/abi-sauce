# src/abi_sauce/ui/plot_helpers.py
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import plotly.graph_objects as go

from abi_sauce.services import trace_processing as tp
from abi_sauce.ui import theme

Number = int | float

# --- Plotly/Pandas compatibility shim (handles partially-initialized pandas) ---
try:
    import _plotly_utils.basevalidators as _bv  # type: ignore

    _pd = getattr(_bv, "pd", None)
    if _pd is not None and not hasattr(_pd, "Series"):
        _bv.pd = None
except Exception:
    pass


# ---- small utils --------------------------------------------------------------
def _clamp(v: int, lo: int, hi: int) -> int:
    return tp.clamp(v, lo, hi)


def _to_simple_list(seq):
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


def _decimate_max_bins_xy(
    xs: Sequence[Number], ys: Sequence[Number], target_points: int
) -> tuple[list[float], list[float]]:
    xs = _to_simple_list(xs) or []
    ys = _to_simple_list(ys) or []
    n = min(len(xs), len(ys))
    if n <= target_points or target_points <= 0:
        return xs[:n], ys[:n]
    bin_sz = max(1, n // target_points)
    out_x, out_y = [], []
    for i in range(0, n, bin_sz):
        block_y = ys[i : i + bin_sz]
        if not block_y:
            continue
        m = max(block_y)
        j = block_y.index(m)
        out_y.append(float(m))
        out_x.append(float(xs[i + j]))
    return out_x, out_y


# ---- plotting primitives ------------------------------------------------------
def _safe_add_scatter(fig, xs, ys, **scatter_kwargs):
    xs = _to_simple_list(xs)
    ys = _to_simple_list(ys)
    n = len(ys) if ys is not None else 0

    MAX_POINTS = 30000
    if n > MAX_POINTS:
        bin_sz = max(1, n // MAX_POINTS)
        xs = xs[::bin_sz]
        ys = [max(ys[i : i + bin_sz]) for i in range(0, n, bin_sz)]

    trace_cls = go.Scattergl if (len(ys or []) > 5000) else go.Scatter
    fig.add_trace(trace_cls(x=xs, y=ys, **scatter_kwargs))
    return fig


def _safe_add_bar(fig, xs, ys, **bar_kwargs):
    xs = _to_simple_list(xs)
    ys = _to_simple_list(ys)
    trace = go.Bar(x=xs, y=ys, **bar_kwargs)
    fig.add_trace(trace)
    return fig


def _add_slider_overview_trace(
    fig: go.Figure,
    xs: Sequence[Number],
    ys: Sequence[Number],
    *,
    target_points: int = 3000,
    name: str = "overview",
    color: str | None = None,
) -> None:
    dx, dy = _decimate_max_bins_xy(xs, ys, target_points)
    fig.add_trace(
        go.Scatter(
            x=dx,
            y=dy,
            mode="lines",
            name=name,
            hoverinfo="skip",
            showlegend=False,
            xaxis="x",
            yaxis="y",
            line={"width": 1, "color": color},
            opacity=0.25,
        )
    )


def _add_hatched_band(
    fig: go.Figure,
    *,
    x0: float,
    x1: float,
    yref: str = "paper",
    base_alpha: float = theme.HATCH_BASE_ALPHA,
    stripe_alpha: float = theme.HATCH_STRIPE_ALPHA,
    max_lines: int = 300,
    stripe_layer: str = "above",
) -> None:
    if x1 <= x0:
        return
    fig.add_shape(
        type="rect",
        x0=x0,
        x1=x1,
        y0=0,
        y1=1,
        xref="x",
        yref=yref,
        fillcolor=f"rgba(128,128,128,{base_alpha})",
        layer="below",
        line_width=0,
    )
    width = max(1.0, float(x1 - x0))
    n_lines = min(max_lines, max(8, int(width / 10)))
    step = width / n_lines
    cur = float(x0)
    for _ in range(n_lines + 1):
        fig.add_shape(
            type="line",
            x0=cur,
            x1=cur,
            y0=0,
            y1=1,
            xref="x",
            yref=yref,
            line={"width": 1, "color": f"rgba(120,120,120,{stripe_alpha})"},
            layer=stripe_layer,
        )
        cur += step


# ---- helpers for peak/hover text ---------------------------------------------
def _dominant_base_at_peak(
    p: int, channels: Mapping[str, Sequence[Number]], window: int = 2
) -> str | None:
    best_base, best_val = None, float("-inf")
    for b in "ACGT":
        arr = channels.get(b) or []
        if not arr:
            continue
        s = _clamp(p - window, 0, len(arr) - 1)
        e = _clamp(p + window, 0, len(arr) - 1)
        if e < s:
            continue
        val = max(arr[s : e + 1])
        if val > best_val:
            best_val = float(val)
            best_base = b
    return best_base


def _peak_heights_samples(
    ploc: list[int], channels: Mapping[str, Sequence[Number]], window: int = 2
) -> list[float]:
    heights: list[float] = []
    for p in ploc:
        best: float = 0.0
        for arr in channels.values():
            if not arr:
                continue
            s = _clamp(p - window, 0, len(arr) - 1)
            e = _clamp(p + window, 0, len(arr) - 1)
            if e >= s:
                best = max(best, float(max(arr[s : e + 1])))
        heights.append(best)
    return heights


def _peak_hover_texts(
    ploc: list[int],
    seq_len: int,
    seq: str | None,
    channels: Mapping[str, Sequence[Number]],
    quals: Sequence[Number] | None,
) -> list[str]:
    texts: list[str] = []
    for i in range(seq_len):
        p = ploc[i] if i < len(ploc) else (ploc[-1] if ploc else 0)
        peak_base = _dominant_base_at_peak(p, channels, window=2)
        called = None
        if seq and i < len(seq):
            c = seq[i].upper()
            if c in "ACGTN":
                called = c
        main = peak_base or (called if called and called != "N" else "?")
        extra = (
            f"<br>call: {called}"
            if (called and called != "N" and peak_base and called != peak_base)
            else ""
        )
        qtxt = f"<br>Q{int(quals[i])}" if (quals and i < len(quals)) else ""
        texts.append(f"Base {i+1}: {main}{extra}<br>sample {p}{qtxt}")
    return texts


# ---- public API (split) ------------------------------------------------------
@dataclass
class TraceData:
    channels: Mapping[str, Sequence[Number]]
    ploc: list[int]
    seq_len: int
    quals: Sequence[Number] | None = None
    seq: str | None = None


@dataclass
class TracePlotConfig:
    mode: Literal["samples", "bases"] = "samples"
    resample_method: Literal["max", "mean", "median"] = "max"
    show_trim: bool = True
    trim_range: tuple[int, int] | None = None  # 1-based inclusive
    uirevision: str = "chrom_trim_base_ticks_v1"
    height: int = 340
    show_grid: bool = True
    peak_only_hover: bool = True
    show_rangeslider: bool = True
    initial_base_span: int = 20
    base_tick_every: int = 10
    # NEW: control range application
    x_range: tuple[float, float] | None = None
    use_initial_range: bool = False


def prepare_trace_layers(data: TraceData, cfg: TracePlotConfig) -> dict[str, object]:
    channels = data.channels
    ploc = list(data.ploc or [])
    seq_len = int(data.seq_len if data.seq_len is not None else len(ploc))
    quals = list(data.quals) if data.quals is not None else None
    seq = data.seq

    total_samples = max((len(arr) for arr in channels.values() if arr), default=0)
    if len(ploc) < seq_len:
        ploc = (ploc[:] + [total_samples - 1] * (seq_len - len(ploc)))[:seq_len]
    elif len(ploc) > seq_len:
        ploc = ploc[:seq_len]

    windows = tp.compute_base_windows(ploc, total_samples)
    base_span = max(1, min(cfg.initial_base_span, max(seq_len, 1)))
    max_signal = max((max(arr) for arr in channels.values() if arr), default=0.0)
    y_max = float(max_signal) * 1.08 if max_signal > 0 else 1.0
    peak_htexts = _peak_hover_texts(ploc, seq_len, seq, channels, quals)

    out: dict[str, object] = {
        "seq_len": seq_len,
        "total_samples": total_samples,
        "windows": windows,
        "y_max": y_max,
        "peak_htexts": peak_htexts,
        "quals": None,
    }

    if quals is not None:
        q_list = list(quals)
        if len(q_list) < seq_len:
            q_list.extend([0] * (seq_len - len(q_list)))
        elif len(q_list) > seq_len:
            q_list = q_list[:seq_len]
        out["quals"] = q_list

    if cfg.mode == "samples":
        xs = list(range(total_samples))
        peak_ys = _peak_heights_samples(ploc, channels, window=2)
        if windows:
            x0 = windows[0][0]
            x1 = windows[min(base_span - 1, len(windows) - 1)][1]
        else:
            sppb = total_samples / max(seq_len, 1)
            x0, x1 = 0, int(round(base_span * sppb))
        tick_bases = list(range(cfg.base_tick_every, seq_len + 1, cfg.base_tick_every))
        if ploc:
            tickvals = [ploc[b - 1] for b in tick_bases if (b - 1) < len(ploc)]
        else:
            sppb = total_samples / max(seq_len, 1)
            tickvals = [int(round(b * sppb)) for b in tick_bases]
        out.update(
            {
                "xs_samples": xs,
                "peak_ys": peak_ys,
                "initial_x_range": [x0, x1],
                "tickvals": tickvals,
                "ticktext": [str(b) for b in tick_bases],
                "q_centers": [(s + e) / 2.0 for (s, e) in windows][:seq_len],
                "q_widths": [(e - s + 1) for (s, e) in windows][:seq_len],
            }
        )
    else:
        xs = list(range(1, seq_len + 1))
        if windows:
            resampled = tp.per_base_aggregate(
                channels, windows, method=cfg.resample_method
            )
            peak_ys = [
                max(
                    resampled.get("A", [0] * seq_len)[i],
                    resampled.get("C", [0] * seq_len)[i],
                    resampled.get("G", [0] * seq_len)[i],
                    resampled.get("T", [0] * seq_len)[i],
                )
                for i in range(seq_len)
            ]
        else:
            resampled = {k: [0.0] * seq_len for k in "ACGT"}
            peak_ys = [0.0] * seq_len
        out.update(
            {
                "xs_bases": xs,
                "resampled": resampled,
                "peak_ys": peak_ys,
                "initial_x_range": [0.5, base_span + 0.5],
            }
        )

    return out


def render_trace_layers(
    data: TraceData, cfg: TracePlotConfig, layers: dict[str, object]
) -> go.Figure:
    """
    IMPORTANT: we only set x-axis 'range' when cfg.x_range is provided OR
    cfg.use_initial_range is True. Otherwise, uirevision preserves the user's view.
    """
    fig = go.Figure()
    channels = data.channels

    channel_hoverinfo = "skip" if cfg.peak_only_hover else "x+y+name"
    bar_hoverinfo = "skip" if cfg.peak_only_hover else "x+y"
    invisible_marker = {
        "size": 16,
        "color": "rgba(0,0,0,0)",
        "line": {"width": 0},
        "symbol": "circle",
    }

    # Layout
    fig.update_layout(
        height=cfg.height,
        uirevision=cfg.uirevision,
        dragmode="pan",
        margin=theme.MARGIN_DEFAULT,
        legend=theme.LEGEND_DEFAULT,
        hovermode="closest",
        modebar={"remove": ["lasso2d", "select2d"]},
    )

    # PHRED bars
    q = layers.get("quals")
    if q is not None:
        if cfg.mode == "samples":
            _safe_add_bar(
                fig,
                layers["q_centers"],
                q[: len(layers["q_centers"])],
                name="quality",
                width=layers["q_widths"],
                marker={"opacity": theme.QUALITY_BAR_OPACITY},
                yaxis="y2",
                hoverinfo=bar_hoverinfo,
                showlegend=True,
            )
        else:
            _safe_add_bar(
                fig,
                layers["xs_bases"],
                q,
                name="quality",
                width=0.98,
                marker={"opacity": theme.QUALITY_BAR_OPACITY},
                yaxis="y2",
                hoverinfo=bar_hoverinfo,
                showlegend=True,
            )

    # Hatched OUTSIDE the keep region
    if cfg.show_trim and cfg.trim_range:
        left_base, right_base = map(int, cfg.trim_range)
        if cfg.mode == "samples":
            windows = layers["windows"]
            total = int(layers["total_samples"])
            if 1 < left_base <= len(windows):
                _add_hatched_band(fig, x0=0, x1=windows[left_base - 1][0])
            if 1 <= right_base < len(windows):
                _add_hatched_band(fig, x0=windows[right_base - 1][1], x1=total - 1)
        else:
            seq_len = int(layers["seq_len"])
            if left_base > 1:
                _add_hatched_band(fig, x0=0.5, x1=(left_base - 0.5))
            if right_base < seq_len:
                _add_hatched_band(fig, x0=(right_base + 0.5), x1=(seq_len + 0.5))

    # Slider overview
    if cfg.show_rangeslider:
        if cfg.mode == "samples":
            xs_over = layers["xs_samples"]
            for base_letter, arr in channels.items():
                _add_slider_overview_trace(
                    fig,
                    xs_over,
                    arr or [],
                    target_points=3000,
                    name=f"{base_letter} overview",
                    color=theme.base_color(base_letter),
                )
        else:
            xs_over = layers["xs_bases"]
            resampled = layers["resampled"]
            for base_letter in "ACGT":
                _add_slider_overview_trace(
                    fig,
                    xs_over,
                    resampled.get(base_letter, []),
                    target_points=2000,
                    name=f"{base_letter} overview",
                    color=theme.base_color(base_letter),
                )

    # Main lines + invisible peak markers
    if cfg.mode == "samples":
        xs_main = layers["xs_samples"]
        for base_letter, arr in channels.items():
            _safe_add_scatter(
                fig,
                xs_main,
                arr,
                mode="lines",
                name=base_letter,
                hoverinfo=channel_hoverinfo,
                hovertemplate=(
                    None
                    if cfg.peak_only_hover
                    else "%{y}<extra>%{fullData.name}</extra>"
                ),
                line={"width": 1.5, "color": theme.base_color(base_letter)},
                opacity=theme.WAVE_OPACITY,
            )
        fig.add_trace(
            go.Scatter(
                x=data.ploc,
                y=layers["peak_ys"],
                mode="markers",
                marker=invisible_marker,
                name="peaks",
                hoverinfo="text",
                hovertext=layers["peak_htexts"],
                showlegend=False,
            )
        )
        xaxis_title = "Sample index (base ticks shown below)"
        xaxis_tick_kwargs = {
            "tickmode": "linear",
            "dtick": max(1, int(cfg.base_tick_every)),
            "tick0": cfg.base_tick_every,
        }

    else:
        xs_main = layers["xs_bases"]
        resampled = layers["resampled"]
        for base_letter in "ACGT":
            _safe_add_scatter(
                fig,
                xs_main,
                resampled.get(base_letter, [0] * len(xs_main)),
                mode="lines",
                name=base_letter,
                hoverinfo=channel_hoverinfo,
                hovertemplate=(
                    None
                    if cfg.peak_only_hover
                    else "%{y}<extra>%{fullData.name}</extra>"
                ),
                line={"width": 1.5, "color": theme.base_color(base_letter)},
            )
        fig.add_trace(
            go.Scatter(
                x=xs_main[: len(layers["peak_htexts"])],
                y=layers["peak_ys"][: len(layers["peak_htexts"])],
                mode="markers",
                marker=invisible_marker,
                name="peaks",
                hoverinfo="text",
                hovertext=layers["peak_htexts"],
                showlegend=False,
            )
        )
        xaxis_title = "Base number"
        xaxis_tick_kwargs = {
            "tickmode": "linear",
            "dtick": max(1, int(cfg.base_tick_every)),
            "tick0": cfg.base_tick_every,
        }

    # Axes — apply range only when requested
    x_range_to_apply: tuple[float, float] | None = None
    if cfg.x_range is not None:
        x_range_to_apply = cfg.x_range
    elif cfg.use_initial_range:
        x0, x1 = layers["initial_x_range"]
        x_range_to_apply = (float(x0), float(x1))

    fig.update_xaxes(
        title_text=xaxis_title,
        **({"range": x_range_to_apply} if x_range_to_apply is not None else {}),
        showgrid=False,
        zeroline=False,
        showticklabels=True,
        rangeslider={
            "visible": cfg.show_rangeslider,
            "thickness": theme.RANGE_SLIDER_THICKNESS,
        },
        **xaxis_tick_kwargs,
    )
    fig.update_yaxes(
        title_text="Signal intensity",
        range=[0, layers["y_max"]],
        showgrid=cfg.show_grid,
        zeroline=False,
    )

    if layers.get("quals") is not None:
        qmax = max(layers["quals"]) if layers["quals"] else 1
        fig.update_layout(
            yaxis2={
                "title": "PHRED",
                "overlaying": "y",
                "side": "right",
                "range": [0, max(1, float(qmax) * 1.05)],
                "showgrid": False,
            },
        )

    return fig


# ---- Back-compat wrapper ------------------------------------------------------
def build_trace_fig(
    channels: Mapping[str, Sequence[Number]],
    ploc: list[int],
    seq_len: int | None = None,
    quals: Sequence[Number] | None = None,
    seq: str | None = None,
    mode: Literal["samples", "bases"] = "samples",
    resample_method: Literal["max", "mean", "median"] = "max",
    show_trim: bool = True,
    trim_range: tuple[int, int] | None = None,
    uirevision: str = "chrom_trim_base_ticks_v1",
    height: int = 340,
    show_grid: bool = True,
    peak_only_hover: bool = True,
    show_rangeslider: bool = True,
    *,
    initial_base_span: int = 20,
    base_tick_every: int = 10,
    x_range: tuple[float, float] | None = None,
    use_initial_range: bool = False,
) -> go.Figure:
    data = TraceData(
        channels=channels,
        ploc=ploc,
        seq_len=(len(ploc) if seq_len is None else seq_len),
        quals=quals,
        seq=seq,
    )
    cfg = TracePlotConfig(
        mode=mode,
        resample_method=resample_method,
        show_trim=show_trim,
        trim_range=trim_range,
        uirevision=uirevision,
        height=height,
        show_grid=show_grid,
        peak_only_hover=peak_only_hover,
        show_rangeslider=show_rangeslider,
        initial_base_span=initial_base_span,
        base_tick_every=base_tick_every,
        x_range=x_range,
        use_initial_range=use_initial_range,
    )
    layers = prepare_trace_layers(data, cfg)
    return render_trace_layers(data, cfg, layers)
