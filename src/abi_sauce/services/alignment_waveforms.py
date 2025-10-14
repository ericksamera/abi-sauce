#!/usr/bin/env python3
from __future__ import annotations

from collections.abc import Sequence

from abi_sauce.models import TraceAsset
from abi_sauce.services.trace_processing import compute_base_windows
from abi_sauce.types import BaseChar

NAN: float = float("nan")


def _to_float_list(seq: Sequence[float] | Sequence[int] | None) -> list[float]:
    if not seq:
        return []
    return [float(x) for x in seq]


def raw_channels_and_windows(
    trace_asset: TraceAsset,
) -> tuple[
    dict[BaseChar, list[float]], list[tuple[int, int]], int, str, list[int] | None
]:
    """
    Prepare float channel arrays and base windows for a TraceAsset.

    Returns:
      channels: per-base-letter float arrays
      windows: inclusive sample windows per base (via PLOC midpoint partitioning)
      total_samples: length of the raw trace vectors
      abif_order: the ABI-reported channel order (if known)
      ploc: raw PLOC list (peak sample indices) if present
    """
    ch_raw = trace_asset.channels or {}
    A = _to_float_list(ch_raw.get("A"))
    C = _to_float_list(ch_raw.get("C"))
    G = _to_float_list(ch_raw.get("G"))
    T = _to_float_list(ch_raw.get("T"))
    channels: dict[BaseChar, list[float]] = {"A": A, "C": C, "G": G, "T": T}

    total_samples = max((len(A), len(C), len(G), len(T)), default=0)
    ploc = list(trace_asset.base_positions or [])
    abif_order = str((trace_asset.meta or {}).get("abif_order", "ACGT"))

    windows: list[tuple[int, int]] = []
    if ploc and total_samples:
        windows = compute_base_windows(ploc, total_samples)

    return channels, windows, total_samples, abif_order, (ploc or None)


def _linspace(n_points: int, start: float, end: float) -> list[float]:
    """Inclusive linspace-like helper that always returns at least 2 points."""
    if n_points <= 1:
        return [start, end]
    step = (end - start) / float(n_points - 1)
    return [start + i * step for i in range(n_points)]


def _segment_max(values: Sequence[float]) -> float:
    mv = 0.0
    for v in values:
        fv = float(v)
        if fv > mv:
            mv = fv
    return mv


def build_aligned_waveforms(
    *,
    columns: list[tuple[int | None, int | None]],
    windows: list[tuple[int, int]],
    channels: dict[BaseChar, Sequence[float]],
    for_A: bool,
    uniform_samples_per_base: bool = True,
    samples_per_base: int = 21,
) -> tuple[dict[BaseChar, tuple[list[float], list[float]]], list[float]]:
    """
    Map raw sample-space channels into alignment-column space.

    Returns (series, peak_y):
      - series: {base: (xs, ys)} with NaN separators between columns (Scattergl-friendly)
      - peak_y: local max per column (for hover markers)
    """
    out: dict[BaseChar, tuple[list[float], list[float]]] = {
        "A": ([], []),
        "C": ([], []),
        "G": ([], []),
        "T": ([], []),
    }
    peak_y: list[float] = []

    chA = list(channels.get("A", []))
    chC = list(channels.get("C", []))
    chG = list(channels.get("G", []))
    chT = list(channels.get("T", []))
    ch_map: dict[BaseChar, list[float]] = {"A": chA, "C": chC, "G": chG, "T": chT}

    for k, (iQ, iR) in enumerate(columns, start=1):
        base_idx = iQ if for_A else iR
        x0 = float(k) - 0.5
        x1 = float(k) + 0.5

        if base_idx is None or base_idx < 0 or base_idx >= len(windows):
            for b in ("A", "C", "G", "T"):
                xs, ys = out[b]  # type: ignore[index]
                xs.append(NAN)
                ys.append(NAN)
            peak_y.append(0.0)
            continue

        s, e = windows[base_idx]
        if e < s:
            s, e = e, s
        seg_len = max(1, (e - s + 1))

        if uniform_samples_per_base:
            n_points = max(2, int(samples_per_base))
            x_seg = _linspace(n_points, x0, x1)
            local_max = 0.0
            for b in ("A", "C", "G", "T"):
                src = ch_map[b]
                ys_seg: list[float] = []
                if seg_len == 1:
                    v = float(src[s]) if s < len(src) else 0.0
                    ys_seg = [v, v]
                else:
                    for t in range(n_points):
                        pos = s + (t * (seg_len - 1)) / float(n_points - 1)
                        i0 = int(pos)
                        i1 = min(i0 + 1, len(src) - 1) if len(src) > 0 else i0
                        frac = pos - i0
                        v0 = float(src[i0]) if 0 <= i0 < len(src) else 0.0
                        v1 = float(src[i1]) if 0 <= i1 < len(src) else 0.0
                        ys_seg.append((1.0 - frac) * v0 + frac * v1)
                xs, ys = out[b]  # type: ignore[index]
                xs.extend(x_seg)
                ys.extend(ys_seg)
                xs.append(NAN)
                ys.append(NAN)
                local_max = max(local_max, _segment_max(ys_seg))
            peak_y.append(local_max)
        else:
            x_seg = _linspace(seg_len, x0, x1)
            local_max = 0.0
            for b in ("A", "C", "G", "T"):
                src = ch_map[b]
                ys_seg = [
                    float(src[i]) if 0 <= i < len(src) else 0.0 for i in range(s, e + 1)
                ]
                xs, ys = out[b]  # type: ignore[index]
                xs.extend(x_seg)
                ys.extend(ys_seg)
                xs.append(NAN)
                ys.append(NAN)
                local_max = max(local_max, _segment_max(ys_seg))
            peak_y.append(local_max)

    return out, peak_y
