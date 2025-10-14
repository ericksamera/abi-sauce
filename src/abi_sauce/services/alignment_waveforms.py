#!/usr/bin/env python3
from __future__ import annotations
from typing import Dict, List, Optional, Sequence, Tuple

from abi_sauce.models import TraceAsset
from abi_sauce.services.trace_processing import compute_base_windows

NAN = float("nan")


def raw_channels_and_windows(
    trace_asset: TraceAsset,
) -> tuple[
    Dict[str, List[float]], List[Tuple[int, int]], int, str, Optional[List[int]]
]:
    """
    Returns (channels, windows, seq_len, seq_upper, quals_or_none).
      - channels: dict A/C/G/T -> list[float]
      - windows: per-base inclusive (L, R) sample windows
      - seq_len: number of bases (from called seq if present, else len(ploc))
      - seq_upper: uppercase called sequence ("" if none)
      - quals_or_none: PHRED per base if present (len == seq_len)
    """
    seq_upper = (trace_asset.sequence or "").upper()
    ploc = list(trace_asset.base_positions or [])
    channels = {b: list(trace_asset.channels.get(b) or []) for b in "ACGT"}
    quals = (
        list(trace_asset.qualities or []) if trace_asset.qualities is not None else None
    )

    total_samples = max((len(v) for v in channels.values() if v), default=0)
    seq_len = len(seq_upper) if seq_upper else (len(ploc) or 0)
    if seq_len == 0:
        return {b: [] for b in "ACGT"}, [], 0, "", None

    # Normalize PLOC length to seq_len
    if len(ploc) < seq_len:
        last = (total_samples - 1) if total_samples else 0
        ploc = (ploc[:] + [last] * (seq_len - len(ploc)))[:seq_len]
    elif len(ploc) > seq_len:
        ploc = ploc[:seq_len]

    windows = compute_base_windows(ploc, total_samples)

    if quals is not None:
        if len(quals) < seq_len:
            quals = quals + [0] * (seq_len - len(quals))
        elif len(quals) > seq_len:
            quals = quals[:seq_len]

    return channels, windows, seq_len, seq_upper, quals


def _map_window_to_column_x(k: int, npts: int) -> List[float]:
    """Map npts samples for alignment column k into [k-0.5, k+0.5], >=2 points."""
    if npts <= 1:
        return [k - 0.5, k + 0.5]
    start = k - 0.5
    step = 1.0 / (npts - 1)
    return [start + j * step for j in range(npts)]


def _resample_linear(segment: Sequence[float], n_points: int) -> List[float]:
    """Linear resample to exactly n_points (duplicate if singleton)."""
    m = len(segment)
    if n_points <= 1:
        v = float(segment[0] if segment else 0.0)
        return [v, v]
    if m <= 1:
        v = float(segment[0] if segment else 0.0)
        return [v] * max(2, n_points)
    out: List[float] = []
    for j in range(n_points):
        t = j * (m - 1) / (n_points - 1)
        i0 = int(t)
        i1 = min(i0 + 1, m - 1)
        frac = t - i0
        out.append((1.0 - frac) * float(segment[i0]) + frac * float(segment[i1]))
    if len(out) == 1:
        out.append(out[0])
    return out


def build_aligned_waveforms(
    columns: List[Tuple[Optional[int], Optional[int]]],
    windows: List[Tuple[int, int]],
    channels: Dict[str, Sequence[float]],
    for_A: bool,
    *,
    uniform_samples_per_base: bool = True,
    samples_per_base: int = 21,
) -> tuple[Dict[str, Tuple[List[float], List[float]]], List[float]]:
    """
    Build per-row waveforms mapped into equidistant alignment columns.

    Returns (series, peak_y):
      - series: {base: (xs, ys)} with NaN separators between columns (Scattergl-friendly)
      - peak_y: local max per column (for hover markers)
    """
    out: Dict[str, Tuple[List[float], List[float]]] = {b: ([], []) for b in "ACGT"}
    peak_y: List[float] = []
    A, C, G, T = channels["A"], channels["C"], channels["G"], channels["T"]

    for k, (iA, iB) in enumerate(columns, start=1):
        idx = iA if for_A else iB
        if idx is None:
            for b in "ACGT":
                out[b][0].append(NAN)
                out[b][1].append(NAN)
            peak_y.append(0.0)
            continue

        L, R = windows[idx]
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
        peak_y.append(float(local_max))

        if uniform_samples_per_base:
            npts = max(2, int(samples_per_base))
            segs = {
                "A": _resample_linear(sA, npts),
                "C": _resample_linear(sC, npts),
                "G": _resample_linear(sG, npts),
                "T": _resample_linear(sT, npts),
            }
        else:
            npts = max(2, R - L + 1)

            def _pad2(arr: Sequence[float]) -> List[float]:
                return [float(arr[0])] * 2 if len(arr) == 1 else [float(v) for v in arr]

            segs = {"A": _pad2(sA), "C": _pad2(sC), "G": _pad2(sG), "T": _pad2(sT)}

        x_segment = _map_window_to_column_x(k, npts)
        for b in "ACGT":
            xs, ys = out[b]
            ys_seg = segs[b]
            if len(ys_seg) != len(x_segment):
                n = min(len(ys_seg), len(x_segment))
                ys_seg = ys_seg[:n]
                xuse = x_segment[:n]
            else:
                xuse = x_segment
            xs.extend(xuse)
            ys.extend(ys_seg)
            xs.append(NAN)
            ys.append(NAN)

    return out, peak_y
