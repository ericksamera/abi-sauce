# src/abi_sauce/services/trace_processing.py
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

Number = int | float


def clamp(v: int, lo: int, hi: int) -> int:
    return max(lo, min(hi, v))


# ---------- Per-sample transforms ----------


def baseline_min(channel: Sequence[Number]) -> list[float]:
    """Subtract the minimum (simple baseline correction)."""
    if not channel:
        return []
    m = float(min(channel))
    return [float(x) - m for x in channel]


def smooth_moving_average(channel: Sequence[Number], window: int = 5) -> list[float]:
    """Centered moving average (odd window)."""
    n = len(channel)
    if n == 0 or window <= 1:
        return list(map(float, channel))
    w = max(1, int(window))
    w -= 1 if w % 2 == 0 else 0
    r = w // 2
    out: list[float] = []
    pref = [0.0]
    for x in channel:
        pref.append(pref[-1] + float(x))
    for i in range(n):
        s = clamp(i - r, 0, n - 1)
        e = clamp(i + r, 0, n - 1)
        # adjust indices for prefix sum (e inclusive)
        val = (pref[e + 1] - pref[s]) / (e - s + 1)
        out.append(val)
    return out


def normalize_0_1(channel: Sequence[Number]) -> list[float]:
    """Scale to [0,1]."""
    if not channel:
        return []
    lo = float(min(channel))
    hi = float(max(channel))
    span = hi - lo if hi > lo else 1.0
    return [(float(x) - lo) / span for x in channel]


def downsample_max_bins(
    channel: Sequence[Number], target_points: int = 8000
) -> list[float]:
    """Downsample by taking max in equal-width bins to approximately target_points."""
    n = len(channel)
    if n <= target_points or target_points <= 0:
        return list(map(float, channel))
    binsz = max(1, n // target_points)
    out: list[float] = []
    for i in range(0, n, binsz):
        out.append(float(max(channel[i : i + binsz])))
    return out


# ---------- Base-window aggregation ----------


def compute_base_windows(ploc: list[int], total_samples: int) -> list[tuple[int, int]]:
    """Midpoint-partitioning into inclusive sample windows per base."""
    n = len(ploc)
    if n == 0:
        return []
    windows: list[tuple[int, int]] = []
    mids = [(ploc[i] + ploc[i + 1]) / 2.0 for i in range(n - 1)]
    for i in range(n):
        start = 0 if i == 0 else (math.floor(mids[i - 1]) + 1)
        end = (total_samples - 1) if i == n - 1 else math.floor(mids[i])
        start = clamp(start, 0, total_samples - 1)
        end = clamp(end, 0, total_samples - 1)
        if end < start:
            end = start
        windows.append((start, end))
    return windows


def resample_channel_to_bases(
    channel: Sequence[Number],
    windows: list[tuple[int, int]],
    method: str = "max",
) -> list[float]:
    """Aggregate per window."""
    import statistics

    res: list[float] = []
    n = len(channel)
    for s, e in windows:
        s = clamp(s, 0, n - 1)
        e = clamp(e, 0, n - 1)
        if e < s:
            e = s
        block = channel[s : e + 1]
        if not block:
            res.append(0.0)
            continue
        if method == "mean":
            res.append(float(sum(block) / len(block)))
        elif method == "median":
            res.append(float(statistics.median(block)))
        else:
            res.append(float(max(block)))
    return res


def per_base_aggregate(
    channels: Mapping[str, Sequence[Number]],
    windows: list[tuple[int, int]],
    method: str = "max",
) -> dict[str, list[float]]:
    """Aggregate every channel to bases."""
    out: dict[str, list[float]] = {}
    for b, arr in channels.items():
        out[b] = resample_channel_to_bases(arr or [], windows, method=method)
    return out
