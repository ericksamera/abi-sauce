from __future__ import annotations

from collections.abc import Sequence


def clamp_unit(value: float) -> float:
    """Clamp one numeric value into the closed unit interval."""
    return max(0.0, min(value, 1.0))


def linspace(start: float, end: float, count: int) -> tuple[float, ...]:
    """Return one inclusive evenly spaced float sequence."""
    if count < 2:
        raise ValueError("count must be >= 2")
    start_float = float(start)
    step = (float(end) - start_float) / float(count - 1)
    return tuple(start_float + (step * float(index)) for index in range(count))


def interpolate_signal(signal: Sequence[int | float], position: float) -> float:
    """Return one linearly interpolated signal value at a raw sample position."""
    if not signal:
        return 0.0
    if len(signal) == 1:
        return float(signal[0])

    clamped_position = max(0.0, min(float(position), float(len(signal) - 1)))
    left_index = int(clamped_position)
    right_index = min(left_index + 1, len(signal) - 1)
    if right_index == left_index:
        return float(signal[left_index])

    fraction = clamped_position - float(left_index)
    left_value = float(signal[left_index])
    right_value = float(signal[right_index])
    return left_value + ((right_value - left_value) * fraction)


def resample_signal_window(
    signal: Sequence[int | float],
    *,
    raw_left: float,
    raw_right: float,
    cell_left: float,
    cell_right: float,
    sample_count: int,
    signal_scale: float = 1.0,
    clamp_to_unit: bool = False,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Project one raw trace window onto one fixed-width display window."""
    if sample_count < 2:
        raise ValueError("sample_count must be >= 2")
    if cell_right <= cell_left:
        raise ValueError("cell_right must be > cell_left")
    if raw_right <= raw_left:
        raise ValueError("raw_right must be > raw_left")

    x_values = linspace(cell_left, cell_right, sample_count)
    raw_positions = linspace(raw_left, raw_right, sample_count)
    resolved_signal_scale = float(signal_scale) if signal_scale > 0 else 1.0

    sampled_signal = tuple(
        interpolate_signal(signal, raw_position) / resolved_signal_scale
        for raw_position in raw_positions
    )
    if clamp_to_unit:
        sampled_signal = tuple(clamp_unit(value) for value in sampled_signal)

    return x_values, sampled_signal
