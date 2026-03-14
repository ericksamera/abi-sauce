from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

from abi_sauce.models import SequenceOrientation, SequenceRecord, TraceData
from abi_sauce.orientation import complement_base
from abi_sauce.trimming import TrimResult

_TRACE_CHANNEL_KEYS: Final[tuple[str, ...]] = ("DATA9", "DATA10", "DATA11", "DATA12")
_DEFAULT_CHANNEL_ORDER: Final[str] = "GATC"
_CHANNEL_COLORS: Final[dict[str, str]] = {
    "A": "green",
    "C": "blue",
    "G": "black",
    "T": "red",
}


@dataclass(frozen=True, slots=True)
class ChromatogramChannel:
    """One normalized trace channel ready for plotting."""

    data_key: str
    base: str
    color: str
    signal: tuple[int, ...] = ()


@dataclass(frozen=True, slots=True)
class ChromatogramBaseCall:
    """One called base positioned along the raw trace x-axis."""

    base_index: int
    base: str
    position: int
    color: str


@dataclass(frozen=True, slots=True)
class ChromatogramBaseSpan:
    """One called base projected onto a bounded raw trace x-interval."""

    base_index: int
    left: float
    right: float
    center: float
    width: float


@dataclass(frozen=True, slots=True)
class ChromatogramQualitySegment:
    """One per-base quality score rendered across a chromatogram span."""

    base_index: int
    left: float
    right: float
    center: float
    width: float
    quality: int


@dataclass(frozen=True, slots=True)
class ChromatogramTrimBoundaries:
    """Trim-boundary markers expressed in raw trace sample coordinates."""

    left: float | None = None
    right: float | None = None


@dataclass(frozen=True, slots=True)
class ChromatogramView:
    """Pure plotting/view-model state for one chromatogram."""

    is_renderable: bool
    render_failure_reason: str | None = None
    x_values: tuple[int, ...] = ()
    channels: tuple[ChromatogramChannel, ...] = ()
    base_calls: tuple[ChromatogramBaseCall, ...] = ()
    base_spans: tuple[ChromatogramBaseSpan, ...] = ()
    quality_segments: tuple[ChromatogramQualitySegment, ...] = ()
    trim_boundaries: ChromatogramTrimBoundaries = field(
        default_factory=ChromatogramTrimBoundaries
    )
    retained_sample_range: tuple[float, float] | None = None
    has_any_retained_samples: bool = True

    @property
    def trace_length(self) -> int:
        """Return the shared sanitized trace length."""
        return len(self.x_values)

    @property
    def has_quality_overlay(self) -> bool:
        """Return whether quality markers are available for plotting."""
        return bool(self.quality_segments)


def build_chromatogram_view(
    record: SequenceRecord,
    trim_result: TrimResult | None = None,
) -> ChromatogramView:
    """Build a normalized chromatogram view model from one sequence record."""
    trace_data = record.trace_data
    if trace_data is None:
        return ChromatogramView(
            is_renderable=False,
            render_failure_reason="missing_trace_data",
        )

    channels = _build_channels(trace_data)
    if not channels:
        return ChromatogramView(
            is_renderable=False,
            render_failure_reason="missing_trace_channels",
        )

    trace_length = min(len(channel.signal) for channel in channels)
    if trace_length <= 0:
        return ChromatogramView(
            is_renderable=False,
            render_failure_reason="empty_trace_channels",
        )

    x_values = tuple(range(trace_length))
    sanitized_channels = tuple(
        ChromatogramChannel(
            data_key=channel.data_key,
            base=channel.base,
            color=channel.color,
            signal=channel.signal[:trace_length],
        )
        for channel in channels
    )
    base_calls = _build_base_calls(
        sequence=record.sequence,
        base_positions=trace_data.base_positions,
        trace_length=trace_length,
    )
    base_spans = _build_base_spans(
        base_calls,
        trace_length=trace_length,
    )
    quality_segments = _build_quality_segments(
        record.qualities,
        base_spans,
    )
    trim_boundaries = _resolve_trim_boundaries(
        trim_result=trim_result,
        base_calls=base_calls,
        trace_length=trace_length,
    )
    retained_sample_range, has_any_retained_samples = _resolve_retained_sample_range(
        trim_result=trim_result,
        base_calls=base_calls,
        trace_length=trace_length,
    )

    raw_view = ChromatogramView(
        is_renderable=True,
        x_values=x_values,
        channels=sanitized_channels,
        base_calls=base_calls,
        base_spans=base_spans,
        quality_segments=quality_segments,
        trim_boundaries=trim_boundaries,
        retained_sample_range=retained_sample_range,
        has_any_retained_samples=has_any_retained_samples,
    )
    return orient_chromatogram_view(raw_view, record.orientation)


def orient_chromatogram_view(
    view: ChromatogramView,
    orientation: SequenceOrientation,
) -> ChromatogramView:
    """Return a chromatogram view in the requested orientation."""
    return _apply_orientation_to_view(view, orientation)


def reverse_complement_chromatogram_view(
    view: ChromatogramView,
) -> ChromatogramView:
    """Return a reverse-complemented chromatogram view."""
    return orient_chromatogram_view(view, "reverse_complement")


def _build_channels(trace_data: TraceData) -> tuple[ChromatogramChannel, ...]:
    channel_order = _normalize_channel_order(trace_data.channel_order)
    base_by_data_key = dict(zip(_TRACE_CHANNEL_KEYS, channel_order, strict=True))

    channels: list[ChromatogramChannel] = []
    for data_key in _TRACE_CHANNEL_KEYS:
        signal = trace_data.channels.get(data_key)
        if signal is None:
            continue

        signal_tuple = tuple(int(value) for value in signal)
        if not signal_tuple:
            continue

        base = base_by_data_key[data_key]
        channels.append(
            ChromatogramChannel(
                data_key=data_key,
                base=base,
                color=_CHANNEL_COLORS[base],
                signal=signal_tuple,
            )
        )

    return tuple(channels)


def _normalize_channel_order(channel_order: str | None) -> str:
    if channel_order is None:
        return _DEFAULT_CHANNEL_ORDER

    normalized = channel_order.strip().upper()
    if len(normalized) != len(_TRACE_CHANNEL_KEYS):
        return _DEFAULT_CHANNEL_ORDER
    if set(normalized) != {"A", "C", "G", "T"}:
        return _DEFAULT_CHANNEL_ORDER
    return normalized


def _build_base_calls(
    *,
    sequence: str,
    base_positions: list[int],
    trace_length: int,
) -> tuple[ChromatogramBaseCall, ...]:
    base_calls: list[ChromatogramBaseCall] = []
    usable_base_count = min(len(sequence), len(base_positions))

    for base_index in range(usable_base_count):
        position = int(base_positions[base_index])
        if position < 0 or position >= trace_length:
            continue

        base = sequence[base_index].upper()
        base_calls.append(
            ChromatogramBaseCall(
                base_index=base_index,
                base=base,
                position=position,
                color=_CHANNEL_COLORS.get(base, "magenta"),
            )
        )

    return tuple(base_calls)


def _build_base_spans(
    base_calls: tuple[ChromatogramBaseCall, ...],
    *,
    trace_length: int,
) -> tuple[ChromatogramBaseSpan, ...]:
    if not base_calls:
        return ()

    positions = tuple(base_call.position for base_call in base_calls)
    base_spans: list[ChromatogramBaseSpan] = []
    for position_index, base_call in enumerate(base_calls):
        left = _quality_left_edge(
            positions,
            position_index=position_index,
            trace_length=trace_length,
        )
        right = _quality_right_edge(
            positions,
            position_index=position_index,
            trace_length=trace_length,
        )
        width = max(right - left, 0.0)
        if width <= 0:
            continue

        base_spans.append(
            ChromatogramBaseSpan(
                base_index=base_call.base_index,
                left=left,
                right=right,
                center=float(base_call.position),
                width=width,
            )
        )

    return tuple(base_spans)


def _build_quality_segments(
    qualities: list[int] | None,
    base_spans: tuple[ChromatogramBaseSpan, ...],
) -> tuple[ChromatogramQualitySegment, ...]:
    if qualities is None or not base_spans:
        return ()

    quality_segments: list[ChromatogramQualitySegment] = []
    for base_span in base_spans:
        if base_span.base_index >= len(qualities):
            break

        quality_segments.append(
            ChromatogramQualitySegment(
                base_index=base_span.base_index,
                left=base_span.left,
                right=base_span.right,
                center=base_span.center,
                width=base_span.width,
                quality=int(qualities[base_span.base_index]),
            )
        )

    return tuple(quality_segments)


def _quality_left_edge(
    positions: tuple[int, ...],
    *,
    position_index: int,
    trace_length: int,
) -> float:
    if position_index <= 0:
        return _clamp_boundary(
            _extrapolated_left_edge(positions),
            trace_length=trace_length,
        )
    return _clamp_boundary(
        _midpoint(positions[position_index - 1], positions[position_index]),
        trace_length=trace_length,
    )


def _quality_right_edge(
    positions: tuple[int, ...],
    *,
    position_index: int,
    trace_length: int,
) -> float:
    if position_index >= len(positions) - 1:
        return _clamp_boundary(
            _extrapolated_right_edge(positions),
            trace_length=trace_length,
        )
    return _clamp_boundary(
        _midpoint(positions[position_index], positions[position_index + 1]),
        trace_length=trace_length,
    )


def _resolve_trim_boundaries(
    *,
    trim_result: TrimResult | None,
    base_calls: tuple[ChromatogramBaseCall, ...],
    trace_length: int,
) -> ChromatogramTrimBoundaries:
    if trim_result is None or not base_calls:
        return ChromatogramTrimBoundaries()

    positions = tuple(base_call.position for base_call in base_calls)
    base_count = len(positions)
    left_removed = min(max(trim_result.bases_removed_left, 0), base_count)
    right_removed = min(max(trim_result.bases_removed_right, 0), base_count)
    kept_end = max(left_removed, base_count - right_removed)

    left_boundary = _left_trim_boundary(
        positions=positions,
        bases_removed_left=left_removed,
        trace_length=trace_length,
    )
    right_boundary = _right_trim_boundary(
        positions=positions,
        kept_end=kept_end,
        base_count=base_count,
        trace_length=trace_length,
    )

    return ChromatogramTrimBoundaries(
        left=left_boundary,
        right=right_boundary,
    )


def _resolve_retained_sample_range(
    *,
    trim_result: TrimResult | None,
    base_calls: tuple[ChromatogramBaseCall, ...],
    trace_length: int,
) -> tuple[tuple[float, float] | None, bool]:
    if trim_result is None or not base_calls:
        return (None, True)

    positions = tuple(base_call.position for base_call in base_calls)
    base_count = len(positions)
    left_removed = min(max(trim_result.bases_removed_left, 0), base_count)
    right_removed = min(max(trim_result.bases_removed_right, 0), base_count)
    kept_end = max(left_removed, base_count - right_removed)

    if trim_result.trimmed_length <= 0 or kept_end <= left_removed:
        return (None, False)

    retained_left = _quality_left_edge(
        positions,
        position_index=left_removed,
        trace_length=trace_length,
    )
    retained_right = _quality_right_edge(
        positions,
        position_index=kept_end - 1,
        trace_length=trace_length,
    )
    if retained_right <= retained_left:
        return (None, False)

    return ((retained_left, retained_right), True)


def _apply_orientation_to_view(
    view: ChromatogramView,
    orientation: SequenceOrientation,
) -> ChromatogramView:
    if orientation == "forward" or not view.is_renderable:
        return view

    trace_length = view.trace_length
    return ChromatogramView(
        is_renderable=True,
        x_values=view.x_values,
        channels=_reverse_complement_channels(view.channels),
        base_calls=_reverse_complement_base_calls(
            view.base_calls,
            trace_length=trace_length,
        ),
        base_spans=_reverse_complement_base_spans(
            view.base_spans,
            trace_length=trace_length,
        ),
        quality_segments=_reverse_complement_quality_segments(
            view.quality_segments,
            trace_length=trace_length,
        ),
        trim_boundaries=_reverse_complement_trim_boundaries(
            view.trim_boundaries,
            trace_length=trace_length,
        ),
        retained_sample_range=_reverse_complement_retained_sample_range(
            view.retained_sample_range,
            trace_length=trace_length,
        ),
        has_any_retained_samples=view.has_any_retained_samples,
    )


def _reverse_complement_channels(
    channels: tuple[ChromatogramChannel, ...],
) -> tuple[ChromatogramChannel, ...]:
    return tuple(
        ChromatogramChannel(
            data_key=channel.data_key,
            base=_complement_display_base(channel.base),
            color=_display_color_for_base(_complement_display_base(channel.base)),
            signal=tuple(reversed(channel.signal)),
        )
        for channel in channels
    )


def _reverse_complement_base_calls(
    base_calls: tuple[ChromatogramBaseCall, ...],
    *,
    trace_length: int,
) -> tuple[ChromatogramBaseCall, ...]:
    oriented_base_calls: list[ChromatogramBaseCall] = []
    for display_index, base_call in enumerate(reversed(base_calls)):
        base = _complement_display_base(base_call.base)
        oriented_base_calls.append(
            ChromatogramBaseCall(
                base_index=display_index,
                base=base,
                position=int(
                    _mirror_coordinate(
                        base_call.position,
                        trace_length=trace_length,
                    )
                ),
                color=_display_color_for_base(base),
            )
        )
    return tuple(oriented_base_calls)


def _reverse_complement_base_spans(
    base_spans: tuple[ChromatogramBaseSpan, ...],
    *,
    trace_length: int,
) -> tuple[ChromatogramBaseSpan, ...]:
    oriented_base_spans: list[ChromatogramBaseSpan] = []
    for display_index, base_span in enumerate(reversed(base_spans)):
        left = _mirror_coordinate(base_span.right, trace_length=trace_length)
        right = _mirror_coordinate(base_span.left, trace_length=trace_length)
        oriented_base_spans.append(
            ChromatogramBaseSpan(
                base_index=display_index,
                left=left,
                right=right,
                center=_mirror_coordinate(
                    base_span.center,
                    trace_length=trace_length,
                ),
                width=max(right - left, 0.0),
            )
        )
    return tuple(oriented_base_spans)


def _reverse_complement_quality_segments(
    quality_segments: tuple[ChromatogramQualitySegment, ...],
    *,
    trace_length: int,
) -> tuple[ChromatogramQualitySegment, ...]:
    oriented_quality_segments: list[ChromatogramQualitySegment] = []
    for display_index, quality_segment in enumerate(reversed(quality_segments)):
        left = _mirror_coordinate(quality_segment.right, trace_length=trace_length)
        right = _mirror_coordinate(quality_segment.left, trace_length=trace_length)
        oriented_quality_segments.append(
            ChromatogramQualitySegment(
                base_index=display_index,
                left=left,
                right=right,
                center=_midpoint(left, right),
                width=max(right - left, 0.0),
                quality=quality_segment.quality,
            )
        )
    return tuple(oriented_quality_segments)


def _reverse_complement_trim_boundaries(
    trim_boundaries: ChromatogramTrimBoundaries,
    *,
    trace_length: int,
) -> ChromatogramTrimBoundaries:
    return ChromatogramTrimBoundaries(
        left=_mirrored_boundary(trim_boundaries.right, trace_length=trace_length),
        right=_mirrored_boundary(trim_boundaries.left, trace_length=trace_length),
    )


def _reverse_complement_retained_sample_range(
    retained_sample_range: tuple[float, float] | None,
    *,
    trace_length: int,
) -> tuple[float, float] | None:
    if retained_sample_range is None:
        return None

    retained_left, retained_right = retained_sample_range
    return (
        _mirror_coordinate(retained_right, trace_length=trace_length),
        _mirror_coordinate(retained_left, trace_length=trace_length),
    )


def _complement_display_base(base: str) -> str:
    return complement_base(base).upper()


def _display_color_for_base(base: str) -> str:
    return _CHANNEL_COLORS.get(base, "magenta")


def _mirrored_boundary(
    boundary: float | None,
    *,
    trace_length: int,
) -> float | None:
    if boundary is None:
        return None
    return _mirror_coordinate(boundary, trace_length=trace_length)


def _mirror_coordinate(value: float, *, trace_length: int) -> float:
    return float(trace_length - 1) - value


def _left_trim_boundary(
    *,
    positions: tuple[int, ...],
    bases_removed_left: int,
    trace_length: int,
) -> float | None:
    if bases_removed_left <= 0:
        return None
    if bases_removed_left >= len(positions):
        return _clamp_boundary(
            _extrapolated_right_edge(positions),
            trace_length=trace_length,
        )
    return _clamp_boundary(
        _midpoint(positions[bases_removed_left - 1], positions[bases_removed_left]),
        trace_length=trace_length,
    )


def _right_trim_boundary(
    *,
    positions: tuple[int, ...],
    kept_end: int,
    base_count: int,
    trace_length: int,
) -> float | None:
    if kept_end >= base_count:
        return None
    if kept_end <= 0:
        return _clamp_boundary(
            _extrapolated_left_edge(positions),
            trace_length=trace_length,
        )
    return _clamp_boundary(
        _midpoint(positions[kept_end - 1], positions[kept_end]),
        trace_length=trace_length,
    )


def _midpoint(left: float, right: float) -> float:
    return (left + right) / 2.0


def _extrapolated_left_edge(positions: tuple[int, ...]) -> float:
    if len(positions) == 1:
        return positions[0] / 2.0
    step = positions[1] - positions[0]
    return positions[0] - (step / 2.0)


def _extrapolated_right_edge(positions: tuple[int, ...]) -> float:
    if len(positions) == 1:
        return positions[0] + 0.5
    step = positions[-1] - positions[-2]
    return positions[-1] + (step / 2.0)


def _clamp_boundary(value: float, *, trace_length: int) -> float:
    if trace_length <= 0:
        return value
    return max(0.0, min(value, float(trace_length - 1)))
