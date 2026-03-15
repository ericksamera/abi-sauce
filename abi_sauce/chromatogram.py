from __future__ import annotations

from dataclasses import dataclass, field
import math
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


@dataclass(frozen=True, slots=True)
class ChromatogramColumnChannel:
    """One resampled trace-channel segment inside a fixed-width base column."""

    base: str
    color: str
    x_values: tuple[float, ...] = ()
    signal: tuple[float, ...] = ()


@dataclass(frozen=True, slots=True)
class ChromatogramColumn:
    """One called base projected onto a fixed-width base-index column."""

    column_index: int
    base_index: int
    query_pos: int
    base: str
    color: str
    quality: int | None
    trace_x: int | None
    peak_height: int | None
    is_retained: bool
    cell_left: float
    cell_right: float
    cell_center: float
    raw_left: float | None = None
    raw_right: float | None = None
    raw_center: float | None = None
    channels: tuple[ChromatogramColumnChannel, ...] = ()

    @property
    def has_trace_signal(self) -> bool:
        """Return whether the column carries any resampled trace segments."""
        return bool(self.channels)


@dataclass(frozen=True, slots=True)
class ChromatogramColumnView:
    """Pure fixed-width base-column view for one chromatogram."""

    is_renderable: bool
    render_failure_reason: str | None = None
    columns: tuple[ChromatogramColumn, ...] = ()
    trim_boundaries: ChromatogramTrimBoundaries = field(
        default_factory=ChromatogramTrimBoundaries
    )
    cell_width: float = 1.0
    samples_per_base: int = 16

    @property
    def base_count(self) -> int:
        """Return the number of fixed-width base columns in the view."""
        return len(self.columns)

    @property
    def has_quality_overlay(self) -> bool:
        """Return whether any base column carries a quality value."""
        return any(column.quality is not None for column in self.columns)

    @property
    def has_any_retained_bases(self) -> bool:
        """Return whether any base column falls inside the retained trim window."""
        return any(column.is_retained for column in self.columns)

    @property
    def x_range(self) -> tuple[float, float]:
        """Return the x-extent in base-column coordinates."""
        return (0.0, float(self.base_count) * self.cell_width)


def build_chromatogram_column_view(
    record: SequenceRecord,
    trim_result: TrimResult | None = None,
    *,
    cell_width: float = 1.0,
    samples_per_base: int = 16,
) -> ChromatogramColumnView:
    """Project one chromatogram onto a fixed-width base-column grid."""
    if cell_width <= 0:
        raise ValueError("cell_width must be > 0")
    if samples_per_base < 2:
        raise ValueError("samples_per_base must be >= 2")

    raw_view = build_chromatogram_view(record, trim_result)
    if not raw_view.is_renderable:
        return ChromatogramColumnView(
            is_renderable=False,
            render_failure_reason=raw_view.render_failure_reason,
            cell_width=cell_width,
            samples_per_base=samples_per_base,
        )

    base_spans_by_index = {
        base_span.base_index: base_span for base_span in raw_view.base_spans
    }
    quality_segments_by_index = {
        quality_segment.base_index: quality_segment
        for quality_segment in raw_view.quality_segments
    }

    columns: list[ChromatogramColumn] = []
    for column_index, base_call in enumerate(raw_view.base_calls, start=1):
        cell_left = float(column_index - 1) * cell_width
        cell_right = cell_left + cell_width
        cell_center = _midpoint(cell_left, cell_right)

        base_span = base_spans_by_index.get(base_call.base_index)
        quality_segment = quality_segments_by_index.get(base_call.base_index)
        raw_left = None if base_span is None else base_span.left
        raw_right = None if base_span is None else base_span.right
        raw_center = (
            float(base_call.position) if base_span is None else float(base_span.center)
        )

        columns.append(
            ChromatogramColumn(
                column_index=column_index,
                base_index=base_call.base_index,
                query_pos=base_call.base_index + 1,
                base=base_call.base,
                color=base_call.color,
                quality=(
                    None if quality_segment is None else int(quality_segment.quality)
                ),
                trace_x=base_call.position,
                peak_height=_peak_height_for_base_call(raw_view, base_call),
                is_retained=_is_base_index_retained(
                    base_index=base_call.base_index,
                    trim_result=trim_result,
                    display_orientation=record.orientation,
                ),
                cell_left=cell_left,
                cell_right=cell_right,
                cell_center=cell_center,
                raw_left=raw_left,
                raw_right=raw_right,
                raw_center=raw_center,
                channels=(
                    ()
                    if base_span is None
                    else _build_column_channels(
                        channels=raw_view.channels,
                        raw_left=base_span.left,
                        raw_right=base_span.right,
                        cell_left=cell_left,
                        cell_right=cell_right,
                        samples_per_base=samples_per_base,
                    )
                ),
            )
        )

    return ChromatogramColumnView(
        is_renderable=True,
        columns=tuple(columns),
        trim_boundaries=_resolve_column_trim_boundaries(
            trim_result=trim_result,
            base_count=len(columns),
            cell_width=cell_width,
            display_orientation=record.orientation,
        ),
        cell_width=cell_width,
        samples_per_base=samples_per_base,
    )


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


def _peak_height_for_base_call(
    view: ChromatogramView,
    base_call: ChromatogramBaseCall,
) -> int | None:
    for channel in view.channels:
        if channel.base == base_call.base and base_call.position < len(channel.signal):
            return channel.signal[base_call.position]

    fallback_heights = [
        channel.signal[base_call.position]
        for channel in view.channels
        if base_call.position < len(channel.signal)
    ]
    if not fallback_heights:
        return None
    return max(fallback_heights)


def _build_column_channels(
    *,
    channels: tuple[ChromatogramChannel, ...],
    raw_left: float,
    raw_right: float,
    cell_left: float,
    cell_right: float,
    samples_per_base: int,
) -> tuple[ChromatogramColumnChannel, ...]:
    if not channels or raw_right <= raw_left:
        return ()

    x_values = _linspace(cell_left, cell_right, samples_per_base)
    raw_x_values = _linspace(raw_left, raw_right, samples_per_base)
    return tuple(
        ChromatogramColumnChannel(
            base=channel.base,
            color=channel.color,
            x_values=x_values,
            signal=tuple(
                _interpolate_signal(channel.signal, raw_x_value)
                for raw_x_value in raw_x_values
            ),
        )
        for channel in channels
    )


def _interpolate_signal(signal: tuple[int, ...], x_value: float) -> float:
    if not signal:
        return 0.0
    if x_value <= 0:
        return float(signal[0])

    max_index = len(signal) - 1
    if x_value >= max_index:
        return float(signal[max_index])

    left_index = int(math.floor(x_value))
    right_index = int(math.ceil(x_value))
    if left_index == right_index:
        return float(signal[left_index])

    right_weight = x_value - float(left_index)
    left_weight = 1.0 - right_weight
    return (
        float(signal[left_index]) * left_weight
        + float(signal[right_index]) * right_weight
    )


def _linspace(start: float, end: float, count: int) -> tuple[float, ...]:
    if count <= 1:
        return ((_midpoint(start, end)),)
    step = (end - start) / float(count - 1)
    return tuple(start + (step * index) for index in range(count))


def _display_trim_start(
    trim_result: TrimResult,
    *,
    display_orientation: SequenceOrientation,
) -> int:
    if display_orientation == "forward":
        return trim_result.bases_removed_left
    return trim_result.bases_removed_right


def _is_base_index_retained(
    *,
    base_index: int,
    trim_result: TrimResult | None,
    display_orientation: SequenceOrientation,
) -> bool:
    if trim_result is None:
        return True
    if trim_result.trimmed_length <= 0:
        return False

    display_trim_start = _display_trim_start(
        trim_result,
        display_orientation=display_orientation,
    )
    display_trim_end = display_trim_start + trim_result.trimmed_length
    return display_trim_start <= base_index < display_trim_end


def _resolve_column_trim_boundaries(
    *,
    trim_result: TrimResult | None,
    base_count: int,
    cell_width: float,
    display_orientation: SequenceOrientation,
) -> ChromatogramTrimBoundaries:
    if trim_result is None or base_count <= 0:
        return ChromatogramTrimBoundaries()

    display_trim_start = min(
        max(
            _display_trim_start(trim_result, display_orientation=display_orientation), 0
        ),
        base_count,
    )
    display_trim_end = min(
        max(display_trim_start + trim_result.trimmed_length, display_trim_start),
        base_count,
    )
    return ChromatogramTrimBoundaries(
        left=(
            None if display_trim_start <= 0 else float(display_trim_start) * cell_width
        ),
        right=(
            None
            if display_trim_end >= base_count
            else float(display_trim_end) * cell_width
        ),
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
