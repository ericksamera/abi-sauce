from __future__ import annotations

from dataclasses import dataclass
import math

from abi_sauce.assembly import (
    AssemblyColumn,
    AssemblyResult,
    AssemblyStrand,
    ConflictResolution,
    MultiAssemblyColumn,
    MultiAssemblyResult,
)
from abi_sauce.chromatogram import (
    ChromatogramBaseCall,
    ChromatogramBaseSpan,
    ChromatogramChannel,
    ChromatogramView,
    build_chromatogram_view,
    reverse_complement_chromatogram_view,
)
from abi_sauce.models import SequenceOrientation, SequenceRecord
from abi_sauce.trimming import TrimResult


@dataclass(frozen=True, slots=True)
class AssemblyTraceChannelSegment:
    """One resampled trace-channel segment for a single alignment cell."""

    base: str
    color: str
    x_values: tuple[float, ...] = ()
    normalized_signal: tuple[float, ...] = ()


@dataclass(frozen=True, slots=True)
class AssemblyTraceCell:
    """One alignment-column cell carrying normalized electropherogram state."""

    column_index: int
    base: str
    consensus_base: str
    resolution: ConflictResolution
    is_gap: bool
    is_overlap: bool
    is_match: bool
    query_index: int | None
    query_pos: int | None
    quality: int | None
    trace_x: int | None
    cell_left: float
    cell_right: float
    cell_center: float
    raw_left: float | None = None
    raw_right: float | None = None
    raw_center: float | None = None
    channels: tuple[AssemblyTraceChannelSegment, ...] = ()

    @property
    def has_trace_signal(self) -> bool:
        """Return whether the cell carries any resampled trace segments."""
        return bool(self.channels)


@dataclass(frozen=True, slots=True)
class AssemblyTraceColumn:
    """One generic alignment column summary for the stacked trace view."""

    column_index: int
    consensus_base: str
    resolution: ConflictResolution
    hover_text: str


@dataclass(frozen=True, slots=True)
class AssemblyTraceRow:
    """One stacked aligned-trace row for an assembly member."""

    label: str
    source_filename: str
    display_name: str
    strand: AssemblyStrand
    y_bottom: float
    y_top: float
    signal_scale: float
    has_trace_signal: bool
    cells: tuple[AssemblyTraceCell, ...] = ()

    @property
    def aligned_sequence(self) -> str:
        """Return the row's gapped aligned sequence."""
        return "".join(cell.base for cell in self.cells)


@dataclass(frozen=True, slots=True)
class AssemblyTraceView:
    """Pure aligned-trace view state for pairwise or multi-read assembly."""

    columns: tuple[AssemblyTraceColumn, ...] = ()
    rows: tuple[AssemblyTraceRow, ...] = ()
    cell_width: float = 1.0
    samples_per_cell: int = 16
    trace_row_height: float = 3.0

    @property
    def alignment_length(self) -> int:
        """Return the number of alignment columns represented in the view."""
        return len(self.columns)

    @property
    def total_height(self) -> float:
        """Return the full y-extent occupied by all trace rows."""
        return float(len(self.rows)) * self.trace_row_height

    @property
    def x_range(self) -> tuple[float, float]:
        """Return the full x-extent in alignment-column space."""
        return (0.0, float(self.alignment_length) * self.cell_width)


@dataclass(frozen=True, slots=True)
class _AssemblyTraceRowSource:
    view: ChromatogramView
    trimmed_start: int
    trimmed_length: int
    channels: tuple[ChromatogramChannel, ...] = ()
    base_calls_by_full_index: dict[int, ChromatogramBaseCall] | None = None
    base_spans_by_full_index: dict[int, ChromatogramBaseSpan] | None = None
    signal_scale: float = 1.0


@dataclass(frozen=True, slots=True)
class _AssemblyColumnProjection:
    base: str
    query_index: int | None
    query_pos: int | None
    quality: int | None
    trace_x: int | None
    is_overlap: bool
    is_match: bool


def build_pairwise_assembly_trace_view(
    *,
    result: AssemblyResult,
    left_source_filename: str,
    left_raw_record: SequenceRecord,
    left_trim_result: TrimResult,
    right_source_filename: str,
    right_raw_record: SequenceRecord,
    right_trim_result: TrimResult,
    cell_width: float = 1.0,
    samples_per_cell: int = 16,
    trace_row_height: float = 3.0,
) -> AssemblyTraceView:
    """Build a stacked alignment-column trace view for one pairwise assembly."""
    _validate_trace_view_parameters(
        cell_width=cell_width,
        samples_per_cell=samples_per_cell,
        trace_row_height=trace_row_height,
    )

    columns = _pairwise_trace_columns(result.columns)
    left_row_source = _build_row_source(
        raw_record=left_raw_record,
        trim_result=left_trim_result,
        strand="forward",
    )
    right_row_source = _build_row_source(
        raw_record=right_raw_record,
        trim_result=right_trim_result,
        strand=result.chosen_right_orientation,
    )

    row_specs: tuple[
        tuple[
            str,
            str,
            str,
            AssemblyStrand,
            _AssemblyTraceRowSource,
            tuple[_AssemblyColumnProjection, ...],
        ],
        ...,
    ] = (
        (
            result.left_display_name,
            left_source_filename,
            result.left_display_name,
            "forward",
            left_row_source,
            tuple(_left_column_projection(column) for column in result.columns),
        ),
        (
            result.right_display_name,
            right_source_filename,
            result.right_display_name,
            result.chosen_right_orientation,
            right_row_source,
            tuple(_right_column_projection(column) for column in result.columns),
        ),
    )

    return _build_trace_view(
        columns=columns,
        row_specs=row_specs,
        cell_width=cell_width,
        samples_per_cell=samples_per_cell,
        trace_row_height=trace_row_height,
    )


def build_multi_assembly_trace_view(
    *,
    result: MultiAssemblyResult,
    raw_records_by_source_filename: dict[str, SequenceRecord],
    trim_results_by_source_filename: dict[str, TrimResult],
    cell_width: float = 1.0,
    samples_per_cell: int = 16,
    trace_row_height: float = 3.0,
) -> AssemblyTraceView:
    """Build a stacked alignment-column trace view for one multi-read assembly."""
    _validate_trace_view_parameters(
        cell_width=cell_width,
        samples_per_cell=samples_per_cell,
        trace_row_height=trace_row_height,
    )

    columns = _multi_trace_columns(result.columns)
    members_by_index = {member.member_index: member for member in result.members}
    row_specs: list[
        tuple[
            str,
            str,
            str,
            AssemblyStrand,
            _AssemblyTraceRowSource,
            tuple[_AssemblyColumnProjection, ...],
        ]
    ] = []
    member_cell_position_by_index = {
        member_index: position
        for position, member_index in enumerate(result.included_member_indices)
    }

    for member_index in result.included_member_indices:
        member = members_by_index[member_index]
        row_source = _build_row_source(
            raw_record=raw_records_by_source_filename[member.source_filename],
            trim_result=trim_results_by_source_filename[member.source_filename],
            strand=member.chosen_orientation,
        )
        row_specs.append(
            (
                member.display_name,
                member.source_filename,
                member.display_name,
                member.chosen_orientation,
                row_source,
                tuple(
                    _multi_member_projection(
                        column,
                        member_cell_position=member_cell_position_by_index[
                            member_index
                        ],
                    )
                    for column in result.columns
                ),
            )
        )

    return _build_trace_view(
        columns=columns,
        row_specs=tuple(row_specs),
        cell_width=cell_width,
        samples_per_cell=samples_per_cell,
        trace_row_height=trace_row_height,
    )


def _validate_trace_view_parameters(
    *,
    cell_width: float,
    samples_per_cell: int,
    trace_row_height: float,
) -> None:
    if cell_width <= 0:
        raise ValueError("cell_width must be > 0")
    if samples_per_cell < 2:
        raise ValueError("samples_per_cell must be >= 2")
    if trace_row_height <= 0:
        raise ValueError("trace_row_height must be > 0")


def _build_trace_view(
    *,
    columns: tuple[AssemblyTraceColumn, ...],
    row_specs: tuple[
        tuple[
            str,
            str,
            str,
            AssemblyStrand,
            _AssemblyTraceRowSource,
            tuple[_AssemblyColumnProjection, ...],
        ],
        ...,
    ],
    cell_width: float,
    samples_per_cell: int,
    trace_row_height: float,
) -> AssemblyTraceView:
    total_rows = len(row_specs)
    rows = tuple(
        _build_trace_row(
            columns=columns,
            label=label,
            source_filename=source_filename,
            display_name=display_name,
            strand=strand,
            row_source=row_source,
            projections=projections,
            row_order_index=row_order_index,
            total_rows=total_rows,
            cell_width=cell_width,
            samples_per_cell=samples_per_cell,
            trace_row_height=trace_row_height,
        )
        for row_order_index, (
            label,
            source_filename,
            display_name,
            strand,
            row_source,
            projections,
        ) in enumerate(row_specs)
    )

    return AssemblyTraceView(
        columns=columns,
        rows=rows,
        cell_width=cell_width,
        samples_per_cell=samples_per_cell,
        trace_row_height=trace_row_height,
    )


def _build_trace_row(
    *,
    columns: tuple[AssemblyTraceColumn, ...],
    label: str,
    source_filename: str,
    display_name: str,
    strand: AssemblyStrand,
    row_source: _AssemblyTraceRowSource,
    projections: tuple[_AssemblyColumnProjection, ...],
    row_order_index: int,
    total_rows: int,
    cell_width: float,
    samples_per_cell: int,
    trace_row_height: float,
) -> AssemblyTraceRow:
    y_bottom = float(total_rows - row_order_index - 1) * trace_row_height
    y_top = y_bottom + trace_row_height
    cells = tuple(
        _build_trace_cell(
            column=column,
            row_source=row_source,
            projection=projection,
            cell_width=cell_width,
            samples_per_cell=samples_per_cell,
        )
        for column, projection in zip(columns, projections, strict=True)
    )
    return AssemblyTraceRow(
        label=label,
        source_filename=source_filename,
        display_name=display_name,
        strand=strand,
        y_bottom=y_bottom,
        y_top=y_top,
        signal_scale=row_source.signal_scale,
        has_trace_signal=any(cell.has_trace_signal for cell in cells),
        cells=cells,
    )


def _build_trace_cell(
    *,
    column: AssemblyTraceColumn,
    row_source: _AssemblyTraceRowSource,
    projection: _AssemblyColumnProjection,
    cell_width: float,
    samples_per_cell: int,
) -> AssemblyTraceCell:
    cell_left = float(column.column_index - 1) * cell_width
    cell_right = cell_left + cell_width
    cell_center = (cell_left + cell_right) / 2.0

    raw_left: float | None = None
    raw_right: float | None = None
    raw_center: float | None = None
    channels: tuple[AssemblyTraceChannelSegment, ...] = ()

    if projection.query_index is not None and projection.base != "-":
        full_base_index = row_source.trimmed_start + projection.query_index
        base_call = _mapping_lookup(
            row_source.base_calls_by_full_index, full_base_index
        )
        base_span = _mapping_lookup(
            row_source.base_spans_by_full_index, full_base_index
        )
        if base_call is not None:
            raw_center = float(base_call.position)
        if base_span is not None:
            raw_left = base_span.left
            raw_right = base_span.right
            raw_center = base_span.center if raw_center is None else raw_center
            channels = _build_cell_channels(
                channels=row_source.channels,
                raw_left=raw_left,
                raw_right=raw_right,
                cell_left=cell_left,
                cell_right=cell_right,
                samples_per_cell=samples_per_cell,
                signal_scale=row_source.signal_scale,
            )

    return AssemblyTraceCell(
        column_index=column.column_index,
        base=projection.base,
        consensus_base=column.consensus_base,
        resolution=column.resolution,
        is_gap=projection.query_index is None or projection.base == "-",
        is_overlap=projection.is_overlap,
        is_match=projection.is_match,
        query_index=projection.query_index,
        query_pos=projection.query_pos,
        quality=projection.quality,
        trace_x=projection.trace_x,
        cell_left=cell_left,
        cell_right=cell_right,
        cell_center=cell_center,
        raw_left=raw_left,
        raw_right=raw_right,
        raw_center=raw_center,
        channels=channels,
    )


def _build_row_source(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    strand: AssemblyStrand,
) -> _AssemblyTraceRowSource:
    view = build_chromatogram_view(raw_record, trim_result)
    oriented_view = (
        view if strand == "forward" else reverse_complement_chromatogram_view(view)
    )
    trimmed_start, _trimmed_end = _assembly_oriented_trim_interval(
        raw_record=raw_record,
        trim_result=trim_result,
        strand=strand,
    )
    base_calls_by_full_index = {
        base_call.base_index: base_call for base_call in oriented_view.base_calls
    }
    base_spans_by_full_index = {
        base_span.base_index: base_span for base_span in oriented_view.base_spans
    }
    signal_scale = _row_signal_scale(
        view=oriented_view,
        base_spans_by_full_index=base_spans_by_full_index,
        trimmed_start=trimmed_start,
        trimmed_length=trim_result.trimmed_length,
    )
    return _AssemblyTraceRowSource(
        view=oriented_view,
        trimmed_start=trimmed_start,
        trimmed_length=trim_result.trimmed_length,
        channels=oriented_view.channels,
        base_calls_by_full_index=base_calls_by_full_index,
        base_spans_by_full_index=base_spans_by_full_index,
        signal_scale=signal_scale,
    )


def _assembly_oriented_trim_interval(
    *,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    strand: AssemblyStrand,
) -> tuple[int, int]:
    display_start = _display_trim_start(
        trim_result,
        display_orientation=raw_record.orientation,
    )
    display_end = display_start + trim_result.trimmed_length
    if strand == "forward":
        return (display_start, display_end)

    full_length = len(raw_record.sequence)
    return (
        full_length - display_end,
        full_length - display_start,
    )


def _display_trim_start(
    trim_result: TrimResult,
    *,
    display_orientation: SequenceOrientation,
) -> int:
    if display_orientation == "forward":
        return trim_result.bases_removed_left
    return trim_result.bases_removed_right


def _row_signal_scale(
    *,
    view: ChromatogramView,
    base_spans_by_full_index: dict[int, ChromatogramBaseSpan],
    trimmed_start: int,
    trimmed_length: int,
) -> float:
    if not view.is_renderable or not view.channels or trimmed_length <= 0:
        return 1.0

    trimmed_spans = [
        base_spans_by_full_index[full_index]
        for full_index in range(trimmed_start, trimmed_start + trimmed_length)
        if full_index in base_spans_by_full_index
    ]
    if not trimmed_spans:
        return 1.0

    left = max(0, int(math.floor(min(span.left for span in trimmed_spans))))
    right = min(
        max(view.trace_length - 1, 0),
        int(math.ceil(max(span.right for span in trimmed_spans))),
    )
    if right < left:
        return 1.0

    max_signal = max(
        channel.signal[sample_index]
        for channel in view.channels
        for sample_index in range(left, right + 1)
    )
    return float(max(max_signal, 1))


def _build_cell_channels(
    *,
    channels: tuple[ChromatogramChannel, ...],
    raw_left: float,
    raw_right: float,
    cell_left: float,
    cell_right: float,
    samples_per_cell: int,
    signal_scale: float,
) -> tuple[AssemblyTraceChannelSegment, ...]:
    if not channels or raw_right <= raw_left:
        return ()

    x_values = _linspace(cell_left, cell_right, samples_per_cell)
    raw_x_values = _linspace(raw_left, raw_right, samples_per_cell)
    scale = signal_scale if signal_scale > 0 else 1.0

    return tuple(
        AssemblyTraceChannelSegment(
            base=channel.base,
            color=channel.color,
            x_values=x_values,
            normalized_signal=tuple(
                _clamp_unit(_interpolate_signal(channel.signal, raw_x) / scale)
                for raw_x in raw_x_values
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
        return ((start + end) / 2.0,)
    step = (end - start) / float(count - 1)
    return tuple(start + (step * index) for index in range(count))


def _clamp_unit(value: float) -> float:
    return max(0.0, min(value, 1.0))


def _pairwise_trace_columns(
    columns: tuple[AssemblyColumn, ...],
) -> tuple[AssemblyTraceColumn, ...]:
    return tuple(
        AssemblyTraceColumn(
            column_index=column.column_index,
            consensus_base=column.consensus_base,
            resolution=column.resolution,
            hover_text=(
                f"column={column.column_index}"
                f"<br>left={column.left_base}"
                f"<br>right={column.right_base}"
                f"<br>consensus={column.consensus_base}"
                f"<br>resolution={column.resolution}"
            ),
        )
        for column in columns
    )


def _multi_trace_columns(
    columns: tuple[MultiAssemblyColumn, ...],
) -> tuple[AssemblyTraceColumn, ...]:
    return tuple(
        AssemblyTraceColumn(
            column_index=column.column_index,
            consensus_base=column.consensus_base,
            resolution=column.resolution,
            hover_text=(
                f"column={column.column_index}"
                f"<br>consensus={column.consensus_base}"
                f"<br>resolution={column.resolution}"
                f"<br>support={_multi_support_summary(column)}"
                f"<br>non_gap_members={column.non_gap_member_count}"
                f"<br>gap_members={column.gap_member_count}"
            ),
        )
        for column in columns
    )


def _multi_support_summary(column: MultiAssemblyColumn) -> str:
    if not column.support_counts:
        return "NA"
    return ", ".join(f"{base}:{count}" for base, count in column.support_counts)


def _multi_member_projection(
    column: MultiAssemblyColumn,
    *,
    member_cell_position: int,
) -> _AssemblyColumnProjection:
    member_cell = column.member_cells[member_cell_position]
    return _AssemblyColumnProjection(
        base=member_cell.base,
        query_index=member_cell.query_index,
        query_pos=member_cell.query_pos,
        quality=member_cell.quality,
        trace_x=member_cell.trace_x,
        is_overlap=column.non_gap_member_count > 1,
        is_match=(
            not member_cell.is_gap
            and member_cell.base == column.consensus_base
            and column.consensus_base != "N"
        ),
    )


def _left_column_projection(column: AssemblyColumn) -> _AssemblyColumnProjection:
    return _AssemblyColumnProjection(
        base=column.left_base,
        query_index=column.left_query_index,
        query_pos=column.left_query_pos,
        quality=column.left_quality,
        trace_x=column.left_trace_x,
        is_overlap=column.is_overlap,
        is_match=column.is_match,
    )


def _right_column_projection(column: AssemblyColumn) -> _AssemblyColumnProjection:
    return _AssemblyColumnProjection(
        base=column.right_base,
        query_index=column.right_query_index,
        query_pos=column.right_query_pos,
        quality=column.right_quality,
        trace_x=column.right_trace_x,
        is_overlap=column.is_overlap,
        is_match=column.is_match,
    )


def _mapping_lookup(mapping, key: int):
    if mapping is None:
        return None
    return mapping.get(key)
