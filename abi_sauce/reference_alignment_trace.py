from __future__ import annotations

from dataclasses import dataclass

from abi_sauce.assembly_trace import (
    AssemblyTraceRowSource,
    build_assembly_trace_row_source,
    resolve_assembly_trace_samples_per_cell,
)
from abi_sauce.models import SequenceRecord
from abi_sauce.reference_alignment_types import (
    AlignmentResult,
    ChosenStrand,
    ReferenceAlignmentColumn,
)
from abi_sauce.signal_sampling import resample_signal_window
from abi_sauce.trimming import TrimResult


@dataclass(frozen=True, slots=True)
class ReferenceAlignmentTraceChannelSegment:
    """One resampled trace-channel segment for a reference-alignment cell."""

    base: str
    color: str
    x_values: tuple[float, ...] = ()
    normalized_signal: tuple[float, ...] = ()


@dataclass(frozen=True, slots=True)
class ReferenceAlignmentTraceCell:
    """One fixed-width alignment cell carrying query-trace state."""

    column_index: int
    ref_base: str
    query_base: str
    event_type: str
    is_gap: bool
    is_match: bool
    ref_pos: int | None
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
    channels: tuple[ReferenceAlignmentTraceChannelSegment, ...] = ()

    @property
    def has_trace_signal(self) -> bool:
        """Return whether the cell carries any resampled trace segments."""
        return bool(self.channels)


@dataclass(frozen=True, slots=True)
class ReferenceAlignmentTraceRow:
    """The aligned query row for the reference-comparison view."""

    label: str
    source_filename: str
    display_name: str
    strand: ChosenStrand
    y_bottom: float
    y_top: float
    signal_scale: float
    has_trace_signal: bool
    cells: tuple[ReferenceAlignmentTraceCell, ...] = ()

    @property
    def aligned_sequence(self) -> str:
        """Return the row's gapped aligned query sequence."""
        return "".join(cell.query_base for cell in self.cells)


@dataclass(frozen=True, slots=True)
class ReferenceAlignmentTraceView:
    """Pure aligned reference-vs-query electropherogram view state."""

    reference_name: str
    source_filename: str
    sample_name: str
    strand: ChosenStrand
    columns: tuple[ReferenceAlignmentColumn, ...] = ()
    rows: tuple[ReferenceAlignmentTraceRow, ...] = ()
    cell_width: float = 1.0
    samples_per_cell: int = 16
    trace_row_height: float = 3.0

    @property
    def alignment_length(self) -> int:
        """Return the number of alignment columns represented in the view."""
        return len(self.columns)

    @property
    def total_height(self) -> float:
        """Return the total y-extent occupied by all query rows."""
        return float(len(self.rows)) * self.trace_row_height

    @property
    def x_range(self) -> tuple[float, float]:
        """Return the full x-extent in alignment-column space."""
        return (0.0, float(self.alignment_length) * self.cell_width)


def build_reference_alignment_trace_view(
    *,
    result: AlignmentResult,
    source_filename: str,
    raw_record: SequenceRecord,
    trim_result: TrimResult,
    cell_width: float = 1.0,
    samples_per_cell: int | None = None,
    trace_row_height: float = 3.0,
    row_source: AssemblyTraceRowSource | None = None,
) -> ReferenceAlignmentTraceView:
    """Build one reference-vs-query aligned electropherogram view."""
    resolved_samples_per_cell = resolve_assembly_trace_samples_per_cell(
        alignment_length=len(result.columns),
        row_count=1,
        requested_samples_per_cell=samples_per_cell,
    )
    _validate_trace_view_parameters(
        cell_width=cell_width,
        samples_per_cell=resolved_samples_per_cell,
        trace_row_height=trace_row_height,
    )

    resolved_row_source = (
        build_assembly_trace_row_source(
            raw_record=raw_record,
            trim_result=trim_result,
            strand=result.strand,
        )
        if row_source is None
        else row_source
    )
    query_row = _build_trace_row(
        result=result,
        source_filename=source_filename,
        row_source=resolved_row_source,
        cell_width=cell_width,
        samples_per_cell=resolved_samples_per_cell,
        trace_row_height=trace_row_height,
    )
    return ReferenceAlignmentTraceView(
        reference_name=result.reference_name,
        source_filename=source_filename,
        sample_name=result.sample_name,
        strand=result.strand,
        columns=result.columns,
        rows=(query_row,),
        cell_width=cell_width,
        samples_per_cell=resolved_samples_per_cell,
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


def _build_trace_row(
    *,
    result: AlignmentResult,
    source_filename: str,
    row_source: AssemblyTraceRowSource,
    cell_width: float,
    samples_per_cell: int,
    trace_row_height: float,
) -> ReferenceAlignmentTraceRow:
    y_bottom = 0.0
    y_top = trace_row_height
    cells = tuple(
        _build_trace_cell(
            column=column,
            row_source=row_source,
            cell_width=cell_width,
            samples_per_cell=samples_per_cell,
        )
        for column in result.columns
    )
    return ReferenceAlignmentTraceRow(
        label=result.sample_name,
        source_filename=source_filename,
        display_name=result.sample_name,
        strand=result.strand,
        y_bottom=y_bottom,
        y_top=y_top,
        signal_scale=row_source.signal_scale,
        has_trace_signal=any(cell.has_trace_signal for cell in cells),
        cells=cells,
    )


def _build_trace_cell(
    *,
    column: ReferenceAlignmentColumn,
    row_source: AssemblyTraceRowSource,
    cell_width: float,
    samples_per_cell: int,
) -> ReferenceAlignmentTraceCell:
    cell_left = float(column.column_index - 1) * cell_width
    cell_right = cell_left + cell_width
    cell_center = (cell_left + cell_right) / 2.0

    raw_left: float | None = None
    raw_right: float | None = None
    raw_center: float | None = None
    channels: tuple[ReferenceAlignmentTraceChannelSegment, ...] = ()

    if column.query_index is not None and column.query_base != "-":
        full_base_index = row_source.trimmed_start + column.query_index
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
                row_source=row_source,
                raw_left=raw_left,
                raw_right=raw_right,
                cell_left=cell_left,
                cell_right=cell_right,
                samples_per_cell=samples_per_cell,
            )

    return ReferenceAlignmentTraceCell(
        column_index=column.column_index,
        ref_base=column.ref_base,
        query_base=column.query_base,
        event_type=column.event_type,
        is_gap=column.query_index is None or column.query_base == "-",
        is_match=column.is_match,
        ref_pos=column.ref_pos,
        query_index=column.query_index,
        query_pos=column.query_pos,
        quality=column.qscore,
        trace_x=column.trace_x,
        cell_left=cell_left,
        cell_right=cell_right,
        cell_center=cell_center,
        raw_left=raw_left,
        raw_right=raw_right,
        raw_center=raw_center,
        channels=channels,
    )


def _build_cell_channels(
    *,
    row_source: AssemblyTraceRowSource,
    raw_left: float,
    raw_right: float,
    cell_left: float,
    cell_right: float,
    samples_per_cell: int,
) -> tuple[ReferenceAlignmentTraceChannelSegment, ...]:
    if not row_source.channels or raw_right <= raw_left:
        return ()

    channel_segments: list[ReferenceAlignmentTraceChannelSegment] = []
    for channel in row_source.channels:
        x_values, sampled_signal = resample_signal_window(
            channel.signal,
            raw_left=raw_left,
            raw_right=raw_right,
            cell_left=cell_left,
            cell_right=cell_right,
            sample_count=samples_per_cell,
            signal_scale=row_source.signal_scale,
            clamp_to_unit=True,
        )
        channel_segments.append(
            ReferenceAlignmentTraceChannelSegment(
                base=channel.base,
                color=channel.color,
                x_values=x_values,
                normalized_signal=sampled_signal,
            )
        )
    return tuple(channel_segments)


def _mapping_lookup(mapping, key: int):
    if mapping is None:
        return None
    return mapping.get(key)


__all__ = [
    "ReferenceAlignmentTraceChannelSegment",
    "ReferenceAlignmentTraceCell",
    "ReferenceAlignmentTraceRow",
    "ReferenceAlignmentTraceView",
    "build_reference_alignment_trace_view",
]
