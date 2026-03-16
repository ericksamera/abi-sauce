from __future__ import annotations

from dataclasses import dataclass

from abi_sauce.assembly_trace import (
    AssemblyTraceRowSource,
    build_assembly_trace_row_source,
    resolve_assembly_trace_samples_per_cell,
)
from abi_sauce.models import SequenceRecord
from abi_sauce.reference_alignment_types import (
    ReferenceMultiAlignmentColumn,
    ReferenceMultiAlignmentResult,
)
from abi_sauce.signal_sampling import resample_signal_window
from abi_sauce.trimming import TrimResult


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentTraceChannelSegment:
    """One resampled trace-channel segment for one shared-reference cell."""

    base: str
    color: str
    x_values: tuple[float, ...] = ()
    normalized_signal: tuple[float, ...] = ()


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentTraceCell:
    """One shared-reference alignment cell carrying row-specific trace state."""

    column_index: int
    anchor_kind: str
    ref_base: str
    query_base: str
    consensus_base: str
    resolution: str
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
    channels: tuple[ReferenceMultiAlignmentTraceChannelSegment, ...] = ()

    @property
    def has_trace_signal(self) -> bool:
        return bool(self.channels)


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentTraceRow:
    """One aligned read row in a shared-reference electropherogram view."""

    label: str
    source_filename: str
    display_name: str
    strand: str
    y_bottom: float
    y_top: float
    signal_scale: float
    has_trace_signal: bool
    cells: tuple[ReferenceMultiAlignmentTraceCell, ...] = ()

    @property
    def aligned_sequence(self) -> str:
        return "".join(cell.query_base for cell in self.cells)


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentTraceView:
    """Pure stacked shared-reference aligned electropherogram view state."""

    reference_name: str
    columns: tuple[ReferenceMultiAlignmentColumn, ...] = ()
    rows: tuple[ReferenceMultiAlignmentTraceRow, ...] = ()
    cell_width: float = 1.0
    samples_per_cell: int = 16
    trace_row_height: float = 3.0

    @property
    def alignment_length(self) -> int:
        return len(self.columns)

    @property
    def total_height(self) -> float:
        return float(len(self.rows)) * self.trace_row_height

    @property
    def x_range(self) -> tuple[float, float]:
        return (0.0, float(self.alignment_length) * self.cell_width)


def build_reference_multi_alignment_trace_view(
    *,
    result: ReferenceMultiAlignmentResult,
    raw_records_by_source_filename: dict[str, SequenceRecord],
    trim_results_by_source_filename: dict[str, TrimResult],
    cell_width: float = 1.0,
    samples_per_cell: int | None = None,
    trace_row_height: float = 3.0,
    row_sources_by_member_index: dict[int, AssemblyTraceRowSource] | None = None,
) -> ReferenceMultiAlignmentTraceView:
    """Build one stacked shared-reference electropherogram view."""
    resolved_samples_per_cell = resolve_assembly_trace_samples_per_cell(
        alignment_length=len(result.columns),
        row_count=len(result.included_member_indices),
        requested_samples_per_cell=samples_per_cell,
    )
    _validate_trace_view_parameters(
        cell_width=cell_width,
        samples_per_cell=resolved_samples_per_cell,
        trace_row_height=trace_row_height,
    )

    members_by_index = {member.member_index: member for member in result.members}
    member_cell_position_by_index = {
        member_index: position
        for position, member_index in enumerate(result.included_member_indices)
    }

    rows: list[ReferenceMultiAlignmentTraceRow] = []
    total_rows = len(result.included_member_indices)
    for row_order_index, member_index in enumerate(result.included_member_indices):
        member = members_by_index[member_index]
        row_source = (
            row_sources_by_member_index[member_index]
            if row_sources_by_member_index is not None
            and member_index in row_sources_by_member_index
            else build_assembly_trace_row_source(
                raw_record=raw_records_by_source_filename[member.source_filename],
                trim_result=trim_results_by_source_filename[member.source_filename],
                strand=member.chosen_strand,
            )
        )
        rows.append(
            _build_trace_row(
                columns=result.columns,
                member_cell_position=member_cell_position_by_index[member_index],
                label=member.display_name,
                source_filename=member.source_filename,
                display_name=member.display_name,
                strand=member.chosen_strand,
                row_source=row_source,
                row_order_index=row_order_index,
                total_rows=total_rows,
                cell_width=cell_width,
                samples_per_cell=resolved_samples_per_cell,
                trace_row_height=trace_row_height,
            )
        )

    return ReferenceMultiAlignmentTraceView(
        reference_name=result.reference_name,
        columns=result.columns,
        rows=tuple(rows),
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
    columns: tuple[ReferenceMultiAlignmentColumn, ...],
    member_cell_position: int,
    label: str,
    source_filename: str,
    display_name: str,
    strand: str,
    row_source: AssemblyTraceRowSource,
    row_order_index: int,
    total_rows: int,
    cell_width: float,
    samples_per_cell: int,
    trace_row_height: float,
) -> ReferenceMultiAlignmentTraceRow:
    y_bottom = float(total_rows - row_order_index - 1) * trace_row_height
    y_top = y_bottom + trace_row_height
    cells = tuple(
        _build_trace_cell(
            column=column,
            member_cell_position=member_cell_position,
            row_source=row_source,
            cell_width=cell_width,
            samples_per_cell=samples_per_cell,
        )
        for column in columns
    )
    return ReferenceMultiAlignmentTraceRow(
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
    column: ReferenceMultiAlignmentColumn,
    member_cell_position: int,
    row_source: AssemblyTraceRowSource,
    cell_width: float,
    samples_per_cell: int,
) -> ReferenceMultiAlignmentTraceCell:
    member_cell = column.member_cells[member_cell_position]
    cell_left = float(column.column_index - 1) * cell_width
    cell_right = cell_left + cell_width
    cell_center = (cell_left + cell_right) / 2.0

    raw_left: float | None = None
    raw_right: float | None = None
    raw_center: float | None = None
    channels: tuple[ReferenceMultiAlignmentTraceChannelSegment, ...] = ()

    if member_cell.query_index is not None and member_cell.base != "-":
        full_base_index = row_source.trimmed_start + member_cell.query_index
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

    return ReferenceMultiAlignmentTraceCell(
        column_index=column.column_index,
        anchor_kind=column.anchor_kind,
        ref_base=column.ref_base,
        query_base=member_cell.base,
        consensus_base=column.consensus_base,
        resolution=column.resolution,
        is_gap=member_cell.is_gap,
        is_match=(
            not member_cell.is_gap
            and member_cell.base == column.consensus_base
            and column.consensus_base not in {"N", "-"}
        ),
        ref_pos=column.ref_pos,
        query_index=member_cell.query_index,
        query_pos=member_cell.query_pos,
        quality=member_cell.qscore,
        trace_x=member_cell.trace_x,
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
) -> tuple[ReferenceMultiAlignmentTraceChannelSegment, ...]:
    if not row_source.channels or raw_right <= raw_left:
        return ()

    channel_segments: list[ReferenceMultiAlignmentTraceChannelSegment] = []
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
            ReferenceMultiAlignmentTraceChannelSegment(
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
    "ReferenceMultiAlignmentTraceChannelSegment",
    "ReferenceMultiAlignmentTraceCell",
    "ReferenceMultiAlignmentTraceRow",
    "ReferenceMultiAlignmentTraceView",
    "build_reference_multi_alignment_trace_view",
]
