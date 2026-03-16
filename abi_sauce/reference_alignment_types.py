from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeAlias

from abi_sauce.alignment_policy import AlignmentStrand, AlignmentStrandPolicy

StrandPolicy: TypeAlias = AlignmentStrandPolicy
ChosenStrand: TypeAlias = AlignmentStrand
ReferenceConsensusResolution: TypeAlias = Literal[
    "concordant",
    "single_read",
    "quality_resolved",
    "majority_resolved",
    "ambiguous",
    "deleted",
]
ReferenceMultiAnchorKind: TypeAlias = Literal["reference", "insertion"]


@dataclass(frozen=True, slots=True)
class AlignmentEvent:
    ref_pos: int | None
    query_pos: int | None
    event_type: str
    ref_base: str
    query_base: str
    qscore: int | None
    flank_q_left: int | None
    flank_q_right: int | None
    trace_x: int | None
    context_ref: str
    context_query: str
    column_index: int | None = None


@dataclass(frozen=True, slots=True)
class ReferenceAlignmentColumn:
    column_index: int
    ref_base: str
    query_base: str
    event_type: str
    ref_index: int | None
    query_index: int | None
    ref_pos: int | None
    query_pos: int | None
    qscore: int | None
    trace_x: int | None
    is_match: bool


@dataclass(frozen=True, slots=True)
class AlignmentResult:
    sample_name: str
    reference_name: str
    strand: ChosenStrand
    score: float
    reference_start: int | None
    reference_end: int | None
    query_start: int | None
    query_end: int | None
    percent_identity: float
    mismatch_count: int
    insertion_count: int
    deletion_count: int
    aligned_reference: str
    match_line: str
    aligned_query: str
    events: tuple[AlignmentEvent, ...]
    columns: tuple[ReferenceAlignmentColumn, ...] = ()


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentMember:
    member_index: int
    source_filename: str
    display_name: str
    chosen_strand: ChosenStrand
    included: bool
    inclusion_reason: str | None
    trimmed_length: int
    has_trace_data: bool
    has_qualities: bool
    alignment_score: float | None = None
    overlap_length: int | None = None
    percent_identity: float | None = None
    mismatch_count: int | None = None
    insertion_count: int | None = None
    deletion_count: int | None = None


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentMemberCell:
    member_index: int
    base: str
    query_index: int | None
    query_pos: int | None
    qscore: int | None
    trace_x: int | None
    is_gap: bool


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentColumn:
    column_index: int
    anchor_kind: ReferenceMultiAnchorKind
    anchor_index: int
    ref_index: int | None
    ref_pos: int | None
    ref_base: str
    consensus_base: str
    resolution: ReferenceConsensusResolution
    member_cells: tuple[ReferenceMultiAlignmentMemberCell, ...]
    support_counts: tuple[tuple[str, int], ...] = ()
    quality_sums: tuple[tuple[str, int], ...] = ()
    non_gap_member_count: int = 0
    gap_member_count: int = 0
    ambiguous: bool = False
    matches_reference: bool | None = None


@dataclass(frozen=True, slots=True)
class ReferenceMultiAlignmentResult:
    reference_name: str
    reference_sequence: str
    accepted: bool
    rejection_reason: str | None
    members: tuple[ReferenceMultiAlignmentMember, ...]
    included_member_indices: tuple[int, ...]
    aligned_member_sequences: tuple[str, ...]
    gapped_reference: str
    gapped_consensus: str
    consensus_sequence: str
    columns: tuple[ReferenceMultiAlignmentColumn, ...] = ()
    included_member_count: int = 0
    excluded_member_count: int = 0
    ambiguous_column_count: int = 0


__all__ = [
    "StrandPolicy",
    "ChosenStrand",
    "ReferenceConsensusResolution",
    "ReferenceMultiAnchorKind",
    "AlignmentEvent",
    "ReferenceAlignmentColumn",
    "AlignmentResult",
    "ReferenceMultiAlignmentMember",
    "ReferenceMultiAlignmentMemberCell",
    "ReferenceMultiAlignmentColumn",
    "ReferenceMultiAlignmentResult",
]
