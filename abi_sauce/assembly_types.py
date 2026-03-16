from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeAlias

from abi_sauce.alignment_policy import AlignmentStrand

AssemblyStrand: TypeAlias = AlignmentStrand
ConflictResolution: TypeAlias = Literal[
    "concordant",
    "single_read",
    "quality_resolved",
    "majority_resolved",
    "ambiguous",
]


@dataclass(frozen=True, slots=True)
class AssemblyConfig:
    """Scoring and acceptance thresholds for pairwise read assembly."""

    match_score: float = 1.0
    mismatch_score: float = -1.0
    open_internal_gap_score: float = -3.0
    extend_internal_gap_score: float = -1.0
    min_overlap_length: int = 25
    min_percent_identity: float = 70.0
    quality_margin: int = 3


@dataclass(frozen=True, slots=True)
class AssemblyConflict:
    """One consensus-support column worth surfacing in the UI."""

    column_index: int
    left_base: str
    right_base: str
    consensus_base: str
    resolution: ConflictResolution
    left_query_pos: int | None
    right_query_pos: int | None
    left_quality: int | None
    right_quality: int | None
    left_trace_x: int | None
    right_trace_x: int | None


@dataclass(frozen=True, slots=True)
class AssemblyColumn:
    """One gapped alignment column carrying enough state for downstream views."""

    column_index: int
    left_base: str
    right_base: str
    consensus_base: str
    resolution: ConflictResolution
    left_query_index: int | None
    right_query_index: int | None
    left_query_pos: int | None
    right_query_pos: int | None
    left_quality: int | None
    right_quality: int | None
    left_trace_x: int | None
    right_trace_x: int | None
    is_overlap: bool
    is_match: bool


@dataclass(frozen=True, slots=True)
class AssemblyResult:
    """Derived pairwise assembly result for two trimmed reads."""

    left_source_filename: str
    right_source_filename: str
    left_display_name: str
    right_display_name: str
    chosen_right_orientation: AssemblyStrand
    score: float
    accepted: bool
    rejection_reason: str | None
    overlap_length: int
    percent_identity: float
    mismatch_count: int
    aligned_left: str
    match_line: str
    aligned_right: str
    gapped_consensus: str
    consensus_sequence: str
    columns: tuple[AssemblyColumn, ...] = ()
    conflicts: tuple[AssemblyConflict, ...] = ()

    @property
    def conflict_count(self) -> int:
        return len(self.conflicts)


@dataclass(frozen=True, slots=True)
class MultiAssemblyMember:
    """One member selected for a multi-read assembly attempt."""

    member_index: int
    source_filename: str
    display_name: str
    chosen_orientation: AssemblyStrand
    included: bool
    inclusion_reason: str | None
    trimmed_length: int
    has_trace_data: bool
    has_qualities: bool
    is_seed: bool = False
    alignment_score: float | None = None
    overlap_length: int | None = None
    percent_identity: float | None = None


@dataclass(frozen=True, slots=True)
class MultiAssemblyMemberCell:
    """One included-member cell projected onto the multi-read alignment grid."""

    member_index: int
    base: str
    query_index: int | None
    query_pos: int | None
    quality: int | None
    trace_x: int | None
    is_gap: bool


@dataclass(frozen=True, slots=True)
class MultiAssemblyColumn:
    """One global alignment column for a multi-read assembly."""

    column_index: int
    consensus_base: str
    resolution: ConflictResolution
    member_cells: tuple[MultiAssemblyMemberCell, ...]
    support_counts: tuple[tuple[str, int], ...] = ()
    quality_sums: tuple[tuple[str, int], ...] = ()
    non_gap_member_count: int = 0
    gap_member_count: int = 0
    ambiguous: bool = False


@dataclass(frozen=True, slots=True)
class MultiAssemblyResult:
    """Derived anchored multi-read assembly result."""

    members: tuple[MultiAssemblyMember, ...]
    seed_member_index: int
    accepted: bool
    rejection_reason: str | None
    included_member_indices: tuple[int, ...]
    aligned_member_sequences: tuple[str, ...]
    gapped_consensus: str
    consensus_sequence: str
    columns: tuple[MultiAssemblyColumn, ...] = ()
    included_member_count: int = 0
    excluded_member_count: int = 0
    ambiguous_column_count: int = 0

    @property
    def conflict_count(self) -> int:
        return self.ambiguous_column_count


AssemblyComputationResult: TypeAlias = AssemblyResult | MultiAssemblyResult
