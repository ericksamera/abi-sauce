from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

from abi_sauce.alignment_policy import AlignmentStrand, AlignmentStrandPolicy

StrandPolicy: TypeAlias = AlignmentStrandPolicy
ChosenStrand: TypeAlias = AlignmentStrand


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


__all__ = [
    "StrandPolicy",
    "ChosenStrand",
    "AlignmentEvent",
    "AlignmentResult",
]
