from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

from Bio import Align

from abi_sauce.orientation import reverse_complement_sequence

AlignmentStrand = Literal["forward", "reverse_complement"]
AlignmentStrandPolicy = Literal["auto", "forward", "reverse_complement"]


@dataclass(frozen=True, slots=True)
class SelectedOrientedAlignment:
    """One chosen alignment plus the oriented query sequence it used."""

    strand: AlignmentStrand
    sequence: str
    qualities: list[int] | None
    alignment: Any
    score: float


def build_semiglobal_aligner(
    *,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    open_internal_gap_score: float = -3.0,
    extend_internal_gap_score: float = -1.0,
) -> Align.PairwiseAligner:
    """Build a global aligner with free terminal gaps and penalized internal gaps."""
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score

    aligner.open_internal_insertion_score = open_internal_gap_score
    aligner.extend_internal_insertion_score = extend_internal_gap_score
    aligner.open_internal_deletion_score = open_internal_gap_score
    aligner.extend_internal_deletion_score = extend_internal_gap_score

    aligner.open_left_insertion_score = 0.0
    aligner.extend_left_insertion_score = 0.0
    aligner.open_right_insertion_score = 0.0
    aligner.extend_right_insertion_score = 0.0
    aligner.open_left_deletion_score = 0.0
    aligner.extend_left_deletion_score = 0.0
    aligner.open_right_deletion_score = 0.0
    aligner.extend_right_deletion_score = 0.0

    return aligner


def oriented_qualities(
    display_qualities: list[int] | None,
    *,
    strand: AlignmentStrand,
) -> list[int] | None:
    """Return display-oriented qualities in the requested strand orientation."""
    if display_qualities is None:
        return None
    if strand == "forward":
        return [int(value) for value in display_qualities]
    return [int(value) for value in reversed(display_qualities)]


def alignment_strands_for_policy(
    strand_policy: AlignmentStrandPolicy,
) -> tuple[AlignmentStrand, ...]:
    """Return the strand search order for one alignment policy."""
    if strand_policy == "auto":
        return ("forward", "reverse_complement")
    if strand_policy == "forward":
        return ("forward",)
    if strand_policy == "reverse_complement":
        return ("reverse_complement",)
    raise ValueError(f"Unsupported strand policy: {strand_policy}")


def select_best_oriented_alignment(
    *,
    target_sequence: str,
    display_query_sequence: str,
    display_query_qualities: list[int] | None,
    aligner: Align.PairwiseAligner,
    allowed_strands: tuple[AlignmentStrand, ...] | None = None,
) -> SelectedOrientedAlignment | None:
    """Align forward and/or reverse-complement query candidates and keep the best."""
    candidate_strands = (
        allowed_strands
        if allowed_strands is not None
        else ("forward", "reverse_complement")
    )

    best_alignment: SelectedOrientedAlignment | None = None
    for strand in candidate_strands:
        oriented_sequence = (
            display_query_sequence
            if strand == "forward"
            else reverse_complement_sequence(display_query_sequence)
        )
        alignments = aligner.align(target_sequence, oriented_sequence)
        if len(alignments) == 0:
            continue
        alignment = alignments[0]
        candidate = SelectedOrientedAlignment(
            strand=strand,
            sequence=oriented_sequence,
            qualities=oriented_qualities(display_query_qualities, strand=strand),
            alignment=alignment,
            score=float(alignment.score),
        )
        if best_alignment is None or candidate.score > best_alignment.score:
            best_alignment = candidate

    return best_alignment


def alignment_overlap_metrics(
    alignment: Any,
    *,
    target_sequence: str,
    query_sequence: str,
) -> tuple[int, float]:
    """Return overlap length and percent identity across non-gap aligned columns."""
    overlap_length = 0
    match_count = 0
    for target_index, query_index in zip(
        alignment.indices[0],
        alignment.indices[1],
        strict=True,
    ):
        resolved_target_index = int(target_index)
        resolved_query_index = int(query_index)
        if resolved_target_index < 0 or resolved_query_index < 0:
            continue
        overlap_length += 1
        if (
            target_sequence[resolved_target_index]
            == query_sequence[resolved_query_index]
        ):
            match_count += 1

    percent_identity = (match_count / overlap_length) * 100.0 if overlap_length else 0.0
    return overlap_length, percent_identity
