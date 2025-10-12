from __future__ import annotations
from dataclasses import dataclass
from typing import List, Literal
from Bio.Align import PairwiseAligner

# Keep the same external mode names for your UI.
Mode = Literal["globalms", "localms", "globalxx", "localxx"]


@dataclass
class PairwiseResult:
    a_aln: str
    b_aln: str
    score: float
    start: int
    end: int
    identity: float  # 0..1


def _build_gapped_from_alignment(aln, seqA: str, seqB: str) -> tuple[str, str]:
    """Reconstruct gapped strings from PairwiseAligner Alignment."""
    a_blocks = aln.aligned[0]  # Nx2
    b_blocks = aln.aligned[1]
    outA, outB = [], []
    ia = ib = 0
    for (sa, ea), (sb, eb) in zip(a_blocks, b_blocks):
        # gaps before the next aligned block
        if ia < sa:
            outA.append(seqA[ia:sa])  # insertion in A
            outB.append("-" * (sa - ia))  # gap in B
        if ib < sb:
            outA.append("-" * (sb - ib))  # gap in A
            outB.append(seqB[ib:sb])  # insertion in B
        # aligned block
        outA.append(seqA[sa:ea])
        outB.append(seqB[sb:eb])
        ia, ib = ea, eb
    # trailing
    if ia < len(seqA):
        outA.append(seqA[ia:])
        outB.append("-" * (len(seqA) - ia))
    if ib < len(seqB):
        outA.append("-" * (len(seqB) - ib))
        outB.append(seqB[ib:])
    return "".join(outA), "".join(outB)


def _identity(a_aln: str, b_aln: str) -> float:
    pairs = [(a, b) for a, b in zip(a_aln, b_aln) if a != "-" and b != "-"]
    if not pairs:
        return 0.0
    return sum(1 for a, b in pairs if a == b) / len(pairs)


def pairwise_align(
    seqA: str,
    seqB: str,
    mode: Mode = "globalms",
    match: float = 2.0,
    mismatch: float = -1.0,
    gap_open: float = -10.0,
    gap_extend: float = -0.5,
    top_n: int = 1,
) -> List[PairwiseResult]:
    aligner = PairwiseAligner()
    aligner.mode = "local" if mode.startswith("local") else "global"

    if mode.endswith("xx"):
        # identity-like scoring; no penalties
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = 0.0
        aligner.extend_gap_score = 0.0
    else:
        aligner.match_score = match
        aligner.mismatch_score = mismatch
        aligner.open_gap_score = gap_open
        aligner.extend_gap_score = gap_extend

    results: List[PairwiseResult] = []
    for i, aln in enumerate(aligner.align(seqA, seqB)):
        a_aln, b_aln = _build_gapped_from_alignment(aln, seqA, seqB)
        # start/end based on A's aligned coords (best-effort for display)
        start = aln.aligned[0][0][0] if len(aln.aligned[0]) else 0
        end = aln.aligned[0][-1][1] if len(aln.aligned[0]) else len(seqA)
        results.append(
            PairwiseResult(
                a_aln=a_aln,
                b_aln=b_aln,
                score=float(aln.score),
                start=start,
                end=end,
                identity=_identity(a_aln, b_aln),
            )
        )
        if len(results) >= top_n:
            break
    return results
