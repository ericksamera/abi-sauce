from __future__ import annotations
from typing import List, Optional, Tuple

# --- Core Mott implementation (Phred -> error prob -> Kadane) ---


def mott_trim_quals(
    phred: List[int],
    error_limit: float = 0.05,  # typical default
    min_len: int = 20,  # keep at least this many bases
) -> Tuple[int, int, float]:
    """
    Return (left_idx, right_idx, score) for the max-scoring subarray using
    Mott's approach: score_i = error_limit - 10**(-Q/10).
    Indexes are 0-based and inclusive.
    """
    if not phred:
        return 0, -1, 0.0

    # Per-base scores
    scores = []
    for q in phred:
        p_err = 10 ** (-q / 10.0)
        scores.append(error_limit - p_err)

    # Kadane's algorithm (max subarray)
    best = float("-inf")
    best_l = 0
    best_r = -1
    cur = 0.0
    cur_l = 0
    for i, s in enumerate(scores):
        if cur <= 0:
            cur = s
            cur_l = i
        else:
            cur += s
        if cur > best:
            best = cur
            best_l = cur_l
            best_r = i

    if best_r < best_l:
        return 0, -1, best

    # Enforce minimum length by symmetric expansion as available
    span = best_r - best_l + 1
    if span < min_len:
        need = min_len - span
        extend_left = min(best_l, need // 2)
        extend_right = min(len(phred) - 1 - best_r, need - extend_left)
        best_l -= extend_left
        best_r += extend_right

    return best_l, best_r, best


def trim_trace_asset_mott(
    sequence: str | None,
    phred: List[int] | None,
    *,
    error_limit: float = 0.05,
    min_len: int = 20,
) -> Optional[Tuple[int, int, str]]:
    """
    Convenience wrapper for TraceAsset fields.
    Returns (left_idx, right_idx, trimmed_sequence) or None if not possible.
    """
    if not sequence or not phred:
        return None
    left_trim, right_trim, _ = mott_trim_quals(
        phred, error_limit=error_limit, min_len=min_len
    )
    if right_trim < left_trim:
        return None
    return left_trim, right_trim, sequence[left_trim : right_trim + 1]
