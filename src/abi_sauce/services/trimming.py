from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Tuple

# ------------------------------
# Domain objects
# ------------------------------


@dataclass(frozen=True)
class TrimConfig:
    error_limit: float = 0.05  # typical
    min_len: int = 20  # typical


@dataclass(frozen=True)
class TrimResult:
    left: int  # 0-based inclusive
    right: int  # 0-based inclusive
    score: float  # Kadane/Mott score
    trimmed_seq: str  # sequence[left:right+1]

    @property
    def length(self) -> int:
        return (self.right - self.left + 1) if self.right >= self.left else 0


# ------------------------------
# Core Mott implementation
# ------------------------------


def _mott_window(phred: List[int], cfg: TrimConfig) -> Tuple[int, int, float]:
    """Return (L, R, score) for the max-scoring subarray using Mott:
    score_i = cfg.error_limit - 10**(-Q/10)."""
    if not phred:
        return 0, -1, 0.0

    # per-base scores
    scores = [cfg.error_limit - (10 ** (-q / 10.0)) for q in phred]

    # Kadane
    best = float("-inf")
    best_l, best_r = 0, -1
    cur = 0.0
    cur_l = 0
    for i, s in enumerate(scores):
        if cur <= 0:
            cur, cur_l = s, i
        else:
            cur += s
        if cur > best:
            best, best_l, best_r = cur, cur_l, i

    if best_r < best_l:
        return 0, -1, best

    # Enforce min length by symmetric expansion
    span = best_r - best_l + 1
    if span < cfg.min_len:
        need = cfg.min_len - span
        extend_left = min(best_l, need // 2)
        extend_right = min(len(phred) - 1 - best_r, need - extend_left)
        best_l -= extend_left
        best_r += extend_right

    return best_l, best_r, best


def trim_with_mott(
    sequence: str | None, phred: List[int] | None, cfg: TrimConfig
) -> Optional[TrimResult]:
    """Typed API that returns a TrimResult or None."""
    if not sequence or not phred:
        return None
    L, R, score = _mott_window(phred, cfg)
    if R < L:
        return None
    return TrimResult(left=L, right=R, score=score, trimmed_seq=sequence[L : R + 1])


# ------------------------------
# Back-compat functions (existing callers keep working)
# ------------------------------


def mott_trim_quals(
    phred: List[int],
    error_limit: float = 0.05,
    min_len: int = 20,
) -> Tuple[int, int, float]:
    cfg = TrimConfig(error_limit=error_limit, min_len=min_len)
    L, R, score = _mott_window(phred, cfg)
    return L, R, score


def trim_trace_asset_mott(
    sequence: str | None,
    phred: List[int] | None,
    *,
    error_limit: float = 0.05,
    min_len: int = 20,
) -> Optional[Tuple[int, int, str]]:
    cfg = TrimConfig(error_limit=error_limit, min_len=min_len)
    res = trim_with_mott(sequence, phred, cfg)
    if not res:
        return None
    return res.left, res.right, res.trimmed_seq
