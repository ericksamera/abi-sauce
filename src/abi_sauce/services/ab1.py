# src/abi_sauce/services/ab1.py
from __future__ import annotations
from io import BytesIO
from typing import Dict, Any, List, Optional, Tuple

from Bio import SeqIO


def _moving_avg(xs: List[float], win: int = 5) -> List[float]:
    """Simple moving average with small window; win should be odd."""
    n = len(xs)
    if n == 0 or win <= 1:
        return xs[:]
    win = max(3, win | 1)  # ensure odd >=3
    half = win // 2
    csum = [0.0]
    for v in xs:
        csum.append(csum[-1] + float(v))
    out = [0.0] * n
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        out[i] = (csum[hi] - csum[lo]) / (hi - lo)
    return out


def _composite_signal(channels: Dict[str, Optional[List[int]]]) -> List[float]:
    """
    Build a composite intensity curve as the pointwise max across A/C/G/T.
    Handles channels of different lengths.
    """
    lengths = [len(v) for v in channels.values() if v]
    N = max(lengths) if lengths else 0
    if N == 0:
        return []
    out = [0.0] * N
    for base in "ACGT":
        arr = channels.get(base)
        if not arr:
            continue
        m = min(N, len(arr))
        for i in range(m):
            vi = float(arr[i])
            if vi > out[i]:
                out[i] = vi
    return out


def estimate_ploc(
    channels: Dict[str, Optional[List[int]]], seq_len: Optional[int] = None
) -> Optional[List[int]]:
    """
    Heuristic fallback to estimate base peak positions (PLOC) from the raw channels.

    Strategy:
      1) Build composite signal = max(A,C,G,T).
      2) Smooth lightly with a moving average to reduce wiggles.
      3) Pick local maxima above a dynamic threshold (~12% of max).
      4) If seq_len is known: split the sample axis into seq_len bins and choose
         the best peak per bin (fallback: argmax in the bin). This yields exactly
         seq_len monotonically increasing positions.
         If seq_len is unknown: just return all candidates (viewer will use len(ploc)).

    Returns None if we can't compute anything useful.
    """
    comp = _composite_signal(channels)
    if not comp:
        return None

    sm = _moving_avg(comp, win=5)
    vmax = max(sm) if sm else 0.0
    if vmax <= 0:
        return None
    thresh = max(5.0, 0.12 * vmax)  # 12% of global max intensity

    # 3) find local maxima (>= right neighbor, > left) above threshold
    N = len(sm)
    cands: List[int] = []
    for i in range(1, N - 1):
        if sm[i] > sm[i - 1] and sm[i] >= sm[i + 1] and sm[i] >= thresh:
            cands.append(i)

    if not cands:
        # fallback: just take the single global maximum
        imax = max(range(N), key=lambda t: sm[t])
        return [imax] if seq_len in (None, 1) else [imax] * (seq_len or 1)

    if not seq_len or seq_len <= 0:
        # Let the viewer use len(ploc) as seq_len if we don’t have a sequence.
        return cands

    # 4) bin the axis and choose best peak per bin (or local argmax if empty)
    peaks: List[int] = []
    for k in range(seq_len):
        start = (k * N) // seq_len
        end = ((k + 1) * N) // seq_len - 1
        if end < start:
            end = start

        # pick best candidate inside the bin
        best_i = None
        best_v = -1.0
        # fast scan since cands is sorted
        # binary search left bound
        # (manual scan is fine; arrays are small enough for AB1)
        for idx in cands:
            if idx < start:
                continue
            if idx > end:
                break
            v = sm[idx]
            if v > best_v:
                best_v = v
                best_i = idx
        if best_i is None:
            # no candidate in bin -> choose local argmax in [start,end]
            best_i = max(range(start, end + 1), key=lambda t: sm[t])
        peaks.append(best_i)

    # ensure strictly increasing (no overlaps)
    for i in range(1, len(peaks)):
        if peaks[i] <= peaks[i - 1]:
            peaks[i] = min(N - 1, max(peaks[i - 1] + 1, peaks[i]))
    return peaks


def parse_ab1(
    raw: bytes,
) -> Tuple[
    Dict[str, Optional[List[int]]],  # channels
    Optional[str],  # sequence
    Optional[List[int]],  # qualities
    Optional[List[int]],  # ploc (original or estimated)
    Dict[str, Any],  # meta flags
]:
    """
    Parse AB1 bytes with Biopython, returning (channels, sequence, phred, ploc, meta).

    If the ABIF lacks PLOC2, we estimate PLOC positions from the channel traces.
    Meta includes:
      - abif_keys: first ~200 raw ABIF keys (debug)
      - abif_order: FWO_ base order used to map DATA9..12 to A/C/G/T
      - has_phred: whether PHRED qualities were present (letter_annotations/PCON2)
      - has_ploc2: whether PLOC2 existed in the file
      - ploc_estimated: True if we synthesized PLOC positions
    """
    handle = BytesIO(raw)
    rec = SeqIO.read(handle, "abi")  # may raise on invalid file

    abif_raw: Dict[str, Any] = dict(rec.annotations.get("abif_raw", {}))

    # Base order (FWO_ tag), e.g. "GATC" -> map onto DATA9..DATA12
    fwo = abif_raw.get("FWO_", b"GATC")
    order = (
        fwo.decode(errors="ignore") if isinstance(fwo, (bytes, bytearray)) else str(fwo)
    )
    data_keys = ["DATA9", "DATA10", "DATA11", "DATA12"]

    channels: Dict[str, Optional[List[int]]] = {b: None for b in "ACGT"}
    for base, key in zip(order, data_keys):
        arr = abif_raw.get(key)
        if arr is None:
            continue
        try:
            channels[base] = list(map(int, arr))
        except Exception:
            channels[base] = None

    seq = str(rec.seq) if getattr(rec, "seq", None) else None

    # Qualities (prefer Biopython letter_annotations; fallback to PCON2)
    quals: Optional[List[int]] = None
    if hasattr(rec, "letter_annotations") and "phred_quality" in rec.letter_annotations:
        try:
            quals = list(map(int, rec.letter_annotations["phred_quality"]))
        except Exception:
            quals = None
    if quals is None:
        pcon = abif_raw.get("PCON2")
        if pcon is not None:
            try:
                quals = list(map(int, pcon))
            except Exception:
                quals = None

    # Base positions (PLOC2) or estimate if missing
    had_ploc2 = "PLOC2" in abif_raw
    ploc: Optional[List[int]] = None
    if had_ploc2:
        try:
            ploc = list(map(int, abif_raw["PLOC2"]))
        except Exception:
            ploc = None

    ploc_estimated = False
    if not ploc:
        # estimate using channels; prefer aligning count to called-seq length if present
        seq_len = len(seq) if seq else None
        est = estimate_ploc(channels, seq_len=seq_len)
        if est:
            ploc = est
            ploc_estimated = True

    meta = {
        "abif_keys": list(abif_raw.keys())[:200],
        "abif_order": order,
        "has_phred": bool(quals),
        "has_ploc2": bool(had_ploc2),  # whether file contained PLOC2 tag
        "ploc_estimated": bool(ploc_estimated),
    }
    return channels, seq, quals, ploc, meta
